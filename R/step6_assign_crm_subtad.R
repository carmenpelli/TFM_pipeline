################################################################################
##  PASO 6 (revisado) — ANOTACIÓN CRM -> SUB-TAD MÁS ESPECÍFICO
##  --------------------------------------------------------------------------
##  Cambia el criterio anterior ("máximo solape en bases", que sesgaba hacia
##  TADs gigantes) por el criterio biológicamente respaldado:
##
##    Asignar cada CRM al TAD MÁS PEQUEÑO que lo contiene (el sub-dominio más
##    específico). La literatura sitúa el aislamiento de super-enhancers en
##    sub-TADs / dominios de super-enhancer, no en los TADs grandes.
##
##  Modo de pertenencia:
##    - "contained" (preferido): el TAD cubre el CRM por completo.
##    - si ningún TAD lo contiene del todo -> respaldo "overlap": entre los TADs
##      que solapan el CRM, se elige el más pequeño (para no perder CRMs).
##    Los CRMs sin ningún TAD (ni contenedor ni solapante) se descartan.
##
##  Coordenadas: repr_start/repr_end en CRMs y TADs reducidos.
##
##  Dependencias: data.table, GenomicRanges
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

#' Anota cada CRM con el TAD más pequeño que lo contiene (o que lo solapa).
#'
#' @param crm_reduced data.table de CRMs reducidos (cluster_id, chr, repr_*).
#' @param tad_reduced data.table de TADs reducidos (cluster_id, chr, repr_*).
#' @param crm_id_col,tad_id_col columnas id.
#' @param start_col,end_col coordenadas (def. repr_start/repr_end).
#' @return lista:
#'   crm_tad    : crm_id, tad_id, tad_length, assignment_mode ("contained"/"overlap")
#'   unassigned : CRMs sin ningún TAD
#'   summary    : resumen
annotate_crm_subtad <- function(crm_reduced, tad_reduced,
                                crm_id_col   = "cluster_id",
                                tad_id_col   = "cluster_id",
                                crm_start_col = "repr_start",
                                crm_end_col   = "repr_end",
                                tad_start_col = "repr_start",
                                tad_end_col   = "repr_end") {

  crm <- as.data.table(copy(crm_reduced))
  tad <- as.data.table(copy(tad_reduced))

  # Uso de [[ ]] (no get()) para evitar que un nombre como "start" se resuelva
  # a la función start() de IRanges cuando la columna no existe.
  crm[, .crm_id := crm[[crm_id_col]]]
  crm[, .s := crm[[crm_start_col]]]; crm[, .e := crm[[crm_end_col]]]
  tad[, .tad_id := tad[[tad_id_col]]]
  tad[, .s := tad[[tad_start_col]]]; tad[, .e := tad[[tad_end_col]]]
  tad[, .tad_len := .e - .s + 1L]

  # Comprobación explícita: si alguna coordenada salió no-numérica, avisar claro
  if (!is.numeric(crm$.s) || !is.numeric(crm$.e))
    stop("Columnas de coordenadas de CRM no numéricas: revisa crm_start_col/crm_end_col (",
         crm_start_col, "/", crm_end_col, ") en crm_reduced.")
  if (!is.numeric(tad$.s) || !is.numeric(tad$.e))
    stop("Columnas de coordenadas de TAD no numéricas: revisa tad_start_col/tad_end_col (",
         tad_start_col, "/", tad_end_col, ") en tad_reduced.")

  n_crm_in <- nrow(crm)
  .msg("PASO 6 (sub-TAD): ", n_crm_in, " CRMs vs ", nrow(tad), " TADs. ",
       "Criterio: TAD más pequeño que contiene (respaldo: solapa).")

  gr_crm <- GRanges(crm$chr, IRanges(crm$.s, crm$.e))
  gr_tad <- GRanges(tad$chr, IRanges(tad$.s, tad$.e))

  ## --- (1) Contención estricta: TAD que cubre el CRM por completo ----------
  ## findOverlaps con type="within": query (CRM) dentro de subject (TAD).
  hits_in <- findOverlaps(gr_crm, gr_tad, type = "within", ignore.strand = TRUE)
  contained <- data.table(
    crm_row = queryHits(hits_in),
    crm_id  = crm$.crm_id[queryHits(hits_in)],
    tad_id  = tad$.tad_id[subjectHits(hits_in)],
    tad_len = tad$.tad_len[subjectHits(hits_in)]
  )
  # Entre los TADs que CONTIENEN el CRM, el más pequeño (desempate por tad_id)
  setorder(contained, crm_id, tad_len, tad_id)
  contained_best <- contained[, .SD[1L], by = crm_id]
  contained_best[, assignment_mode := "contained"]

  ## --- (2) Respaldo por solape para los CRMs no contenidos -----------------
  not_contained <- setdiff(crm$.crm_id, contained_best$crm_id)
  overlap_best <- NULL
  if (length(not_contained) > 0L) {
    idx_nc <- which(crm$.crm_id %in% not_contained)
    hits_ov <- findOverlaps(gr_crm[idx_nc], gr_tad, ignore.strand = TRUE)
    if (length(hits_ov) > 0L) {
      ov <- data.table(
        crm_id  = crm$.crm_id[idx_nc][queryHits(hits_ov)],
        tad_id  = tad$.tad_id[subjectHits(hits_ov)],
        tad_len = tad$.tad_len[subjectHits(hits_ov)]
      )
      # Entre los TADs que SOLAPAN, el más pequeño
      setorder(ov, crm_id, tad_len, tad_id)
      overlap_best <- ov[, .SD[1L], by = crm_id]
      overlap_best[, assignment_mode := "overlap"]
    }
  }

  ## --- Combinar -------------------------------------------------------------
  crm_tad <- rbind(
    contained_best[, .(crm_id, tad_id, tad_length = tad_len, assignment_mode)],
    if (!is.null(overlap_best))
      overlap_best[, .(crm_id, tad_id, tad_length = tad_len, assignment_mode)]
  )

  assigned_ids <- unique(crm_tad$crm_id)
  unassigned <- crm[!(.crm_id %in% assigned_ids), .(crm_id = .crm_id)]

  summary_dt <- data.table(
    n_crm_in        = n_crm_in,
    n_assigned      = nrow(crm_tad),
    n_contained     = crm_tad[assignment_mode == "contained", .N],
    n_overlap_only  = crm_tad[assignment_mode == "overlap", .N],
    n_unassigned    = nrow(unassigned),
    frac_assigned   = round(nrow(crm_tad) / n_crm_in, 4),
    n_tads_used     = uniqueN(crm_tad$tad_id),
    median_tad_len_used = as.numeric(median(crm_tad$tad_length))
  )
  .msg("PASO 6 OK: ", nrow(crm_tad), " asignados (",
       summary_dt$n_contained, " contenidos, ", summary_dt$n_overlap_only,
       " solo solape), ", nrow(unassigned), " descartados. TADs usados: ",
       summary_dt$n_tads_used, ".")

  list(crm_tad = crm_tad[], unassigned = unassigned[], summary = summary_dt[])
}

## ============================================================================
## EJEMPLO DE USO
## ----------------------------------------------------------------------------
## source("step6_assign_crm_subtad.R")
##
## # CASO A — P6 sobre CRMs colapsados por ID (enh_unique: ID, chr, start, end):
## a6 <- annotate_crm_subtad(
##   crm_reduced  = enh_unique,            # los de P1 (1.343.498)
##   tad_reduced  = red_tad$tad_reduced,
##   crm_id_col   = "ID",
##   crm_start_col = "start", crm_end_col = "end",   # enh_unique usa start/end
##   tad_start_col = "repr_start", tad_end_col = "repr_end"  # TAD usa repr_*
## )
##
## # CASO B — P6 sobre CRMs reducidos (repr_* en ambas tablas): valores por defecto
## # a6 <- annotate_crm_subtad(red_crm$crm_reduced, red_tad$tad_reduced)
##
## print(a6$summary)
## a6$crm_tad[, .N, by = assignment_mode]
## ============================================================================
