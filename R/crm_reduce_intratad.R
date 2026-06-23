################################################################################
##  crm_reduce_intratad.R — Reducción intra-TAD de CRMs
##  --------------------------------------------------------------------------
##  Este script aplica la reducción de redundancia de CRMs dentro de cada TAD
##  reducido. Cada CRM se asigna previamente a un TAD y se compara únicamente
##  con CRMs pertenecientes al mismo dominio topológico.
##
##  Este enfoque mantiene la coherencia topológica del análisis y divide el
##  cálculo global de solapamientos en subproblemas por TAD.
##
##  Dependencias:
##    data.table y funciones definidas en crm_explore.R,
##    crm_antichaining_fast.R y crm_reduce_final.R.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

#' Reduce la redundancia de CRMs dentro de cada TAD.
#'
#' @param crm_unique CRMs colapsados por ID (ID, chr, start, end, ...).
#' @param crm_tad    tabla crm_id → tad_id (de annotate_crm_subtad: a6$crm_tad).
#' @param crm_id_in_unique nombre de la columna ID en crm_unique (por defecto "ID").
#' @param crm_id_in_tad    nombre de la columna id en crm_tad (por defecto "crm_id").
#' @param thr_recip,jaccard_thr,simpson_thr criterio (idéntico al global).
#' @param min_per_tad número mínimo de CRMs en un TAD para intentar reducir
#'        (TADs con < min_per_tad CRMs se mantienen sin reducción, sin comparar).
#' @param ... se pasa a reduce_redundancy_crms (meta_cols, sum_cols, chunk_size...).
#' @return lista:
#'   crm_reduced : todos los CRMs reducidos (con tad_id y cluster_id global único).
#'   per_tad     : data.table con n_in / n_reduced por TAD.
#'   summary     : resumen global.
reduce_redundancy_crms_intratad <- function(crm_unique, crm_tad,
                                           crm_id_in_unique = "ID",
                                           crm_id_in_tad    = "crm_id",
                                           thr_recip   = 0.50,
                                           jaccard_thr = 0.70,
                                           simpson_thr = 0.99,
                                           min_per_tad = 2L,
                                           ...) {

  if (!exists("reduce_redundancy_crms"))
    stop("Carga crm_reduce_final.R (y sus dependencias) antes.")

  crm <- as.data.table(copy(crm_unique))
  ct  <- as.data.table(copy(crm_tad))
  stopifnot(crm_id_in_unique %in% names(crm))
  stopifnot(all(c(crm_id_in_tad, "tad_id") %in% names(ct)))

  # Unir el tad_id a cada CRM
  setnames(ct, crm_id_in_tad, ".cid")
  crm[, .cid := get(crm_id_in_unique)]
  crm <- merge(crm, ct[, .(.cid, tad_id)], by = ".cid", all.x = TRUE)

  n_total   <- nrow(crm)
  n_no_tad  <- crm[is.na(tad_id), .N]
  .msg("Reducción intra-TAD: ", n_total, " CRMs, ",
       crm[!is.na(tad_id), uniqueN(tad_id)], " TADs. ",
       n_no_tad, " CRMs sin TAD (se omiten).")

  work <- crm[!is.na(tad_id)]
  tad_ids <- work[, unique(tad_id)]

  out_list   <- vector("list", length(tad_ids))
  per_tad    <- vector("list", length(tad_ids))
  map_list   <- vector("list", length(tad_ids))  # trazabilidad original_id->cluster
  cluster_off <- 0L   # offset para hacer cluster_id únicos entre TADs

  for (k in seq_along(tad_ids)) {
    tid <- tad_ids[k]
    sub <- work[tad_id == tid]
    n_in <- nrow(sub)

    if (n_in < min_per_tad) {
      # TAD con pocos CRMs: se mantiene cada CRM como clúster independiente.
      red_sub <- copy(sub)
      red_sub[, cluster_id := cluster_off + seq_len(n_in)]
      red_sub[, n_entities := 1L]
      out_list[[k]] <- red_sub[, .(cluster_id, tad_id,
                                   repr_start = start, repr_end = end,
                                   chr, n_entities)]
      map_list[[k]] <- red_sub[, .(original_id = ID, cluster_id, tad_id)]
      per_tad[[k]] <- data.table(tad_id = tid, n_in = n_in, n_reduced = n_in)
      cluster_off <- cluster_off + n_in
      next
    }

    # Reducción global PERo sobre el subconjunto = intra-TAD
    red <- tryCatch(
      reduce_redundancy_crms(sub, thr_recip = thr_recip,
                             jaccard_thr = jaccard_thr,
                             simpson_thr = simpson_thr,
                             verbose_ac = FALSE, ...),
      error = function(e) { .msg("  TAD ", tid, " error: ", conditionMessage(e)); NULL }
    )
    if (is.null(red)) {  # fallback: sin reducir
      red_sub <- copy(sub)
      red_sub[, cluster_id := cluster_off + seq_len(n_in)]
      red_sub[, n_entities := 1L]
      out_list[[k]] <- red_sub[, .(cluster_id, tad_id,
                                   repr_start = start, repr_end = end,
                                   chr, n_entities)]
      map_list[[k]] <- red_sub[, .(original_id = ID, cluster_id, tad_id)]
      per_tad[[k]] <- data.table(tad_id = tid, n_in = n_in, n_reduced = n_in)
      cluster_off <- cluster_off + n_in
      next
    }

    rr <- as.data.table(red$crm_reduced)
    # cluster_id local → global único con offset
    rr[, cluster_id := cluster_id + cluster_off]
    rr[, tad_id := tid]
    # Coordenada representativa: usar recommended_coord si existe; si no, union.
    if (!"repr_start" %in% names(rr)) {
      rr[, repr_start := union_start]
      rr[, repr_end   := union_end]
    }
    out_list[[k]] <- rr[, .(cluster_id, tad_id, repr_start, repr_end,
                            chr, n_entities)]
    # Mapping de trazabilidad: original_id → cluster_id (con offset y tad_id)
    if (!is.null(red$mapping)) {
      mp <- as.data.table(red$mapping)[, .(original_id, cluster_id)]
      mp[, cluster_id := cluster_id + cluster_off]
      mp[, tad_id := tid]
      map_list[[k]] <- mp
    }
    per_tad[[k]] <- data.table(tad_id = tid, n_in = n_in, n_reduced = nrow(rr))
    cluster_off <- cluster_off + nrow(rr)

    if (k %% 100 == 0L) .msg("  ", k, "/", length(tad_ids), " TADs procesados.")
  }

  crm_reduced <- rbindlist(out_list, use.names = TRUE, fill = TRUE)
  per_tad_dt  <- rbindlist(per_tad)
  mapping     <- rbindlist(map_list, use.names = TRUE, fill = TRUE)

  summary_dt <- data.table(
    n_crm_in        = n_total,
    n_crm_no_tad    = n_no_tad,
    n_crm_processed = nrow(work),
    n_reduced       = nrow(crm_reduced),
    n_tads          = length(tad_ids),
    reduction_ratio = round(1 - nrow(crm_reduced) / nrow(work), 4)
  )
  .msg("Intra-TAD OK: ", nrow(work), " -> ", nrow(crm_reduced),
       " (reducción ", round(100 * summary_dt$reduction_ratio, 1), "%).")

  list(crm_reduced = crm_reduced[], mapping = mapping[],
       per_tad = per_tad_dt[], summary = summary_dt[])
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("crm_explore.R"); source("crm_antichaining_fast.R")
## source("crm_reduce_final.R"); source("crm_reduce_intratad.R")
##
## # P6 primero, sobre los CRMs colapsados por ID (enh_unique, 1.343.498):
## a6 <- annotate_crm_subtad(crm_reduced = enh_unique,   # entrada procedente de P1
##                           tad_reduced = red_tad$tad_reduced,
##                           crm_id_col = "ID", start_col = "start", end_col = "end")
## # (annotate_crm_subtad usa repr_* por defecto; pásale start/end de enh_unique)
##
## red_crm_it <- reduce_redundancy_crms_intratad(
##   crm_unique = enh_unique, crm_tad = a6$crm_tad,
##   crm_id_in_unique = "ID", crm_id_in_tad = "crm_id"
## )
## print(red_crm_it$summary)
## red_crm_it$per_tad[order(-n_in)][1:10]   # los TADs con más CRMs
## ============================================================================
