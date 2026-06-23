################################################################################
##  drr_propagate_subtad.R — Propagación de anotación topológica a DRRs
##  --------------------------------------------------------------------------
##  Este script propaga la asignación de TAD reducido desde los CRMs consenso
##  hacia las DRRs. Para cada DRR se identifica el conjunto de TADs reducidos
##  representado por sus CRMs y se resume su coherencia topológica.
##
##  Dependencias:
##    data.table y GenomicRanges.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

#' Propaga el sub-TAD de cada CRM a las DRRs que lo contienen.
#'
#' @param drr         tabla de DRRs (cmp$stacking). Usa drr_id, chr, drr_start,
#'                    drr_end, candidate_class. Se excluye Extensive_overlap.
#' @param crm_reduced CRMs reducidos (cluster_id, chr, repr_start, repr_end).
#' @param crm_tad     a6$crm_tad (crm_id, tad_id, tad_length, assignment_mode).
#' @return lista:
#'   drr_tad_long : una fila por (DRR, sub-TAD) con número de CRMs que lo soportan.
#'   drr_tad_sum  : una fila por DRR: número de sub-TADs distintos, sub-TAD dominante,
#'                  fracción de CRMs en el dominante, flag single_subtad.
#'   by_class     : resumen por clase: % de DRRs en un único sub-TAD.
propagate_subtad_to_drr <- function(drr, crm_reduced, crm_tad = NULL) {
  cand <- drr[candidate_class != "Extensive_overlap"]

  # (1) Pertenencia CRM → DRR por solape posicional
  drr_gr <- GRanges(cand$chr, IRanges(cand$drr_start, cand$drr_end),
                    drr_id = cand$drr_id, candidate_class = cand$candidate_class)
  # crm_reduced ya trae tad_id (la reducción intra-TAD lo añade): lo propagamos
  # directamente, sin cruzar con crm_tad (que usa IDs originales, no cluster_id).
  crm_gr <- GRanges(crm_reduced$chr,
                    IRanges(crm_reduced$repr_start, crm_reduced$repr_end),
                    cluster_id = crm_reduced$cluster_id,
                    tad_id     = crm_reduced$tad_id)
  ov <- findOverlaps(drr_gr, crm_gr, ignore.strand = TRUE)
  member_tad <- data.table(
    drr_id          = mcols(drr_gr)$drr_id[queryHits(ov)],
    candidate_class = mcols(drr_gr)$candidate_class[queryHits(ov)],
    cluster_id      = mcols(crm_gr)$cluster_id[subjectHits(ov)],
    tad_id          = mcols(crm_gr)$tad_id[subjectHits(ov)]
  )
  .msg("Pertenencia CRM-DRR: ", nrow(member_tad), " pares (",
       uniqueN(member_tad$drr_id), " DRRs).")

  # tad_length desde crm_reduced (un valor por tad_id, si está disponible)
  if ("tad_length" %in% names(crm_reduced)) {
    tl <- unique(as.data.table(crm_reduced)[, .(tad_id, tad_length)])
    member_tad <- merge(member_tad, tl, by = "tad_id", all.x = TRUE)
  } else {
    member_tad[, tad_length := NA_integer_]
  }

  # Conteo de CRMs por (DRR, sub-TAD)
  drr_tad_long <- member_tad[!is.na(tad_id), .(
    n_crms_in_tad = .N,
    tad_length    = tad_length[1L]
  ), by = .(drr_id, candidate_class, tad_id)][order(drr_id, -n_crms_in_tad)]

  # (3) Resumen por DRR: número de sub-TADs, dominante, fracción en dominante
  drr_tad_sum <- drr_tad_long[, {
    tot <- sum(n_crms_in_tad)
    dom_n <- n_crms_in_tad[1L]            # ya ordenado desc
    .(n_subtads          = .N,
      n_crms_assigned    = tot,
      dominant_tad_id    = tad_id[1L],
      dominant_tad_len   = tad_length[1L],
      frac_in_dominant   = round(dom_n / tot, 3),
      single_subtad      = .N == 1L)
  }, by = .(drr_id, candidate_class)]

  # Resumen por clase
  by_class <- drr_tad_sum[, .(
    n_drr                 = .N,
    pct_single_subtad     = round(100 * mean(single_subtad), 1),
    median_n_subtads      = as.numeric(median(n_subtads)),
    median_frac_dominant  = as.numeric(median(frac_in_dominant))
  ), by = candidate_class][order(-pct_single_subtad)]

  .msg("DRRs en un único sub-TAD: ",
       drr_tad_sum[single_subtad == TRUE, .N], " de ", nrow(drr_tad_sum),
       " (", round(100*mean(drr_tad_sum$single_subtad),1), "%).")

  list(drr_tad_long = drr_tad_long[], drr_tad_sum = drr_tad_sum[],
       by_class = by_class[])
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("step6_assign_crm_subtad.R")   # a6 <- annotate_crm_subtad(...)
## source("drr_propagate_subtad.R")
##
## prop <- propagate_subtad_to_drr(
##   drr         = cmp$stacking,
##   crm_reduced = red_crm$crm_reduced,
##   crm_tad     = a6$crm_tad
## )
##
## print(prop$by_class)
## # Lectura: pct_single_subtad alto en Dense_complex apoyaría la predicción de
## # la literatura (super-enhancers aislados en un único sub-dominio).
##
## head(prop$drr_tad_sum[order(-frac_in_dominant)])
## # DRRs multi-subTAD (posibles cruces de frontera):
## prop$drr_tad_sum[single_subtad == FALSE][order(-n_subtads)][1:20]
## ============================================================================
