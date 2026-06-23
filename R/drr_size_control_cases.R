################################################################################
##  drr_size_control_cases.R — Control por tamaño y casos destacados DRR-SEdb
##  --------------------------------------------------------------------------
##  Este script evalúa la concordancia entre DRRs y SEdb controlando por longitud
##  de las regiones, y genera tablas de casos destacados según criterios de
##  solapamiento, recurrencia en contextos celulares y genes asociados.
##
##  El análisis permite separar parcialmente el efecto de la longitud de la DRR
##  del efecto de la densidad estructural.
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

## ============================================================================
## (A) CONTROL EMPAREJADo POR TAMAÑO
## ============================================================================
#' Compara solapamiento DRR–SEdb por clase, dentro de bins de longitud.
#'
#' @param drr     tabla DRRs (cmp$stacking); se excluye Extensive_overlap.
#' @param overlap val$positional$overlap (drr_id, candidate_class, se_id).
#' @param breaks  cortes de longitud (pb) para los bins.
#' @return lista: per_drr (con bin y número SE), by_bin_class (mediana SE por bin×clase),
#'         contrast (Dense_complex vs Simple dentro de cada bin).
sedb_size_matched <- function(drr, overlap,
                              breaks = c(0, 1000, 2500, 5000, 10000, 25000,
                                         50000, 100000)) {
  cand <- drr[candidate_class != "Extensive_overlap", .(drr_id, candidate_class, drr_length)]

  # número de SE distintos por DRR (0 si no solapa)
  se_count <- overlap[, .(n_SE = uniqueN(se_id)), by = drr_id]
  cand <- merge(cand, se_count, by = "drr_id", all.x = TRUE)
  cand[is.na(n_SE), n_SE := 0L]

  # Bin de longitud
  cand[, size_bin := cut(drr_length, breaks = breaks, include.lowest = TRUE,
                         dig.lab = 8)]

  # Mediana y media de SE por bin × clase
  # as.numeric() en las medianas evita el error de data.table de tipos
  # inconsistentes entre grupos (integer en unos, double en otros).
  by_bin_class <- cand[, .(
    n_drr      = .N,
    median_SE  = as.numeric(median(n_SE)),
    mean_SE    = as.numeric(mean(n_SE)),
    frac_with_SE = mean(n_SE > 0)
  ), by = .(size_bin, candidate_class)][order(size_bin, candidate_class)]

  # Contraste directo dentro de cada bin de tamaño.
  # Se utiliza Dense_complex frente a Extended_complex, ya que a tamaños grandes
  # Simple_DRR está escasamente representada, mientras que Dense y Extended coexisten en los
  # bins y se diferencian principalmente en densidad. Este contraste ayuda a aislar
  # el efecto de densidad frente al efecto de tamaño.
  wide <- dcast(by_bin_class, size_bin ~ candidate_class,
                value.var = "median_SE")
  contrast <- wide[, .(size_bin)]
  if ("Dense_complex_DRR" %in% names(wide) && "Extended_complex_DRR" %in% names(wide)) {
    contrast[, dense_median    := wide$Dense_complex_DRR]
    contrast[, extended_median := wide$Extended_complex_DRR]
    contrast[, dense_vs_extended_median :=
               wide$Dense_complex_DRR / wide$Extended_complex_DRR]
  }
  # También la fracción que solapa con algún SE (robusta a outliers de conteo)
  wfrac <- dcast(by_bin_class, size_bin ~ candidate_class,
                 value.var = "frac_with_SE")
  if ("Dense_complex_DRR" %in% names(wfrac) && "Extended_complex_DRR" %in% names(wfrac)) {
    contrast[, dense_frac_with_SE    := wfrac$Dense_complex_DRR]
    contrast[, extended_frac_with_SE := wfrac$Extended_complex_DRR]
  }

  list(per_drr = cand[], by_bin_class = by_bin_class[], contrast = contrast[])
}

## ============================================================================
## (B) CASOS DESTACADOS DRR–SEdb
## ============================================================================
#' Construye una tabla de coincidencias DRR–SE destacadas (3 criterios).
#'
#' @param drr        tabla DRRs (cmp$stacking).
#' @param se_chr8    SEdb chr8 (de load_sedb_chr8): se_id, se_chr, se_start, se_end,
#'                   cell_id, se_rank, se_gene_* ...
#' @param genes_real opcional: salida de genes_in_drr() para añadir genes reales.
#' @param top_n      número de casos por criterio.
#' @return lista de data.tables: reciprocal, multicell, gene_centric.
highlight_drr_se_cases <- function(drr, se_chr8, genes_real = NULL, top_n = 20) {
  cand <- drr[candidate_class != "Extensive_overlap"]

  drr_gr <- GRanges(cand$chr, IRanges(cand$drr_start, cand$drr_end),
                    drr_id = cand$drr_id, candidate_class = cand$candidate_class,
                    drr_length = cand$drr_length)
  se_gr  <- GRanges(se_chr8$se_chr, IRanges(se_chr8$se_start, se_chr8$se_end))

  hits <- findOverlaps(drr_gr, se_gr)
  qi <- queryHits(hits); si <- subjectHits(hits)

  # Solape en bases y recíproco
  ov <- pmin(end(drr_gr)[qi], end(se_gr)[si]) -
        pmax(start(drr_gr)[qi], start(se_gr)[si]) + 1L
  ov[ov < 0L] <- 0L
  len_drr <- width(drr_gr)[qi]
  len_se  <- width(se_gr)[si]

  pairs <- data.table(
    drr_id          = mcols(drr_gr)$drr_id[qi],
    candidate_class = mcols(drr_gr)$candidate_class[qi],
    drr_length      = mcols(drr_gr)$drr_length[qi],
    se_id           = se_chr8$se_id[si],
    se_start        = se_chr8$se_start[si],
    se_end          = se_chr8$se_end[si],
    cell_id         = if ("cell_id" %in% names(se_chr8)) se_chr8$cell_id[si] else NA,
    se_rank         = if ("se_rank" %in% names(se_chr8)) se_chr8$se_rank[si] else NA,
    se_gene_closest = if ("se_gene_closest" %in% names(se_chr8)) se_chr8$se_gene_closest[si] else NA,
    overlap_bp      = ov,
    reciprocal      = pmin(ov / len_drr, ov / len_se)
  )

  ## --- Criterio 1: solape recíproco alto (coincidencia ajustada) -----------
  reciprocal <- pairs[order(-reciprocal)][
    , head(.SD, 1), by = drr_id   # mejor SE por DRR
  ][order(-reciprocal)][seq_len(min(top_n, .N))]

  ## --- Criterio 2: SE en muchos contextos celulares ------------------------
  ## Nº de contextos celulares distintos del SE region (mismo se_*coords) que
  ## solapan cada DRR. Aproximamos por número de cell_id distintos entre los SE
  ## que solapan cada DRR.
  multicell <- pairs[, .(
    candidate_class = candidate_class[1L],
    n_cell_contexts = uniqueN(cell_id),
    n_SE            = uniqueN(se_id),
    example_gene    = na.omit(se_gene_closest)[1L]
  ), by = drr_id][order(-n_cell_contexts, -n_SE)][seq_len(min(top_n, .N))]

  ## --- Criterio 3: centrados en genes conocidos ----------------------------
  ## DRRs cuyo SE solapante está anotado a un gen (se_gene_closest no vacío),
  ## priorizando los de mayor recíproco.
  gene_centric <- pairs[!is.na(se_gene_closest) & se_gene_closest != "" &
                          se_gene_closest != "-"][
    order(-reciprocal)][, head(.SD, 1), by = drr_id][
    order(-reciprocal)][seq_len(min(top_n, .N))]

  # Añadir genes reales (human_genes) si se aportan
  if (!is.null(genes_real)) {
    gr <- as.data.table(genes_real)[, .(drr_id, genes_real)]
    reciprocal   <- merge(reciprocal,   gr, by = "drr_id", all.x = TRUE)
    gene_centric <- merge(gene_centric, gr, by = "drr_id", all.x = TRUE)
  }

  list(reciprocal = reciprocal[], multicell = multicell[],
       gene_centric = gene_centric[])
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("drr_sedb_validation.R")     # load_sedb_chr8, validate_drr_vs_sedb
## source("drr_gene_analysis.R")       # genes_in_drr
## source("drr_size_control_cases.R")  # este archivo
##
## se_chr8 <- load_sedb_chr8("~/TFM/SEdb_Human_SE.bed")
## drr     <- cmp$stacking
##
## ## (A) Control por tamaño:
## sm <- sedb_size_matched(drr, val$positional$overlap)
## print(sm$by_bin_class)   # mediana de SE por bin de longitud × clase
## print(sm$contrast)       # Dense_complex / Simple dentro de cada bin
## # Interpretación: valores superiores a 1 indican mayor concordancia en la
## # aporta señal independiente del tamaño. Si ~1 o <1, el efecto era tamaño.
##
## ## (B) Casos destacados:
## g_in  <- genes_in_drr(drr, load_human_genes_chr8("~/TFM/human_genes.tsv"))
## cases <- highlight_drr_se_cases(drr, se_chr8, genes_real = g_in, top_n = 20)
## print(cases$reciprocal)     # coincidencias ajustadas
## print(cases$multicell)      # reproducibles (muchos contextos)
## print(cases$gene_centric)   # centradas en genes
## ============================================================================
