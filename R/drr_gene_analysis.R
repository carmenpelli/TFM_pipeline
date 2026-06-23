################################################################################
##  ANÁLISIS GÉNICO Y NORMALIZACIÓN DE DRRs
##  --------------------------------------------------------------------------
##  Usa human_genes.tsv (chr,start,end,symbol,score,strand) para:
##    (1) Validación posicional: qué genes REALES caen en/cerca de cada DRR.
##    (2) Marcado de promotores estilo ROSE (TSS ± tss_window, def. 2500 pb).
##    (3) Anotación del gen más cercano (contexto).
##  Y añade:
##    (4) Normalización por tamaño de la validación SEdb (solapes SEdb / kb),
##        para separar el efecto "densidad" del efecto "tamaño".
##
##  RECORDATORIO: candidatos estructurales, NO super-enhancers confirmados.
##  Dependencias: data.table, GenomicRanges
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

## ----------------------------------------------------------------------------
## Cargar human_genes y calcular TSS según strand.
## ----------------------------------------------------------------------------
load_human_genes <- function(path = "human_genes.tsv", chr = "chr8") {
  g <- fread(path)
  setnames(g, tolower(names(g)))
  stopifnot(all(c("chr","start","end","symbol","strand") %in% names(g)))
  g <- g[chr == ..chr & !is.na(start) & !is.na(end)]
  g[, start := as.integer(start)]
  g[, end   := as.integer(end)]
  # TSS = start si strand +, end si strand - (fin 5')
  g[, tss := fifelse(strand == "+", start, end)]
  .msg("human_genes ", chr, ": ", nrow(g), " genes.")
  g[]
}

## Alias retrocompatible.
load_human_genes_chr8 <- function(path = "human_genes.tsv", chr = "chr8") {
  load_human_genes(path = path, chr = chr)
}

## ============================================================================
## (1) GENES REALES QUE CAEN EN / CERCA DE CADA DRR
## ============================================================================
#' Asigna genes reales (por posición) a cada DRR.
#'
#' @param drr        tabla de DRRs (cmp$stacking), se excluye Extensive_overlap.
#' @param genes      human_genes chr8 (de load_human_genes_chr8).
#' @param flank      pb de margen alrededor de la DRR para considerar "cerca".
#' @return data.table por DRR con genes solapantes y nº.
genes_in_drr <- function(drr, genes, flank = 0L) {
  cand <- drr[candidate_class != "Extensive_overlap"]
  drr_gr <- GRanges(cand$chr,
                    IRanges(pmax(cand$drr_start - flank, 1L), cand$drr_end + flank),
                    drr_id = cand$drr_id, candidate_class = cand$candidate_class)
  g_gr   <- GRanges(genes$chr, IRanges(genes$start, genes$end),
                    symbol = genes$symbol)
  ov <- findOverlaps(drr_gr, g_gr)
  dt <- data.table(
    drr_id          = mcols(drr_gr)$drr_id[queryHits(ov)],
    candidate_class = mcols(drr_gr)$candidate_class[queryHits(ov)],
    symbol          = mcols(g_gr)$symbol[subjectHits(ov)]
  )
  per_drr <- dt[, .(
    n_genes_real = uniqueN(symbol),
    genes_real   = paste(sort(unique(symbol)), collapse = ";")
  ), by = .(drr_id, candidate_class)]
  per_drr[]
}

## ============================================================================
## (2) MARCADO DE PROMOTORES (estilo ROSE: TSS ± tss_window)
## ============================================================================
#' Marca DRRs que solapan una ventana de promotor alrededor de algún TSS.
#'
#' @param drr        tabla de DRRs.
#' @param genes      human_genes con columna tss.
#' @param tss_window pb a cada lado del TSS (def. 2500, como ROSE).
#' @return data.table drr_id, overlaps_promoter (lógico), n_promoters.
flag_promoter_drr <- function(drr, genes, tss_window = 2500L) {
  cand <- drr[candidate_class != "Extensive_overlap"]
  drr_gr <- GRanges(cand$chr, IRanges(cand$drr_start, cand$drr_end),
                    drr_id = cand$drr_id)
  prom_gr <- GRanges(genes$chr,
                     IRanges(pmax(genes$tss - tss_window, 1L),
                             genes$tss + tss_window),
                     symbol = genes$symbol)
  ov <- findOverlaps(drr_gr, prom_gr)
  prom <- data.table(
    drr_id = mcols(drr_gr)$drr_id[queryHits(ov)],
    symbol = mcols(prom_gr)$symbol[subjectHits(ov)]
  )
  per_drr <- prom[, .(n_promoters = uniqueN(symbol),
                      promoter_genes = paste(sort(unique(symbol)), collapse=";")),
                  by = drr_id]
  # Todas las candidatas, con flag
  out <- merge(data.table(drr_id = cand$drr_id), per_drr, by = "drr_id", all.x = TRUE)
  out[is.na(n_promoters), n_promoters := 0L]
  out[is.na(promoter_genes), promoter_genes := ""]
  out[, overlaps_promoter := n_promoters > 0L]
  .msg("DRRs que solapan promotor (TSS±", tss_window, "): ",
       out[overlaps_promoter == TRUE, .N], " de ", nrow(out), ".")
  out[]
}

## ============================================================================
## (3) GEN MÁS CERCANO A CADA DRR (contexto)
## ============================================================================
#' Asigna a cada DRR el gen con TSS más cercano (distancia con signo no; absoluta).
nearest_gene_drr <- function(drr, genes) {
  cand <- drr[candidate_class != "Extensive_overlap"]
  drr_gr <- GRanges(cand$chr,
                    IRanges((cand$drr_start + cand$drr_end) %/% 2L,
                            (cand$drr_start + cand$drr_end) %/% 2L),
                    drr_id = cand$drr_id)
  tss_gr <- GRanges(genes$chr, IRanges(genes$tss, genes$tss), symbol = genes$symbol)
  nr <- distanceToNearest(drr_gr, tss_gr)
  data.table(
    drr_id        = mcols(drr_gr)$drr_id[queryHits(nr)],
    nearest_gene  = mcols(tss_gr)$symbol[subjectHits(nr)],
    distance_bp   = mcols(nr)$distance
  )[]
}

## ============================================================================
## (4) NORMALIZACIÓN POR TAMAÑO DE LA VALIDACIÓN SEdb
## ============================================================================
#' Recalcula el solapamiento SEdb normalizado por longitud de DRR (solapes/kb).
#'
#' @param overlap salida de validate_drr_vs_sedb()$overlap (drr_id, candidate_class, se_id).
#' @param drr     tabla de DRRs (para obtener drr_length por drr_id).
#' @return lista: per_drr (con se_per_kb) y summary por clase.
sedb_overlap_normalized <- function(overlap, drr) {
  cand <- drr[candidate_class != "Extensive_overlap",
              .(drr_id, drr_length)]
  # nº de SE distintos por DRR
  per_drr <- overlap[, .(n_SE = uniqueN(se_id)), by = .(drr_id, candidate_class)]
  per_drr <- merge(per_drr, cand, by = "drr_id", all.x = TRUE)
  per_drr[, se_per_kb := n_SE / (drr_length / 1000)]

  summ <- per_drr[, .(
    n_candidates       = .N,
    median_SE          = median(n_SE),
    median_SE_per_kb   = median(se_per_kb),
    mean_SE_per_kb     = mean(se_per_kb)
  ), by = candidate_class]

  base <- summ[candidate_class == "Simple_DRR", median_SE_per_kb]
  if (length(base) == 1L && base > 0) {
    summ[, enrichment_per_kb_vs_simple := round(median_SE_per_kb / base, 2)]
  } else {
    summ[, enrichment_per_kb_vs_simple := NA_real_]
  }
  setorder(summ, -enrichment_per_kb_vs_simple)
  list(per_drr = per_drr[], summary = summ[])
}

## ============================================================================
## EJEMPLO DE USO
## ----------------------------------------------------------------------------
## source("drr_gene_analysis.R")
## genes <- load_human_genes_chr8("human_genes.tsv")
##
## drr <- cmp$stacking
##
## # (1) Genes reales en cada DRR:
## g_in   <- genes_in_drr(drr, genes, flank = 0L)
## # resumen por clase:
## g_in[, .(median_genes = median(n_genes_real), .N), by = candidate_class]
##
## # (2) Marcado de promotores:
## prom   <- flag_promoter_drr(drr, genes, tss_window = 2500L)
## prom[, mean(overlaps_promoter), by = drr_id][, mean(V1)]   # fracción global
##
## # (3) Gen más cercano:
## near   <- nearest_gene_drr(drr, genes)
##
## # (4) Normalización por tamaño de la validación SEdb:
## #     (usa el overlap de la validación previa: val$positional$overlap)
## norm   <- sedb_overlap_normalized(val$positional$overlap, drr)
## print(norm$summary)   # ¿Dense_complex sigue enriquecido tras normalizar por kb?
## ============================================================================
