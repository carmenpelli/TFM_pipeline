################################################################################
##  VALIDACIÓN DE DRRs (criterio B, apilamiento) CONTRA SEdb
##  --------------------------------------------------------------------------
##  Adapta la validación SEdb antigua a la salida del pipeline nuevo:
##    - DRRs: cmp$stacking (criterio apilamiento), excluyendo Extensive_overlap.
##    - CRMs reducidos: red_crm$crm_reduced  (cluster_id, repr_*, n_entities...)
##    - Trazabilidad CRM original->cluster: red_crm$mapping (original_id, cluster_id)
##    - enh2gene: relación crm original -> gen
##
##  Valida dos cosas:
##    (1) Solapamiento posicional DRR vs super-enhancers de SEdb, por clase.
##        Hipótesis: Dense_complex se enriquece frente a Simple.
##    (2) Recuperación génica: genes ligados a nuestras DRRs vs genes SEdb.
##
##  RECORDATORIO (documento): son candidatos estructurales / regiones densas,
##  NUNCA super-enhancers confirmados. SEdb sirve como validación cruzada
##  posicional, no como verdad biológica absoluta.
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

## ----------------------------------------------------------------------------
## 1) Cargar y limpiar SEdb (chr8). Reutiliza tu lógica de limpieza de comillas.
## ----------------------------------------------------------------------------
load_sedb <- function(bed_path = "SEdb_Human_SE.bed", chr = "chr8") {
  se <- fread(bed_path, sep = "\t", header = TRUE, quote = "", fill = TRUE)
  setnames(se, gsub('"', '', names(se)))
  char_cols <- names(se)[sapply(se, is.character)]
  se[, (char_cols) := lapply(.SD, function(x) gsub('"', '', x)), .SDcols = char_cols]

  num_cols <- intersect(c("se_start","se_end","se_rank"), names(se))
  se[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

  se_chr <- se[se_chr == chr]
  if (nrow(se_chr) == 0L) {
    .msg("AVISO: 0 super-enhancers de SEdb para '", chr, "'. ",
         "Cromosomas disponibles en el BED: ",
         paste(head(sort(unique(se$se_chr)), 30), collapse = ", "),
         ". Revisa la nomenclatura (¿'chr8' vs '8'?).")
  }
  .msg("SEdb ", chr, ": ", nrow(se_chr), " super-enhancers cargados.")
  se_chr[]
}

## Alias retrocompatible: el código antiguo llamaba load_sedb_chr8(bed_path).
## Por defecto mantiene chr8, pero ahora acepta 'chr' para no quedar fijado.
load_sedb_chr8 <- function(bed_path = "SEdb_Human_SE.bed", chr = "chr8") {
  load_sedb(bed_path = bed_path, chr = chr)
}

## ----------------------------------------------------------------------------
## 2) Validación POSICIONAL: solapamiento DRR vs SEdb, enriquecimiento por clase.
## ----------------------------------------------------------------------------
validate_drr_vs_sedb <- function(drr, se_chr8) {
  # Sólo candidatas (excluye Extensive_overlap)
  cand <- drr[candidate_class != "Extensive_overlap"]

  drr_gr <- GRanges(cand$chr, IRanges(cand$drr_start, cand$drr_end),
                    drr_id = cand$drr_id, candidate_class = cand$candidate_class)
  se_gr  <- GRanges(se_chr8$se_chr, IRanges(se_chr8$se_start, se_chr8$se_end),
                    se_id = se_chr8$se_id,
                    cell_id = if ("cell_id" %in% names(se_chr8)) se_chr8$cell_id else NA,
                    se_rank = if ("se_rank" %in% names(se_chr8)) se_chr8$se_rank else NA)

  hits <- findOverlaps(drr_gr, se_gr)
  overlap <- data.table(
    drr_id          = mcols(drr_gr)$drr_id[queryHits(hits)],
    candidate_class = mcols(drr_gr)$candidate_class[queryHits(hits)],
    se_id           = mcols(se_gr)$se_id[subjectHits(hits)],
    cell_id         = mcols(se_gr)$cell_id[subjectHits(hits)],
    se_rank         = mcols(se_gr)$se_rank[subjectHits(hits)]
  )

  # Totales por clase (denominador)
  class_totals <- cand[, .(n_total = .N), by = candidate_class]

  # DRRs con al menos un solape SEdb
  with_ovl <- overlap[, .(n_with_ovl = uniqueN(drr_id)), by = candidate_class]

  by_class <- merge(class_totals, with_ovl, by = "candidate_class", all.x = TRUE)
  by_class[is.na(n_with_ovl), n_with_ovl := 0L]
  by_class[, fraction_with_SE := n_with_ovl / n_total]

  # Intensidad: nº de SE distintos y contextos celulares por DRR
  intensity <- overlap[, .(
    n_SE_overlaps   = uniqueN(se_id),
    n_cell_contexts = uniqueN(cell_id),
    median_se_rank  = median(se_rank, na.rm = TRUE)
  ), by = .(drr_id, candidate_class)]

  intensity_summary <- intensity[, .(
    n_candidates               = .N,
    median_SE_overlaps_per_DRR = median(n_SE_overlaps),
    mean_SE_overlaps_per_DRR   = mean(n_SE_overlaps),
    median_cell_contexts       = median(n_cell_contexts)
  ), by = candidate_class]

  # Enriquecimiento vs Simple_DRR (igual que tu script antiguo)
  base_simple <- intensity_summary[candidate_class == "Simple_DRR",
                                   median_SE_overlaps_per_DRR]
  if (length(base_simple) == 1L && base_simple > 0) {
    intensity_summary[, enrichment_vs_simple :=
                        round(median_SE_overlaps_per_DRR / base_simple, 2)]
  } else {
    intensity_summary[, enrichment_vs_simple := NA_real_]
  }
  setorder(intensity_summary, -enrichment_vs_simple)

  list(overlap = overlap[], by_class = by_class[],
       intensity = intensity[], intensity_summary = intensity_summary[])
}

## ----------------------------------------------------------------------------
## 3) Construir el enlace DRR -> clusters -> CRMs originales -> genes.
## ----------------------------------------------------------------------------
## Para asignar genes a cada DRR necesitamos saber qué clusters (CRMs reducidos)
## caen en cada DRR. Lo hacemos por solapamiento posicional entre la coordenada
## del cluster (repr_*) y la DRR, ya que ambos viven en el mismo espacio.
link_drr_to_genes <- function(drr, crm_reduced, mapping, enh2gene,
                              crm_id_col = "crm_ID",
                              gene_col   = "hgnc_symbol_target_genes") {

  cand <- drr[candidate_class != "Extensive_overlap"]

  # a) DRR <-> cluster por solapamiento de coordenadas
  drr_gr <- GRanges(cand$chr, IRanges(cand$drr_start, cand$drr_end),
                    drr_id = cand$drr_id, candidate_class = cand$candidate_class)
  cl_gr  <- GRanges(crm_reduced$chr,
                    IRanges(crm_reduced$repr_start, crm_reduced$repr_end),
                    cluster_id = crm_reduced$cluster_id)
  ov <- findOverlaps(drr_gr, cl_gr)
  drr_cluster <- data.table(
    drr_id          = mcols(drr_gr)$drr_id[queryHits(ov)],
    candidate_class = mcols(drr_gr)$candidate_class[queryHits(ov)],
    cluster_id      = mcols(cl_gr)$cluster_id[subjectHits(ov)]
  )

  # b) cluster <-> CRM original (mapping: original_id, cluster_id)
  map <- as.data.table(mapping)[, .(cluster_id, original_id)]

  # c) CRM original <-> gen (enh2gene)
  e2g <- as.data.table(copy(enh2gene))
  setnames(e2g, old = c(crm_id_col, gene_col),
           new = c("original_id", "gene_symbol"), skip_absent = TRUE)
  e2g <- e2g[!is.na(original_id) & !is.na(gene_symbol) & gene_symbol != "",
             .(original_id, gene_symbol)]

  # Encadenar: drr -> cluster -> crm -> gen
  dcg <- merge(drr_cluster, map, by = "cluster_id", allow.cartesian = TRUE)
  dcg <- merge(dcg, e2g, by = "original_id", allow.cartesian = TRUE)
  dcg <- unique(dcg[, .(drr_id, candidate_class, gene_symbol)])

  # Genes por DRR
  genes_per_drr <- dcg[, .(
    n_genes_ours = uniqueN(gene_symbol),
    genes_ours   = paste(sort(unique(gene_symbol)), collapse = ";")
  ), by = .(drr_id, candidate_class)]

  genes_per_drr[]
}

## ----------------------------------------------------------------------------
## 4) Validación GÉNICA: genes nuestros vs genes SEdb por DRR.
## ----------------------------------------------------------------------------
.split_genes <- function(x) {
  x <- unique(unlist(strsplit(as.character(x), "[,;]")))
  x <- trimws(x); x <- x[!is.na(x) & x != ""]; unique(x)
}

validate_genes_vs_sedb <- function(genes_per_drr, overlap, se_chr8,
                                   gene_cols = c("se_gene_overlap",
                                                 "se_gene_proximal",
                                                 "se_gene_closest")) {
  gene_cols <- intersect(gene_cols, names(se_chr8))

  # Genes SEdb por se_id
  se_genes <- se_chr8[, c("se_id", gene_cols), with = FALSE]

  # Unir genes SEdb a cada DRR vía la tabla de solapamiento posicional
  ov_genes <- merge(overlap[, .(drr_id, candidate_class, se_id)],
                    se_genes, by = "se_id", allow.cartesian = TRUE)
  sedb_per_drr <- ov_genes[, .(
    genes_SEdb = paste(sort(.split_genes(do.call(paste, c(.SD, sep=";")))),
                       collapse = ";")
  ), by = .(drr_id, candidate_class), .SDcols = gene_cols]

  # Comparar genes nuestros vs SEdb
  gv <- merge(genes_per_drr, sedb_per_drr,
              by = c("drr_id", "candidate_class"), all.x = TRUE)
  gv[is.na(genes_SEdb), genes_SEdb := ""]

  gv[, c("n_ours","n_sedb","n_shared","shared") := {
    o <- .split_genes(genes_ours); s <- .split_genes(genes_SEdb)
    sh <- intersect(o, s)
    list(length(o), length(s), length(sh), paste(sort(sh), collapse=";"))
  }, by = drr_id]

  gv[, gene_match := n_shared > 0]

  gv_summary <- gv[, .(
    n_candidates       = .N,
    n_with_gene_match  = sum(gene_match),
    fraction_match     = mean(gene_match),
    median_shared      = median(n_shared)
  ), by = candidate_class][order(-fraction_match)]

  list(gene_validation = gv[], summary = gv_summary[])
}

## ----------------------------------------------------------------------------
## 5) Orquestador
## ----------------------------------------------------------------------------
run_sedb_validation <- function(drr, crm_reduced, mapping, enh2gene,
                                bed_path = "SEdb_Human_SE.bed",
                                chr = "chr8",
                                crm_id_col = "crm_ID",
                                gene_col   = "hgnc_symbol_target_genes") {
  se_chr8 <- load_sedb(bed_path, chr = chr)

  .msg("Validación posicional DRR vs SEdb...")
  pos <- validate_drr_vs_sedb(drr, se_chr8)

  if (is.null(mapping) || is.null(enh2gene)) {
    motivo <- if (is.null(mapping)) "mapping = NULL" else "enh2gene = NULL"
    .msg(motivo, ": se omite la parte génica (solo validación posicional).")
    return(list(positional = pos, genes_per_drr = NULL, genic = NULL))
  }

  # Comprobación de columnas de enh2gene (error claro si no cuadran)
  if (!all(c(crm_id_col, gene_col) %in% names(enh2gene))) {
    .msg("AVISO: enh2gene no tiene las columnas '", crm_id_col, "' / '", gene_col,
         "'. Columnas disponibles: ", paste(names(enh2gene), collapse = ", "),
         ". Se omite la parte génica.")
    return(list(positional = pos, genes_per_drr = NULL, genic = NULL))
  }

  .msg("Enlazando DRR -> clusters -> CRMs -> genes...")
  genes_per_drr <- link_drr_to_genes(drr, crm_reduced, mapping, enh2gene,
                                     crm_id_col = crm_id_col, gene_col = gene_col)

  .msg("Validación génica vs SEdb...")
  gen <- validate_genes_vs_sedb(genes_per_drr, pos$overlap, se_chr8)

  list(
    positional = pos,
    genes_per_drr = genes_per_drr,
    genic = gen
  )
}

## ============================================================================
## EJEMPLO DE USO
## ----------------------------------------------------------------------------
## source("drr_compare_criteria.R")
## source("drr_sedb_validation.R")
##
## cmp <- compare_drr_criteria(red_crm$crm_reduced, max_gap = 12500)
##
## val <- run_sedb_validation(
##   drr         = cmp$stacking,
##   crm_reduced = red_crm$crm_reduced,
##   mapping     = red_crm$mapping,
##   enh2gene    = enh2gene,
##   bed_path    = "SEdb_Human_SE.bed",
##   crm_id_col  = "crm_ID",                  # ajusta al nombre real
##   gene_col    = "hgnc_symbol_target_genes" # ajusta al nombre real
## )
##
## print(val$positional$by_class)            # fracción con solape SEdb por clase
## print(val$positional$intensity_summary)   # enriquecimiento vs Simple
## print(val$genic$summary)                  # recuperación génica por clase
## ============================================================================
