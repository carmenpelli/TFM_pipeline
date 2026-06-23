################################################################################
##  crm_reduce_final.R — Reducción de redundancia de CRMs
##  --------------------------------------------------------------------------
##  Este script reduce la redundancia de CRMs mediante un criterio compuesto de
##  similitud y un procedimiento de anti-chaining basado en representantes.
##
##  Criterio de arista:
##    reciprocal_overlap >= 0.50
##    y (jaccard >= 0.70 o simpson >= 0.99)
##
##  La reducción no define nuevas entidades funcionales, sino que consolida
##  anotaciones estructuralmente redundantes procedentes de múltiples fuentes,
##  conservando la trazabilidad con los identificadores originales.
##
##  Dependencias:
##    data.table, GenomicRanges e igraph.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(igraph)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

## ----------------------------------------------------------------------------
## Helper de concatenación de valores únicos (idéntico criterio que en TADs):
## desdobla multivaluados con ";", quita placeholder "-", deduplica.
## ----------------------------------------------------------------------------
.collapse_unique_crm <- function(x, sep = ";") {
  parts <- unlist(strsplit(as.character(x), sep, fixed = TRUE))
  parts <- trimws(parts)
  parts <- parts[parts != "" & !is.na(parts)]
  real  <- parts[parts != "-"]
  if (length(real) > 0L) parts <- real
  uq <- unique(parts)
  if (length(uq) == 0L) return(NA_character_)
  paste(uq, collapse = sep)
}

## ============================================================================
## Función principal: reduce_redundancy_crms()
## ============================================================================

#' Reduce la redundancia de CRMs por anti-chaining + criterio compuesto.
#'
#' @param crm_unique  data.table de CRMs colapsados por ID (chr,start,end,ID +
#'        metadatos).
#' @param thr_recip   umbral recíproco (por defecto 0.50, definido para el análisis).
#' @param jaccard_thr,simpson_thr umbrales del criterio compuesto (0.70/0.99).
#' @param meta_cols   columnas de trazabilidad a concatenar (valores únicos).
#' @param sum_cols    columnas numéricas a sumar por cluster (por defecto n_rows_collapsed).
#' @param chunk_size  tamaño de bloque para el cálculo de aristas.
#' @param id_sep      separador para IDs y metadatos.
#' @param verbose_ac  progreso del anti-chaining.
#' @return lista: crm_reduced, mapping, summary.
reduce_redundancy_crms <- function(crm_unique,
                                   thr_recip   = 0.50,
                                   jaccard_thr = 0.70,
                                   simpson_thr = 0.99,
                                   meta_cols   = c("biosample_name",
                                                   "biological_sample_type",
                                                   "cell_line_CLO",
                                                   "cell_type_CL",
                                                   "anatomical_structures_UBERON",
                                                   "BTO"),
                                   sum_cols    = c("n_rows_collapsed"),
                                   chunk_size  = 5e6,
                                   id_sep      = ";",
                                   verbose_ac  = FALSE) {

  if (!exists("compute_crm_edges_chunked"))
    stop("Carga crm_explore.R (compute_crm_edges_chunked).")
  if (!exists("anti_chaining_fast"))
    stop("Carga crm_antichaining_fast.R (anti_chaining_fast).")

  dt <- as.data.table(copy(crm_unique))
  stopifnot(all(c("chr", "start", "end", "ID") %in% names(dt)))
  dt[, start := as.integer(start)]
  dt[, end   := as.integer(end)]
  dt[, length := end - start + 1L]

  meta_present <- intersect(meta_cols, names(dt))
  sum_present  <- intersect(sum_cols,  names(dt))
  for (cc in sum_present) dt[, (cc) := as.numeric(get(cc))]

  n_input <- nrow(dt)
  .msg("Reducción CRMs: ", n_input, " CRMs de entrada. Criterio: recíproco>=",
       thr_recip, " Y (J>=", jaccard_thr, " O S>=", simpson_thr, ").")

  ## --- 1) Aristas (criterio compuesto, por bloques) ------------------------
  edges <- compute_crm_edges_chunked(
    dt, thr_recip = thr_recip, jaccard_thr = jaccard_thr,
    simpson_thr = simpson_thr, criterion = "composite",
    keep_metrics = FALSE, chunk_size = chunk_size
  )

  ## --- 2) Clustering anti-chaining (optimizado) --------------------------------
  .msg("Clustering anti-chaining...")
  membership <- anti_chaining_fast(dt, edges, verbose = verbose_ac)
  n_clusters <- membership[, uniqueN(cluster_id)]
  .msg("Entidades no redundantes: ", n_clusters,
       " (reducción ", round(100 * (1 - n_clusters / n_input), 1), "%).")

  ## --- 3) Consenso / union por cluster -------------------------------------
  d <- merge(dt, membership[, .(ID, cluster_id, representative_id)], by = "ID")

  coords <- d[, .(
    chr             = chr[1L],
    consensus_start = max(start),
    consensus_end   = min(end),
    union_start     = min(start),
    union_end       = max(end),
    n_entities      = .N
  ), by = cluster_id]

  coords[, has_core := consensus_end >= consensus_start]
  coords[, consensus_length := pmax(consensus_end - consensus_start + 1L, 0L)]
  coords[, union_length     := union_end - union_start + 1L]
  coords[, boundary_variability := union_length - consensus_length]
  coords[, recommended_coord := ifelse(has_core, "consensus", "union")]
  coords[, repr_start := ifelse(has_core, consensus_start, union_start)]
  coords[, repr_end   := ifelse(has_core, consensus_end,   union_end)]

  ## --- 4) representative_id (el del anti-chaining: la semilla del cluster) --
  ## El anti-chaining ya define un representante (la semilla = CRM más grande).
  ## Lo conservamos como representative_id del cluster.
  rep_dt <- unique(d[, .(cluster_id, representative_id)])

  ## --- 5) Lista de IDs originales ------------------------------------------
  ids_dt <- d[, .(original_ids = paste(ID, collapse = id_sep),
                  n_original_ids = .N), by = cluster_id]

  ## --- 6) Metadatos concatenados (valores únicos) --------------------------
  if (length(meta_present) > 0L) {
    meta_dt <- d[, lapply(.SD, .collapse_unique_crm),
                 by = cluster_id, .SDcols = meta_present]
  } else {
    meta_dt <- unique(d[, .(cluster_id)])
  }

  ## --- 6b) Columnas numéricas sumadas --------------------------------------
  if (length(sum_present) > 0L) {
    sum_dt <- d[, lapply(.SD, sum, na.rm = TRUE),
                by = cluster_id, .SDcols = sum_present]
    setnames(sum_dt, sum_present, paste0("sum_", sum_present))
  } else {
    sum_dt <- NULL
  }

  ## --- 7) Ensamblar crm_reduced --------------------------------------------
  parts <- list(coords, rep_dt, ids_dt, meta_dt)
  if (!is.null(sum_dt)) parts <- c(parts, list(sum_dt))
  crm_reduced <- Reduce(function(a, b) merge(a, b, by = "cluster_id"), parts)

  front <- c("cluster_id", "representative_id", "chr",
             "repr_start", "repr_end", "recommended_coord", "has_core",
             "consensus_start", "consensus_end",
             "union_start", "union_end",
             "consensus_length", "union_length", "boundary_variability",
             "n_entities", "n_original_ids",
             paste0("sum_", sum_present), "original_ids")
  front <- intersect(front, names(crm_reduced))
  setcolorder(crm_reduced, c(front, setdiff(names(crm_reduced), front)))
  setkey(crm_reduced, cluster_id)

  ## --- 8) Tabla de mapeo (trazabilidad pura) -------------------------------
  map_cols <- c("ID", "cluster_id", "chr", "start", "end", "length",
                meta_present, sum_present)
  map_cols <- intersect(map_cols, names(d))
  mapping <- d[, ..map_cols]
  setnames(mapping, "ID", "original_id")
  setkey(mapping, cluster_id, original_id)

  ## --- 9) Resumen -----------------------------------------------------------
  summary_dt <- data.table(
    chr              = dt$chr[1L],
    thr_recip        = thr_recip,
    jaccard_thr      = jaccard_thr,
    simpson_thr      = simpson_thr,
    clustering       = "anti_chaining_greedy",
    n_input          = n_input,
    n_clusters       = n_clusters,
    n_singletons     = coords[n_entities == 1L, .N],
    n_merged         = coords[n_entities > 1L, .N],
    max_cluster_size = max(coords$n_entities),
    n_without_core   = coords[has_core == FALSE, .N],
    reduction_ratio  = 1 - n_clusters / n_input
  )
  .msg("Clusters sin núcleo común (usan union): ", summary_dt$n_without_core, ".")

  list(
    crm_reduced = crm_reduced[],
    mapping     = mapping[],
    summary     = summary_dt[]
  )
}

## ----------------------------------------------------------------------------
## Guardado (coherente con el estilo de save_tad_reduction).
## ----------------------------------------------------------------------------
save_crm_reduction <- function(reduction, chr,
                               output_dir = "results/intermediate") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  f1 <- file.path(output_dir, paste0("crm_reduced_", chr, ".tsv.gz"))
  f2 <- file.path(output_dir, paste0("crm_mapping_", chr, ".tsv.gz"))
  f3 <- file.path(output_dir, paste0("crm_reduction_summary_", chr, ".tsv"))
  fwrite(reduction$crm_reduced, f1, sep = "\t")
  fwrite(reduction$mapping,     f2, sep = "\t")
  fwrite(reduction$summary,     f3, sep = "\t")
  .msg("Guardado: ", f1); .msg("Guardado: ", f2); .msg("Guardado: ", f3)
  invisible(c(f1, f2, f3))
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("crm_explore.R")
## source("crm_antichaining_fast.R")
## source("crm_reduce_final.R")
##
## # Sobre una región (prueba):
## reg <- subset_region(enh_unique, "chr8", 125e6, 130e6)
## red <- reduce_redundancy_crms(reg, thr_recip = 0.50)
## print(red$summary)
## head(red$crm_reduced[, .(cluster_id, n_entities, sum_n_rows_collapsed,
##                          repr_start, repr_end, recommended_coord)])
##
## # Sobre el cromosoma completo (cuando proceda):
## # red_chr8 <- reduce_redundancy_crms(enh_unique, thr_recip = 0.50)
## # save_crm_reduction(red_chr8, chr = "chr8")
## ============================================================================
