################################################################################
##  crm_antichaining_cohesion.R — Evaluación de cohesión en clústeres de CRMs
##  --------------------------------------------------------------------------
##  Este script inspecciona la cohesión interna de los clústeres generados por
##  el procedimiento anti-chaining aplicado a CRMs.
##
##  El objetivo es comprobar si los clústeres presentan una región común de
##  intersección entre sus miembros, compatible con redundancia estructural, o
##  si persisten patrones de encadenamiento transitivo.
##
##  Métricas principales:
##    - has_common_core: existencia de una región común compartida por todos
##      los CRMs del clúster.
##    - core_fraction: fracción del intervalo total cubierta por dicho núcleo.
##    - len_ratio: relación entre las longitudes máxima y mínima dentro del
##      clúster.
##
##  Reutiliza anti_chaining_from_edges() y compute_crm_edges_chunked().
##
##  Dependencias:
##    data.table, GenomicRanges, igraph.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(igraph)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

#' Inspecciona la cohesión de los clústeres del anti-chaining sobre una región.
#'
#' @param dt          data.table de CRMs (región): chr,start,end,ID + metadatos.
#' @param thr_recip   umbral recíproco (def. 0.50).
#' @param jaccard_thr,simpson_thr umbrales del criterio compuesto.
#' @param top_n       número de clústeres mayores a examinar en detalle.
#' @param sample_rows número de CRMs a mostrar del clúster mayor.
#' @param meta_cols   metadatos a incluir en la muestra.
#' @param chunk_size  tamaño de bloque.
#' @return lista con sizes, cohesion, biggest, global_summary y membership.
inspect_antichaining_cohesion <- function(dt,
                                          thr_recip   = 0.50,
                                          jaccard_thr = 0.70,
                                          simpson_thr = 0.99,
                                          top_n       = 10,
                                          sample_rows = 40,
                                          meta_cols   = c("biosample_name",
                                                          "biological_sample_type",
                                                          "interval_length"),
                                          chunk_size  = 5e6) {

  if (!exists("compute_crm_edges_chunked"))
    stop("Carga crm_explore.R (compute_crm_edges_chunked).")
  if (!exists("anti_chaining_from_edges"))
    stop("Carga crm_compare_strategies.R (anti_chaining_from_edges).")

  dt <- as.data.table(dt)
  dt[, length := end - start + 1L]
  n_input <- nrow(dt)

  ## Cálculo de aristas mediante el criterio compuesto y aplicación de anti-chaining.
  edges <- compute_crm_edges_chunked(
    dt, thr_recip = thr_recip, jaccard_thr = jaccard_thr,
    simpson_thr = simpson_thr, criterion = "composite",
    keep_metrics = FALSE, chunk_size = chunk_size
  )
  membership <- anti_chaining_from_edges(dt, edges)

  d <- merge(dt, membership, by = "ID")
  sizes <- membership[, .(n = .N), by = cluster_id][order(-n)]
  .msg("Anti-chaining: ", nrow(sizes), " clústeres | mayor: ", sizes$n[1L], ".")

  ## Cohesión de los clústeres de mayor tamaño.
  top_ids <- sizes[n > 1L][seq_len(min(top_n, .N))]$cluster_id
  cohesion <- d[cluster_id %in% top_ids, {
    sp   <- max(end) - min(start) + 1L
    core <- min(end) - max(start) + 1L
    .(n               = .N,
      span            = sp,
      core_len        = core,
      has_common_core = core > 0L,
      core_fraction   = ifelse(sp > 0, core / sp, NA_real_),
      len_min         = min(length),
      len_max         = max(length),
      len_ratio       = max(length) / min(length))
  }, by = cluster_id][order(-n)]

  ## Muestra del clúster de mayor tamaño.
  big_id <- sizes$cluster_id[1L]
  big <- d[cluster_id == big_id][order(start)]
  show_cols <- intersect(c("ID","chr","start","end","length", meta_cols),
                         names(big))
  ns <- nrow(big)
  biggest <- if (ns > sample_rows)
    big[unique(round(seq(1, ns, length.out = sample_rows))), ..show_cols] else
    big[, ..show_cols]

  ## Diagnóstico global sobre los clústeres con más de un miembro.
  multi_ids <- sizes[n > 1L]$cluster_id
  gc <- d[cluster_id %in% multi_ids, {
    core <- min(end) - max(start) + 1L
    sp   <- max(end) - min(start) + 1L
    .(has_common_core = core > 0L,
      core_fraction   = ifelse(sp > 0, core / sp, NA_real_),
      len_ratio       = max(length) / min(length))
  }, by = cluster_id]
  global_summary <- gc[, .(
    n_multi_clusters  = .N,
    frac_con_nucleo   = mean(has_common_core),
    core_frac_mediana = median(core_fraction),
    len_ratio_mediana = median(len_ratio)
  )]
  .msg("Clústeres multi-miembro con núcleo común: ",
       round(100 * global_summary$frac_con_nucleo, 1), "% | ",
       "core_fraction mediana: ", round(global_summary$core_frac_mediana, 3), ".")

  list(
    sizes          = sizes[],
    cohesion       = cohesion[],
    biggest        = biggest[],
    big_cluster_id = big_id,
    global_summary = global_summary[],
    membership     = membership[]
  )
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("crm_explore.R")
## source("crm_compare_strategies.R")
## source("crm_antichaining_cohesion.R")
##
## reg <- subset_region(enh_unique, "chr8", 125e6, 130e6)
## ac <- inspect_antichaining_cohesion(reg, thr_recip = 0.50)
##
## print(ac$global_summary)   # ¿los clusters tienen núcleo común?
## print(ac$cohesion)         # los 10 mayores en detalle
## print(ac$biggest)          # muestra del clúster de mayor tamaño
##
## # Interpretación:
## #   valores altos de frac_con_nucleo y core_fraction > 0 sugieren clústeres cohesionados.
## #   has_common_core = FALSE en clústeres grandes sugiere encadenamiento residual,
## #   que puede evaluarse aumentando thr_recip y repitiendo la inspección.
## ============================================================================
