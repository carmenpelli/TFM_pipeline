################################################################################
##  crm_inspect.R — Inspección de clústeres de CRMs
##  --------------------------------------------------------------------------
##  Este script construye clústeres de CRMs mediante componentes conexas bajo
##  un criterio de similitud y evalúa la cohesión de los clústeres de mayor
##  tamaño.
##
##  El análisis permite distinguir clústeres con núcleo común, compatibles con
##  redundancia estructural o acumulación local de CRMs, de clústeres generados
##  por encadenamiento transitivo.
##
##  Reutiliza compute_crm_edges_chunked() de crm_explore.R.
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

## ============================================================================
## Inspección de clústeres de CRMs
## ============================================================================

#' Construye clústeres de CRMs y evalúa la cohesión de los de mayor tamaño.
#'
#' @param dt          data.table de CRMs (región): chr,start,end,ID + metadatos.
#' @param thr_recip   umbral recíproco (def. 0.50, criterio del documento).
#' @param jaccard_thr,simpson_thr umbrales del criterio compuesto.
#' @param criterion   "composite" (del documento) o "reciprocal_only".
#' @param top_n       número de clústeres mayores a examinar en detalle.
#' @param sample_rows nº de CRMs a mostrar del cluster más grande (muestra).
#' @param meta_cols   columnas de metadatos a incluir en la muestra.
#' @param chunk_size  tamaño de bloque para compute_crm_edges_chunked.
#' @return lista:
#'   sizes      : data.table(cluster_id, n) ordenada de forma descendente.
#'   cohesion   : diagnóstico por clúster de los top_n mayores.
#'   biggest    : muestra de CRMs del clúster más grande.
#'   membership : data.table(ID, cluster_id) completa
inspect_crm_clusters <- function(dt,
                                 thr_recip    = 0.50,
                                 jaccard_thr  = 0.70,
                                 simpson_thr  = 0.99,
                                 criterion    = "composite",
                                 top_n        = 10,
                                 sample_rows  = 40,
                                 meta_cols    = c("biosample_name",
                                                  "biological_sample_type",
                                                  "interval_length"),
                                 chunk_size   = 5e6) {

  dt <- as.data.table(dt)
  stopifnot(all(c("chr", "start", "end", "ID") %in% names(dt)))
  all_ids <- dt$ID

  ## --- Cálculo de aristas mediante la función por bloques --------------------
  if (!exists("compute_crm_edges_chunked"))
    stop("Falta compute_crm_edges_chunked(). Carga primero crm_explore.R.")

  edges <- compute_crm_edges_chunked(
    dt, thr_recip = thr_recip, jaccard_thr = jaccard_thr,
    simpson_thr = simpson_thr, criterion = criterion,
    keep_metrics = FALSE, chunk_size = chunk_size
  )

  ## --- Construcción del grafo y componentes conexas --------------------------
  edf <- if (nrow(edges) > 0L)
    as.data.frame(edges[, .(id_i, id_j)], stringsAsFactors = FALSE) else
    data.frame(id_i = character(), id_j = character(), stringsAsFactors = FALSE)
  g <- graph_from_data_frame(edf, directed = FALSE,
                             vertices = data.frame(name = all_ids))
  comp <- igraph::components(g)
  membership <- data.table(ID = names(comp$membership),
                           cluster_id = as.integer(comp$membership))

  sizes <- membership[, .(n = .N), by = cluster_id][order(-n)]
  .msg("Clústeres: ", nrow(sizes), " | clúster mayor: ", sizes$n[1L], " CRMs.")

  ## --- Diagnóstico de cohesión de los clústeres de mayor tamaño --------------
  top_ids <- sizes[n > 1L][seq_len(min(top_n, .N))]$cluster_id
  d <- merge(dt, membership, by = "ID")
  d[, length := end - start + 1L]

  cohesion <- d[cluster_id %in% top_ids, {
    sp   <- max(end) - min(start) + 1L         # extensión total del clúster
    core <- min(end) - max(start) + 1L         # núcleo común compartido por los miembros
    .(n               = .N,
      span            = sp,
      core_len        = core,
      has_common_core = core > 0L,
      core_fraction   = ifelse(sp > 0, core / sp, NA_real_),
      len_min         = min(length),
      len_max         = max(length),
      len_ratio       = max(length) / min(length))
  }, by = cluster_id][order(-n)]

  ## --- Muestra del clúster de mayor tamaño ----------------------------------
  big_id <- sizes$cluster_id[1L]
  big <- d[cluster_id == big_id][order(start)]
  show_cols <- intersect(c("ID", "chr", "start", "end", "length",
                           meta_cols), names(big))
  # Muestreo distribuido a lo largo de las coordenadas del clúster.
  ns <- nrow(big)
  if (ns > sample_rows) {
    idx <- unique(round(seq(1, ns, length.out = sample_rows)))
    biggest <- big[idx, ..show_cols]
  } else {
    biggest <- big[, ..show_cols]
  }

  ## --- Diagnóstico global de cohesión de clústeres multi-miembro -------------
  global_coh <- d[cluster_id %in% sizes[n > 1L]$cluster_id, {
    core <- min(end) - max(start) + 1L
    sp   <- max(end) - min(start) + 1L
    .(has_common_core = core > 0L,
      core_fraction   = ifelse(sp > 0, core / sp, NA_real_),
      len_ratio       = max(length) / min(length))
  }, by = cluster_id]
  global_summary <- global_coh[, .(
    n_multi_clusters     = .N,
    frac_con_nucleo      = mean(has_common_core),
    core_frac_mediana    = median(core_fraction),
    len_ratio_mediana    = median(len_ratio)
  )]
  .msg("Clústeres multi-miembro con núcleo común: ",
       round(100 * global_summary$frac_con_nucleo, 1), "% | ",
       "core_fraction mediana: ", round(global_summary$core_frac_mediana, 3), ".")

  list(
    sizes          = sizes[],
    cohesion       = cohesion[],
    biggest        = biggest[],
    big_cluster_id = big_id,
    big_n          = sizes$n[1L],
    global_summary = global_summary[],
    membership     = membership[]
  )
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("crm_explore.R")        # define compute_crm_edges_chunked
## source("crm_inspect.R")        # este archivo
##
## reg <- subset_region(enh_unique, "chr8", 125e6, 130e6)
##
## insp <- inspect_crm_clusters(reg, thr_recip = 0.50, top_n = 10)
##
## # Diagnóstico global de existencia de núcleo común.
## print(insp$global_summary)
## #   frac_con_nucleo cercano a 1 y core_frac_mediana > 0 sugieren clústeres cohesionados.
## #   valores bajos sugieren posible encadenamiento transitivo.
##
## # Clústeres de mayor tamaño en detalle:
## print(insp$cohesion)
## #   has_common_core = FALSE o core_fraction <= 0 sugieren ausencia de núcleo común.
##
## # Muestra del clúster más grande:
## print(insp$biggest)
## #   La concentración de coordenadas en una región estrecha sugiere cohesión.
## #   La dispersión progresiva a lo largo del span sugiere encadenamiento.
## ============================================================================
