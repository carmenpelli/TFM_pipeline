################################################################################
##  tad_threshold_sweep.R — Exploración de umbrales de redundancia para TADs
##  --------------------------------------------------------------------------
##  Este script realiza un análisis exploratorio de umbrales de solapamiento
##  recíproco para evaluar su efecto sobre la reducción de redundancia de TADs.
##
##  Para una rejilla de umbrales, calcula el número de pares solapantes,
##  aristas, clústeres y entidades fusionadas que se obtendrían bajo distintos
##  criterios de similitud.
##
##  El objetivo es apoyar la selección del umbral mediante resúmenes
##  cuantitativos, sin modificar ni reducir definitivamente las anotaciones.
##
##  Métricas calculadas:
##    overlap_length, Jaccard, Simpson y reciprocal_overlap.
##
##  Dependencias:
##    data.table, GenomicRanges, igraph.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(igraph)
})

## ----------------------------------------------------------------------------
## 1) Cálculo de pares solapantes y métricas de similitud.
##    Este paso se realiza una sola vez; posteriormente se filtra por umbral.
## ----------------------------------------------------------------------------

#' Calcula pares solapados (i<j) y sus métricas de similitud para un conjunto.
#'
#' @param dt data.table con chr,start,end,ID (un registro por ID).
#' @return data.table: id_i,id_j,overlap_length,len_i,len_j,
#'         jaccard,simpson,reciprocal_overlap
compute_overlap_pairs_metrics <- function(dt) {
  dt <- as.data.table(dt)
  stopifnot(all(c("chr", "start", "end", "ID") %in% names(dt)))

  gr <- GRanges(
    seqnames = dt$chr,
    ranges   = IRanges(start = dt$start, end = dt$end),
    ID       = dt$ID
  )

  hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
  qh <- queryHits(hits); sh <- subjectHits(hits)

  keep <- qh < sh           # pares i<j, sin autopares ni duplicados simétricos
  qh <- qh[keep]; sh <- sh[keep]

  if (length(qh) == 0L) {
    return(data.table(
      id_i = character(), id_j = character(),
      overlap_length = integer(), len_i = integer(), len_j = integer(),
      jaccard = numeric(), simpson = numeric(), reciprocal_overlap = numeric()
    ))
  }

  st <- start(gr); en <- end(gr); len <- width(gr)
  s_i <- st[qh]; e_i <- en[qh]; L_i <- len[qh]
  s_j <- st[sh]; e_j <- en[sh]; L_j <- len[sh]

  ov <- pmin(e_i, e_j) - pmax(s_i, s_j) + 1L
  ov[ov < 0L] <- 0L
  union_len <- L_i + L_j - ov

  pairs <- data.table(
    id_i = mcols(gr)$ID[qh],
    id_j = mcols(gr)$ID[sh],
    overlap_length = as.integer(ov),
    len_i = as.integer(L_i),
    len_j = as.integer(L_j)
  )
  pairs[, jaccard := overlap_length / union_len]
  pairs[, simpson := overlap_length / pmin(len_i, len_j)]
  pairs[, reciprocal_overlap := pmin(overlap_length / len_i,
                                     overlap_length / len_j)]
  pairs[]
}

## ----------------------------------------------------------------------------
## 2) Conteo de clústeres a partir de un conjunto de aristas.
## ----------------------------------------------------------------------------

#' Cuenta componentes conexas a partir del universo de IDs y las aristas.
#'
#' @param all_ids vector con todos los IDs, incluidos singletons.
#' @param edges   data.table con columnas id_i,id_j (puede estar vacía).
#' @return lista: n_clusters, sizes (vector de tamaños de cada cluster)
count_clusters <- function(all_ids, edges) {
  if (nrow(edges) > 0L) {
    edf <- as.data.frame(edges[, .(id_i, id_j)], stringsAsFactors = FALSE)
  } else {
    edf <- data.frame(id_i = character(), id_j = character(),
                      stringsAsFactors = FALSE)
  }
  g <- graph_from_data_frame(
    d = edf, directed = FALSE,
    vertices = data.frame(name = all_ids, stringsAsFactors = FALSE)
  )
  comp <- igraph::components(g)
  list(n_clusters = comp$no, sizes = as.integer(comp$csize))
}

## ----------------------------------------------------------------------------
## 3) Barrido de umbrales.
## ----------------------------------------------------------------------------

#' Evalúa una rejilla de umbrales de reciprocal_overlap y similitud sobre TADs.
#'
#' Para cada umbral de recíproco evalúa DOS escenarios de "arista":
#'   (A) "reciprocal_only": arista si reciprocal_overlap >= thr_recip
#'       -> evalúa la redundancia mediante solapamiento recíproco puro.
#'   (B) "recip_plus_sim": arista si reciprocal_overlap >= thr_recip Y
#'          (jaccard >= jaccard_thr O simpson >= simpson_thr)
#'       -> reproduce el criterio compuesto usado en CRMs con fines comparativos.
#'
#' @param tad_unique   data.table TADs colapsados (chr,start,end,ID).
#' @param recip_grid   vector de umbrales de reciprocal_overlap a probar.
#' @param jaccard_thr  umbral de jaccard para el escenario (B).
#' @param simpson_thr  umbral de simpson para el escenario (B).
#' @return data.table con una fila por (escenario, umbral) y columnas:
#'   scenario, thr_recip, n_input, n_pairs_overlap, n_edges,
#'   n_clusters, n_merged_clusters, max_cluster_size,
#'   n_entities_in_merges, reduction_ratio
sweep_tad_thresholds <- function(tad_unique,
                                 recip_grid  = c(0.50, 0.60, 0.70, 0.80, 0.90),
                                 jaccard_thr = 0.70,
                                 simpson_thr = 0.99) {

  dt <- as.data.table(tad_unique)
  all_ids <- dt$ID
  n_input <- length(all_ids)

  message("[sweep] Calculando pares solapados y métricas...")
  pairs <- compute_overlap_pairs_metrics(dt)
  message("[sweep] Pares solapados (i<j): ", nrow(pairs))

  # Función interna para evaluar un escenario y un umbral.
  eval_one <- function(scenario, thr) {
    if (nrow(pairs) == 0L) {
      edges <- pairs
    } else if (scenario == "reciprocal_only") {
      edges <- pairs[reciprocal_overlap >= thr]
    } else { # recip_plus_sim
      edges <- pairs[reciprocal_overlap >= thr &
                       (jaccard >= jaccard_thr | simpson >= simpson_thr)]
    }

    cl <- count_clusters(all_ids, edges)
    sizes <- cl$sizes
    merged <- sizes[sizes > 1L]   # clústeres que fusionan más de un TAD

    data.table(
      scenario             = scenario,
      thr_recip            = thr,
      n_input              = n_input,
      n_pairs_overlap      = nrow(pairs),
      n_edges              = nrow(edges),
      n_clusters           = cl$n_clusters,
      n_merged_clusters    = length(merged),               # número de grupos fusionados
      max_cluster_size     = if (length(sizes)) max(sizes) else 0L,
      n_entities_in_merges = if (length(merged)) sum(merged) else 0L,
      # fracción de TADs eliminados bajo este umbral
      reduction_ratio      = 1 - cl$n_clusters / n_input
    )
  }

  res <- rbindlist(lapply(recip_grid, function(t) {
    rbind(
      eval_one("reciprocal_only", t),
      eval_one("recip_plus_sim", t)
    )
  }))

  setorder(res, scenario, thr_recip)
  res[]
}

## ----------------------------------------------------------------------------
## 4) Distribución de reciprocal_overlap para orientar la selección de umbrales
## ----------------------------------------------------------------------------

#' Devuelve la distribución del reciprocal_overlap de los pares solapados.
#' Resume la concentración de valores antes de seleccionar una rejilla de umbrales.
#'
#' @param tad_unique data.table TADs (chr,start,end,ID).
#' @param probs      cuantiles a reportar.
#' @return lista: summary (quantiles) y pairs (data.table completa de métricas)
describe_tad_overlap <- function(tad_unique,
                                 probs = c(0, .1, .25, .5, .75, .9, .95, .99, 1)) {
  pairs <- compute_overlap_pairs_metrics(as.data.table(tad_unique))
  if (nrow(pairs) == 0L) {
    message("No hay pares solapados.")
    return(list(summary = NULL, pairs = pairs))
  }
  q_recip   <- quantile(pairs$reciprocal_overlap, probs = probs)
  q_jaccard <- quantile(pairs$jaccard,            probs = probs)
  q_simpson <- quantile(pairs$simpson,            probs = probs)
  summ <- data.table(
    quantile           = names(q_recip),
    reciprocal_overlap = as.numeric(q_recip),
    jaccard            = as.numeric(q_jaccard),
    simpson            = as.numeric(q_simpson)
  )
  list(summary = summ[], pairs = pairs)
}

## ============================================================================
## Ejemplo de uso exploratorio
## ----------------------------------------------------------------------------
## # tad_unique debe estar cargado en memoria.
##
## # (a) Resumen de la distribución de solapamiento recíproco:
## desc <- describe_tad_overlap(tad_unique)
## print(desc$summary)
##
## # (b) Barrido de umbrales y comparación de escenarios:
## sweep <- sweep_tad_thresholds(
##   tad_unique,
##   recip_grid = c(0.50, 0.60, 0.70, 0.80, 0.90)
## )
## print(sweep)
##
## # Interpretación de columnas:
## #   n_clusters           -> nº de TADs que quedarían tras reducir
## #   n_merged_clusters    -> cuántos grupos fusionan >1 TAD
## #   max_cluster_size     -> tamaño del grupo más grande.
## #   reduction_ratio      -> fracción de TADs eliminados
## #
## # Seleccionar un umbral en el que la reducción se estabilice sin un aumento
## # excesivo de max_cluster_size.
## ============================================================================
