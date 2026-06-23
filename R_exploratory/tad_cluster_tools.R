################################################################################
##  tad_cluster_tools.R — Inspección exploratoria de clústeres de TADs
##  --------------------------------------------------------------------------
##  Este script reúne utilidades para la fase exploratoria de reducción de
##  redundancia de TADs.
##
##  Funciones principales:
##    (1) inspect_tad_clusters()
##        Construye clústeres a partir de un umbral de solapamiento recíproco y
##        permite examinar los clústeres de mayor tamaño, sus coordenadas,
##        longitudes y métricas de cohesión.
##
##    (2) cluster_tads_anti_chaining()
##        Implementa una variante alternativa de agrupamiento greedy basada en
##        similitud respecto al representante del clúster.
##
##  El script incluye una versión autocontenida de compute_overlap_pairs_metrics()
##  para poder ejecutarse de forma independiente.
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
## Versión mínima para uso independiente.
## Si tad_threshold_sweep.R ya está cargado, esta definición mantiene la misma
## interfaz y cálculo.
## ----------------------------------------------------------------------------
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
  keep <- qh < sh
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

## ============================================================================
## (1) Inspección de clústeres mediante componentes conexas
## ============================================================================

#' Construye clústeres por componentes conexas e inspecciona los de mayor tamaño.
#'
#' @param tad_unique data.table TADs (chr,start,end,ID).
#' @param thr_recip  umbral de reciprocal_overlap para crear arista.
#' @param use_similarity si TRUE, replica el criterio estilo CRM:
#'        arista si recíproco>=thr_recip Y (jaccard>=jaccard_thr O simpson>=simpson_thr).
#'        Si FALSE, se utiliza únicamente el solapamiento recíproco.
#' @param jaccard_thr,simpson_thr umbrales para use_similarity=TRUE.
#' @param top_n      nº de clusters más grandes a devolver en detalle.
#' @return lista con:
#'   membership : data.table(ID, cluster_id)
#'   sizes      : data.table(cluster_id, n) ordenada desc
#'   top        : data.table con los TADs de los top_n clusters más grandes,
#'                con coordenadas, longitud y métricas de cohesión por cluster
#'   cohesion   : data.table por clúster con métricas de cohesión.
inspect_tad_clusters <- function(tad_unique,
                                 thr_recip      = 0.80,
                                 use_similarity = FALSE,
                                 jaccard_thr    = 0.70,
                                 simpson_thr    = 0.99,
                                 top_n          = 5) {

  dt <- as.data.table(tad_unique)
  all_ids <- dt$ID

  pairs <- compute_overlap_pairs_metrics(dt)

  # Selección de aristas según criterio
  if (nrow(pairs) == 0L) {
    edges <- pairs
  } else if (use_similarity) {
    edges <- pairs[reciprocal_overlap >= thr_recip &
                     (jaccard >= jaccard_thr | simpson >= simpson_thr)]
  } else {
    edges <- pairs[reciprocal_overlap >= thr_recip]
  }

  # Grafo y componentes conexas.
  edf <- if (nrow(edges) > 0L)
    as.data.frame(edges[, .(id_i, id_j)], stringsAsFactors = FALSE) else
    data.frame(id_i = character(), id_j = character(), stringsAsFactors = FALSE)

  g <- graph_from_data_frame(
    d = edf, directed = FALSE,
    vertices = data.frame(name = all_ids, stringsAsFactors = FALSE)
  )
  comp <- igraph::components(g)

  membership <- data.table(
    ID         = names(comp$membership),
    cluster_id = as.integer(comp$membership)
  )

  sizes <- membership[, .(n = .N), by = cluster_id][order(-n)]

  # Identificadores de los clústeres de mayor tamaño.
  top_ids <- sizes[n > 1L][seq_len(min(top_n, .N))]$cluster_id

  # Detalle de los TADs incluidos en esos clústeres.
  d_annot <- merge(dt, membership, by = "ID")
  d_annot[, length := end - start + 1L]

  top <- d_annot[cluster_id %in% top_ids][order(cluster_id, start)]

  ## --- Métricas de cohesión por clúster -------------------------------------
  ## Para cada clúster se calcula:
  ##   n              : número de TADs.
  ##   span           : extensión total del clúster.
  ##   core_len       : longitud del núcleo común.
  ##   len_min/len_max: longitudes extremas de los TADs del clúster.
  ##   len_ratio      : relación entre la longitud máxima y mínima.
  ##   has_common_core: TRUE si todos los TADs comparten una región común.
  ## Estas métricas permiten detectar clústeres con núcleo común y clústeres
  ## potencialmente afectados por encadenamiento transitivo.
  cohesion <- d_annot[, {
    sp   <- max(end) - min(start) + 1L
    core <- min(end) - max(start) + 1L
    lmin <- min(length); lmax <- max(length)
    .(n               = .N,
      span            = sp,
      core_len        = core,
      has_common_core = core > 0L,
      len_min         = lmin,
      len_max         = lmax,
      len_ratio       = lmax / lmin,
      # fracción del span cubierta por el núcleo común
      core_fraction   = ifelse(sp > 0, core / sp, NA_real_))
  }, by = cluster_id][order(-n)]

  list(
    membership = membership[],
    sizes      = sizes[],
    top        = top[],
    cohesion   = cohesion[]
  )
}

#' Resumen rápido de un cluster concreto (para mirar uno en detalle).
#'
#' @param inspection salida de inspect_tad_clusters().
#' @param cluster    cluster_id a examinar.
#' @return data.table con los TADs de ese cluster (coordenadas, longitud).
show_cluster <- function(inspection, cluster) {
  memb <- inspection$membership
  ids  <- memb[cluster_id == cluster]$ID
  top  <- inspection$top
  out  <- top[cluster_id == cluster]
  if (nrow(out) == 0L) {
    message("El clúster ", cluster, " no está entre los clústeres mostrados; ",
            "contiene ", length(ids), " TADs. Aumente top_n para visualizarlo.")
  }
  out[]
}

## ============================================================================
## (2) Agrupamiento anti-chaining como variante metodológica alternativa
## ----------------------------------------------------------------------------
## Esta función implementa una alternativa al agrupamiento por componentes
## conexas. Su uso debe documentarse explícitamente si se adopta en el análisis
## principal.
##
## Estrategia (greedy por representante):
##   1. Se calculan los pares que superan el umbral recíproco (aristas válidas).
##   2. Se ordenan los TADs por longitud (de mayor a menor) como semillas.
##   3. Cada TAD no asignado inicia un nuevo clúster y se define como representante.
##   4. Se incorporan al clúster los TADs no asignados que superan el umbral
##      de solapamiento recíproco respecto al representante.
##   5. Se repite hasta asignar todos.
##
## La diferencia principal respecto a componentes conexas es que la asignación
## exige similitud directa con el representante del clúster.
## ============================================================================

#' Agrupamiento greedy anti-chaining basado en similitud al representante.
#'
#' @param tad_unique data.table TADs (chr,start,end,ID).
#' @param thr_recip  umbral de reciprocal_overlap contra el representante.
#' @param seed_by    criterio de semilla: "length" (más largo primero, por
#'                   defecto) o "degree" (más conexiones primero).
#' @return lista:
#'   membership : data.table(ID, cluster_id, representative_id, is_representative)
#'   clusters   : data.table por cluster (cluster_id, representative_id, n,
#'                consensus/union coords, longitudes)
#'   summary    : data.table con n_input, n_clusters, max_cluster_size, reduction_ratio
cluster_tads_anti_chaining <- function(tad_unique,
                                       thr_recip = 0.80,
                                       seed_by   = c("length", "degree")) {

  seed_by <- match.arg(seed_by)
  dt <- as.data.table(tad_unique)
  dt[, length := end - start + 1L]

  pairs <- compute_overlap_pairs_metrics(dt)
  edges <- if (nrow(pairs) > 0L) pairs[reciprocal_overlap >= thr_recip] else pairs

  # Mapa de adyacencia válida (sólo pares que superan el umbral recíproco).
  # Lo hacemos bidireccional para poder consultar vecinos de cualquier nodo.
  if (nrow(edges) > 0L) {
    adj <- rbind(
      edges[, .(from = id_i, to = id_j)],
      edges[, .(from = id_j, to = id_i)]
    )
    setkey(adj, from)
  } else {
    adj <- data.table(from = character(), to = character())
    setkey(adj, from)
  }

  # Orden de semillas
  if (seed_by == "length") {
    seed_order <- dt[order(-length)]$ID
  } else { # degree
    deg <- if (nrow(adj) > 0L) adj[, .(d = .N), by = from] else
      data.table(from = character(), d = integer())
    deg_map <- setNames(deg$d, deg$from)
    dt[, deg := ifelse(ID %in% names(deg_map), deg_map[ID], 0L)]
    seed_order <- dt[order(-deg, -length)]$ID
  }

  assigned <- setNames(rep(FALSE, nrow(dt)), dt$ID)
  cluster_of <- setNames(rep(NA_integer_, nrow(dt)), dt$ID)
  rep_of     <- setNames(rep(NA_character_, nrow(dt)), dt$ID)

  cid <- 0L
  for (seed in seed_order) {
    if (assigned[[seed]]) next   # ya asignado a un cluster previo
    cid <- cid + 1L

    # La semilla actúa como representante del nuevo clúster.
    assigned[[seed]]   <- TRUE
    cluster_of[[seed]] <- cid
    rep_of[[seed]]     <- seed

    # Candidatos: vecinos del representante que superan el umbral recíproco.
    nb <- adj[.(seed), to, nomatch = 0L]
    if (length(nb)) {
      nb <- unique(nb)
      nb <- nb[!assigned[nb]]    # únicamente los elementos aún no asignados
      if (length(nb)) {
        assigned[nb]   <- TRUE
        cluster_of[nb] <- cid
        rep_of[nb]     <- seed
      }
    }
  }

  membership <- data.table(
    ID                = names(cluster_of),
    cluster_id        = as.integer(cluster_of),
    representative_id = unname(rep_of[names(cluster_of)])
  )
  membership[, is_representative := (ID == representative_id)]

  # Estadísticas por clúster.
  d <- merge(dt, membership, by = "ID")
  clusters <- d[, .(
    representative_id = representative_id[1L],
    n                 = .N,
    consensus_start   = max(start),
    consensus_end     = min(end),
    union_start       = min(start),
    union_end         = max(end),
    len_min           = min(length),
    len_max           = max(length)
  ), by = cluster_id][order(cluster_id)]
  clusters[, consensus_length := pmax(consensus_end - consensus_start, 0L)]
  clusters[, union_length     := union_end - union_start]

  n_input <- nrow(dt)
  summary_dt <- data.table(
    method           = "anti_chaining_greedy",
    seed_by          = seed_by,
    thr_recip        = thr_recip,
    n_input          = n_input,
    n_clusters       = nrow(clusters),
    max_cluster_size = max(clusters$n),
    reduction_ratio  = 1 - nrow(clusters) / n_input
  )

  list(
    membership = membership[],
    clusters   = clusters[],
    summary    = summary_dt[]
  )
}

#' Compara componentes conexas vs anti-chaining a un mismo umbral.
#'
#' @param tad_unique data.table TADs.
#' @param thr_recip  umbral recíproco común.
#' @return data.table con métricas lado a lado de ambos métodos.
compare_clustering_methods <- function(tad_unique, thr_recip = 0.80) {
  # Componentes conexas reutilizando la lógica de inspección.
  cc <- inspect_tad_clusters(tad_unique, thr_recip = thr_recip,
                             use_similarity = FALSE, top_n = 0)
  cc_sizes <- cc$sizes
  cc_row <- data.table(
    method           = "connected_components",
    thr_recip        = thr_recip,
    n_clusters       = nrow(cc_sizes),
    max_cluster_size = max(cc_sizes$n),
    reduction_ratio  = 1 - nrow(cc_sizes) / nrow(as.data.table(tad_unique))
  )

  # Agrupamiento anti-chaining.
  ac <- cluster_tads_anti_chaining(tad_unique, thr_recip = thr_recip,
                                   seed_by = "length")
  ac_row <- data.table(
    method           = "anti_chaining_greedy",
    thr_recip        = thr_recip,
    n_clusters       = ac$summary$n_clusters,
    max_cluster_size = ac$summary$max_cluster_size,
    reduction_ratio  = ac$summary$reduction_ratio
  )

  rbind(cc_row, ac_row)
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("tad_threshold_sweep.R")          # opcional (métricas)
## source("tad_cluster_tools.R")            # este archivo
##
## ## (1) INSPECCIÓN a umbral 0.80 (recíproco puro)
## insp <- inspect_tad_clusters(tad_unique, thr_recip = 0.80, top_n = 5)
## print(insp$sizes[1:10])      # tamaños de los 10 clusters mayores
## print(insp$cohesion[1:10])   # diagnóstico de chaining por cluster
## # Inspección del clúster de mayor tamaño:
## big_id <- insp$sizes$cluster_id[1]
## print(show_cluster(insp, big_id))
## # En 'cohesion', valores bajos de core_fraction o len_ratio elevado sugieren
## # posible encadenamiento transitivo.
##
## ## (2) Agrupamiento anti-chaining a umbral 0.80
## ac <- cluster_tads_anti_chaining(tad_unique, thr_recip = 0.80)
## print(ac$summary)
##
## ## Comparar ambos métodos lado a lado:
## print(compare_clustering_methods(tad_unique, thr_recip = 0.80))
## ============================================================================
