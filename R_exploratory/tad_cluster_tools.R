################################################################################
##  INSPECCIÓN DE CLUSTERS Y CLUSTERING ANTI-CHAINING PARA TADs
##  --------------------------------------------------------------------------
##  Dos utilidades para la fase exploratoria de reducción de redundancia de TADs:
##
##    (1) inspect_tad_clusters()  -> DIAGNÓSTICO (no cambia metodología).
##        Construye los clusters a un umbral dado (componentes conexas, igual
##        que tu pipeline) y te deja examinar los más grandes: qué TADs los
##        forman, sus coordenadas, longitudes y dispersión. Sirve para ver si
##        un cluster grande son anotaciones equivalentes o una CADENA espuria.
##
##    (2) cluster_tads_anti_chaining() -> VARIANTE METODOLÓGICA ALTERNATIVA.
##        Rompe el encadenamiento transitivo exigiendo que cada miembro supere
##        el umbral recíproco CONTRA EL REPRESENTANTE del cluster, no sólo
##        contra algún vecino. Es DISTINTA del clustering por componentes
##        conexas de tu pipeline de CRMs: úsala sólo como opción consciente.
##
##  Depende de compute_overlap_pairs_metrics() definida en tad_threshold_sweep.R
##  (o equivalente). Aquí se incluye una copia mínima por si se usa suelto.
##
##  Dependencias: data.table, GenomicRanges, igraph
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(igraph)
})

## ----------------------------------------------------------------------------
## (copia mínima por si este archivo se usa de forma independiente)
## Si ya cargaste tad_threshold_sweep.R, esta definición simplemente la
## sobrescribe con una idéntica.
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
## (1) INSPECCIÓN DE CLUSTERS (DIAGNÓSTICO, mismo método que el pipeline)
## ============================================================================

#' Construye clusters por componentes conexas a un umbral recíproco dado
#' y permite inspeccionar los más grandes en detalle.
#'
#' @param tad_unique data.table TADs (chr,start,end,ID).
#' @param thr_recip  umbral de reciprocal_overlap para crear arista.
#' @param use_similarity si TRUE, replica el criterio estilo CRM:
#'        arista si recíproco>=thr_recip Y (jaccard>=jaccard_thr O simpson>=simpson_thr).
#'        Si FALSE (por defecto), arista por recíproco puro (recomendado en TADs).
#' @param jaccard_thr,simpson_thr umbrales para use_similarity=TRUE.
#' @param top_n      nº de clusters más grandes a devolver en detalle.
#' @return lista con:
#'   membership : data.table(ID, cluster_id)
#'   sizes      : data.table(cluster_id, n) ordenada desc
#'   top        : data.table con los TADs de los top_n clusters más grandes,
#'                con coordenadas, longitud y métricas de cohesión por cluster
#'   cohesion   : data.table por cluster con diagnóstico de chaining
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

  # Grafo + componentes conexas (idéntico al pipeline)
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

  # IDs de los top_n clusters más grandes
  top_ids <- sizes[n > 1L][seq_len(min(top_n, .N))]$cluster_id

  # Detalle de los TADs en esos clusters (coordenadas + longitud)
  d_annot <- merge(dt, membership, by = "ID")
  d_annot[, length := end - start + 1L]

  top <- d_annot[cluster_id %in% top_ids][order(cluster_id, start)]

  ## --- Cohesión / diagnóstico de chaining por cluster ----------------------
  ## Para cada cluster medimos:
  ##   n              : nº de TADs
  ##   span           : extensión total (max(end)-min(start)) del cluster
  ##   core_len       : "núcleo" común (min(end)-max(start)); <=0 => SIN núcleo
  ##   len_min/len_max: longitudes extremas de los TADs del cluster
  ##   len_ratio      : len_max/len_min (escalas mezcladas => jerarquía)
  ##   has_common_core: TRUE si todos los TADs comparten una base común
  ## Un cluster "sano" (anotaciones equivalentes) tiene core_len>0 y len_ratio
  ## moderado. Un cluster por CHAINING suele tener core_len<=0 (sin base común)
  ## y/o len_ratio grande (mezcla de escalas).
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
      # fracción del span cubierta por el núcleo común (1=idénticos, <=0=chaining)
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
    message("El cluster ", cluster, " no está entre los top mostrados; ",
            "contiene ", length(ids), " TADs. Aumenta top_n para verlo.")
  }
  out[]
}

## ============================================================================
## (2) CLUSTERING ANTI-CHAINING (VARIANTE METODOLÓGICA ALTERNATIVA)
## ----------------------------------------------------------------------------
## ¡ATENCIÓN! Esto NO es el clustering por componentes conexas de tu pipeline.
## Es una alternativa que ROMPE el encadenamiento transitivo. Úsala de forma
## consciente y decláralo en la memoria si la adoptas.
##
## Estrategia (greedy por representante):
##   1. Se calculan los pares que superan el umbral recíproco (aristas válidas).
##   2. Se ordenan los TADs por longitud (de mayor a menor) como semillas.
##   3. Se recorre la lista: cada TAD aún no asignado abre un cluster nuevo y
##      se convierte en REPRESENTANTE.
##   4. Se asignan a ese cluster todos los TADs no asignados que superen el
##      umbral recíproco CONTRA EL REPRESENTANTE (no contra cualquier vecino).
##   5. Se repite hasta asignar todos.
##
## Diferencia clave vs. componentes conexas:
##   - Componentes conexas: A-B-C-D se funden aunque A y D no se parezcan.
##   - Anti-chaining: D sólo entra si se parece al REPRESENTANTE del cluster,
##     evitando cadenas de TADs de escalas distintas.
## ============================================================================

#' Clustering greedy anti-chaining basado en similitud al representante.
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

    # El seed es el representante del nuevo cluster
    assigned[[seed]]   <- TRUE
    cluster_of[[seed]] <- cid
    rep_of[[seed]]     <- seed

    # Candidatos: vecinos del representante que superan el umbral recíproco.
    # (adj ya sólo contiene pares >= thr_recip, así que todos sus vecinos valen.)
    nb <- adj[.(seed), to, nomatch = 0L]
    if (length(nb)) {
      nb <- unique(nb)
      nb <- nb[!assigned[nb]]    # sólo los aún no asignados
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

  # Estadísticas por cluster
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
  # Componentes conexas (vía inspect, reutilizando su lógica)
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

  # Anti-chaining
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
## EJEMPLO DE USO
## ----------------------------------------------------------------------------
## source("tad_threshold_sweep.R")          # opcional (métricas)
## source("tad_cluster_tools.R")            # este archivo
##
## ## (1) INSPECCIÓN a umbral 0.80 (recíproco puro)
## insp <- inspect_tad_clusters(tad_unique, thr_recip = 0.80, top_n = 5)
## print(insp$sizes[1:10])      # tamaños de los 10 clusters mayores
## print(insp$cohesion[1:10])   # diagnóstico de chaining por cluster
## # Mirar el cluster más grande en detalle:
## big_id <- insp$sizes$cluster_id[1]
## print(show_cluster(insp, big_id))
## # En 'cohesion': has_common_core=FALSE o core_fraction<=0 o len_ratio grande
## # => el cluster es probablemente una CADENA, no anotaciones equivalentes.
##
## ## (2) ANTI-CHAINING a umbral 0.80
## ac <- cluster_tads_anti_chaining(tad_unique, thr_recip = 0.80)
## print(ac$summary)
##
## ## Comparar ambos métodos lado a lado:
## print(compare_clustering_methods(tad_unique, thr_recip = 0.80))
## ============================================================================
