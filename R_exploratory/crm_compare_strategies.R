################################################################################
##  crm_compare_strategies.R — Comparación de estrategias de agrupamiento de CRMs
##  --------------------------------------------------------------------------
##  Este script compara, sobre una región genómica, distintas estrategias para
##  limitar el encadenamiento transitivo durante el agrupamiento de CRMs.
##
##  Estrategia A — Separación por escala:
##    Los CRMs de mayor tamaño se analizan como un conjunto diferenciado, y el
##    agrupamiento se aplica al subconjunto restante. El umbral de tamaño por
##    defecto es 12.500 pb.
##
##  Estrategia B — Anti-chaining global:
##    Se aplica un agrupamiento greedy sobre el conjunto completo, exigiendo
##    similitud suficiente respecto al representante de cada clúster.
##
##  Ambas estrategias reutilizan compute_crm_edges_chunked() para calcular las
##  aristas filtradas sin materializar todos los pares simultáneamente.
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

## ----------------------------------------------------------------------------
## Utilidad para resumir tamaños de clúster a partir de una tabla de pertenencia.
## ----------------------------------------------------------------------------
.cluster_size_stats <- function(membership, n_input) {
  sizes <- membership[, .(n = .N), by = cluster_id]$n
  data.table(
    n_clusters       = length(sizes),
    n_merged         = sum(sizes > 1L),
    max_cluster_size = max(sizes),
    reduction_ratio  = 1 - length(sizes) / n_input
  )
}

## ============================================================================
## Agrupamiento anti-chaining mediante estrategia greedy por representante.
## Recibe aristas ya filtradas y evita recalcular pares. Los CRMs se recorren
## por tamaño decreciente; cada CRM no asignado inicia un clúster y agrupa sus
## vecinos no asignados.
## ============================================================================

#' Anti-chaining greedy a partir de aristas ya filtradas.
#'
#' @param dt    data.table de CRMs (chr,start,end,ID).
#' @param edges data.table de aristas (id_i,id_j) que pasan el criterio.
#' @return data.table membership con ID, cluster_id y representative_id.
anti_chaining_from_edges <- function(dt, edges) {
  dt <- as.data.table(dt)
  dt[, length := end - start + 1L]

  # Adyacencia bidireccional basada en aristas válidas.
  if (nrow(edges) > 0L) {
    adj <- rbind(
      edges[, .(from = id_i, to = id_j)],
      edges[, .(from = id_j, to = id_i)]
    )
    setkey(adj, from)
  } else {
    adj <- data.table(from = character(), to = character()); setkey(adj, from)
  }

  # Semillas ordenadas de mayor a menor tamaño.
  seed_order <- dt[order(-length)]$ID

  assigned   <- setNames(rep(FALSE, nrow(dt)), dt$ID)
  cluster_of <- setNames(rep(NA_integer_, nrow(dt)), dt$ID)
  rep_of     <- setNames(rep(NA_character_, nrow(dt)), dt$ID)

  cid <- 0L
  for (seed in seed_order) {
    if (assigned[[seed]]) next
    cid <- cid + 1L
    assigned[[seed]]   <- TRUE
    cluster_of[[seed]] <- cid
    rep_of[[seed]]     <- seed

    nb <- adj[.(seed), to, nomatch = 0L]
    if (length(nb)) {
      nb <- unique(nb)
      nb <- nb[!assigned[nb]]
      if (length(nb)) {
        assigned[nb]   <- TRUE
        cluster_of[nb] <- cid
        rep_of[nb]     <- seed
      }
    }
  }

  data.table(
    ID                = names(cluster_of),
    cluster_id        = as.integer(cluster_of),
    representative_id = unname(rep_of[names(cluster_of)])
  )
}

## ----------------------------------------------------------------------------
## Componentes conexas a partir de aristas previamente filtradas.
## ----------------------------------------------------------------------------
connected_from_edges <- function(dt, edges) {
  all_ids <- as.data.table(dt)$ID
  edf <- if (nrow(edges) > 0L)
    as.data.frame(edges[, .(id_i, id_j)], stringsAsFactors = FALSE) else
    data.frame(id_i = character(), id_j = character(), stringsAsFactors = FALSE)
  g <- graph_from_data_frame(edf, directed = FALSE,
                             vertices = data.frame(name = all_ids))
  comp <- igraph::components(g)
  data.table(ID = names(comp$membership),
             cluster_id = as.integer(comp$membership))
}

## ============================================================================
## COMPARACIÓN DE ESTRATEGIAS
## ============================================================================

#' Compara estrategias de clustering de CRMs sobre una región.
#'
#' @param dt           data.table de CRMs (región).
#' @param thr_recip    umbral recíproco (def. 0.50, criterio del documento).
#' @param jaccard_thr,simpson_thr umbrales del criterio compuesto.
#' @param size_thr     umbral de tamaño para apartar CRM grandes (def. 12500).
#' @param chunk_size   tamaño de bloque.
#' @return lista con:
#'   comparison : data.table con una fila por estrategia y sus métricas.
#'   memberships: lista de los membership de cada estrategia (para inspección).
compare_crm_strategies <- function(dt,
                                   thr_recip   = 0.50,
                                   jaccard_thr = 0.70,
                                   simpson_thr = 0.99,
                                   size_thr    = 12500,
                                   chunk_size  = 5e6) {

  dt <- as.data.table(copy(dt))
  dt[, length := end - start + 1L]
  n_input <- nrow(dt)

  ## --- Aristas del conjunto completo, calculadas una sola vez ----------------
  .msg("Calculando aristas del conjunto completo...")
  edges_all <- compute_crm_edges_chunked(
    dt, thr_recip = thr_recip, jaccard_thr = jaccard_thr,
    simpson_thr = simpson_thr, criterion = "composite",
    keep_metrics = FALSE, chunk_size = chunk_size
  )

  rows <- list()
  memberships <- list()

  ## === Escenario de referencia: componentes conexas sobre el conjunto completo ===
  m_cc <- connected_from_edges(dt, edges_all)
  st <- .cluster_size_stats(m_cc, n_input)
  rows[["baseline_connected"]] <- cbind(
    data.table(strategy = "baseline_connected", n_input = n_input,
               n_large_aside = 0L), st)
  memberships[["baseline_connected"]] <- m_cc

  ## === Estrategia B: anti-chaining sobre el conjunto completo ===============
  .msg("Estrategia B: anti-chaining sobre el conjunto completo...")
  m_ac <- anti_chaining_from_edges(dt, edges_all)
  st <- .cluster_size_stats(m_ac[, .(cluster_id)], n_input)
  rows[["B_antichaining_all"]] <- cbind(
    data.table(strategy = "B_antichaining_all", n_input = n_input,
               n_large_aside = 0L), st)
  memberships[["B_antichaining_all"]] <- m_ac

  ## === Estrategia A: separación por escala =================================
  ## Se separan los CRMs de mayor tamaño y se agrupa el subconjunto restante.
  .msg("Estrategia A: separar por escala (size_thr=", size_thr, ")...")
  large_ids  <- dt[length > size_thr]$ID
  normal_dt  <- dt[length <= size_thr]
  n_large    <- length(large_ids)
  .msg("  CRM grandes apartados: ", n_large,
       " (", round(100 * n_large / n_input, 2), "%). ",
       "Normales: ", nrow(normal_dt), ".")

  # Aristas únicamente entre CRMs del subconjunto restante.
  edges_norm <- compute_crm_edges_chunked(
    normal_dt, thr_recip = thr_recip, jaccard_thr = jaccard_thr,
    simpson_thr = simpson_thr, criterion = "composite",
    keep_metrics = FALSE, chunk_size = chunk_size
  )

  # A.1: separación por escala y componentes conexas en el subconjunto restante.
  m_a_cc <- connected_from_edges(normal_dt, edges_norm)
  # Los CRMs separados por tamaño se tratan como singletons.
  st <- .cluster_size_stats(m_a_cc[, .(cluster_id)], nrow(normal_dt))
  rows[["A_separate_connected"]] <- cbind(
    data.table(strategy = "A_separate_connected", n_input = n_input,
               n_large_aside = n_large), st)
  memberships[["A_separate_connected_normals"]] <- m_a_cc

  # A.2: separación por escala y anti-chaining en el subconjunto restante.
  m_a_ac <- anti_chaining_from_edges(normal_dt, edges_norm)
  st <- .cluster_size_stats(m_a_ac[, .(cluster_id)], nrow(normal_dt))
  rows[["A_separate_antichaining"]] <- cbind(
    data.table(strategy = "A_separate_antichaining", n_input = n_input,
               n_large_aside = n_large), st)
  memberships[["A_separate_antichaining_normals"]] <- m_a_ac

  comparison <- rbindlist(rows, use.names = TRUE)

  list(comparison = comparison[], memberships = memberships,
       large_ids = large_ids)
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("crm_explore.R")
## source("crm_compare_strategies.R")
##
## reg <- subset_region(enh_unique, "chr8", 125e6, 130e6)
##
## cmp <- compare_crm_strategies(reg, thr_recip = 0.50, size_thr = 12500)
## print(cmp$comparison)
## # Columnas principales:
## #   max_cluster_size permite detectar agrupamientos excesivamente grandes.
## #   n_large_aside indica cuántos CRMs se separan por tamaño en la estrategia A.
## #
## # Interpretación:
## #   - Una reducción marcada de max_cluster_size tras separar por escala
## #     sugiere que los CRMs de mayor tamaño contribuían al encadenamiento.
## #   - Si la estrategia B controla el tamaño de los clústeres, el anti-chaining
## #     puede ser suficiente sin separación previa.
## #   - n_clusters y reduction_ratio permiten comparar la intensidad de reducción.
## ============================================================================
