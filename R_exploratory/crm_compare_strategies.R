################################################################################
##  COMPARACIÓN DE ESTRATEGIAS ANTI-CHAINING PARA CRMs
##  --------------------------------------------------------------------------
##  Compara, sobre una región, dos formas de evitar el chaining detectado:
##
##    Estrategia A — SEPARAR POR ESCALA:
##      Aparta los CRM grandes (> size_thr) como categoría "regiones densas /
##      candidatos a super-enhancer" y aplica el clustering SOLO a los CRM
##      normales. size_thr por defecto = 12500 pb (escala del documento del TFM,
##      DOI 10.1080/15592294.2018.1514231; reutilizada como frontera de tamaño).
##      Dentro de los normales puede usar componentes conexas o anti-chaining.
##
##    Estrategia B — ANTI-CHAINING SOBRE TODO:
##      Clustering greedy anti-chaining sobre el conjunto completo, sin separar.
##      Rompe cadenas exigiendo similitud recíproca CONTRA EL REPRESENTANTE.
##
##  Ambas reutilizan compute_crm_edges_chunked() (crm_explore.R) para no
##  materializar los millones de pares.
##
##  Dependencias: data.table, GenomicRanges, igraph
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
## Utilidad: tamaños de cluster a partir de un membership.
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
## CLUSTERING ANTI-CHAINING (greedy por representante) — versión escalable.
## Recibe ARISTAS ya filtradas (no recalcula pares). Cada CRM no asignado, en
## orden de mayor a menor tamaño, abre un cluster y absorbe a sus vecinos no
## asignados (vecinos = pares que ya pasaron el criterio recíproco contra él).
## ============================================================================

#' Anti-chaining greedy a partir de aristas ya filtradas.
#'
#' @param dt    data.table de CRMs (chr,start,end,ID).
#' @param edges data.table de aristas (id_i,id_j) que pasan el criterio.
#' @return data.table membership (ID, cluster_id, representative_id).
anti_chaining_from_edges <- function(dt, edges) {
  dt <- as.data.table(dt)
  dt[, length := end - start + 1L]

  # Adyacencia bidireccional (sólo aristas válidas)
  if (nrow(edges) > 0L) {
    adj <- rbind(
      edges[, .(from = id_i, to = id_j)],
      edges[, .(from = id_j, to = id_i)]
    )
    setkey(adj, from)
  } else {
    adj <- data.table(from = character(), to = character()); setkey(adj, from)
  }

  # Semillas: de mayor a menor tamaño (el grande define el cluster)
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
## Componentes conexas a partir de aristas ya filtradas.
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

  ## --- Aristas del conjunto COMPLETO (una vez, criterio compuesto) ----------
  .msg("Calculando aristas del conjunto completo...")
  edges_all <- compute_crm_edges_chunked(
    dt, thr_recip = thr_recip, jaccard_thr = jaccard_thr,
    simpson_thr = simpson_thr, criterion = "composite",
    keep_metrics = FALSE, chunk_size = chunk_size
  )

  rows <- list()
  memberships <- list()

  ## === BASELINE: componentes conexas sobre todo (lo que ya vimos) ==========
  m_cc <- connected_from_edges(dt, edges_all)
  st <- .cluster_size_stats(m_cc, n_input)
  rows[["baseline_connected"]] <- cbind(
    data.table(strategy = "baseline_connected", n_input = n_input,
               n_large_aside = 0L), st)
  memberships[["baseline_connected"]] <- m_cc

  ## === ESTRATEGIA B: anti-chaining sobre todo ==============================
  .msg("Estrategia B: anti-chaining sobre el conjunto completo...")
  m_ac <- anti_chaining_from_edges(dt, edges_all)
  st <- .cluster_size_stats(m_ac[, .(cluster_id)], n_input)
  rows[["B_antichaining_all"]] <- cbind(
    data.table(strategy = "B_antichaining_all", n_input = n_input,
               n_large_aside = 0L), st)
  memberships[["B_antichaining_all"]] <- m_ac

  ## === ESTRATEGIA A: separar por escala ====================================
  ## Aparta CRM grandes; clustering (anti-chaining) sólo sobre los normales.
  .msg("Estrategia A: separar por escala (size_thr=", size_thr, ")...")
  large_ids  <- dt[length > size_thr]$ID
  normal_dt  <- dt[length <= size_thr]
  n_large    <- length(large_ids)
  .msg("  CRM grandes apartados: ", n_large,
       " (", round(100 * n_large / n_input, 2), "%). ",
       "Normales: ", nrow(normal_dt), ".")

  # Aristas SÓLO entre CRM normales (recalcular sobre el subconjunto, para no
  # arrastrar los puentes grandes que causaban el chaining)
  edges_norm <- compute_crm_edges_chunked(
    normal_dt, thr_recip = thr_recip, jaccard_thr = jaccard_thr,
    simpson_thr = simpson_thr, criterion = "composite",
    keep_metrics = FALSE, chunk_size = chunk_size
  )

  # A.1: separar + componentes conexas sobre normales
  m_a_cc <- connected_from_edges(normal_dt, edges_norm)
  # los grandes van como singletons (cada uno su propio "cluster"/categoría)
  st <- .cluster_size_stats(m_a_cc[, .(cluster_id)], nrow(normal_dt))
  rows[["A_separate_connected"]] <- cbind(
    data.table(strategy = "A_separate_connected", n_input = n_input,
               n_large_aside = n_large), st)
  memberships[["A_separate_connected_normals"]] <- m_a_cc

  # A.2: separar + anti-chaining sobre normales
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
## EJEMPLO DE USO
## ----------------------------------------------------------------------------
## source("crm_explore.R")
## source("crm_compare_strategies.R")
##
## reg <- subset_region(enh_unique, "chr8", 125e6, 130e6)
##
## cmp <- compare_crm_strategies(reg, thr_recip = 0.50, size_thr = 12500)
## print(cmp$comparison)
## # Columnas clave:
## #   max_cluster_size -> el detector de chaining. Bajará mucho en A y B
## #                       respecto al baseline (17670).
## #   n_large_aside    -> cuántos CRM grandes se apartaron en la estrategia A.
## #
## # Lectura:
## #   - Si A_separate_* reduce max_cluster_size pero B_antichaining no tanto,
## #     el problema eran los PUENTES grandes -> separar por escala es la clave.
## #   - Si B ya controla el tamaño sin separar, el anti-chaining basta solo.
## #   - Compara n_clusters y reduction_ratio para ver cuánto reduce cada una.
## ============================================================================
