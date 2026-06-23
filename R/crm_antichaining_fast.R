################################################################################
##  crm_antichaining_fast.R — Implementación optimizada de anti-chaining
##  --------------------------------------------------------------------------
##  Este script implementa un algoritmo greedy de anti-chaining a partir de
##  aristas de redundancia previamente filtradas. La implementación utiliza
##  índices enteros, listas de adyacencia y vectores planos para reducir el
##  coste computacional en conjuntos grandes de CRMs.
##
##  La lógica de agrupamiento se mantiene equivalente a la versión basada en
##  representantes; la diferencia principal está en la eficiencia de la
##  implementación.
##
##  Dependencias:
##    data.table.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

#' Anti-chaining greedy optimizado a partir de aristas ya filtradas.
#'
#' @param dt    data.table de CRMs (chr,start,end,ID). Debe contener ID.
#' @param edges data.table de aristas (id_i,id_j) que pasan el criterio.
#' @param seed_metric vector OPCIONAL de prioridad de semilla por fila de dt
#'        (mayor = antes). Por defecto la longitud (end-start+1): el elemento
#'        elemento de mayor longitud se utiliza como semilla. Debe tener length == nrow(dt).
#' @param verbose si TRUE, informa progreso cada cierto número de semillas.
#' @return data.table membership (ID, cluster_id, representative_id).
anti_chaining_fast <- function(dt, edges, seed_metric = NULL, verbose = TRUE) {

  dt <- as.data.table(dt)
  ids <- dt$ID
  n   <- length(ids)

  ## --- 1) Mapear IDs a índices enteros 1..n --------------------------------
  ## match() una sola vez; a partir de aquí todo opera con enteros.
  if (nrow(edges) > 0L) {
    ei <- match(edges$id_i, ids)
    ej <- match(edges$id_j, ids)
    # Descartar aristas con IDs no presentes, si existieran.
    ok <- !is.na(ei) & !is.na(ej)
    ei <- ei[ok]; ej <- ej[ok]
  } else {
    ei <- integer(0); ej <- integer(0)
  }

  ## --- 2) Adyacencia como lista indexada por entero ------------------------
  ## Duplicamos cada arista en ambos sentidos y hacemos UN solo split.
  ## adj[[v]] = vector de vecinos de v. Construcción O(E), sin consultas
  ## repetidas dentro del bucle.
  if (length(ei) > 0L) {
    from <- c(ei, ej)
    to   <- c(ej, ei)
    # split por 'from'; factor con niveles 1..n para incluir nodos sin aristas
    adj <- split(to, factor(from, levels = seq_len(n)))
  } else {
    adj <- vector("list", n)
  }

  ## --- 3) Orden de semillas (mayor métrica primero) ------------------------
  if (is.null(seed_metric)) {
    seed_metric <- dt$end - dt$start + 1L      # longitud
  }
  stopifnot(length(seed_metric) == n)
  seed_order <- order(seed_metric, decreasing = TRUE)

  ## --- 4) Bucle greedy con estado plano ------------------------------------
  assigned   <- logical(n)                      # FALSE por defecto
  cluster_of <- integer(n)                      # 0 = sin asignar
  rep_of     <- integer(n)                      # índice del representante

  cid <- 0L
  report_every <- max(1L, as.integer(n / 20))   # ~20 mensajes de progreso

  for (k in seq_len(n)) {
    seed <- seed_order[k]
    if (assigned[seed]) next

    cid <- cid + 1L
    assigned[seed]   <- TRUE
    cluster_of[seed] <- cid
    rep_of[seed]     <- seed

    nb <- adj[[seed]]
    if (length(nb)) {
      # Vecinos aún no asignados mediante indexación lógica por posición.
      nb <- nb[!assigned[nb]]
      if (length(nb)) {
        assigned[nb]   <- TRUE
        cluster_of[nb] <- cid
        rep_of[nb]     <- seed
      }
    }

    if (verbose && (k %% report_every == 0L))
      .msg("  anti-chaining: ", k, "/", n, " semillas procesadas, ",
           cid, " clusters.")
  }

  data.table(
    ID                = ids,
    cluster_id        = cluster_of,
    representative_id = ids[rep_of]
  )
}

## ============================================================================
## VERIFICACIÓN DE EQUIVALENCIA (opcional, para validar contra la versión de referencia)
## ----------------------------------------------------------------------------
## Comprueba que la versión optimizada produce la misma partición que la original.
## La comparación se basa en la partición obtenida, no en los
## cluster_id concretos (que pueden numerarse distinto si hay empates de orden).
## ============================================================================

#' Comprueba si dos tablas de pertenencia representan la misma partición.
#' @return TRUE si agrupan idénticamente los IDs.
same_partition <- function(m1, m2) {
  setkey(m1, ID); setkey(m2, ID)
  m <- merge(m1[, .(ID, c1 = cluster_id)],
             m2[, .(ID, c2 = cluster_id)], by = "ID")
  # Dos particiones son iguales si el mapeo c1<→c2 es biyectivo y consistente.
  ok1 <- m[, uniqueN(c2) == 1L, by = c1][, all(V1)]
  ok2 <- m[, uniqueN(c1) == 1L, by = c2][, all(V1)]
  ok1 && ok2
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("crm_explore.R")
## source("crm_compare_strategies.R")     # versión de referencia (para validar)
## source("crm_antichaining_fast.R")      # este archivo
##
## reg <- subset_region(enh_unique, "chr8", 125e6, 130e6)
## edges <- compute_crm_edges_chunked(reg, thr_recip = 0.50, criterion="composite")
##
## # Versión optimizada:
## system.time(m_fast <- anti_chaining_fast(reg, edges))
##
## # (Opcional) validar equivalencia contra la versión de referencia:
## # m_slow <- anti_chaining_from_edges(reg, edges)
## # same_partition(m_fast, m_slow)   # debe ser TRUE
##
## # Tamaños de cluster:
## m_fast[, .N, by = cluster_id][order(-N)][1:10]
## ============================================================================
