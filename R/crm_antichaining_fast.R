################################################################################
##  ANTI-CHAINING OPTIMIZADO (R PURO) — escalable a ~1.3M CRMs
##  --------------------------------------------------------------------------
##  Misma lógica greedy por representante que anti_chaining_from_edges(), pero
##  reescrita para eliminar los cuellos de botella que la hacían lenta:
##
##    1. IDs -> índices ENTEROS (match una sola vez). Búsquedas O(1) por
##       posición en lugar de por nombre de carácter.
##    2. Adyacencia como LISTA indexada por entero (split una sola vez), en
##       lugar de consultar un data.table con clave dentro del bucle.
##    3. Estado en vectores lógicos/enteros planos; sin nombres.
##    4. Semillas en orden de mayor a menor tamaño (idéntico al original).
##
##  El resultado (asignación a clusters) es EQUIVALENTE al anti-chaining previo;
##  sólo cambia la implementación, no el método.
##
##  Dependencias: data.table  (igraph/GenomicRanges no se usan aquí)
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
#'        más grande siembra el cluster. Debe tener length == nrow(dt).
#' @param verbose si TRUE, informa progreso cada cierto nº de semillas.
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
    # Por seguridad, descartar aristas con IDs no presentes (no debería ocurrir)
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
      # vecinos aún no asignados (indexación lógica por posición: rápida)
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
## VERIFICACIÓN DE EQUIVALENCIA (opcional, para validar contra la versión lenta)
## ----------------------------------------------------------------------------
## Comprueba que la versión rápida produce la MISMA partición que la original.
## La comparación es por PARTICIÓN (qué elementos caen juntos), no por los
## cluster_id concretos (que pueden numerarse distinto si hay empates de orden).
## ============================================================================

#' ¿Dos memberships representan la misma partición?
#' @return TRUE si agrupan idénticamente los IDs.
same_partition <- function(m1, m2) {
  setkey(m1, ID); setkey(m2, ID)
  m <- merge(m1[, .(ID, c1 = cluster_id)],
             m2[, .(ID, c2 = cluster_id)], by = "ID")
  # Dos particiones son iguales si el mapeo c1<->c2 es biyectivo y consistente.
  ok1 <- m[, uniqueN(c2) == 1L, by = c1][, all(V1)]
  ok2 <- m[, uniqueN(c1) == 1L, by = c2][, all(V1)]
  ok1 && ok2
}

## ============================================================================
## EJEMPLO DE USO
## ----------------------------------------------------------------------------
## source("crm_explore.R")
## source("crm_compare_strategies.R")     # versión lenta (para validar)
## source("crm_antichaining_fast.R")      # este archivo
##
## reg <- subset_region(enh_unique, "chr8", 125e6, 130e6)
## edges <- compute_crm_edges_chunked(reg, thr_recip = 0.50, criterion="composite")
##
## # Versión rápida:
## system.time(m_fast <- anti_chaining_fast(reg, edges))
##
## # (Opcional) validar equivalencia contra la versión lenta:
## # m_slow <- anti_chaining_from_edges(reg, edges)
## # same_partition(m_fast, m_slow)   # debe ser TRUE
##
## # Tamaños de cluster:
## m_fast[, .N, by = cluster_id][order(-N)][1:10]
## ============================================================================
