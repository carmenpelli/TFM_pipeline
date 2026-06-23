################################################################################
##  crm_dense_regions.R — Detección exploratoria de regiones reguladoras densas
##  --------------------------------------------------------------------------
##  Este script agrupa CRMs no redundantes mediante un criterio de proximidad
##  lineal a lo largo del cromosoma. Fue utilizado como aproximación exploratoria
##  para comparar métricas de distancia y evaluar la formación de regiones
##  reguladoras densas.
##
##  La agrupación se realiza sobre CRMs consecutivos ordenados por posición. Se
##  inicia una nueva región cuando la distancia entre elementos supera el umbral
##  definido por dist_thr.
##
##  Se implementan dos métricas de distancia:
##    - gap: distancia entre bordes de intervalos consecutivos.
##    - center: distancia entre puntos medios consecutivos.
##
##  Las regiones obtenidas se interpretan como agrupaciones estructurales de
##  CRMs, no como super-enhancers funcionalmente validados.
##
##  Coordenadas utilizadas:
##    repr_start y repr_end de los CRMs reducidos.
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

## ============================================================================
## Núcleo: agrupamiento por proximidad sobre elementos ordenados por posición.
## ============================================================================

#' Agrupa CRMs consecutivos cuya distancia al elemento anterior es inferior a dist_thr.
#'
#' @param dt       data.table con al menos: ID(o cluster_id), chr, s, e
#'                 (s,e = coordenadas a usar; aquí repr_start/repr_end).
#' @param dist_thr umbral de distancia (pb). Def. 12500.
#' @param metric   "gap" (borde-borde) o "center" (centro-centro).
#' @return vector entero region_id, en el mismo orden que dt tras ordenar por chr y s.
#'         Cada cambio de cromosoma o distancia >= dist_thr inicia una nueva región.
.chain_by_proximity <- function(dt, dist_thr = 12500, metric = c("gap","center")) {
  metric <- match.arg(metric)
  n <- nrow(dt)
  if (n == 0L) return(integer(0))

  s <- dt$s; e <- dt$e; chr <- dt$chr
  region_id <- integer(n)
  region_id[1L] <- 1L
  rid <- 1L

  if (n >= 2L) {
    if (metric == "gap") {
      # Distancia respecto al elemento anterior, con los intervalos ya ordenados por s.
      # Se utiliza el máximo end alcanzado en la región para manejar
      # adecuadamente intervalos solapados o anidados.
      cur_max_end <- e[1L]
      for (i in 2L:n) {
        if (chr[i] != chr[i-1L]) {
          rid <- rid + 1L
          cur_max_end <- e[i]
        } else {
          gap <- s[i] - cur_max_end           # distancia al límite derecho acumulado de la región
          if (gap >= dist_thr) rid <- rid + 1L
          if (e[i] > cur_max_end) cur_max_end <- e[i]
        }
        region_id[i] <- rid
      }
    } else { # center
      ctr <- (s + e) / 2
      for (i in 2L:n) {
        if (chr[i] != chr[i-1L]) {
          rid <- rid + 1L
        } else {
          d <- ctr[i] - ctr[i-1L]             # distancia entre puntos medios consecutivos
          if (d >= dist_thr) rid <- rid + 1L
        }
        region_id[i] <- rid
      }
    }
  }
  region_id
}

## ============================================================================
## Función principal: detect_dense_regions()
## ============================================================================

#' Detecta regiones reguladoras densas mediante proximidad entre CRMs.
#'
#' @param crm_reduced data.table de CRMs no redundantes (salida de
#'        reduce_redundancy_crms). Debe contener cluster_id, chr, repr_start y repr_end.
#' @param dist_thr    umbral de proximidad (pb). Def. 12500 (criterio del documento).
#' @param metric      "gap" o "center".
#' @param min_members nº mínimo de CRMs para considerar "región densa"
#'        (def. 2: al menos dos elementos próximos).
#' @return lista:
#'   marked  : crm_reduced + columnas region_id, region_size, is_dense_candidate
#'   regions : data.table por región densa con coordenadas y estadísticas
#'   summary : resumen global
detect_dense_regions <- function(crm_reduced,
                                 dist_thr    = 12500,
                                 metric      = c("gap","center"),
                                 min_members = 2L) {
  metric <- match.arg(metric)
  dt <- as.data.table(copy(crm_reduced))
  req <- c("cluster_id", "chr", "repr_start", "repr_end")
  stopifnot(all(req %in% names(dt)))

  # Coordenadas de trabajo: repr_start y repr_end.
  dt[, s := repr_start]
  dt[, e := repr_end]

  # Orden por cromosoma y posición.
  setorder(dt, chr, s, e)

  .msg("Detección de regiones densas: ", nrow(dt), " CRMs no redundantes | ",
       "métrica=", metric, " | umbral=", dist_thr, " pb.")

  dt[, region_id := .chain_by_proximity(.SD, dist_thr = dist_thr, metric = metric)]

  # Tamaño de cada región y marca de región candidata.
  dt[, region_size := .N, by = region_id]
  dt[, is_dense_candidate := region_size >= min_members]

  ## --- Estadísticas y coordenadas por región densa ---------------------------
  regions <- dt[is_dense_candidate == TRUE, {
    n   <- .N
    rs  <- min(s); re <- max(e)               # extensión total de la región
    span <- re - rs + 1L
    .(chr            = chr[1L],
      region_start   = rs,
      region_end     = re,
      region_span    = span,
      n_crms         = n,
      # densidad expresada como CRMs por kb de región
      density_per_kb = n / (span / 1000),
      # longitud media de los CRMs incluidos
      mean_crm_len   = mean(e - s + 1L),
      # identificadores de los CRMs reducidos incluidos
      member_cluster_ids = paste(cluster_id, collapse = ";"))
  }, by = region_id][order(chr, region_start)]

  ## --- Resumen global --------------------------------------------------------
  summary_dt <- data.table(
    metric             = metric,
    dist_thr           = dist_thr,
    min_members        = min_members,
    n_crms_input       = nrow(dt),
    n_regions_total    = dt[, uniqueN(region_id)],
    n_dense_candidates = nrow(regions),
    n_crms_in_dense    = dt[is_dense_candidate == TRUE, .N],
    max_region_size    = if (nrow(regions)) max(regions$n_crms) else 0L,
    max_region_span    = if (nrow(regions)) max(regions$region_span) else 0L
  )
  .msg("Regiones reguladoras densas detectadas: ", nrow(regions),
       " (", summary_dt$n_crms_in_dense, " CRMs agrupados).")

  # Limpiar columnas auxiliares en el marcado de salida
  dt[, c("s","e") := NULL]

  list(
    marked  = dt[],
    regions = regions[],
    summary = summary_dt[]
  )
}

## ============================================================================
## Comparación de las dos métricas de distancia
## ============================================================================

#' Ejecuta detect_dense_regions con ambas métricas y compara los resúmenes.
#'
#' @param crm_reduced data.table de CRMs no redundantes.
#' @param dist_thr,min_members parámetros comunes.
#' @return lista: gap, center (resultados completos) y comparison (resúmenes).
compare_distance_metrics <- function(crm_reduced, dist_thr = 12500,
                                     min_members = 2L) {
  r_gap <- detect_dense_regions(crm_reduced, dist_thr = dist_thr,
                                metric = "gap", min_members = min_members)
  r_ctr <- detect_dense_regions(crm_reduced, dist_thr = dist_thr,
                                metric = "center", min_members = min_members)
  comparison <- rbind(r_gap$summary, r_ctr$summary)
  list(gap = r_gap, center = r_ctr, comparison = comparison[])
}

## ============================================================================
## EJEMPLO DE USO
## ----------------------------------------------------------------------------
## source("crm_dense_regions.R")
##
## # red_crm$crm_reduced contiene los CRMs reducidos del cromosoma analizado.
##
## # Comparar ambas métricas de distancia:
## cmp <- compare_distance_metrics(red_crm$crm_reduced, dist_thr = 12500)
## print(cmp$comparison)
## #   Compara n_dense_candidates, n_crms_in_dense, max_region_size/span entre
## #   "gap" y "center". La métrica gap tiende a ser más robusta a diferencias de tamaño;
## #   center puede fragmentar regiones cuando existen CRMs de mayor longitud.
##
## # Resultados de la métrica seleccionada, por ejemplo gap:
## head(cmp$gap$regions)         # regiones densas con coords y estadísticas
## cmp$gap$marked[is_dense_candidate == TRUE]   # CRMs marcados
##
## # Exportación opcional:
## # fwrite(cmp$gap$regions, "results/dense_regions_chr8_gap.tsv.gz", sep="\t")
## ============================================================================
