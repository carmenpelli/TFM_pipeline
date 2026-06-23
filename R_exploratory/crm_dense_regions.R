################################################################################
##  FASE 2 — DETECCIÓN DE REGIONES REGULATORIAS DENSAS / CANDIDATOS A SUPER-ENHANCER
##  --------------------------------------------------------------------------
##  Sobre los CRMs YA NO REDUNDANTES (salida de reduce_redundancy_crms),
##  agrupa los que están a < dist_thr (def. 12.5 kb) entre elementos
##  CONSECUTIVOS, según el criterio del documento del TFM
##  (DOI 10.1080/15592294.2018.1514231).
##
##  IMPORTANTE — terminología (del documento):
##    Los grupos resultantes se etiquetan como
##      "structural super-enhancer candidates" / "dense regulatory regions".
##    NUNCA como super-enhancers confirmados (no hay H3K27ac/MED1/BRD4/ChIP-seq).
##
##  A diferencia de la reducción de redundancia, aquí el ENCADENAMIENTO lineal
##  es DELIBERADO: un candidato a super-enhancer es, por definición, una cadena
##  de enhancers próximos consecutivos a lo largo del cromosoma. Se ordena por
##  posición y se corta donde el gap supera el umbral.
##
##  Dos métricas de distancia (se calculan AMBAS para comparar):
##    - gap   : hueco entre bordes = start(siguiente) - end(anterior).
##              Negativo si se solapan (=> distancia efectiva 0). Fiel al paper.
##    - center: distancia entre puntos medios (centro a centro).
##
##  Coordenada usada: repr_start/repr_end (consensus/recommended) del CRM reducido.
##
##  Dependencias: data.table
################################################################################

suppressPackageStartupMessages({
  library(data.table)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

## ============================================================================
## Núcleo: encadenamiento por proximidad sobre elementos ordenados por posición.
## ============================================================================

#' Agrupa CRMs consecutivos cuya distancia al anterior < dist_thr.
#'
#' @param dt       data.table con al menos: ID(o cluster_id), chr, s, e
#'                 (s,e = coordenadas a usar; aquí repr_start/repr_end).
#' @param dist_thr umbral de distancia (pb). Def. 12500.
#' @param metric   "gap" (borde-borde) o "center" (centro-centro).
#' @return vector entero region_id (mismo orden que dt YA ORDENADO por chr,s).
#'         Cada salto de cromosoma o gap>=dist_thr inicia una región nueva.
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
      # gap respecto al elemento anterior (ya ordenados por s).
      # Usamos el MÁXIMO end visto hasta ahora en la región para no romper
      # cadenas por un elemento previo más corto (manejo de solapamientos
      # anidados): el "frente" de la región es el mayor end alcanzado.
      cur_max_end <- e[1L]
      for (i in 2L:n) {
        if (chr[i] != chr[i-1L]) {
          rid <- rid + 1L
          cur_max_end <- e[i]
        } else {
          gap <- s[i] - cur_max_end           # hueco al frente de la región
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
          d <- ctr[i] - ctr[i-1L]             # centro a centro consecutivos
          if (d >= dist_thr) rid <- rid + 1L
        }
        region_id[i] <- rid
      }
    }
  }
  region_id
}

## ============================================================================
## FUNCIÓN PRINCIPAL: detect_dense_regions()
## ============================================================================

#' Detecta regiones densas / candidatos a super-enhancer por proximidad.
#'
#' @param crm_reduced data.table de CRMs no redundantes (salida de
#'        reduce_redundancy_crms). Debe tener cluster_id, chr, repr_start, repr_end.
#' @param dist_thr    umbral de proximidad (pb). Def. 12500 (criterio del documento).
#' @param metric      "gap" o "center".
#' @param min_members nº mínimo de CRMs para considerar "región densa"
#'        (def. 2: al menos dos elementos próximos. ROSE original suele exigir
#'        un mínimo de constituyentes; ajústalo si tu criterio difiere).
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

  # Coordenada a usar = consensus/recommended (repr_*)
  dt[, s := repr_start]
  dt[, e := repr_end]

  # Orden por cromosoma y posición (imprescindible para el encadenamiento)
  setorder(dt, chr, s, e)

  .msg("Detección de regiones densas: ", nrow(dt), " CRMs no redundantes | ",
       "métrica=", metric, " | umbral=", dist_thr, " pb.")

  dt[, region_id := .chain_by_proximity(.SD, dist_thr = dist_thr, metric = metric)]

  # Tamaño de cada región y marca de candidato
  dt[, region_size := .N, by = region_id]
  dt[, is_dense_candidate := region_size >= min_members]

  ## --- Estadísticas + coordenadas por región densa -------------------------
  regions <- dt[is_dense_candidate == TRUE, {
    n   <- .N
    rs  <- min(s); re <- max(e)               # extensión total de la región
    span <- re - rs + 1L
    .(chr            = chr[1L],
      region_start   = rs,
      region_end     = re,
      region_span    = span,
      n_crms         = n,
      # densidad: CRMs por kb de región
      density_per_kb = n / (span / 1000),
      # tamaño medio de gap interno (sólo informativo)
      mean_crm_len   = mean(e - s + 1L),
      # IDs de los CRMs reducidos que la forman (trazabilidad)
      member_cluster_ids = paste(cluster_id, collapse = ";"))
  }, by = region_id][order(chr, region_start)]

  ## --- Resumen global ------------------------------------------------------
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
  .msg("Candidatos a super-enhancer / regiones densas: ", nrow(regions),
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
## COMPARACIÓN DE LAS DOS MÉTRICAS (gap vs center)
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
## # red_crm$crm_reduced son los 530.974 CRMs no redundantes de chr8.
##
## # Comparar ambas métricas de distancia:
## cmp <- compare_distance_metrics(red_crm$crm_reduced, dist_thr = 12500)
## print(cmp$comparison)
## #   Compara n_dense_candidates, n_crms_in_dense, max_region_size/span entre
## #   "gap" y "center". Recuerda: gap es fiel al paper y robusto a tamaños;
## #   center tiende a fragmentar regiones con CRMs grandes.
##
## # Resultados de la métrica elegida (p.ej. gap):
## head(cmp$gap$regions)         # regiones densas con coords y estadísticas
## cmp$gap$marked[is_dense_candidate == TRUE]   # CRMs marcados
##
## # Guardar:
## # fwrite(cmp$gap$regions, "results/dense_regions_chr8_gap.tsv.gz", sep="\t")
## ============================================================================
