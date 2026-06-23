################################################################################
##  crm_explore.R — Exploración de redundancia de CRMs y cálculo de aristas
##  --------------------------------------------------------------------------
##  Este script contiene funciones auxiliares para explorar la redundancia de
##  CRMs a gran escala y calcular pares de regiones solapantes mediante
##  procesamiento por bloques.
##
##  Incluye:
##    - selección de ventanas genómicas contiguas para análisis exploratorios;
##    - cálculo eficiente de aristas de redundancia entre CRMs;
##    - muestreo de pares solapantes para resumir distribuciones de métricas;
##    - análisis de sensibilidad de umbrales.
##
##  Dependencias:
##    data.table, GenomicRanges e igraph.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(igraph)
})

.msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

## ============================================================================
## 1) SELECTOR DE REGIÓN CONTIGUA (muestreo que preserva solapamientos)
## ----------------------------------------------------------------------------
## Se evita el muestreo aleatorio de CRMs, ya que alteraría los solapamientos. Se extrae
## una ventana genómica contigua [from, to], conservando todos los CRMs que
## la solapan, de modo que la estructura de redundancia queda intacta.
## ============================================================================

#' Extrae los CRMs de una ventana genómica contigua.
#'
#' @param dt   data.table con chr,start,end,ID (+ metadatos).
#' @param chr  cromosoma (por ejemplo "chr8").
#' @param from inicio de la ventana (pb).
#' @param to   fin de la ventana (pb).
#' @param mode "overlap" conserva CRMs que solapan la ventana (recomendado);
#'             "contained" sólo los totalmente contenidos.
#' @return data.table subconjunto.
subset_region <- function(dt, chr, from, to, mode = c("overlap", "contained")) {
  mode <- match.arg(mode)
  d <- as.data.table(dt)
  # Renombramos el argumento para evitar colisión con la columna 'chr'
  target_chr <- chr
  if (mode == "overlap") {
    out <- d[chr == target_chr & end >= from & start <= to]
  } else {
    out <- d[chr == target_chr & start >= from & end <= to]
  }
  .msg("Región ", target_chr, ":", from, "-", to, " (", mode, "): ",
       nrow(out), " CRMs.")
  out[]
}

## ============================================================================
## 2) CÁLCULo DE PARES POR BLOQUES CON FILTRADo AL VUELO
## ----------------------------------------------------------------------------
## Procesa los hits de findOverlaps en trozos. Para cada trozo calcula las
## métricas y descarta de inmediato los pares que no cumplen el criterio,
## acumulando sólo las aristas supervivientes. Así nunca se materializan
## los cientos de millones de pares completos.
## ============================================================================

#' Calcula aristas de redundancia de CRMs por bloques (filtrado al vuelo).
#'
#' @param dt          data.table no redundante (chr,start,end,ID).
#' @param thr_recip   umbral de reciprocal_overlap (por defecto 0.50, definido para el análisis).
#' @param jaccard_thr umbral de jaccard (por defecto 0.70).
#' @param simpson_thr umbral de simpson (por defecto 0.99).
#' @param criterion   "composite" (recíproco y (J o S); criterio definido para el análisis)
#'                    o "reciprocal_only" (sólo recíproco; para comparar).
#' @param keep_metrics si TRUE, conserva las métricas en las aristas (útil para
#'                    inspección); si FALSE, sólo id_i,id_j (más ligero).
#' @param chunk_size  número de hits a procesar por bloque (por defecto 5e6).
#' @return data.table de aristas que pasan el criterio:
#'         id_i,id_j (+ métricas si keep_metrics=TRUE).
compute_crm_edges_chunked <- function(dt,
                                      thr_recip   = 0.50,
                                      jaccard_thr = 0.70,
                                      simpson_thr = 0.99,
                                      criterion   = c("composite", "reciprocal_only"),
                                      keep_metrics = FALSE,
                                      chunk_size  = 5e6) {

  criterion <- match.arg(criterion)
  dt <- as.data.table(dt)
  gr <- GRanges(seqnames = dt$chr,
                ranges   = IRanges(start = dt$start, end = dt$end),
                ID       = dt$ID)

  st <- start(gr); en <- end(gr); len <- width(gr); ids <- mcols(gr)$ID

  .msg("findOverlaps sobre ", length(gr), " CRMs...")
  hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
  qh_all <- queryHits(hits); sh_all <- subjectHits(hits)

  # Sólo pares i<j (evita auto-pares y simétricos)
  keep <- qh_all < sh_all
  qh_all <- qh_all[keep]; sh_all <- sh_all[keep]
  n_pairs <- length(qh_all)
  .msg("Pares solapados (i<j): ", n_pairs, ". Procesando por bloques de ",
       format(chunk_size, scientific = FALSE), "...")

  if (n_pairs == 0L) {
    return(data.table(id_i = character(), id_j = character()))
  }

  # Procesamiento por bloques
  n_chunks <- ceiling(n_pairs / chunk_size)
  edge_list <- vector("list", n_chunks)

  for (k in seq_len(n_chunks)) {
    i0 <- (k - 1L) * chunk_size + 1L
    i1 <- min(k * chunk_size, n_pairs)
    qi <- qh_all[i0:i1]; sj <- sh_all[i0:i1]

    ov <- pmin(en[qi], en[sj]) - pmax(st[qi], st[sj]) + 1L
    ov[ov < 0L] <- 0L
    Li <- len[qi]; Lj <- len[sj]

    recip <- pmin(ov / Li, ov / Lj)

    # Filtro recíproco primero (barato y muy selectivo)
    pass <- recip >= thr_recip
    if (criterion == "composite") {
      # Sólo calcular J y S sobre los que ya pasaron recíproco (eficiencia)
      if (any(pass)) {
        jac <- ov[pass] / (Li[pass] + Lj[pass] - ov[pass])
        sim <- ov[pass] / pmin(Li[pass], Lj[pass])
        pass2 <- (jac >= jaccard_thr) | (sim >= simpson_thr)
        idx <- which(pass)[pass2]
      } else {
        idx <- integer(0)
      }
    } else {
      idx <- which(pass)
    }

    if (length(idx)) {
      if (keep_metrics) {
        ovk <- ov[idx]; Lik <- Li[idx]; Ljk <- Lj[idx]
        edge_list[[k]] <- data.table(
          id_i = ids[qi[idx]], id_j = ids[sj[idx]],
          overlap_length = as.integer(ovk),
          len_i = as.integer(Lik), len_j = as.integer(Ljk),
          jaccard = ovk / (Lik + Ljk - ovk),
          simpson = ovk / pmin(Lik, Ljk),
          reciprocal_overlap = pmin(ovk / Lik, ovk / Ljk)
        )
      } else {
        edge_list[[k]] <- data.table(id_i = ids[qi[idx]], id_j = ids[sj[idx]])
      }
    }

    if (k %% 5L == 0L || k == n_chunks)
      .msg("  bloque ", k, "/", n_chunks, " procesado.")
  }

  edges <- rbindlist(edge_list, use.names = TRUE, fill = TRUE)
  .msg("Aristas que pasan el criterio (", criterion, "): ", nrow(edges),
       " de ", n_pairs, " pares (", round(100 * nrow(edges) / n_pairs, 3), "%).")
  edges[]
}

## ============================================================================
## 3) Distribución de métricas sobre una muestra de pares
## ----------------------------------------------------------------------------
## Para describir la distribución sin materializar 212M de pares, se toma una
## muestra aleatoria de pares, no de CRMs, y se calculan sus métricas. Este enfoque
## es estadísticamente válido para describir la distribución de solapamientos.
## ============================================================================

#' Distribución de métricas sobre una muestra aleatoria de pares solapados.
#'
#' @param dt        data.table de CRMs (idealmente una región).
#' @param n_sample  número de pares a muestrear para la distribución (por defecto 2e6).
#' @param probs     cuantiles a reportar.
#' @return data.table con cuantiles de recíproco, jaccard, simpson.
describe_crm_overlap <- function(dt, n_sample = 2e6,
                                 probs = c(0,.1,.25,.5,.75,.9,.95,.99,1)) {
  dt <- as.data.table(dt)
  gr <- GRanges(dt$chr, IRanges(dt$start, dt$end), ID = dt$ID)
  st <- start(gr); en <- end(gr); len <- width(gr)

  hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
  qh <- queryHits(hits); sh <- subjectHits(hits)
  keep <- qh < sh; qh <- qh[keep]; sh <- sh[keep]
  n_pairs <- length(qh)
  if (n_pairs == 0L) { .msg("Sin pares solapados."); return(NULL) }

  # Muestreo de pares si hay demasiados
  if (n_pairs > n_sample) {
    sel <- sample.int(n_pairs, n_sample)
    qh <- qh[sel]; sh <- sh[sel]
    .msg("Muestreando ", n_sample, " de ", n_pairs, " pares para la distribución.")
  }

  ov <- pmin(en[qh], en[sh]) - pmax(st[qh], st[sh]) + 1L
  ov[ov < 0L] <- 0L
  Li <- len[qh]; Lj <- len[sh]
  recip <- pmin(ov / Li, ov / Lj)
  jac   <- ov / (Li + Lj - ov)
  sim   <- ov / pmin(Li, Lj)

  data.table(
    quantile           = names(quantile(recip, probs)),
    reciprocal_overlap = as.numeric(quantile(recip, probs)),
    jaccard            = as.numeric(quantile(jac, probs)),
    simpson            = as.numeric(quantile(sim, probs))
  )[]
}

## ============================================================================
## 4) BARRIDo DE UMBRALES (sobre una región, criterio compuesto de CRMs)
## ----------------------------------------------------------------------------
## Para cada umbral recíproco evalúa el criterio compuesto y cuenta clusters
## (componentes conexas), reportando el diagnóstico de chaining clave:
## max_cluster_size.
## ============================================================================

#' Barre umbrales recíprocos con el criterio compuesto de CRMs sobre una región.
#'
#' @param dt          data.table de CRMs (región).
#' @param recip_grid  umbrales de reciprocal_overlap a probar.
#' @param jaccard_thr,simpson_thr umbrales del criterio compuesto.
#' @param chunk_size  tamaño de bloque para el cálculo.
#' @return data.table: por umbral, n_edges, n_clusters, max_cluster_size, etc.
sweep_crm_thresholds <- function(dt,
                                 recip_grid  = c(0.50, 0.60, 0.70, 0.80, 0.90),
                                 jaccard_thr = 0.70,
                                 simpson_thr = 0.99,
                                 chunk_size  = 5e6) {
  dt <- as.data.table(dt)
  all_ids <- dt$ID
  n_input <- length(all_ids)

  res <- rbindlist(lapply(recip_grid, function(thr) {
    edges <- compute_crm_edges_chunked(
      dt, thr_recip = thr, jaccard_thr = jaccard_thr,
      simpson_thr = simpson_thr, criterion = "composite",
      keep_metrics = FALSE, chunk_size = chunk_size
    )
    edf <- if (nrow(edges) > 0L)
      as.data.frame(edges[, .(id_i, id_j)], stringsAsFactors = FALSE) else
      data.frame(id_i = character(), id_j = character(), stringsAsFactors = FALSE)
    g <- graph_from_data_frame(edf, directed = FALSE,
                               vertices = data.frame(name = all_ids))
    comp <- igraph::components(g)
    sizes <- as.integer(comp$csize)
    data.table(
      thr_recip        = thr,
      n_input          = n_input,
      n_edges          = nrow(edges),
      n_clusters       = comp$no,
      n_merged         = sum(sizes > 1L),
      max_cluster_size = max(sizes),
      reduction_ratio  = 1 - comp$no / n_input
    )
  }))
  res[order(thr_recip)][]
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("crm_explore.R")
##
## # 1) Extraer una región de alta densidad para exploración:
## reg <- subset_region(enh_unique, "chr8", 125e6, 130e6)
##
## # 2) Distribución de métricas en la región:
## print(describe_crm_overlap(reg))
##
## # 3) Barrido de umbrales con el criterio compuesto definido para el análisis:
## sw <- sweep_crm_thresholds(reg, recip_grid = c(0.50,0.60,0.70,0.80,0.90))
## print(sw)
## # Vigilar max_cluster_size: si a 0.50 se dispara, hay chaining en zonas densas.
##
## # 4) Si un umbral parece bien, obtener las aristas CON métricas para inspección:
## edges <- compute_crm_edges_chunked(reg, thr_recip = 0.50, keep_metrics = TRUE)
## # (y reutilizar inspect_tad_clusters / cohesión, que son genéricas, sobre 'reg')
## ============================================================================
