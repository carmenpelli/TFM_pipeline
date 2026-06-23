################################################################################
##  exploracion_parametros.R — JUSTIFICACIÓN ESTADÍSTICA DE LOS PARÁMETROS
##  ==========================================================================
##  Reproduce y formaliza la exploración que llevó a elegir los umbrales del
##  criterio compuesto de reducción de CRMs:
##      arista si  reciprocal_overlap >= 0.50
##                 Y (jaccard >= 0.70  O  simpson >= 0.99)
##
##  Produce tres bloques de evidencia, cada uno con su test:
##    (A) Distribución de las métricas de solapamiento (recíproco, Jaccard,
##        Simpson, Dice) sobre los pares solapados reales -> justifica los
##        umbrales a partir de la forma de las distribuciones.
##    (B) Redundancia entre métricas: por qué se DESCARTA Dice (es una función
##        monótona de Jaccard, Dice = 2J/(1+J)) y qué aporta cada métrica
##        retenida. Tests: correlación de Spearman + verificación de la relación
##        analítica Dice–Jaccard.
##    (C) Análisis de sensibilidad: barrido de umbrales y efecto sobre la
##        reducción (nº de clústeres, tamaño máximo de clúster = diagnóstico de
##        chaining). Justifica la elección final.
##
##  Reproducible: trabaja sobre una REGIÓN contigua de CRMs (preserva la
##  estructura de solapamiento). Por defecto usa una ventana de chr8, pero
##  acepta cualquier cromosoma/ventana o el cromosoma completo.
##
##  DEPENDENCIAS: data.table, GenomicRanges, IRanges, ggplot2, patchwork
##  (Reutiliza funciones de R/crm_explore.R si está disponible en el path.)
##
##  USO:
##    source("exploracion_parametros.R")
##    # cargar CRMs crudos de un cromosoma (colapsados por ID, como en P1):
##    crm <- cargar_crms_region("data/enh_per_chr/chr8.tsv.gz",
##                              chr="chr8", from=38e6, to=42e6)
##    expl <- justificar_parametros(crm, out_dir="figs_param")
##    print(expl$tabla_distribucion)
##    print(expl$tabla_correlacion)
##    print(expl$tabla_sensibilidad)
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
})

if (!exists(".msg")) .msg <- function(...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

## Parámetros finales elegidos (para marcarlos en las figuras)
P_RECIP   <- 0.50
P_JACCARD <- 0.70
P_SIMPSON <- 0.99

## ---------------------------------------------------------------------------
## 0. Carga de CRMs crudos de una ventana (colapso por ID, como en P1)
## ---------------------------------------------------------------------------
#' Lee un fichero crudo de CRMs y devuelve una ventana contigua colapsada por ID.
#' Preservar una ventana contigua (no muestreo aleatorio) mantiene intacta la
#' estructura de solapamiento, condición necesaria para describir las métricas.
cargar_crms_region <- function(path, chr, from = NULL, to = NULL) {
  chr_sel <- chr   # evita colisión con la columna 'chr' en data.table
  raw <- fread(path, sep = "\t", header = TRUE, quote = "", fill = TRUE)
  setnames(raw, gsub('"', '', names(raw)))
  idc <- intersect(c("ID", "id", "name"), names(raw))[1]
  if (is.na(idc)) idc <- names(raw)[4]
  setnames(raw, idc, "ID")
  cc <- intersect(c("chr", "chrom", "seqnames"), names(raw))[1]
  if (!is.na(cc) && cc != "chr") setnames(raw, cc, "chr")
  raw[, start := as.integer(start)][, end := as.integer(end)]
  raw <- raw[!is.na(start) & !is.na(end) & end >= start]
  # colapso por ID (P1): min(start), max(end)
  col <- raw[, .(chr = chr[1L], start = min(start), end = max(end)), by = ID]
  col <- col[chr == chr_sel]
  if (!is.null(from)) col <- col[end >= from]
  if (!is.null(to))   col <- col[start <= to]
  .msg("Región ", chr, if (!is.null(from)) paste0(":", from, "-", to) else "",
       ": ", nrow(col), " CRMs colapsados.")
  col[]
}

## ---------------------------------------------------------------------------
## 1. Métricas de TODOS los pares solapados de una región (con Dice incluido)
## ---------------------------------------------------------------------------
#' Calcula recíproco, Jaccard, Simpson y Dice de cada par solapado.
#' @return data.table con una fila por par solapado y sus 4 métricas.
metricas_pares <- function(dt, n_sample = 2e6, seed = 1) {
  dt <- as.data.table(dt)
  gr <- GRanges(dt$chr, IRanges(dt$start, dt$end))
  st <- start(gr); en <- end(gr); len <- width(gr)
  hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
  qh <- queryHits(hits); sh <- subjectHits(hits)
  keep <- qh < sh; qh <- qh[keep]; sh <- sh[keep]
  n_pairs <- length(qh)
  if (n_pairs == 0L) { .msg("Sin pares solapados."); return(NULL) }
  if (n_pairs > n_sample) {
    set.seed(seed); sel <- sample.int(n_pairs, n_sample)
    qh <- qh[sel]; sh <- sh[sel]
    .msg("Muestreando ", n_sample, " de ", n_pairs, " pares solapados.")
  }
  ov <- pmin(en[qh], en[sh]) - pmax(st[qh], st[sh]) + 1L
  ov[ov < 0L] <- 0L
  Li <- len[qh]; Lj <- len[sh]
  data.table(
    reciprocal = pmin(ov / Li, ov / Lj),
    jaccard    = ov / (Li + Lj - ov),
    simpson    = ov / pmin(Li, Lj),
    dice       = 2 * ov / (Li + Lj)          # Dice–Sørensen
  )[]
}

## ---------------------------------------------------------------------------
## 2. (A) Distribución de métricas: cuantiles + figura
## ---------------------------------------------------------------------------
tabla_cuantiles <- function(M, probs = c(0,.1,.25,.5,.75,.9,.95,.99,1)) {
  data.table(
    cuantil    = paste0(probs * 100, "%"),
    reciproco  = round(as.numeric(quantile(M$reciprocal, probs)), 3),
    jaccard    = round(as.numeric(quantile(M$jaccard,    probs)), 3),
    simpson    = round(as.numeric(quantile(M$simpson,    probs)), 3),
    dice       = round(as.numeric(quantile(M$dice,       probs)), 3))
}

fig_distribucion_metricas <- function(M) {
  long <- melt(M, measure.vars = c("reciprocal","jaccard","simpson","dice"),
               variable.name = "metrica", value.name = "valor")
  long[, metrica := factor(metrica,
                           levels = c("reciprocal","jaccard","simpson","dice"),
                           labels = c("Recíproco","Jaccard","Simpson","Dice"))]
  # líneas de umbral elegido
  vlíneas <- data.table(
    metrica = factor(c("Recíproco","Jaccard","Simpson"),
                     levels = levels(long$metrica)),
    umbral = c(P_RECIP, P_JACCARD, P_SIMPSON))
  ggplot(long, aes(valor)) +
    geom_histogram(aes(fill = metrica), bins = 50, color = "white",
                   linewidth = 0.1, show.legend = FALSE) +
    geom_vline(data = vlíneas, aes(xintercept = umbral),
               linetype = "dashed", color = "#C44E52", linewidth = 0.7) +
    geom_text(data = vlíneas, aes(x = umbral, y = Inf, label = umbral),
              vjust = 1.5, hjust = -0.15, color = "#C44E52", size = 3) +
    facet_wrap(~ metrica, scales = "free", ncol = 2) +
    scale_fill_manual(values = c("#5A7A99","#6BAE8E","#E0A458","#9B7FB5")) +
    labs(title = "Distribución de las métricas de solapamiento entre CRMs",
         subtitle = "Pares solapados de la región; línea roja = umbral elegido",
         x = "valor de la métrica", y = "nº de pares") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))
}

## ---------------------------------------------------------------------------
## 3. (B) Redundancia entre métricas: por qué se descarta Dice
## ---------------------------------------------------------------------------
#' Correlaciones de Spearman entre las 4 métricas + verificación de que
#' Dice = 2J/(1+J) (relación analítica que hace a Dice redundante con Jaccard).
analisis_redundancia <- function(M) {
  cols <- c("reciprocal","jaccard","simpson","dice")
  cmat <- cor(M[, ..cols], method = "spearman")
  
  # Verificación de la relación analítica Dice–Jaccard
  dice_pred <- 2 * M$jaccard / (1 + M$jaccard)
  err <- max(abs(M$dice - dice_pred))
  rho_dice_jac <- cor(M$dice, M$jaccard, method = "spearman")
  
  tabla <- as.data.table(round(cmat, 4), keep.rownames = "metrica")
  list(
    correlaciones = tabla,
    dice_jaccard_rho = rho_dice_jac,
    dice_jaccard_maxerr = err,   # ~0 confirma Dice = 2J/(1+J)
    nota = sprintf(
      "Dice y Jaccard: rho de Spearman = %.4f; error máx. frente a 2J/(1+J) = %.2e. Dice es una transformación monótona de Jaccard, luego es redundante y se descarta.",
      rho_dice_jac, err))
}

fig_redundancia <- function(M, n_pts = 30000, seed = 1) {
  set.seed(seed)
  S <- if (nrow(M) > n_pts) M[sample(.N, n_pts)] else M
  # curva analítica Dice = 2J/(1+J)
  jx <- seq(0, 1, length.out = 200); dy <- 2 * jx / (1 + jx)
  curva <- data.table(jaccard = jx, dice = dy)
  p1 <- ggplot(S, aes(jaccard, dice)) +
    geom_point(alpha = 0.08, size = 0.4, color = "#5A7A99") +
    geom_line(data = curva, color = "#C44E52", linewidth = 0.9) +
    annotate("text", x = 0.25, y = 0.9, label = "Dice = 2J/(1+J)",
             color = "#C44E52", size = 3.4, hjust = 0) +
    labs(title = "Dice es redundante con Jaccard",
         subtitle = "Los puntos (pares reales) caen sobre la curva analítica",
         x = "Jaccard", y = "Dice") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  # Jaccard vs Simpson: NO redundantes (capturan cosas distintas)
  p2 <- ggplot(S, aes(jaccard, simpson)) +
    geom_point(alpha = 0.08, size = 0.4, color = "#6BAE8E") +
    labs(title = "Jaccard y Simpson son complementarias",
         subtitle = "Simpson alto con Jaccard bajo = inclusión (un CRM dentro de otro)",
         x = "Jaccard", y = "Simpson") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork); p1 + p2
  } else list(dice = p1, simpson = p2)
}

## ---------------------------------------------------------------------------
## 4. (C) Análisis de sensibilidad: barrido de umbrales
## ---------------------------------------------------------------------------
#' Para cada umbral recíproco (con el criterio compuesto), cuenta clústeres y
#' el tamaño máximo de clúster (diagnóstico de chaining: si explota, el umbral
#' es demasiado laxo). Requiere igraph para componentes conexas.
#'
#' Versión autocontenida (no depende de crm_explore.R): calcula aristas con
#' findOverlaps y el criterio compuesto, y agrupa con igraph.
sensibilidad_umbrales <- function(dt,
                                  recip_grid   = c(0.30,0.40,0.50,0.60,0.70,0.80,0.90),
                                  jaccard_grid = P_JACCARD,
                                  simpson_thr  = P_SIMPSON) {
  stopifnot(requireNamespace("igraph", quietly = TRUE))
  dt <- as.data.table(dt)
  gr <- GRanges(dt$chr, IRanges(dt$start, dt$end))
  st <- start(gr); en <- end(gr); len <- width(gr); ids <- dt$ID; n <- length(ids)
  hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
  qh <- queryHits(hits); sh <- subjectHits(hits)
  keep <- qh < sh; qh <- qh[keep]; sh <- sh[keep]
  ov <- pmin(en[qh], en[sh]) - pmax(st[qh], st[sh]) + 1L; ov[ov < 0L] <- 0L
  Li <- len[qh]; Lj <- len[sh]
  recip <- pmin(ov/Li, ov/Lj); jac <- ov/(Li+Lj-ov); sim <- ov/pmin(Li,Lj)
  
  grid <- CJ(recip_thr = recip_grid, jaccard_thr = jaccard_grid)
  res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
    rt <- grid$recip_thr[i]; jt <- grid$jaccard_thr[i]
    pass <- recip >= rt & (jac >= jt | sim >= simpson_thr)
    ei <- qh[pass]; ej <- sh[pass]
    g <- igraph::graph_from_edgelist(cbind(ei, ej), directed = FALSE)
    g <- igraph::add_vertices(g, max(0, n - igraph::vcount(g)))
    comp <- igraph::components(g)
    data.table(recip_thr = rt, jaccard_thr = jt,
               n_input = n, n_edges = sum(pass),
               n_clusters = comp$no,
               reduction_pct = round(100 * (1 - comp$no / n), 2),
               max_cluster = max(comp$csize))
  }))
  res[]
}

fig_sensibilidad <- function(S) {
  m <- melt(S, id.vars = c("recip_thr","jaccard_thr"),
            measure.vars = c("reduction_pct","max_cluster"),
            variable.name = "indicador", value.name = "valor")
  m[, indicador := factor(indicador, levels = c("reduction_pct","max_cluster"),
                          labels = c("% de reducción", "tamaño máx. de clúster (chaining)"))]
  ggplot(m, aes(recip_thr, valor)) +
    geom_line(color = "#5A7A99", linewidth = 0.8) +
    geom_point(color = "#5A7A99", size = 2) +
    geom_vline(xintercept = P_RECIP, linetype = "dashed", color = "#C44E52") +
    annotate("text", x = P_RECIP, y = Inf, label = paste0("  elegido: ", P_RECIP),
             color = "#C44E52", hjust = 0, vjust = 1.5, size = 3) +
    facet_wrap(~ indicador, scales = "free_y", ncol = 1) +
    labs(title = "Análisis de sensibilidad del umbral recíproco",
         subtitle = sprintf("Criterio compuesto con Jaccard=%.2f, Simpson=%.2f",
                            P_JACCARD, P_SIMPSON),
         x = "umbral de solapamiento recíproco", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

## ---------------------------------------------------------------------------
## 5. Runner: ejecuta los tres bloques y guarda figuras + tablas
## ---------------------------------------------------------------------------
justificar_parametros <- function(crm, out_dir = "figs_param",
                                  n_sample = 2e6, recip_grid = NULL) {
  try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE)
  dir.create(out_dir, showWarnings = FALSE)
  save1 <- function(p, f, w, h) ggsave(file.path(out_dir, f), p, width = w,
                                       height = h, dpi = 300,
                                       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL)
  
  .msg("(A) Métricas de los pares solapados...")
  M <- metricas_pares(crm, n_sample = n_sample)
  q <- tabla_cuantiles(M)
  save1(fig_distribucion_metricas(M), "param_A_distribucion.png", 9, 6)
  
  .msg("(B) Redundancia entre métricas (Dice vs Jaccard)...")
  red <- analisis_redundancia(M)
  save1(fig_redundancia(M), "param_B_redundancia.png", 11, 4.5)
  .msg(red$nota)
  
  .msg("(C) Análisis de sensibilidad...")
  S <- sensibilidad_umbrales(crm,
                             recip_grid = if (is.null(recip_grid))
                               c(0.30,0.40,0.50,0.60,0.70,0.80,0.90) else recip_grid)
  save1(fig_sensibilidad(S), "param_C_sensibilidad.png", 7, 6)
  
  .msg("Hecho. Figuras y tablas en ", out_dir, "/")
  list(tabla_distribucion = q,
       tabla_correlacion = red$correlaciones,
       dice_redundancia = red,
       tabla_sensibilidad = S,
       metricas = M)
}