################################################################################
##  figuras_tfm.R — Generación reproducible de figuras del TFM
##  ==========================================================================
##  Este script genera las figuras de la memoria a partir de los objetos RDS
##  producidos por el pipeline:
##
##    results/lineB/<chr>/pipeline_objects_<chr>.rds
##
##  y, cuando está disponible, del agregado global almacenado en:
##
##    results/global/
##
##  Cada figura se implementa como una función independiente fig_XX(...), que
##  devuelve un objeto ggplot o, en el caso de la visualización genómica con
##  Gviz, dibuja directamente sobre el dispositivo gráfico activo. La función
##  main() ejecuta el conjunto de figuras y las exporta al directorio figs_out/.
##
##  Reproducibilidad:
##    Las figuras se derivan de los objetos generados por el pipeline, evitando
##    depender de tablas externas preagregadas. El directorio de resultados puede
##    definirse mediante la variable de entorno TFM_RESULTS o modificando
##    RESULTS_DIR.
##
##  Uso en terminal:
##    Rscript figuras_tfm.R
##
##  Uso desde R:
##    source("figuras_tfm.R")
##    main()
##
##  Ejecución de una figura concreta:
##    source("figuras_tfm.R")
##    D <- load_all()
##    print(fig01_reduccion(D))
##
##  Dependencias:
##    CRAN: ggplot2, data.table, scales, patchwork, ggrepel
##    Bioconductor: Gviz, GenomicRanges
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

## ---------------------------------------------------------------------------
## 0. Configuración general
## ---------------------------------------------------------------------------
## Instalación de dependencias, si fuese necesario:
##
##   install.packages(c("ggplot2", "data.table", "scales", "patchwork", "ragg"))
##
##   if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
##
##   BiocManager::install(c("Gviz", "GenomicRanges"))
##
## Nota sobre codificación:
##   Si los caracteres acentuados no se renderizan correctamente, es recomendable
##   utilizar un locale UTF-8. La función main() intenta establecer C.UTF-8 de
##   forma automática.
## ---------------------------------------------------------------------------

RESULTS_DIR <- Sys.getenv("TFM_RESULTS", "results")
OUT_DIR     <- "figs_out"
dir.create(OUT_DIR, showWarnings = FALSE)

## Orden y etiquetas de las clases de DRR.
CLASS_ORDER <- c("Simple_DRR", "Compact_DRR",
                 "Extended_complex_DRR", "Dense_complex_DRR")
CLASS_LABEL <- c(Simple_DRR = "Simple", Compact_DRR = "Compact",
                 Extended_complex_DRR = "Extended",
                 Dense_complex_DRR = "Dense complex",
                 Extensive_overlap = "Extensive overlap")
CLASS_COLOR <- c(Simple_DRR = "#7BA7C7", Compact_DRR = "#6BAE8E",
                 Extended_complex_DRR = "#E0A458",
                 Dense_complex_DRR = "#C44E52",
                 Extensive_overlap = "#9B7FB5")

## Familia tipográfica con soporte para caracteres acentuados.
TFM_FONT <- "DejaVu Sans"

## Tema común para las figuras del TFM.
theme_tfm <- function(base = 12) {
  theme_minimal(base_size = base, base_family = TFM_FONT) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", size = base + 1),
          plot.subtitle = element_text(color = "grey35"),
          axis.title = element_text(color = "grey25"),
          legend.position = "right",
          plot.margin = margin(10, 14, 10, 10))
}

.msg <- function(...) cat(sprintf("[%s] %s\n",
                                  format(Sys.time(), "%H:%M:%S"), paste0(...)))

## Guardado de figuras. Se usa ragg si está disponible; en caso contrario,
## se emplea el dispositivo PNG de grDevices.
save_fig <- function(plot, file, width = 8, height = 5, dpi = 300) {
  dev_ragg <- requireNamespace("ragg", quietly = TRUE)
  if (dev_ragg) {
    ggsave(file, plot, width = width, height = height, dpi = dpi,
           device = ragg::agg_png)
  } else {
    ggsave(file, plot, width = width, height = height, dpi = dpi,
           device = grDevices::png, type = "cairo")
  }
  invisible(file)
}

## ---------------------------------------------------------------------------
## 1. Carga de datos desde los objetos RDS generados por el pipeline
## ---------------------------------------------------------------------------

#' Localiza y carga los objetos pipeline_objects_<chr>.rds.
#'
#' Devuelve una lista de tablas combinadas, preparadas para la generación de
#' figuras y resúmenes gráficos.
load_all <- function(results_dir = RESULTS_DIR) {
  rds <- list.files(file.path(results_dir, "lineB"),
                    pattern = "^pipeline_objects_.*\\.rds$",
                    recursive = TRUE, full.names = TRUE)
  if (!length(rds)) stop("No encuentro pipeline_objects_*.rds en ",
                         file.path(results_dir, "lineB"))
  chr_of <- sub(".*pipeline_objects_(.*)\\.rds$", "\\1", rds)
  ord <- order(suppressWarnings(as.integer(gsub("chr", "", chr_of))),
               chr_of)
  rds <- rds[ord]; chr_of <- chr_of[ord]
  .msg("Cargando ", length(rds), " cromosomas: ", paste(chr_of, collapse = ", "))
  
  red    <- list()
  tadcl  <- list()
  crmcl  <- list()
  per_tad_l <- list()
  drr    <- list()
  pos    <- list()
  gen    <- list()
  topo   <- list()
  modtab <- list()
  perdrr <- list()
  
  for (k in seq_along(rds)) {
    ch <- chr_of[k]; o <- readRDS(rds[k])
    
    red[[ch]] <- data.table(
      chr = ch,
      tad_in  = o$red_tad$summary$n_input,
      tad_out = o$red_tad$summary$n_clusters,
      crm_in  = o$red_crm$summary$n_crm_in,
      crm_assigned = o$red_crm$summary$n_crm_processed,
      crm_out = o$red_crm$summary$n_reduced)
    
    tadcl[[ch]] <- data.table(chr = ch, level = "TAD",
                              n = o$red_tad$tad_reduced$n_entities)
    crmcl[[ch]] <- data.table(chr = ch, level = "CRM",
                              n = o$red_crm$crm_reduced$n_entities)
    
    # CRMs por TAD reducido antes y después de la reducción.
    if (!is.null(o$red_crm$per_tad)) {
      pt <- as.data.table(o$red_crm$per_tad); pt[, chr := ch]
      per_tad_l[[ch]] <- pt
    }
    
    d <- as.data.table(o$cmp$stacking); d[, chr := ch]
    drr[[ch]] <- d
    
    if (!is.null(o$val$positional$by_class)) {
      p <- as.data.table(o$val$positional$by_class); p[, chr := ch]; pos[[ch]] <- p
    }
    if (!is.null(o$val$genic$summary)) {
      g <- as.data.table(o$val$genic$summary); g[, chr := ch]; gen[[ch]] <- g
    }
    if (!is.null(o$prop$by_class)) {
      tp <- as.data.table(o$prop$by_class); tp[, chr := ch]; topo[[ch]] <- tp
    }
    if (!is.null(o$mod$coef_table) && !isTRUE(o$mod$degenerate)) {
      m <- as.data.table(o$mod$coef_table); m[, chr := ch]; modtab[[ch]] <- m
    }
    if (!is.null(o$sm$per_drr)) {
      pd <- as.data.table(o$sm$per_drr); pd[, chr := ch]; perdrr[[ch]] <- pd
    }
  }
  
  ## Modelo global, si existe.
  gfile <- file.path(results_dir, "global", "model_coef_global.tsv")
  mod_global <- if (file.exists(gfile)) fread(gfile) else NULL
  
  list(
    chrs       = chr_of,
    reduction  = rbindlist(red),
    tad_clust  = rbindlist(tadcl),
    crm_clust  = rbindlist(crmcl),
    per_tad    = rbindlist(per_tad_l, use.names = TRUE, fill = TRUE),
    drr        = rbindlist(drr, use.names = TRUE, fill = TRUE),
    positional = rbindlist(pos, use.names = TRUE, fill = TRUE),
    genic      = rbindlist(gen, use.names = TRUE, fill = TRUE),
    topology   = rbindlist(topo, use.names = TRUE, fill = TRUE),
    model_chr  = rbindlist(modtab, use.names = TRUE, fill = TRUE),
    per_drr    = rbindlist(perdrr, use.names = TRUE, fill = TRUE),
    model_global = mod_global
  )
}

## Funciones auxiliares para ordenar clases y cromosomas.
.class_factor <- function(x, include_ext = FALSE) {
  lv <- if (include_ext) c(CLASS_ORDER, "Extensive_overlap") else CLASS_ORDER
  factor(x, levels = lv, labels = unname(CLASS_LABEL[lv]))
}
.chr_factor <- function(x) {
  u <- unique(x); n <- suppressWarnings(as.integer(gsub("chr", "", u)))
  factor(x, levels = u[order(n, u)])
}

## ===========================================================================
## Bloque 1. Reducción de redundancia de TADs y CRMs
## ===========================================================================

#' F1. Reducción de TADs y CRMs por cromosoma.
#'
#' Representa el porcentaje de reducción obtenido en las etapas de reducción
#' de TADs y CRMs.
fig01_reduccion <- function(D) {
  r <- copy(D$reduction)
  r[, `:=`(tad_pct = 100 * (1 - tad_out / tad_in),
           crm_pct = 100 * (1 - crm_out / crm_in))]
  m <- melt(r[, .(chr, TAD = tad_pct, CRM = crm_pct)],
            id.vars = "chr", variable.name = "nivel", value.name = "pct")
  m[, chr := .chr_factor(chr)]
  
  ggplot(m, aes(chr, pct, fill = nivel)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72) +
    scale_fill_manual(values = c(TAD = "#5A7A99", CRM = "#8896A6"), name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Reducción de redundancia por cromosoma",
         subtitle = "Porcentaje de entidades eliminadas tras colapso y reducción",
         x = NULL, y = "% de reducción") +
    theme_tfm() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
}

#' F1b. Conteos absolutos antes y después de la reducción.
#'
#' Muestra la magnitud absoluta de la reducción en TADs y CRMs.
fig01b_reduccion_abs <- function(D) {
  r <- copy(D$reduction)
  long <- rbindlist(list(
    r[, .(chr, nivel = "TAD", antes = tad_in, despues = tad_out)],
    r[, .(chr, nivel = "CRM", antes = crm_in, despues = crm_out)]))
  m <- melt(long, id.vars = c("chr", "nivel"),
            variable.name = "fase", value.name = "n")
  m[, chr := .chr_factor(chr)]
  m[, fase := factor(fase, levels = c("antes", "despues"),
                     labels = c("Crudos", "Reducidos"))]
  
  ggplot(m, aes(chr, n, fill = fase)) +
    geom_col(position = "dodge", width = 0.72) +
    facet_wrap(~ nivel, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c(Crudos = "#C9B79C", Reducidos = "#5A7A99"),
                      name = NULL) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       expand = expansion(mult = c(0, 0.06))) +
    labs(title = "Entidades antes y después de la reducción",
         x = NULL, y = "nº de entidades") +
    theme_tfm() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
}

#' F2. Distribución del tamaño de los clústeres de reducción.
#'
#' Resume cuántas entidades originales se agrupan por representante. La cola
#' extrema se recorta al percentil 99 para facilitar la visualización.
fig02_tam_cluster <- function(D) {
  cl <- rbindlist(list(D$tad_clust, D$crm_clust))
  cl[, level := factor(level, levels = c("TAD", "CRM"))]
  # Recorte de la cola extrema para mejorar la legibilidad.
  cl[, cap := quantile(n, 0.99), by = level]
  cl[, n_cap := pmin(n, cap)]
  
  ggplot(cl, aes(n_cap, fill = level)) +
    geom_histogram(bins = 40, alpha = 0.85, color = "white", linewidth = 0.2) +
    facet_wrap(~ level, scales = "free", ncol = 2) +
    scale_fill_manual(values = c(TAD = "#5A7A99", CRM = "#8896A6"),
                      guide = "none") +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    labs(title = "Tamaño de los clústeres de reducción",
         subtitle = "Nº de entidades originales colapsadas por representante (cola al p99)",
         x = "entidades por clúster", y = "frecuencia") +
    theme_tfm() +
    theme(panel.grid.major.x = element_line(color = "grey92"))
}

## ===========================================================================
## Bloque 2. Caracterización de regiones reguladoras densas
## ===========================================================================

#' F3. Número de DRRs por clase en el conjunto analizado.
fig03_conteos <- function(D) {
  d <- D$drr[candidate_class %in% CLASS_ORDER]
  tot <- d[, .N, by = candidate_class]
  tot[, cls := .class_factor(candidate_class)]
  
  ggplot(tot, aes(cls, N, fill = candidate_class)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.4) +
    geom_text(aes(label = label_number(big.mark = ".")(N)),
              vjust = -0.4, size = 3.6, fontface = "bold", color = "grey20") +
    scale_fill_manual(values = CLASS_COLOR, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                       labels = label_number(scale_cut = cut_short_scale())) +
    labs(title = sprintf("DRRs candidatas por clase (total: %s)",
                         label_number(big.mark = ".")(nrow(d))),
         x = NULL, y = "nº de DRRs") +
    theme_tfm()
}

#' F3b. Composición de clases de DRR por cromosoma.
fig03b_conteos_chr <- function(D) {
  d <- D$drr[candidate_class %in% CLASS_ORDER]
  cc <- d[, .N, by = .(chr, candidate_class)]
  cc[, chr := .chr_factor(chr)]
  cc[, cls := .class_factor(candidate_class)]
  
  ggplot(cc, aes(chr, N, fill = cls)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = setNames(CLASS_COLOR[CLASS_ORDER],
                                        CLASS_LABEL[CLASS_ORDER]), name = "Clase") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = label_number(scale_cut = cut_short_scale())) +
    labs(title = "Composición de clases de DRR por cromosoma",
         x = NULL, y = "nº de DRRs") +
    theme_tfm() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' F4. Distribución de tamaños de DRR por clase.
fig04_tamanos <- function(D) {
  d <- D$drr[candidate_class %in% CLASS_ORDER]
  d[, cls := .class_factor(candidate_class)]
  med <- d[, .(m = median(drr_length)), by = cls]
  
  ggplot(d, aes(cls, drr_length, fill = candidate_class)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 8700, ymax = 19000,
             fill = "grey70", alpha = 0.18) +
    geom_jitter(width = 0.28, alpha = 0.06, size = 0.4,
                color = "grey30", show.legend = FALSE) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85,
                 color = "grey20", show.legend = FALSE) +
    geom_text(data = med, aes(cls, m, label = paste0(label_number(big.mark=".")(round(m)), " bp")),
              inherit.aes = FALSE, hjust = -0.25, vjust = -0.6,
              size = 3, fontface = "bold", color = "grey15") +
    scale_fill_manual(values = CLASS_COLOR) +
    scale_y_log10(labels = label_number(big.mark = ".", scale_cut = cut_short_scale()),
                  breaks = c(100, 1000, 10000, 1e5)) +
    labs(title = "Distribución de tamaños por clase de DRR",
         subtitle = "Banda gris = mediana de SE en literatura (8,7–19 kb); escala logarítmica",
         x = NULL, y = "longitud de la DRR (bp)") +
    theme_tfm()
}

## ===========================================================================
## Bloque 3. Comparación con SEdb y análisis topológico
## ===========================================================================

#' F5. Recuperación posicional y génica de anotaciones de SEdb por clase.
fig05_recuperacion <- function(D) {
  pos <- D$positional[candidate_class %in% CLASS_ORDER,
                      .(n_total = sum(n_total), n_ovl = sum(n_with_ovl)),
                      by = candidate_class]
  pos[, `:=`(frac = n_ovl / n_total, tipo = "Posicional (solape SEdb)")]
  gen <- D$genic[candidate_class %in% CLASS_ORDER,
                 .(n_total = sum(n_candidates), n_ovl = sum(n_with_gene_match)),
                 by = candidate_class]
  gen[, `:=`(frac = n_ovl / n_total, tipo = "Génica (gen compartido)")]
  
  m <- rbind(pos, gen)
  m[, cls := .class_factor(candidate_class)]
  m[, tipo := factor(tipo, levels = c("Posicional (solape SEdb)",
                                      "Génica (gen compartido)"))]
  
  ggplot(m, aes(cls, frac, fill = candidate_class)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.4) +
    geom_text(aes(label = percent(frac, accuracy = 1)),
              vjust = -0.4, size = 3.3, fontface = "bold", color = "grey20") +
    facet_wrap(~ tipo) +
    scale_fill_manual(values = CLASS_COLOR, guide = "none") +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.13)),
                       limits = c(0, NA)) +
    labs(title = "Recuperación de super-enhancers anotados (SEdb)",
         subtitle = "Fracción de DRRs por clase con solape posicional / gen compartido",
         x = NULL, y = "% de DRRs") +
    theme_tfm() +
    theme(strip.text = element_text(face = "bold", size = 11),
          axis.text.x = element_text(angle = 20, hjust = 1))
}

#' F6. Forest plot del enriquecimiento de Dense_complex_DRR.
#'
#' Representa las razones de tasas de incidencia (IRR) por cromosoma y, cuando
#' está disponible, la estimación global.
fig06_forest <- function(D, ref_term = "candidate_classDense_complex_DRR") {
  mc <- D$model_chr[term == ref_term]
  mc[, `:=`(lo = exp(log(IRR) - 1.96 * std_error),
            hi = exp(log(IRR) + 1.96 * std_error))]
  setorder(mc, IRR)
  mc[, chr := factor(chr, levels = chr)]
  
  g_irr <- g_lo <- g_hi <- NA_real_
  if (!is.null(D$model_global)) {
    gr <- D$model_global[term == ref_term]
    if (nrow(gr)) {
      g_irr <- gr$IRR
      g_lo <- exp(log(gr$IRR) - 1.96 * gr$std_error)
      g_hi <- exp(log(gr$IRR) + 1.96 * gr$std_error)
    }
  }
  
  p <- ggplot(mc, aes(IRR, chr)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60")
  if (!is.na(g_irr)) {
    p <- p +
      annotate("rect", xmin = g_lo, xmax = g_hi, ymin = -Inf, ymax = Inf,
               fill = "#2C5F8A", alpha = 0.10) +
      geom_vline(xintercept = g_irr, color = "#2C5F8A", linewidth = 0.9)
  }
  p <- p +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.35,
                   color = "#C9908F", linewidth = 0.7) +
    geom_point(color = "#C44E52", fill = "white", shape = 21,
               size = 2.6, stroke = 1.1) +
    labs(title = "Enriquecimiento en super-enhancers por cromosoma",
         subtitle = if (!is.na(g_irr))
           sprintf("IRR de Dense_complex vs Extended (modelo negbin). Línea azul = global (IRR=%.1f)", g_irr)
         else "IRR de Dense_complex vs Extended (modelo binomial negativa)",
         x = "IRR (controlado por tamaño)", y = NULL) +
    theme_tfm() +
    theme(panel.grid.major.y = element_blank())
  p
}

#' F7. Coherencia topológica de las DRRs por clase.
fig07_topologia <- function(D) {
  tp <- D$topology[candidate_class %in% CLASS_ORDER]
  # Media ponderada por el número de DRRs entre cromosomas.
  tp[, single := pct_single_subtad / 100 * n_drr]
  ag <- tp[, .(pct = 100 * sum(single) / sum(n_drr)), by = candidate_class]
  ag[, cls := .class_factor(candidate_class)]
  
  ggplot(ag, aes(cls, pct, fill = candidate_class)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.1f%%", pct)),
              vjust = -0.4, size = 3.5, fontface = "bold", color = "grey20") +
    scale_fill_manual(values = CLASS_COLOR, guide = "none") +
    coord_cartesian(ylim = c(70, 102)) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(title = "Coherencia topológica de las DRRs",
         subtitle = "Porcentaje contenido en un único sub-TAD (media ponderada del genoma)",
         x = NULL, y = "% en sub-TAD único") +
    theme_tfm()
}

## ===========================================================================
## Bloque 4. Visualización genómica tipo IGV
## ===========================================================================
## Se incluyen dos implementaciones alternativas:
##
##   (A) igv_gviz()   — visualización basada en Bioconductor/Gviz.
##   (B) igv_ggplot() — visualización basada únicamente en ggplot2.
##
## Ambas funciones leen un objeto pipeline_objects_<chr>.rds y representan una
## región genómica concreta.
## ---------------------------------------------------------------------------

#' Carga las pistas genómicas de una región a partir del RDS de un cromosoma.
#'
#' Devuelve una lista con TADs, CRMs y DRRs dentro de la ventana seleccionada.
load_region <- function(chr, start, end, results_dir = RESULTS_DIR) {
  chr_sel <- chr
  rds <- file.path(results_dir, "lineB", chr,
                   paste0("pipeline_objects_", chr, ".rds"))
  if (!file.exists(rds)) stop("No existe ", rds)
  o <- readRDS(rds)
  
  tads <- as.data.table(o$red_tad$tad_reduced)[
    chr == chr_sel & repr_start < end & repr_end > start,
    .(id = cluster_id, start = repr_start, end = repr_end)]
  
  crms <- as.data.table(o$red_crm$crm_reduced)[
    repr_start < end & repr_end > start,
    .(start = repr_start, end = repr_end)]
  
  drrs <- as.data.table(o$cmp$stacking)[
    drr_start < end & drr_end > start]
  # Incorporar el número de solapamientos con SEdb por DRR.
  if (!is.null(o$sm$per_drr)) {
    se <- as.data.table(o$sm$per_drr)[, .(drr_id, n_SE)]
    drrs <- merge(drrs, se, by = "drr_id", all.x = TRUE)
  } else drrs[, n_SE := NA_integer_]
  drrs[is.na(n_SE), n_SE := 0L]
  
  list(chr = chr, start = start, end = end,
       tads = tads, crms = crms, drrs = drrs)
}

## ---- (A) Visualización con Gviz ------------------------------------------

#' Visualización genómica con Gviz.
#'
#' Requiere Gviz y GenomicRanges. Dibuja sobre el dispositivo gráfico activo o
#' exporta la figura si se proporciona un nombre de archivo.
igv_gviz <- function(region, file = NULL, width = 11, height = 7, dpi = 300) {
  stopifnot(requireNamespace("Gviz", quietly = TRUE),
            requireNamespace("GenomicRanges", quietly = TRUE))
  G <- asNamespace("Gviz"); GR <- GenomicRanges::GRanges
  IR <- IRanges::IRanges
  ch <- region$chr; s <- region$start; e <- region$end
  
  ## Eje de coordenadas.
  axis <- G$GenomeAxisTrack(littleTicks = TRUE, fontcolor = "#444444",
                            col = "#888888")
  
  ## Pista de sub-TADs.
  tad_gr <- GR(ch, IR(region$tads$start, region$tads$end),
               id = paste0("TAD ", region$tads$id))
  tad_tr <- G$AnnotationTrack(tad_gr, name = "sub-TADs",
                              fill = "#D9E2EC", col = "#5A7A99",
                              stacking = "squish", fontcolor.group = "#33506b")
  
  ## Pista de densidad de CRMs.
  if (nrow(region$crms)) {
    mids <- (region$crms$start + region$crms$end) / 2
    mids <- mids[mids >= s & mids <= e]
    bins <- seq(s, e, length.out = 160)
    h <- hist(mids, breaks = bins, plot = FALSE, include.lowest = TRUE)
    dens_gr <- GR(ch, IR(round(head(bins, -1)), round(bins[-1]) - 1))
    dens_tr <- G$DataTrack(dens_gr, data = h$counts, type = "histogram",
                           name = "CRMs/kb", fill.histogram = "#8896A6",
                           col.histogram = "#8896A6", col.axis = "#777777")
  } else dens_tr <- NULL
  
  ## Pista de DRRs coloreadas por clase.
  drr_gr <- GR(ch, IR(region$drrs$drr_start, region$drrs$drr_end))
  drr_tr <- G$AnnotationTrack(range = drr_gr, name = "DRRs",
                              stacking = "squish", col = "white",
                              feature = region$drrs$candidate_class,
                              id = region$drrs$drr_id)
  ## Asignación de colores por clase presente.
  present <- unique(region$drrs$candidate_class)
  for (nm in present)
    Gviz::displayPars(drr_tr)[[nm]] <- unname(CLASS_COLOR[nm])
  
  ## Pista de solapamientos con SEdb.
  se_gr <- GR(ch, IR(region$drrs$drr_start, region$drrs$drr_end))
  se_tr <- G$DataTrack(se_gr, data = region$drrs$n_SE, type = "histogram",
                       name = "SEdb (nº SE)", fill.histogram = "#C44E52",
                       col.histogram = "#C44E52", col.axis = "#777777")
  
  tracks <- Filter(Negate(is.null),
                   list(axis, tad_tr, dens_tr, drr_tr, se_tr))
  
  open_dev <- !is.null(file)
  if (open_dev) {
    if (requireNamespace("ragg", quietly = TRUE))
      ragg::agg_png(file, width = width, height = height, units = "in", res = dpi)
    else png(file, width = width, height = height, units = "in", res = dpi, type = "cairo")
  }
  Gviz::plotTracks(tracks, from = s, to = e, chromosome = ch,
                   sizes = c(0.8, 1, 1.2, 1.6, 1.2)[seq_along(tracks)],
                   main = sprintf("%s:%s-%s", ch, format(s, big.mark = "."),
                                  format(e, big.mark = ".")),
                   cex.main = 1, background.title = "#33506b")
  if (open_dev) { dev.off(); invisible(file) }
}

## ---- (B) Visualización con ggplot2 ---------------------------------------

#' Visualización genómica basada en ggplot2.
#'
#' Devuelve un objeto patchwork con cuatro pistas apiladas cuando patchwork está
#' disponible. En caso contrario, devuelve la pista de DRRs.
igv_ggplot <- function(region, show_title = TRUE) {
  s <- region$start; e <- region$end; ch <- region$chr
  base <- theme_tfm(11) +
    theme(axis.title.y = element_text(size = 9),
          plot.title = element_blank(), plot.subtitle = element_blank(),
          legend.position = "none",
          plot.margin = margin(2, 12, 2, 6))
  xlims <- coord_cartesian(xlim = c(s, e), expand = FALSE)
  xscale <- scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6))
  
  ## Pista de sub-TADs con apilamiento por carriles.
  td <- copy(region$tads); setorder(td, start)
  lane <- integer(nrow(td)); last_end <- numeric(0)
  for (i in seq_len(nrow(td))) {
    placed <- FALSE
    for (l in seq_along(last_end)) if (td$start[i] >= last_end[l]) {
      lane[i] <- l; last_end[l] <- td$end[i]; placed <- TRUE; break }
    if (!placed) { lane[i] <- length(last_end) + 1; last_end[lane[i]] <- td$end[i] }
  }
  td[, lane := lane]
  p_tad <- ggplot(td) +
    geom_rect(aes(xmin = pmax(start, s), xmax = pmin(end, e),
                  ymin = lane - 0.4, ymax = lane + 0.4),
              fill = "#D9E2EC", color = "#5A7A99", linewidth = 0.4) +
    xscale + xlims + labs(y = "sub-TADs", x = NULL) + base +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          panel.grid = element_blank())
  
  ## Pista de densidad de CRMs.
  mids <- (region$crms$start + region$crms$end) / 2
  mids <- mids[mids >= s & mids <= e]
  dens <- data.table(x = mids)
  p_crm <- ggplot(dens, aes(x)) +
    geom_histogram(bins = 160, fill = "#8896A6", color = NA) +
    xscale + coord_cartesian(xlim = c(s, e), expand = FALSE) +
    labs(y = "CRMs/kb", x = NULL) + base +
    theme(axis.text.x = element_blank())
  
  ## Pista de DRRs por clase con apilamiento por carriles.
  dr <- copy(region$drrs); setorder(dr, drr_start)
  lane <- integer(nrow(dr)); last_end <- numeric(0)
  for (i in seq_len(nrow(dr))) {
    placed <- FALSE
    for (l in seq_along(last_end)) if (dr$drr_start[i] >= last_end[l]) {
      lane[i] <- l; last_end[l] <- dr$drr_end[i]; placed <- TRUE; break }
    if (!placed) { lane[i] <- length(last_end) + 1; last_end[lane[i]] <- dr$drr_end[i] }
  }
  dr[, lane := lane]
  p_drr <- ggplot(dr) +
    geom_rect(aes(xmin = drr_start, xmax = drr_end,
                  ymin = lane - 0.42, ymax = lane + 0.42,
                  fill = candidate_class), color = "white", linewidth = 0.3) +
    scale_fill_manual(values = CLASS_COLOR) +
    xscale + coord_cartesian(xlim = c(s, e), expand = FALSE) +
    labs(y = "DRRs", x = NULL) + base +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          panel.grid = element_blank())
  
  ## Pista de solapamientos con SEdb.
  p_se <- ggplot(dr) +
    geom_rect(aes(xmin = drr_start, xmax = drr_end, ymin = 0, ymax = n_SE,
                  fill = candidate_class), color = NA) +
    scale_fill_manual(values = CLASS_COLOR) +
    xscale + coord_cartesian(xlim = c(s, e), expand = FALSE) +
    labs(y = "SEdb (nº SE)",
         x = sprintf("Posición en %s (Mb)", ch)) + base +
    theme(legend.position = "none")
  
  ttl <- sprintf("%s:%s-%s", ch, format(s, big.mark = "."), format(e, big.mark = "."))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    pw <- (p_tad / p_crm / p_drr / p_se) +
      plot_layout(heights = c(0.8, 1, 1.6, 1.1))
    if (show_title)
      pw <- pw + plot_annotation(
        title = paste0("Visualización tipo IGV — ", ttl),
        theme = theme(plot.title = element_text(
          face = "bold", size = 13, family = TFM_FONT)))
    pw
  } else {
    .msg("patchwork no disponible; se devuelve únicamente la pista de DRRs.")
    p_drr + labs(title = if (show_title) ttl else NULL)
  }
}

## ---- Leyenda de clases para visualizaciones genómicas ---------------------

#' Devuelve un objeto ggplot con la leyenda de clases de DRR.
igv_legend <- function() {
  d <- data.table(cls = .class_factor(CLASS_ORDER), x = seq_along(CLASS_ORDER))
  ggplot(d, aes(x, 1, fill = cls)) + geom_tile() +
    scale_fill_manual(values = setNames(CLASS_COLOR[CLASS_ORDER],
                                        CLASS_LABEL[CLASS_ORDER]), name = "Clase DRR") +
    theme_void(base_family = TFM_FONT) + theme(legend.position = "bottom")
}

#' Genera un panel de regiones genómicas seleccionadas.
#'
#' Cada región debe proporcionarse como una lista con chr, start, end y,
#' opcionalmente, title.
igv_panel <- function(regions, results_dir = RESULTS_DIR) {
  stopifnot(requireNamespace("patchwork", quietly = TRUE))
  library(patchwork)
  
  # Fila de título para cada región del panel.
  title_strip <- function(txt) {
    ggplot() +
      annotate("text", x = 0, y = 0.5, label = txt, hjust = 0, vjust = 0.5,
               size = 4.2, fontface = "bold", family = TFM_FONT, color = "grey15") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void()
  }
  
  # Para cada región se genera una fila de título y un bloque de pistas genómicas.
  blocks <- lapply(regions, function(rg) {
    reg <- load_region(rg$chr, rg$start, rg$end, results_dir)
    ttl <- if (!is.null(rg$title)) rg$title else
      sprintf("%s:%s-%s", rg$chr, format(rg$start, big.mark="."),
              format(rg$end, big.mark="."))
    body <- wrap_elements(full = igv_ggplot(reg, show_title = FALSE))
    wrap_plots(title_strip(ttl), body, ncol = 1, heights = c(0.06, 1))
  })
  
  wrap_plots(blocks, ncol = 1)
}

## Regiones destacadas por defecto para la visualización genómica.
DEFAULT_REGIONS <- list(
  list(chr = "chr1",  start = 149740000, end = 149960000,
       title = "chr1: DRR más densa del genoma (2.743 SE)"),
  list(chr = "chr10", start = 72210000,  end = 72420000,
       title = "chr10: segundo foco de densidad (2.326 SE)"),
  list(chr = "chr8",  start = 38900000,  end = 39200000,
       title = "chr8: región de referencia (672 SE)")
)

## ===========================================================================
## Ejecución principal: generación de figuras en figs_out/
## ===========================================================================

main <- function(results_dir = RESULTS_DIR, out_dir = OUT_DIR,
                 igv_engine = c("ggplot", "gviz", "both")) {
  igv_engine <- match.arg(igv_engine)
  ## Configuración UTF-8 para caracteres acentuados.
  try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE)
  dir.create(out_dir, showWarnings = FALSE)
  
  .msg("Cargando datos...")
  D <- load_all(results_dir)
  
  figs <- list(
    fig01_reduccion        = list(fig01_reduccion(D),        9, 5),
    fig01b_reduccion_abs   = list(fig01b_reduccion_abs(D),   9, 7),
    fig02_tam_cluster      = list(fig02_tam_cluster(D),      9, 4.5),
    fig03_conteos          = list(fig03_conteos(D),          7, 5),
    fig03b_conteos_chr     = list(fig03b_conteos_chr(D),     9, 5),
    fig04_tamanos          = list(fig04_tamanos(D),          8, 5),
    fig05_recuperacion     = list(fig05_recuperacion(D),     9, 5),
    fig06_forest           = list(fig06_forest(D),           7.5, 7.5),
    fig07_topologia        = list(fig07_topologia(D),        6.5, 4.5),
    fig09_flujo            = list(fig09_flujo(D),            10, 5.5),
    fig10_soporte          = list(fig10_soporte(D),          7, 5),
    fig12_crm_por_tad      = list(fig12_crm_por_tad(D),      8, 5),
    fig12b_crm_por_tad_chr = list(fig12b_crm_por_tad_chr(D), 9, 5)
  )
  for (nm in names(figs)) {
    f <- figs[[nm]]
    save_fig(f[[1]], file.path(out_dir, paste0(nm, ".png")),
             width = f[[2]], height = f[[3]])
    .msg("guardada ", nm)
  }
  
  ## Figura de longitud original frente a consenso. Si los ficheros originales
  ## no están disponibles, se genera una variante por nivel de soporte.
  enh_dir <- file.path(dirname(results_dir), "data", "enh_per_chr")
  if (!dir.exists(enh_dir))
    enh_dir <- file.path(results_dir, "..", "data", "enh_per_chr")
  save_fig(fig11_longitud(D, enh_dir = enh_dir, results_dir = results_dir),
           file.path(out_dir, "fig11_longitud.png"), 8, 5)
  .msg("guardada fig11_longitud")
  
  ## Visualización genómica con el motor seleccionado.
  if (igv_engine %in% c("ggplot", "both")) {
    reg <- load_region("chr8", 38900000, 39200000, results_dir)
    save_fig(igv_ggplot(reg), file.path(out_dir, "fig08_igv_ggplot.png"), 11, 7)
    save_fig(igv_panel(DEFAULT_REGIONS, results_dir),
             file.path(out_dir, "fig08_igv_panel.png"), 11, 14)
    save_fig(fig13_contraste(results_dir),
             file.path(out_dir, "fig13_contraste.png"), 11, 16)
    .msg("guardada IGV ggplot + panel + contraste")
  }
  if (igv_engine %in% c("gviz", "both")) {
    if (requireNamespace("Gviz", quietly = TRUE)) {
      reg <- load_region("chr8", 38900000, 39200000, results_dir)
      igv_gviz(reg, file.path(out_dir, "fig08_igv_gviz.png"), 11, 7)
      .msg("guardada IGV Gviz")
    } else .msg("Gviz no instalado; se omite la versión Gviz de la visualización genómica.")
  }
  
  .msg("===== FIGURAS COMPLETADAS en ", out_dir, " =====")
  invisible(D)
}

## Ejecutar main() cuando el archivo se lanza mediante Rscript.
if (sys.nframe() == 0L && !interactive()) {
  main(igv_engine = "both")
}

## ===========================================================================
## Figuras complementarias sobre reducción TAD/CRM y soporte de consensos
## ===========================================================================

#' F9. Esquema cuantitativo del pipeline.
#'
#' Resume los totales globales de reducción de TADs y CRMs a partir de los
#' campos summary almacenados en los objetos RDS.
fig09_flujo <- function(D, raw_dir = NULL) {
  # Totales globales agregados por cromosoma.
  red <- D$reduction
  tad_collapsed <- sum(red$tad_in)
  tad_reduced   <- sum(red$tad_out)
  crm_collapsed <- sum(red$crm_in)
  crm_reduced   <- sum(red$crm_out)
  
  # CRMs asignados a TAD reducido.
  crm_assigned <- if ("crm_assigned" %in% names(red)) sum(red$crm_assigned) else NA_integer_
  
  # Definición de cajas y posiciones para el diagrama de flujo.
  fmt <- function(x) format(x, big.mark = ".", scientific = FALSE)
  boxes <- data.table(
    row = c(rep("TADs", 2), rep("CRMs", if (is.na(crm_assigned)) 3 else 4)),
    x = c(1, 3,
          if (is.na(crm_assigned)) c(1, 2.5, 4) else c(1, 2.33, 3.66, 5)),
    label = c(
      paste0("TADs colapsados\npor ID\n", fmt(tad_collapsed)),
      paste0("TADs reducidos\n(no redundantes)\n", fmt(tad_reduced)),
      if (is.na(crm_assigned))
        c(paste0("CRMs colapsados\npor ID\n", fmt(crm_collapsed)),
          paste0("CRMs asignados\na sub-TAD\n(intra-TAD)"),
          paste0("CRMs consenso\n", fmt(crm_reduced)))
      else
        c(paste0("CRMs colapsados\npor ID\n", fmt(crm_collapsed)),
          paste0("CRMs asignados\na sub-TAD\n", fmt(crm_assigned)),
          paste0("(reducción\nintra-TAD)"),
          paste0("CRMs consenso\n", fmt(crm_reduced)))
    ))
  boxes[, y := ifelse(row == "TADs", 2, 1)]
  boxes[, fillc := ifelse(row == "TADs", "#5A7A99", "#8896A6")]
  
  # Flechas entre cajas consecutivas dentro de cada fila.
  arrows <- boxes[, .(x0 = head(x, -1), x1 = tail(x, -1),
                      y = y[1]), by = row]
  
  bw <- 0.78; bh <- 0.42
  ggplot() +
    geom_segment(data = arrows,
                 aes(x = x0 + bw/2, xend = x1 - bw/2, y = y, yend = y),
                 arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
                 color = "grey55", linewidth = 0.6) +
    geom_tile(data = boxes, aes(x, y, fill = fillc),
              width = bw, height = bh, color = "white") +
    geom_text(data = boxes, aes(x, y, label = label),
              color = "white", size = 3, fontface = "bold", lineheight = 0.95) +
    scale_fill_identity() +
    scale_y_continuous(limits = c(0.5, 2.5),
                       breaks = c(1, 2), labels = c("CRMs", "TADs")) +
    labs(title = "Resumen cuantitativo de la reducción de redundancia",
         subtitle = sprintf("Genoma completo: %s\u2192%s TADs y %s\u2192%s CRMs",
                            fmt(tad_collapsed), fmt(tad_reduced),
                            fmt(crm_collapsed), fmt(crm_reduced)),
         x = NULL, y = NULL) +
    theme_tfm() +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(face = "bold", size = 11))
}

#' F10. Distribución del soporte de los CRMs consenso.
#'
#' Utiliza el número de CRMs originales que contribuyen a cada consenso.
fig10_soporte <- function(D, kind = c("donut", "bar")) {
  
  kind <- match.arg(kind)
  
  cl <- copy(D$crm_clust)[level == "CRM"]
  
  cl[, bin := cut(
    n,
    breaks = c(0, 1, 5, 20, Inf),
    labels = c("1 (singleton)", "2–5", "6–20", ">20")
  )]
  
  ag <- cl[, .(n = .N), by = bin][order(bin)]
  ag[, pct := n / sum(n)]
  
  pal <- c(
    "1 (singleton)" = "#A6B8C8",
    "2–5" = "#7F9BB3",
    "6–20" = "#557592",
    ">20" = "#33506b"
  )
  
  if (kind == "bar") {
    
    ag[, grp := "CRMs consenso"]
    
    return(
      ggplot(ag, aes(grp, n, fill = bin)) +
        geom_col(
          position = "fill",
          width = 0.5,
          color = "white"
        ) +
        coord_flip() +
        scale_fill_manual(
          values = pal,
          name = "CRMs originales\npor consenso"
        ) +
        scale_y_continuous(
          labels = scales::percent
        ) +
        labs(
          title = "Soporte de los CRMs consenso",
          subtitle = "Fracción según nº de CRMs originales colapsados",
          x = NULL,
          y = NULL
        ) +
        theme_tfm() +
        theme(
          axis.text.y = element_blank(),
          panel.grid = element_blank()
        )
    )
  }
  
  # Representación tipo donut.
  ag[, frac := n / sum(n)]
  ag[, ymax := cumsum(frac)]
  ag[, ymin := ymax - frac]
  ag[, ymid := (ymin + ymax) / 2]
  ag[, lab := paste0(scales::percent(frac, accuracy = 0.1))]
  ag[, small := frac < 0.03]
  
  ggplot(
    ag,
    aes(
      ymin = ymin,
      ymax = ymax,
      xmin = 2.5,
      xmax = 4.5,
      fill = bin
    )
  ) +
    geom_rect(
      color = "white",
      linewidth = 0.6
    ) +
    coord_polar(theta = "y") +
    xlim(c(1.5, 4.6)) +
    scale_fill_manual(
      values = pal,
      name = "CRMs originales\npor consenso"
    ) +
    
    # Etiquetas interiores para las secciones principales.
    geom_text(
      data = ag[small == FALSE],
      aes(
        x = 3.5,
        y = ymid,
        label = lab
      ),
      size = 3.35,
      color = "white",
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    
    # Línea guía para secciones pequeñas.
    geom_segment(
      data = ag[small == TRUE],
      aes(
        x = 4.0,
        xend = 4.18,
        y = ymid,
        yend = ymid
      ),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = "black"
    ) +
    
    # Etiqueta exterior para secciones pequeñas.
    geom_label(
      data = ag[small == TRUE],
      aes(
        x = 4.30,
        y = ymid,
        label = lab
      ),
      inherit.aes = FALSE,
      size = 3.5,
      fontface = "bold",
      fill = "white",
      color = "black",
      label.size = 0.25
    ) +
    
    labs(
      title = "Soporte de los CRMs consenso",
      subtitle = sprintf(
        "%s CRMs consenso; nº de CRMs originales colapsados",
        format(sum(ag$n), big.mark = ".")
      )
    ) +
    theme_void(base_family = TFM_FONT) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey35"),
      legend.position = "right"
    )
}

#' F11. Distribución de longitud: CRMs originales frente a CRMs consenso.
#'
#' Las longitudes originales se leen de los ficheros de entrada indicados en
#' enh_dir. Las longitudes de consenso se recuperan de los objetos RDS. Si los
#' ficheros originales no están disponibles, se genera una variante por nivel
#' de soporte del consenso.
fig11_longitud <- function(D, enh_dir = file.path(RESULTS_DIR, "..", "data", "enh_per_chr"),
                           chrs = NULL, sample_n = NULL, results_dir = RESULTS_DIR) {
  # Longitud de los CRMs consenso recuperada desde los objetos RDS.
  rds <- list.files(file.path(results_dir, "lineB"),
                    pattern = "^pipeline_objects_.*\\.rds$",
                    recursive = TRUE, full.names = TRUE)
  if (!is.null(chrs)) {
    keep <- sub(".*pipeline_objects_(.*)\\.rds$", "\\1", rds) %in% chrs
    rds <- rds[keep]
  }
  cons <- rbindlist(lapply(rds, function(f) {
    o <- readRDS(f); cr <- as.data.table(o$red_crm$crm_reduced)
    d <- data.table(len = cr$repr_end - cr$repr_start + 1L)
    if (nrow(d) > 40000L) d <- d[sample(.N, 40000L)]
    d
  }))
  cons[, tipo := "CRMs consenso"]
  
  # Longitud de CRMs originales tras colapso por identificador.
  orig <- NULL
  if (!is.null(enh_dir) && dir.exists(enh_dir)) {
    chs <- if (!is.null(chrs)) chrs else D$chrs
    parts <- lapply(chs, function(ch) {
      fp <- file.path(enh_dir, paste0(ch, ".tsv.gz"))
      if (!file.exists(fp)) return(NULL)
      raw <- fread(fp)
      # Columnas esperadas: chr, start, end, ID.
      idc <- intersect(c("ID", "id", "name"), names(raw))[1]
      if (is.na(idc)) idc <- names(raw)[4]
      setnames(raw, idc, "ID")
      raw[, start := as.integer(start)][, end := as.integer(end)]
      raw <- raw[!is.na(start) & !is.na(end) & end >= start]
      col <- raw[, .(start = min(start), end = max(end)), by = ID]
      d <- data.table(len = col$end - col$start + 1L)
      if (nrow(d) > 40000L) d <- d[sample(.N, 40000L)]
      d
    })
    orig <- rbindlist(Filter(Negate(is.null), parts))
    if (nrow(orig)) orig[, tipo := "CRMs originales\n(colapsados por ID)"]
  }
  
  if (is.null(orig) || !nrow(orig)) {
    .msg("No se encontraron CRMs originales en ", enh_dir,
         "; se genera la variante por nivel de soporte.")
    cl <- copy(D$crm_clust)[level == "CRM"]
    # Recuperar longitud y soporte de CRMs consenso desde los objetos RDS.
    cons2 <- rbindlist(lapply(rds, function(f) {
      o <- readRDS(f); cr <- as.data.table(o$red_crm$crm_reduced)
      d <- data.table(len = cr$repr_end - cr$repr_start + 1L, n = cr$n_entities)
      # Submuestreo por cromosoma para mantener tiempos y memoria acotados.
      if (nrow(d) > 40000L) d <- d[sample(.N, 40000L)]
      d
    }))
    if (!nrow(cons2)) stop("No pude releer CRMs consenso de los .rds.")
    cons2[, sop := cut(n, c(0,1,5,20,Inf),
                       labels = c("1","2\u20135","6\u201320",">20"))]
    return(
      ggplot(cons2, aes(sop, len, fill = sop)) +
        geom_violin(color = "grey30", alpha = 0.8, scale = "width") +
        geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white") +
        scale_y_log10(labels = label_number(big.mark=".", scale_cut=cut_short_scale())) +
        scale_fill_manual(values = c("#C9D3DD","#8FA9BE","#5A7A99","#33506b"),
                          guide = "none") +
        labs(title = "Longitud de los CRMs consenso por nivel de soporte",
             subtitle = "A mayor soporte, ¿mayor tamaño? (escala log)",
             x = "CRMs originales colapsados", y = "longitud (bp)") +
        theme_tfm()
    )
  }
  
  # Submuestreo opcional para agilizar la representación de densidad.
  both <- rbind(orig, cons)
  if (!is.null(sample_n) && nrow(both) > sample_n)
    both <- both[sample(.N, sample_n)]
  
  ggplot(both, aes(len, fill = tipo, color = tipo)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    scale_x_log10(labels = label_number(big.mark=".", scale_cut=cut_short_scale()),
                  breaks = c(100,1000,10000,1e5,1e6)) +
    scale_fill_manual(values = c("#C9B79C","#5A7A99"), name = NULL) +
    scale_color_manual(values = c("#A8906F","#3d5876"), name = NULL) +
    labs(title = "Longitud de los CRMs: originales vs consenso",
         subtitle = "La reducción mantiene tamaños compatibles con regiones reguladoras (escala log)",
         x = "longitud (bp)", y = "densidad") +
    theme_tfm() + theme(legend.position = "top")
}

#' F12a. Distribución global del número de CRMs consenso por TAD reducido.
fig12_crm_por_tad <- function(D) {
  pt <- D$per_tad[!is.na(n_reduced)]
  med <- median(pt$n_reduced)
  ggplot(pt, aes(n_reduced)) +
    geom_histogram(bins = 50, fill = "#5A7A99", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = med, linetype = "dashed", color = "#C44E52",
               linewidth = 0.7) +
    annotate("text", x = med, y = Inf, label = paste0("  mediana = ", med),
             hjust = 0, vjust = 1.6, color = "#C44E52", size = 3.3) +
    scale_x_log10(labels = label_number(big.mark = ".")) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    labs(title = "Carga reguladora por sub-TAD",
         subtitle = "Nº de CRMs consenso asignados a cada sub-TAD reducido (escala log)",
         x = "CRMs consenso por sub-TAD", y = "nº de sub-TADs") +
    theme_tfm() + theme(panel.grid.major.x = element_line(color = "grey92"))
}

#' F12b. Distribución de CRMs por TAD reducido antes y después de la reducción.
fig12b_crm_por_tad_chr <- function(D) {
  pt <- copy(D$per_tad)
  m <- melt(pt[, .(chr, tad_id, `originales` = n_in, `consenso` = n_reduced)],
            id.vars = c("chr","tad_id"), variable.name = "fase",
            value.name = "n")
  m[, chr := .chr_factor(chr)]
  ggplot(m, aes(chr, n, fill = fase)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.85) +
    scale_fill_manual(values = c(originales = "#C9B79C", consenso = "#5A7A99"),
                      name = "CRMs") +
    scale_y_log10(labels = label_number(big.mark = ".")) +
    coord_cartesian(ylim = c(1, NA)) +
    labs(title = "CRMs por sub-TAD, antes y después de la reducción",
         subtitle = "Por cromosoma (escala log); la reducción descongestiona los sub-TADs",
         x = NULL, y = "CRMs por sub-TAD") +
    theme_tfm() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
}

## ---- F13. Panel de regiones de contraste ---------------------------------

#' Regiones seleccionadas para ilustrar distintos patrones de comportamiento.
CONTRAST_REGIONS <- list(
  list(chr = "chr1",  start = 149770000, end = 149930000,
       title = "(1) Dense complex — recuperación máxima (441 CRMs, 2.743 SE)"),
  list(chr = "chr1",  start = 32040000,  end = 32108000,
       title = "(2) Compact pequeña pero con soporte alto (73 CRMs, 50 SE)"),
  list(chr = "chr21", start = 11840000,  end = 11942000,
       title = "(3) Extended de baja densidad (4 CRMs, 0 SE)"),
  list(chr = "chr14", start = 94010000,  end = 94122000,
       title = "(4) Densa pero sin recuperación SEdb (302 CRMs, 0 SE)")
)

#' Genera un panel de contraste usando las regiones predefinidas.
fig13_contraste <- function(results_dir = RESULTS_DIR) {
  igv_panel(CONTRAST_REGIONS, results_dir)
}
