################################################################################
##  drr_lineB_stats_figs.R — Análisis estadístico y generación de figuras
##  --------------------------------------------------------------------------
##  Este script contiene funciones auxiliares para evaluar la concordancia entre
##  regiones reguladoras densas (DRRs) y anotaciones de super-enhancers de SEdb.
##
##  Incluye:
##    - Ajuste de un modelo de conteo para evaluar si la clase estructural de
##      DRR se asocia con el número de solapamientos con SEdb, controlando por
##      la longitud de la región.
##    - Generación de figuras resumen basadas en ggplot2 para visualizar la
##      relación entre clase de DRR, longitud y concordancia con SEdb.
##
##  Entrada esperada:
##    sm$per_drr, generado por la etapa de comparación con SEdb, con columnas:
##    drr_id, candidate_class, drr_length, n_SE y size_bin.
##
##  Dependencias:
##    data.table, ggplot2; MASS es opcional para el ajuste binomial negativo.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

## ============================================================================
## (1) MODELO ESTADÍSTICO: clase + longitud -> nº de solapes SEdb
## ============================================================================
#' Ajusta un modelo de conteo de solapamientos con SEdb.
#'
#' @param per_drr Tabla por DRR con drr_id, candidate_class, drr_length y n_SE.
#' @param ref_class Clase de referencia. Por defecto se utiliza
#'        "Extended_complex_DRR", correspondiente a regiones complejas extensas
#'        con menor densidad relativa.
#' @return Lista con el modelo ajustado, la tabla de coeficientes y la familia usada.
model_density_controlled <- function(per_drr,
                                     ref_class = "Extended_complex_DRR") {
  d <- as.data.table(copy(per_drr))
  d <- d[!is.na(n_SE) & !is.na(drr_length)]
  d[, candidate_class := relevel(factor(candidate_class), ref = ref_class)]
  # log-longitud (en kb) como covariable de tamaño
  d[, log_len_kb := log10(drr_length / 1000)]

  # --- Diagnóstico previo al ajuste -----------------------------------------
  # El modelo requiere variación suficiente en n_SE. Cuando el número de
  # solapamientos positivos es muy bajo, el ajuste puede resultar inestable,
  # especialmente en cromosomas con pocas DRRs informativas.
  n_obs       <- nrow(d)
  n_nonzero   <- sum(d$n_SE > 0, na.rm = TRUE)
  n_classes   <- d[, uniqueN(candidate_class)]
  signal_ok   <- n_nonzero >= 10L && n_classes >= 2L && var(d$n_SE) > 0

  if (!signal_ok) {
    .msg("AVISO: señal insuficiente para el modelo (n_obs=", n_obs,
         ", n_SE>0=", n_nonzero, ", clases=", n_classes,
         "). Modelo marcado como DEGENERADO (típico en cromosomas pequeños).")
    return(list(model = NULL, coef_table = NULL, type = "degenerate",
                degenerate = TRUE,
                diagnostics = list(n_obs = n_obs, n_nonzero = n_nonzero,
                                   n_classes = n_classes),
                data = d[]))
  }

  # Ajuste binomial negativo para conteos sobredispersos; alternativa quasipoisson.
  type <- "negbin"; converged <- TRUE
  model <- tryCatch({
    if (!requireNamespace("MASS", quietly = TRUE)) stop("no MASS")
    withCallingHandlers(
      MASS::glm.nb(n_SE ~ candidate_class + log_len_kb, data = d),
      warning = function(w) {
        if (grepl("did not converge|iteration limit", conditionMessage(w)))
          converged <<- FALSE
        invokeRestart("muffleWarning")
      })
  }, error = function(e) {
    type <<- "quasipoisson"
    withCallingHandlers(
      glm(n_SE ~ candidate_class + log_len_kb, data = d,
          family = quasipoisson(link = "log")),
      warning = function(w) {
        if (grepl("did not converge|iteration limit", conditionMessage(w)))
          converged <<- FALSE
        invokeRestart("muffleWarning")
      })
  })

  if (!converged) {
    .msg("AVISO: el modelo (", type, ") NO convergió. Marcado como DEGENERADO.")
    return(list(model = model, coef_table = NULL, type = type,
                degenerate = TRUE,
                diagnostics = list(n_obs = n_obs, n_nonzero = n_nonzero,
                                   n_classes = n_classes),
                data = d[]))
  }

  sm <- summary(model)
  ct <- as.data.table(coef(sm), keep.rownames = "term")
  setnames(ct, names(ct)[2:min(5,ncol(ct))],
           c("estimate","std_error","stat","p_value")[seq_len(min(4,ncol(ct)-1))])
  # IRR (incidence rate ratio): exp(estimate), interpretable como factor multiplicativo.
  ct[, IRR := round(exp(estimate), 3)]

  .msg("Modelo (", type, "): clase + log_longitud -> n_SE. ",
       "Referencia = ", ref_class, ".")
  .msg("Interpretación: IRR de cada clase = factor de solapes SEdb vs referencia",
       " A IGUALDAD de longitud.")

  list(model = model, coef_table = ct[], type = type, degenerate = FALSE,
       diagnostics = list(n_obs = n_obs, n_nonzero = n_nonzero,
                          n_classes = n_classes),
       data = d[])
}

## ============================================================================
## (2) FIGURAS
## ============================================================================

#' Genera un boxplot de solapamientos con SEdb por clase dentro del bin 50-100 kb.
#' @return objeto ggplot.
fig_control_bin <- function(per_drr, target_bin = "(50000,100000]",
                            file = NULL) {
  d <- as.data.table(per_drr)
  sub <- d[as.character(size_bin) == target_bin]
  if (nrow(sub) == 0L) {
    # Si el bin objetivo no está disponible, se utiliza el bin de mayor longitud presente.
    sub <- d[size_bin == levels(size_bin)[max(as.integer(size_bin), na.rm=TRUE)]]
  }
  ord <- c("Simple_DRR","Compact_DRR","Extended_complex_DRR","Dense_complex_DRR")
  sub[, candidate_class := factor(candidate_class,
                                  levels = intersect(ord, unique(candidate_class)))]

  p <- ggplot(sub, aes(candidate_class, n_SE, fill = candidate_class)) +
    geom_boxplot(outlier.alpha = 0.3) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(x = "DRR class",
         y = "SEdb super-enhancers solapados (nº)",
         title = "Control por tamaño: solapamiento SEdb en DRRs de 50-100 kb",
         subtitle = "A igualdad de tamaño, las regiones densas recuperan más super-enhancers") +
    theme_bw(base_size = 12)

  if (!is.null(file)) { ggsave(file, p, width = 7, height = 5, dpi = 150)
    .msg("Figura guardada: ", file) }
  p
}

#' Resume los solapamientos con SEdb por clase a lo largo de los bins de longitud.
#' @return objeto ggplot.
fig_bins <- function(per_drr, file = NULL) {
  d <- as.data.table(per_drr)
  summ <- d[, .(median_SE = as.numeric(median(n_SE)),
                frac = mean(n_SE > 0), n = .N),
            by = .(size_bin, candidate_class)]
  ord <- c("Simple_DRR","Compact_DRR","Extended_complex_DRR","Dense_complex_DRR")
  summ[, candidate_class := factor(candidate_class,
                                   levels = intersect(ord, unique(candidate_class)))]

  p <- ggplot(summ, aes(size_bin, median_SE,
                        colour = candidate_class, group = candidate_class)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    scale_colour_brewer(palette = "Set2") +
    labs(x = "Bin de longitud de DRR (pb)",
         y = "Mediana de super-enhancers SEdb solapados",
         colour = "DRR class",
         title = "Solapamiento SEdb por clase y tamaño de DRR") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  if (!is.null(file)) { ggsave(file, p, width = 8, height = 5, dpi = 150)
    .msg("Figura guardada: ", file) }
  p
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("R/drr_lineB_stats_figs.R")
##
## # sm$per_drr debe proceder de la etapa de comparación con SEdb.
## mod <- model_density_controlled(sm$per_drr, ref_class = "Extended_complex_DRR")
## print(mod$coef_table)
##
## # Un IRR > 1 para candidate_classDense_complex_DRR indica un mayor número
## # esperado de solapamientos con SEdb respecto a la clase de referencia,
## # tras controlar por la longitud de la región.
##
## p1 <- fig_control_bin(sm$per_drr, file = "results/figB_control_50_100kb.png")
## p2 <- fig_bins(sm$per_drr,        file = "results/figB_bins.png")
## print(p1); print(p2)
## ============================================================================
