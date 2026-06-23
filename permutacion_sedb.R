################################################################################
##  permutacion_sedb.R — Test de permutación estratificado por TAD reducido
##  ==========================================================================
##  Este script evalúa si la mayor concordancia de determinadas clases de DRR
##  con anotaciones de SEdb puede explicarse por su localización en regiones
##  genómicas con mayor carga reguladora, o si persiste tras controlar por el
##  contexto topológico local.
##
##  Estrategia:
##    - Cada DRR se asigna a un TAD reducido, utilizado como estrato topológico.
##    - Dentro de cada TAD reducido se permutan las etiquetas de clase entre las
##      DRRs, manteniendo fijas sus coordenadas, su número de solapamientos con
##      SEdb y su pertenencia al TAD.
##    - En cada permutación se recalcula el estadístico de interés por clase,
##      como la fracción de DRRs con al menos un solapamiento con SEdb o el
##      número medio de solapamientos.
##    - El valor observado se compara con la distribución nula generada mediante
##      B permutaciones.
##
##  Interpretación:
##    Si una clase mantiene valores observados superiores a los esperados bajo
##    permutación dentro del mismo TAD reducido, la asociación con SEdb no puede
##    atribuirse únicamente a la vecindad genómica local.
##
##  El análisis se reconstruye a partir de los objetos:
##    results/lineB/<chr>/pipeline_objects_<chr>.rds
##
##  Uso:
##    source("permutacion_sedb.R")
##    res <- run_permutation_all(results_dir = "results", B = 2000, seed = 1)
##    print(res$summary)
##    p <- plot_permutation(res)
##    print(p)
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!exists(".msg")) .msg <- function(...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

CLASS_ORDER <- c("Simple_DRR", "Compact_DRR",
                 "Extended_complex_DRR", "Dense_complex_DRR")

## ---------------------------------------------------------------------------
## 1. Construcción de la tabla DRR -> clase, SEdb y TAD reducido
## ---------------------------------------------------------------------------

#' Construye la tabla de entrada para el test de permutación.
#'
#' @return data.table con drr_id, candidate_class, n_SE, has_SE, tad_id y chr.
build_perm_table <- function(o, chr) {
  pd <- as.data.table(o$sm$per_drr)[, .(drr_id, candidate_class, n_SE)]
  dl <- as.data.table(o$prop$drr_tad_long)[, .(drr_id, tad_id, n_crms_in_tad)]

  # Una DRR puede estar asociada a más de un TAD reducido. En ese caso se
  # asigna al TAD con mayor contribución de CRMs; en caso de empate se utiliza
  # el identificador de TAD como criterio determinista.
  setorder(dl, drr_id, -n_crms_in_tad, tad_id)
  dl1 <- dl[, .SD[1L], by = drr_id]

  dt <- merge(pd, dl1[, .(drr_id, tad_id)], by = "drr_id", all.x = TRUE)
  dt <- dt[!is.na(tad_id) & candidate_class %in% CLASS_ORDER]
  dt[, has_SE := as.integer(n_SE > 0)]
  dt[, chr := chr]
  dt[]
}

## ---------------------------------------------------------------------------
## 2. Estadístico observado por clase
## ---------------------------------------------------------------------------

#' Calcula estadísticos de concordancia con SEdb por clase.
#'
#' Se resumen la fracción de DRRs con al menos un solapamiento con SEdb y el
#' número medio de solapamientos por DRR.
class_stats <- function(dt) {
  dt[, .(frac_withSE = mean(has_SE),
         mean_nSE    = mean(n_SE),
         n           = .N), by = candidate_class]
}

## ---------------------------------------------------------------------------
## 3. Permutación de etiquetas dentro de cada TAD reducido
## ---------------------------------------------------------------------------

#' Permuta las etiquetas de clase dentro de cada TAD reducido.
#'
#' Los TADs reducidos que contienen una sola DRR no modifican su etiqueta y no
#' contribuyen información al contraste estratificado.
permute_within_tad <- function(dt) {
  dt[, .(perm_class = sample(candidate_class)), by = tad_id,
     .SDcols = "candidate_class"]$perm_class
}

## ---------------------------------------------------------------------------
## 4. Test de permutación para uno o varios cromosomas
## ---------------------------------------------------------------------------

#' Ejecuta el test de permutación estratificado por TAD reducido.
#'
#' @param dt Tabla generada por build_perm_table. Puede combinar cromosomas; el
#'        barajado se realiza por estratos definidos como chr:tad_id.
#' @param B Número de permutaciones.
#' @param stat Estadístico a evaluar: "frac" para fracción con SEdb o "mean"
#'        para número medio de solapamientos.
#' @param seed Semilla para reproducibilidad.
#'
#' @return Lista con valores observados, matriz nula, p-valores empíricos,
#'         z-scores y medias de la distribución nula.
permutation_test <- function(dt, B = 2000, stat = c("frac", "mean"), seed = 1) {
  stat <- match.arg(stat)
  set.seed(seed)
  dt <- copy(dt)

  # Estrato de permutación: TAD reducido dentro de cada cromosoma.
  dt[, stratum := paste(chr, tad_id, sep = ":")]

  obs_tab <- class_stats(dt)
  obs <- setNames(if (stat == "frac") obs_tab$frac_withSE else obs_tab$mean_nSE,
                  obs_tab$candidate_class)
  classes <- CLASS_ORDER[CLASS_ORDER %in% names(obs)]
  obs <- obs[classes]

  # Respuesta fija y etiquetas de clase a permutar.
  y <- if (stat == "frac") dt$has_SE else dt$n_SE
  grp <- dt$stratum
  cls <- dt$candidate_class

  # Índices por estrato para acelerar las permutaciones.
  idx_by_stratum <- split(seq_len(nrow(dt)), grp)

  null_mat <- matrix(NA_real_, nrow = B, ncol = length(classes),
                     dimnames = list(NULL, classes))

  for (b in seq_len(B)) {
    perm_cls <- cls
    for (ix in idx_by_stratum) if (length(ix) > 1L)
      perm_cls[ix] <- sample(cls[ix])

    # Media del estadístico de respuesta por clase permutada.
    agg <- tapply(y, perm_cls, mean)
    null_mat[b, ] <- agg[classes]
  }

  # p-valor empírico unilateral de cola superior, con corrección +1.
  pval <- sapply(classes, function(c)
    (sum(null_mat[, c] >= obs[c], na.rm = TRUE) + 1) / (B + 1))

  zsc <- sapply(classes, function(c) {
    mu <- mean(null_mat[, c], na.rm = TRUE)
    sdv <- sd(null_mat[, c], na.rm = TRUE)
    if (is.na(sdv) || sdv == 0) NA_real_ else (obs[c] - mu) / sdv
  })

  list(stat = stat, B = B, classes = classes, observed = obs,
       null = null_mat, p_value = pval, z = zsc,
       null_mean = colMeans(null_mat, na.rm = TRUE))
}

## ---------------------------------------------------------------------------
## 5. Ejecución global del test de permutación
## ---------------------------------------------------------------------------

#' Carga los resultados cromosómicos y ejecuta el test global.
#'
#' @param results_dir Directorio raíz de resultados.
#' @param B Número de permutaciones.
#' @param stat Estadístico a evaluar.
#' @param seed Semilla para reproducibilidad.
#' @param chrs Vector opcional de cromosomas a incluir.
run_permutation_all <- function(results_dir = "results", B = 2000,
                                stat = "frac", seed = 1, chrs = NULL) {
  rds <- list.files(file.path(results_dir, "lineB"),
                    pattern = "^pipeline_objects_.*\\.rds$",
                    recursive = TRUE, full.names = TRUE)

  if (!is.null(chrs)) {
    keep <- sub(".*pipeline_objects_(.*)\\.rds$", "\\1", rds) %in% chrs
    rds <- rds[keep]
  }

  .msg("Construyendo tabla DRR\u2192TAD reducido de ", length(rds), " cromosomas...")

  dt <- rbindlist(lapply(rds, function(f) {
    ch <- sub(".*pipeline_objects_(.*)\\.rds$", "\\1", f)
    build_perm_table(readRDS(f), ch)
  }), use.names = TRUE, fill = TRUE)

  .msg("DRRs analizables: ", nrow(dt), " en ",
       length(unique(paste(dt$chr, dt$tad_id))), " TADs reducidos.")

  # Número de TADs reducidos informativos, definidos como aquellos con al menos
  # dos clases de DRR representadas.
  dt[, stratum := paste(chr, tad_id, sep = ":")]
  mix <- dt[, .(k = uniqueN(candidate_class), n = .N), by = stratum]

  .msg("TADs reducidos con \u22652 clases: ",
       sum(mix$k >= 2), " de ", nrow(mix))

  .msg("Ejecutando ", B, " permutaciones estratificadas por TAD reducido...")

  res <- permutation_test(dt, B = B, stat = stat, seed = seed)

  summary_tab <- data.table(
    clase       = res$classes,
    observado   = round(res$observed, 4),
    nula_media  = round(res$null_mean, 4),
    z           = round(res$z, 2),
    p_perm      = signif(res$p_value, 3))

  res$summary <- summary_tab[]
  res$data <- dt

  .msg("Test de permutación finalizado.")

  res
}

## ---------------------------------------------------------------------------
## 6. Representación de la distribución nula y del valor observado
## ---------------------------------------------------------------------------

#' Genera una figura de la distribución nula frente al valor observado.
plot_permutation <- function(res) {
  CLASS_COLOR <- c(Simple_DRR = "#7BA7C7", Compact_DRR = "#6BAE8E",
                   Extended_complex_DRR = "#E0A458", Dense_complex_DRR = "#C44E52")

  CLASS_LABEL <- c(Simple_DRR = "Simple", Compact_DRR = "Compact",
                   Extended_complex_DRR = "Extended", Dense_complex_DRR = "Dense complex")

  nm <- melt(as.data.table(res$null), measure.vars = res$classes,
             variable.name = "clase", value.name = "valor")

  obs <- data.table(clase = res$classes, obs = res$observed,
                    p = res$p_value)

  nm[, clase := factor(clase, levels = res$classes, labels = CLASS_LABEL[res$classes])]
  obs[, clase := factor(clase, levels = res$classes, labels = CLASS_LABEL[res$classes])]

  obs[, etiqueta := ifelse(p < 1 / (res$B + 1) * 1.5,
                           paste0("p < ", signif(1/(res$B+1), 2)),
                           paste0("p = ", signif(p, 2)))]

  ylab <- if (res$stat == "frac")
    "fracción de DRRs con \u22651 super-enhancer SEdb" else "nº medio de SE por DRR"

  ggplot(nm, aes(valor)) +
    geom_histogram(aes(fill = clase), bins = 40, alpha = 0.85,
                   color = "white", linewidth = 0.15) +
    geom_vline(data = obs, aes(xintercept = obs), color = "#222222",
               linewidth = 0.8) +
    geom_text(data = obs, aes(x = obs, y = Inf, label = etiqueta),
              hjust = 1.05, vjust = 1.6, size = 3, color = "#222222") +
    facet_wrap(~ clase, scales = "free", ncol = 2) +
    scale_fill_manual(values = setNames(CLASS_COLOR[res$classes],
                                        CLASS_LABEL[res$classes]), guide = "none") +
    labs(title = "Test de permutación estratificado por TAD reducido",
         subtitle = paste0("Distribución nula (", res$B,
                           " permutaciones dentro de cada TAD reducido) frente al valor observado"),
         x = ylab, y = "frecuencia de permutaciones") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))
}
