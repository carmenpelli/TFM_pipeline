################################################################################
##  permutacion_sedb.R — TEST DE PERMUTACIÓN CONTROLADO POR VECINDAD
##  ==========================================================================
##  Objetivo: descartar que el mayor solapamiento de las DRRs densas con SEdb
##  se deba al CONTEXTO GENÓMICO (vecindad) en lugar de a su densidad intrínseca.
##
##  Estrategia (label shuffling estratificado por TAD reducido):
##    - Cada DRR se asigna a un TAD reducido (su dominio topológico = "vecindad").
##    - Dentro de cada TAD reducido se BARAJAN las etiquetas de clase entre las DRRs,
##      preservando cuántas DRRs hay de cada clase en ese TAD reducido.
##    - Para cada permutación se recalcula el estadístico de interés por clase
##      (p. ej. fracción de DRRs con >=1 super-enhancer SEdb, o nº medio de SE).
##    - Se compara el valor OBSERVADO con la distribución nula de B permutaciones.
##
##  Interpretación: si dentro del mismo TAD reducido (mismo "barrio") las DRRs densas
##  siguen recuperando más SEdb que lo esperado al barajar etiquetas, el
##  enriquecimiento NO se explica por la vecindad genómica.
##
##  Reproducible 100% desde los pipeline_objects_<chr>.rds. No requiere el BED.
##
##  USO:
##    source("permutacion_sedb.R")
##    res <- run_permutation_all(results_dir = "results", B = 2000, seed = 1)
##    print(res$summary)
##    p <- plot_permutation(res); print(p)
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
## 1. Construir la tabla DRR -> (clase, n_SE, TAD reducido) para un cromosoma
## ---------------------------------------------------------------------------
#' @return data.table: drr_id, candidate_class, n_SE, has_SE, tad_id, chr
build_perm_table <- function(o, chr) {
  pd <- as.data.table(o$sm$per_drr)[, .(drr_id, candidate_class, n_SE)]
  dl <- as.data.table(o$prop$drr_tad_long)[, .(drr_id, tad_id, n_crms_in_tad)]

  # Una DRR puede tocar varios TADs reducidos: asignamos al TAD reducido donde tiene más
  # contexto (n_crms_in_tad máximo); desempate por tad_id menor.
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
#' Fracción de DRRs con >=1 SE y nº medio de SE, por clase.
class_stats <- function(dt) {
  dt[, .(frac_withSE = mean(has_SE),
         mean_nSE    = mean(n_SE),
         n           = .N), by = candidate_class]
}

## ---------------------------------------------------------------------------
## 3. Una permutación: barajar etiquetas DENTRO de cada TAD reducido
## ---------------------------------------------------------------------------
#' Devuelve un vector de clases permutadas, barajando dentro de cada tad_id.
#' TADs reducidos con una sola DRR no aportan (su etiqueta no cambia): es correcto,
#' esos casos no informan sobre vecindad.
permute_within_tad <- function(dt) {
  dt[, .(perm_class = sample(candidate_class)), by = tad_id,
     .SDcols = "candidate_class"]$perm_class
}

## ---------------------------------------------------------------------------
## 4. Test de permutación para un conjunto de DRRs (uno o varios cromosomas)
## ---------------------------------------------------------------------------
#' @param dt  tabla build_perm_table (puede combinar cromosomas: barajado por
#'            (chr,tad_id) para no mezclar dominios de distintos cromosomas).
#' @param B   nº de permutaciones.
#' @param stat "frac" (fracción con SE) o "mean" (nº medio de SE).
#' @return lista con observado, nula (matriz B x clases), p-valores y z-scores.
permutation_test <- function(dt, B = 2000, stat = c("frac", "mean"), seed = 1) {
  stat <- match.arg(stat)
  set.seed(seed)
  dt <- copy(dt)
  # clave de estrato: TAD reducido dentro de cromosoma
  dt[, stratum := paste(chr, tad_id, sep = ":")]

  obs_tab <- class_stats(dt)
  obs <- setNames(if (stat == "frac") obs_tab$frac_withSE else obs_tab$mean_nSE,
                  obs_tab$candidate_class)
  classes <- CLASS_ORDER[CLASS_ORDER %in% names(obs)]
  obs <- obs[classes]

  # respuesta a permutar (fija) y agrupación
  y <- if (stat == "frac") dt$has_SE else dt$n_SE
  grp <- dt$stratum
  cls <- dt$candidate_class

  # precomputar índices por estrato para barajar rápido
  idx_by_stratum <- split(seq_len(nrow(dt)), grp)

  null_mat <- matrix(NA_real_, nrow = B, ncol = length(classes),
                     dimnames = list(NULL, classes))
  for (b in seq_len(B)) {
    perm_cls <- cls
    for (ix in idx_by_stratum) if (length(ix) > 1L)
      perm_cls[ix] <- sample(cls[ix])
    # media de y por clase permutada
    agg <- tapply(y, perm_cls, mean)
    null_mat[b, ] <- agg[classes]
  }

  # p-valor de cola superior (enriquecimiento) con corrección +1
  pval <- sapply(classes, function(c)
    (sum(null_mat[, c] >= obs[c], na.rm = TRUE) + 1) / (B + 1))
  zsc <- sapply(classes, function(c) {
    mu <- mean(null_mat[, c], na.rm = TRUE); sdv <- sd(null_mat[, c], na.rm = TRUE)
    if (is.na(sdv) || sdv == 0) NA_real_ else (obs[c] - mu) / sdv
  })

  list(stat = stat, B = B, classes = classes, observed = obs,
       null = null_mat, p_value = pval, z = zsc,
       null_mean = colMeans(null_mat, na.rm = TRUE))
}

## ---------------------------------------------------------------------------
## 5. Runner: carga todos los cromosomas y ejecuta el test global
## ---------------------------------------------------------------------------
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
  # cuántas DRRs están en TADs reducidos con mezcla de clases (las informativas)
  dt[, stratum := paste(chr, tad_id, sep = ":")]
  mix <- dt[, .(k = uniqueN(candidate_class), n = .N), by = stratum]
  .msg("TADs reducidos con \u22652 clases (informativos): ",
       sum(mix$k >= 2), " de ", nrow(mix))

  .msg("Ejecutando ", B, " permutaciones (estratificadas por TAD reducido)...")
  res <- permutation_test(dt, B = B, stat = stat, seed = seed)

  summary_tab <- data.table(
    clase       = res$classes,
    observado   = round(res$observed, 4),
    nula_media  = round(res$null_mean, 4),
    z           = round(res$z, 2),
    p_perm      = signif(res$p_value, 3))
  res$summary <- summary_tab[]
  res$data <- dt
  .msg("Hecho.")
  res
}

## ---------------------------------------------------------------------------
## 6. Figura: distribución nula vs observado por clase
## ---------------------------------------------------------------------------
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
    labs(title = "Test de permutación controlado por vecindad (TAD reducido)",
         subtitle = paste0("Distribución nula (", res$B,
                           " permutaciones de etiquetas dentro de cada TAD reducido) vs valor observado (línea)"),
         x = ylab, y = "frecuencia (permutaciones)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))
}
