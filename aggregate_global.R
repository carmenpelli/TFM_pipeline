################################################################################
##  aggregate_global.R — Agregación global de resultados cromosómicos
##  --------------------------------------------------------------------------
##  Este script integra los resultados generados por main.R o
##  run_all_chromosomes.R y produce un análisis global del genoma.
##
##  Salidas principales:
##    - Tabla global de DRRs por clase y cromosoma.
##    - Modelo global de densidad basado en regresión binomial negativa,
##      ajustado sobre el conjunto completo de DRRs.
##    - Resúmenes globales de número de DRRs por clase, concordancia con
##      SEdb y propiedades topológicas.
##    - Figuras globales asociadas al análisis de clases.
##
##  Entrada esperada:
##    results/lineB/<chr>/pipeline_objects_<chr>.rds
##
##  Uso en R:
##    setwd("~/TFM/TFM_pipeline")
##    source("aggregate_global.R")
##
##  Uso en terminal:
##    Rscript aggregate_global.R
################################################################################

suppressPackageStartupMessages({
  library(data.table)
})
source(file.path("R", "drr_lineB_stats_figs.R"))   # funciones de modelado y figuras

.msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                  paste0(...)))

RESULTS  <- "results"
OUT_DIR  <- file.path(RESULTS, "global")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. Identificación de objetos RDS por cromosoma
## ---------------------------------------------------------------------------
rds_files <- list.files(file.path(RESULTS, "lineB"),
                        pattern = "^pipeline_objects_.*\\.rds$",
                        recursive = TRUE, full.names = TRUE)
if (length(rds_files) == 0)
  stop("No hay resultados por cromosoma en ", file.path(RESULTS, "lineB"),
       ". Ejecuta primero main.R / run_all_chromosomes.R.")

# Extraer los nombres de cromosoma a partir de las rutas.
chr_of <- sub(".*pipeline_objects_(.*)\\.rds$", "\\1", rds_files)
.msg("Cromosomas encontrados (", length(rds_files), "): ",
     paste(chr_of, collapse = ", "))

## ---------------------------------------------------------------------------
## 2. Integración de tablas por DRR y resúmenes cromosómicos
## ---------------------------------------------------------------------------
per_drr_list <- list()
byclass_list <- list()
genic_list   <- list()
topo_list    <- list()

for (k in seq_along(rds_files)) {
  ch  <- chr_of[k]
  obj <- readRDS(rds_files[k])

  # Tabla por DRR utilizada para el modelo global; se añade el cromosoma.
  if (!is.null(obj$sm) && !is.null(obj$sm$per_drr)) {
    pd <- as.data.table(copy(obj$sm$per_drr)); pd[, chr := ch]
    per_drr_list[[ch]] <- pd
  }
  # Número de DRRs por clase.
  if (!is.null(obj$cmp) && !is.null(obj$cmp$by_class_stack)) {
    bc <- as.data.table(copy(obj$cmp$by_class_stack)); bc[, chr := ch]
    byclass_list[[ch]] <- bc
  }
  # Resumen de concordancia génica, cuando está disponible.
  if (!is.null(obj$val) && !is.null(obj$val$genic) &&
      !is.null(obj$val$genic$summary)) {
    g <- as.data.table(copy(obj$val$genic$summary)); g[, chr := ch]
    genic_list[[ch]] <- g
  }
  # Resumen topológico por clase.
  if (!is.null(obj$prop) && !is.null(obj$prop$by_class)) {
    tp <- as.data.table(copy(obj$prop$by_class)); tp[, chr := ch]
    topo_list[[ch]] <- tp
  }
}

per_drr_all <- rbindlist(per_drr_list, use.names = TRUE, fill = TRUE)
.msg("DRRs totales (todas las clases, todos los chr): ", nrow(per_drr_all))

## ---------------------------------------------------------------------------
## 3. Modelo global de densidad mediante regresión binomial negativa
## ---------------------------------------------------------------------------
.msg("Ajustando modelo global de densidad...")
mod_global <- model_density_controlled(per_drr_all,
                                       ref_class = "Extended_complex_DRR")
fwrite(mod_global$coef_table, file.path(OUT_DIR, "model_coef_global.tsv"), sep = "\t")
.msg("Modelo global (", mod_global$type, "):")
print(mod_global$coef_table)

# Generación de figuras globales.
fig_control_bin(per_drr_all, file = file.path(OUT_DIR, "figB_control_global.png"))
fig_bins(per_drr_all,        file = file.path(OUT_DIR, "figB_bins_global.png"))

## ---------------------------------------------------------------------------
## 4. Generación de resúmenes globales
## ---------------------------------------------------------------------------
# (a) Número total de DRRs por clase, agregando todos los cromosomas.
if (length(byclass_list) > 0) {
  byclass_all <- rbindlist(byclass_list, use.names = TRUE, fill = TRUE)
  drr_by_class_global <- byclass_all[, .(
    n_DRRs        = sum(n),
    median_len    = as.numeric(median(rep(median_len, n))),   # mediana ponderada aproximada
    n_chr         = uniqueN(chr)
  ), by = candidate_class][order(-n_DRRs)]
  fwrite(drr_by_class_global, file.path(OUT_DIR, "drr_by_class_global.tsv"), sep = "\t")
  .msg("DRRs por clase (global):")
  print(drr_by_class_global)
}

# (b) Concordancia génica global por clase, agregando candidatos y coincidencias.
if (length(genic_list) > 0) {
  genic_all <- rbindlist(genic_list, use.names = TRUE, fill = TRUE)
  genic_global <- genic_all[, .(
    n_candidates      = sum(n_candidates),
    n_with_gene_match = sum(n_with_gene_match),
    fraction_match    = round(sum(n_with_gene_match) / sum(n_candidates), 4)
  ), by = candidate_class][order(-fraction_match)]
  fwrite(genic_global, file.path(OUT_DIR, "gene_validation_global.tsv"), sep = "\t")
  .msg("Recuperación génica (global):")
  print(genic_global)
}

# (c) Resumen topológico global por clase.
if (length(topo_list) > 0) {
  topo_all <- rbindlist(topo_list, use.names = TRUE, fill = TRUE)
  # Se contemplan nombres alternativos generados por distintas versiones del pipeline.
  num_cols <- intersect(c("n_drr", "n_single", "n_single_subtad"), names(topo_all))
  if (length(num_cols) >= 2) {
    topo_global <- topo_all[, lapply(.SD, sum), by = candidate_class,
                            .SDcols = num_cols]
    fwrite(topo_global, file.path(OUT_DIR, "topology_global.tsv"), sep = "\t")
    .msg("Topología (global):"); print(topo_global)
  } else {
    fwrite(topo_all, file.path(OUT_DIR, "topology_all_chr.tsv"), sep = "\t")
    .msg("Topología: guardada por cromosoma; no se aplicó agregación automática por incompatibilidad de nombres de columna.")
  }
}

# Guardar la tabla maestra de DRRs y el objeto global del análisis.
fwrite(per_drr_all, file.path(OUT_DIR, "per_drr_all_chr.tsv.gz"), sep = "\t")
saveRDS(list(per_drr_all = per_drr_all, mod_global = mod_global,
             chrs = chr_of),
        file.path(OUT_DIR, "global_objects.rds"))

.msg("===== AGREGACIÓN GLOBAL COMPLETADA =====")
.msg("Resultados en: ", OUT_DIR)
