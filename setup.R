################################################################################
##  setup.R — preparación del entorno del pipeline TFM
##  --------------------------------------------------------------------------
##  Crea las carpetas necesarias y comprueba que los paquetes R están instalados.
##  Ejecutar UNA vez antes del primer main.R:   Rscript setup.R
################################################################################

# 1. Estructura de carpetas
dirs <- c("data", "data/tad_per_chr", "data/enh_per_chr",
          "results", "results/intermediate", "results/lineB",
          "results/checkpoints", "results/logs")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
cat("Carpetas creadas:\n  ", paste(dirs, collapse = "\n  "), "\n\n")

# 2. Comprobación de paquetes
required <- c("data.table", "GenomicRanges", "IRanges", "igraph",
              "MASS", "ggplot2")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]

if (length(missing) == 0) {
  cat("Todos los paquetes requeridos están instalados.\n")
} else {
  cat("FALTAN paquetes:\n  ", paste(missing, collapse = ", "), "\n\n")
  cat("Instálalos con:\n")
  bioc <- intersect(missing, c("GenomicRanges", "IRanges"))
  cran <- setdiff(missing, bioc)
  if (length(cran) > 0)
    cat('  install.packages(c(', paste(sprintf('"%s"', cran), collapse=", "), '))\n')
  if (length(bioc) > 0)
    cat('  BiocManager::install(c(', paste(sprintf('"%s"', bioc), collapse=", "), '))\n')
}

# 3. Recordatorio de datos de entrada
cat("\nColoca en data/ tus ficheros de entrada:\n")
cat("  data/tad_per_chr/<chr>.tsv.gz  (TADs crudos por cromosoma)\n")
cat("  data/enh_per_chr/<chr>.tsv.gz  (CRMs crudos por cromosoma)\n")
cat("  data/enh2gene.tsv.gz           (CRM->gen, global)\n")
cat("  data/SEdb_Human_SE.bed         (super-enhancers SEdb, global)\n")
cat("  data/human_genes.tsv           (genes, global)\n")
cat("\nLuego ejecuta un cromosoma:   Rscript main.R chr8\n")
cat("O todos en paralelo:          Rscript run_all_chromosomes.R\n")
