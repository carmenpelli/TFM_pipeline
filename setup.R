################################################################################
##  setup.R — Preparación del entorno del pipeline
##  --------------------------------------------------------------------------
##  Este script crea la estructura básica de directorios y comprueba que las
##  dependencias principales de R están disponibles.
##
##  Uso:
##    Rscript setup.R
################################################################################

## ---------------------------------------------------------------------------
## 1. Creación de la estructura de directorios
## ---------------------------------------------------------------------------

dirs <- c(
  "data",
  "data/tad_per_chr",
  "data/enh_per_chr",
  "results",
  "results/intermediate",
  "results/lineB",
  "results/checkpoints",
  "results/logs"
)

for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("Directorios creados o ya existentes:\n  ",
    paste(dirs, collapse = "\n  "), "\n\n")

## ---------------------------------------------------------------------------
## 2. Comprobación de dependencias de R
## ---------------------------------------------------------------------------

required <- c(
  "data.table",
  "GenomicRanges",
  "IRanges",
  "igraph",
  "MASS",
  "ggplot2"
)

missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]

if (length(missing) == 0) {
  cat("Todas las dependencias requeridas están instaladas.\n")
} else {
  cat("Dependencias no disponibles:\n  ",
      paste(missing, collapse = ", "), "\n\n")

  cat("Comandos sugeridos para la instalación:\n")

  bioc <- intersect(missing, c("GenomicRanges", "IRanges"))
  cran <- setdiff(missing, bioc)

  if (length(cran) > 0)
    cat('  install.packages(c(',
        paste(sprintf('"%s"', cran), collapse = ", "),
        '))\n', sep = "")

  if (length(bioc) > 0)
    cat('  BiocManager::install(c(',
        paste(sprintf('"%s"', bioc), collapse = ", "),
        '))\n', sep = "")
}

## ---------------------------------------------------------------------------
## 3. Estructura esperada de datos de entrada
## ---------------------------------------------------------------------------

cat("\nEstructura esperada de ficheros de entrada:\n")
cat("  data/tad_per_chr/<chr>.tsv.gz  (TADs por cromosoma)\n")
cat("  data/enh_per_chr/<chr>.tsv.gz  (CRMs por cromosoma)\n")
cat("  data/enh2gene.tsv.gz           (relaciones CRM-gen)\n")
cat("  data/SEdb_Human_SE.bed         (anotaciones de super-enhancers de SEdb)\n")
cat("  data/human_genes.tsv           (anotaciones génicas)\n")

cat("\nEjecución de un cromosoma:\n")
cat("  Rscript main.R chr8\n")

cat("\nEjecución cromosómica general:\n")
cat("  Rscript run_all_chromosomes.R\n")
