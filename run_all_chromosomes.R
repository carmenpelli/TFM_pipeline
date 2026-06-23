################################################################################
##  run_all_chromosomes.R — Ejecución cromosómica del pipeline
##  --------------------------------------------------------------------------
##  Este script ejecuta main.R para varios cromosomas mediante procesos R
##  independientes, utilizando Rscript main.R <chr>. Cada cromosoma se procesa
##  en un proceso aislado, lo que reduce interferencias de memoria y estado entre
##  ejecuciones.
##
##  El paralelismo se implementa con makeCluster/parLapply, compatible con
##  RStudio y entornos donde el uso de forking mediante mclapply puede ser
##  problemático.
##
##  Política de fallos:
##    Si algún cromosoma de una tanda falla, no se lanzan nuevas tandas. El
##    script informa del cromosoma afectado y conserva el log correspondiente
##    para su inspección.
##
##  Estrategia de ejecución:
##    Los cromosomas se ordenan según el tamaño de sus ficheros de CRMs y se
##    intercalan cromosomas grandes y pequeños para equilibrar el uso de memoria
##    entre tandas.
##
##  Uso desde RStudio:
##    setwd("~/TFM/TFM_pipeline")
##    source("run_all_chromosomes.R")
##
##    .args <- c("", "chr21,chr22")
##    source("run_all_chromosomes.R")
##
##    .args <- "force"
##    source("run_all_chromosomes.R")
##
##  Uso desde terminal:
##    Rscript run_all_chromosomes.R [force] [chr1,chr2,...]
##
##  Los logs por cromosoma se guardan en:
##    results/logs/<chr>.log
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# Argumentos de ejecución. Si existe .args, se respeta para permitir ejecución
# desde la consola de R; en caso contrario, se leen argumentos de terminal.
if (!exists(".args")) .args <- commandArgs(trailingOnly = TRUE)

FORCE     <- length(.args) >= 1 && tolower(.args[1]) %in% c("force","true","1")
ONLY_CHRS <- if (length(.args) >= 2 && nzchar(.args[2]))
               strsplit(.args[2], ",")[[1]] else NULL

N_PARALLEL <- 2L
ENH_DIR    <- "data/enh_per_chr"
LOG_DIR    <- "results/logs"
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# Ruta al binario Rscript de la instalación de R actualmente activa.
# Esto evita utilizar accidentalmente otra instalación de R disponible en PATH.
RSCRIPT_BIN <- file.path(R.home("bin"), "Rscript")

.msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                  paste0(...)))

## ---------------------------------------------------------------------------
## 1. Identificación y ordenación de cromosomas disponibles
## ---------------------------------------------------------------------------

all_chrs <- paste0("chr", c(1:22, "X"))

if (!is.null(ONLY_CHRS)) all_chrs <- intersect(all_chrs, ONLY_CHRS)

chr_size <- sapply(all_chrs, function(ch) {
  f <- file.path(ENH_DIR, paste0(ch, ".tsv.gz"))
  if (file.exists(f)) file.size(f) else NA_real_
})

present <- all_chrs[!is.na(chr_size)]

if (length(present) < length(all_chrs))
  .msg("Advertencia: faltan ficheros de CRMs para: ",
       paste(setdiff(all_chrs, present), collapse = ", "))

all_chrs <- present
chr_size <- chr_size[present]

if (length(all_chrs) == 0) {
  .msg("No hay cromosomas que procesar. Revisar directorio: ", ENH_DIR)
  quit(status = 1)
}

# Ordenación por tamaño e intercalado de cromosomas grandes y pequeños.
big <- all_chrs[order(chr_size, decreasing = TRUE)]

queue <- character(0)
i <- 1L
j <- length(big)

while (i <= j) {
  queue <- c(queue, big[i])
  if (i != j) queue <- c(queue, big[j])
  i <- i + 1L
  j <- j - 1L
}

.msg("Cromosomas seleccionados (", length(queue), "): ",
     paste(queue, collapse = ", "))

## ---------------------------------------------------------------------------
## 2. Ejecución de un cromosoma en un proceso R independiente
## ---------------------------------------------------------------------------

run_one <- function(chr, force_flag, log_dir, rscript_bin) {
  log_file <- file.path(log_dir, paste0(chr, ".log"))
  cmd_args <- c("main.R", chr)

  if (force_flag) cmd_args <- c(cmd_args, "force")

  # Se usa el mismo Rscript que la sesión actual y se elimina temporalmente
  # R_HOME heredado para evitar conflictos entre instalaciones de R.
  old_rhome <- Sys.getenv("R_HOME", unset = NA)
  Sys.unsetenv("R_HOME")
  on.exit(if (!is.na(old_rhome)) Sys.setenv(R_HOME = old_rhome), add = TRUE)

  status <- system2(rscript_bin, args = cmd_args,
                    stdout = log_file, stderr = log_file)

  list(chr = chr, status = status)
}

## ---------------------------------------------------------------------------
## 3. Procesamiento por tandas con política de parada ante error
## ---------------------------------------------------------------------------

done <- character(0)
failed <- character(0)

batches <- split(queue, ceiling(seq_along(queue) / N_PARALLEL))

for (b in seq_along(batches)) {
  batch <- batches[[b]]

  .msg("=== Tanda ", b, "/", length(batches), ": ",
       paste(batch, collapse = " + "), " ===")

  # Cluster PSOCK compatible con RStudio y entornos sin forking.
  cl <- makeCluster(min(N_PARALLEL, length(batch)), type = "PSOCK")

  # El directorio de trabajo de cada worker debe coincidir con la raíz del
  # proyecto, donde se encuentra main.R.
  wd <- getwd()
  clusterExport(cl, c("wd"), envir = environment())
  clusterEvalQ(cl, setwd(wd))

  res <- tryCatch(
    parLapply(cl, batch, run_one, force_flag = FORCE, log_dir = LOG_DIR,
              rscript_bin = RSCRIPT_BIN),
    error = function(e) {
      .msg("Error en parLapply: ", conditionMessage(e))
      NULL
    }
  )

  stopCluster(cl)

  if (is.null(res)) {
    failed <- c(failed, batch)
    break
  }

  for (r in res) {
    if (is.null(r$status) || r$status != 0) {
      failed <- c(failed, r$chr)
    } else {
      done <- c(done, r$chr)
    }
  }

  if (length(failed) > 0) {
    .msg("Fallo detectado en la tanda ", b,
         ". No se lanzarán nuevas tandas.")
    break
  }
}

## ---------------------------------------------------------------------------
## 4. Resumen final de ejecución
## ---------------------------------------------------------------------------

.msg("=================================================")
.msg("Completados (", length(done), "): ", paste(done, collapse = ", "))

if (length(failed) > 0) {
  .msg("Fallidos (", length(failed), "): ", paste(failed, collapse = ", "))
  for (ch in failed) {
    .msg("  -> revisar log: ", file.path(LOG_DIR, paste0(ch, ".log")))
  }
} else {
  .msg("Todos los cromosomas se completaron correctamente.")
}
