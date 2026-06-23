################################################################################
##  run_all_chromosomes.R — ESCALADO A TODOS LOS CROMOSOMAS (compatible RStudio)
##  --------------------------------------------------------------------------
##  Ejecuta main.R para cada cromosoma, en paralelo de 2 en 2 (2 CPUs), usando
##  procesos R aislados (Rscript main.R chrN). Cada cromosoma corre en su propio
##  proceso, aislado en memoria y estado — lo mas robusto para corridas largas.
##
##  Paralelismo con makeCluster/parLapply (NO mclapply): funciona dentro de
##  RStudio, que no lleva bien el forking de mclapply.
##
##  Politica de fallos (fail-stop): si un cromosoma de una tanda falla, NO se
##  lanza la siguiente tanda; se reporta cual fallo y se detiene la corrida.
##
##  Estrategia de memoria: ordenamos por tamano e intercalamos grande+pequeno,
##  para que cada tanda mezcle un cromosoma pesado con uno ligero.
##
##  USO (desde la consola de RStudio):
##    setwd("~/TFM/TFM_pipeline")
##    source("run_all_chromosomes.R")                          # todos
##    .args <- c("", "chr21,chr22"); source("run_all_chromosomes.R")  # algunos
##    .args <- "force";              source("run_all_chromosomes.R")  # forzado
##
##  USO (desde terminal):
##    Rscript run_all_chromosomes.R [force] [chr1,chr2,...]
##
##  Logs por cromosoma en results/logs/<chr>.log
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# Argumentos: respeta un .args ya definido (consola de R); si no, lee de terminal.
if (!exists(".args")) .args <- commandArgs(trailingOnly = TRUE)
FORCE     <- length(.args) >= 1 && tolower(.args[1]) %in% c("force","true","1")
ONLY_CHRS <- if (length(.args) >= 2 && nzchar(.args[2]))
               strsplit(.args[2], ",")[[1]] else NULL

N_PARALLEL <- 2L
ENH_DIR    <- "data/enh_per_chr"
LOG_DIR    <- "results/logs"
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# Ruta al Rscript de ESTA instalacion de R (evita usar otro R del PATH que
# provoque conflictos de version al cargar paquetes como data.table).
RSCRIPT_BIN <- file.path(R.home("bin"), "Rscript")

.msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                  paste0(...)))

## ---------------------------------------------------------------------------
## 1. Cromosomas presentes, ordenados (intercalando grande + pequeno)
## ---------------------------------------------------------------------------
all_chrs <- paste0("chr", c(1:22, "X"))
if (!is.null(ONLY_CHRS)) all_chrs <- intersect(all_chrs, ONLY_CHRS)

chr_size <- sapply(all_chrs, function(ch) {
  f <- file.path(ENH_DIR, paste0(ch, ".tsv.gz")); if (file.exists(f)) file.size(f) else NA_real_
})
present <- all_chrs[!is.na(chr_size)]
if (length(present) < length(all_chrs))
  .msg("AVISO: faltan ficheros de CRMs para: ",
       paste(setdiff(all_chrs, present), collapse = ", "))
all_chrs <- present; chr_size <- chr_size[present]

if (length(all_chrs) == 0) { .msg("No hay cromosomas que procesar. Revisa ", ENH_DIR); quit(status = 1) }

big <- all_chrs[order(chr_size, decreasing = TRUE)]
queue <- character(0); i <- 1L; j <- length(big)
while (i <= j) {
  queue <- c(queue, big[i]); if (i != j) queue <- c(queue, big[j]); i <- i + 1L; j <- j - 1L
}
.msg("Cromosomas (", length(queue), "): ", paste(queue, collapse = ", "))

## ---------------------------------------------------------------------------
## 2. Funcion que ejecuta UN cromosoma (proceso R aislado via Rscript)
## ---------------------------------------------------------------------------
run_one <- function(chr, force_flag, log_dir, rscript_bin) {
  log_file <- file.path(log_dir, paste0(chr, ".log"))
  cmd_args <- c("main.R", chr)
  if (force_flag) cmd_args <- c(cmd_args, "force")
  # Usar el MISMO Rscript que la sesion actual y limpiar R_HOME heredado
  # (evita el conflicto 'ignoring environment value of R_HOME' que rompe
  #  la carga de data.table cuando el subproceso hereda otro R_HOME).
  old_rhome <- Sys.getenv("R_HOME", unset = NA)
  Sys.unsetenv("R_HOME")
  on.exit(if (!is.na(old_rhome)) Sys.setenv(R_HOME = old_rhome), add = TRUE)
  status <- system2(rscript_bin, args = cmd_args,
                    stdout = log_file, stderr = log_file)
  list(chr = chr, status = status)      # status 0 = OK
}

## ---------------------------------------------------------------------------
## 3. Procesar por tandas de N_PARALLEL con makeCluster, parando ante fallo
## ---------------------------------------------------------------------------
done <- character(0); failed <- character(0)
batches <- split(queue, ceiling(seq_along(queue) / N_PARALLEL))

for (b in seq_along(batches)) {
  batch <- batches[[b]]
  .msg("=== Tanda ", b, "/", length(batches), ": ",
       paste(batch, collapse = " + "), " ===")

  # Cluster de tantos workers como cromosomas en la tanda (max N_PARALLEL).
  # PSOCK funciona en RStudio (no usa forking).
  cl <- makeCluster(min(N_PARALLEL, length(batch)), type = "PSOCK")
  # El working dir del worker debe ser el del proyecto (donde esta main.R)
  wd <- getwd()
  clusterExport(cl, c("wd"), envir = environment())
  clusterEvalQ(cl, setwd(wd))

  res <- tryCatch(
    parLapply(cl, batch, run_one, force_flag = FORCE, log_dir = LOG_DIR,
              rscript_bin = RSCRIPT_BIN),
    error = function(e) { .msg("Error en parLapply: ", conditionMessage(e)); NULL }
  )
  stopCluster(cl)

  if (is.null(res)) { failed <- c(failed, batch); break }

  for (r in res) {
    if (is.null(r$status) || r$status != 0) failed <- c(failed, r$chr)
    else done <- c(done, r$chr)
  }

  if (length(failed) > 0) {
    .msg("Fallo detectado en la tanda ", b, ". Fail-stop: no se lanzan mas tandas.")
    break
  }
}

## ---------------------------------------------------------------------------
## 4. Resumen final
## ---------------------------------------------------------------------------
.msg("=================================================")
.msg("Completados (", length(done), "): ", paste(done, collapse = ", "))
if (length(failed) > 0) {
  .msg("FALLARON (", length(failed), "): ", paste(failed, collapse = ", "))
  for (ch in failed) .msg("  -> revisa: ", file.path(LOG_DIR, paste0(ch, ".log")))
} else {
  .msg("Todos los cromosomas completados correctamente.")
}
