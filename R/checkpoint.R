################################################################################
##  checkpoint.R — EJECUCIÓN RECUPERABLE (CHECKPOINTS)
##  --------------------------------------------------------------------------
##  Permite que el pipeline RETOME donde se quedó: cada paso costoso se envuelve
##  con cache_step(), que guarda su resultado en un .rds y, en ejecuciones
##  posteriores, lo CARGA en vez de recalcular.
##
##  Uso:
##    red_tad <- cache_step("tad_reduced", CHR, force = FORCE_RECOMPUTE, {
##      reduce_redundancy_tads(tad_unique, thr_recip = 0.80)
##    })
##  - Si existe el checkpoint y force=FALSE -> lo carga (rápido).
##  - Si no existe (o force=TRUE) -> evalúa la expresión, guarda y devuelve.
##
##  Dependencias: ninguna (base R).
################################################################################

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# Directorio de checkpoints (se fija desde main.R vía options o variable global)
.checkpoint_dir <- function() {
  d <- getOption("tfm.checkpoint_dir", "results/checkpoints")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Ejecuta o recupera un paso del pipeline.
#'
#' @param name  identificador del paso (p.ej. "tad_reduced", "crm_intratad").
#' @param chr   cromosoma (para nombrar el fichero por cromosoma).
#' @param expr  expresión a evaluar si no hay checkpoint (entre llaves {}).
#' @param force si TRUE, recalcula aunque exista el checkpoint.
#' @return el resultado del paso (calculado o cargado).
cache_step <- function(name, chr, force = FALSE, expr) {
  f <- file.path(.checkpoint_dir(), paste0(name, "_", chr, ".rds"))
  if (!force && file.exists(f)) {
    .msg("[checkpoint] cargando '", name, "' (", chr, ") desde disco — se omite el cálculo.")
    return(readRDS(f))
  }
  .msg("[checkpoint] calculando '", name, "' (", chr, ")...")
  result <- eval.parent(substitute(expr))
  saveRDS(result, f)
  .msg("[checkpoint] guardado: ", f)
  result
}

#' Borra checkpoints de un cromosoma (para forzar recálculo limpio).
#' @param chr cromosoma. @param names vector de pasos, o NULL = todos.
clear_checkpoints <- function(chr, names = NULL) {
  d <- .checkpoint_dir()
  if (is.null(names)) {
    fs <- list.files(d, pattern = paste0("_", chr, "\\.rds$"), full.names = TRUE)
  } else {
    fs <- file.path(d, paste0(names, "_", chr, ".rds"))
    fs <- fs[file.exists(fs)]
  }
  if (length(fs) > 0) { file.remove(fs); .msg("Borrados ", length(fs), " checkpoints de ", chr, ".") }
  else .msg("No había checkpoints que borrar para ", chr, ".")
  invisible(fs)
}
