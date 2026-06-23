################################################################################
##  P1 — CARGA Y COLAPSO POR ID
##  --------------------------------------------------------------------------
##  Punto de entrada del pipeline: carga los ficheros CRUDOS de un cromosoma
##  (TADs y CRMs de múltiples bases de datos, con IDs que pueden repetirse en
##  varias filas) y los colapsa a un registro por ID.
##
##  Colapso por ID: para cada ID, se toma chr (primero observado), min(start),
##  max(end). Esto une las filas que describen la MISMA entidad anotada con el
##  mismo identificador en distintas fuentes/biosamples.
##
##  Saneamiento: descarta filas con coordenadas inválidas (NA o end < start).
##
##  Dependencias: data.table
################################################################################

suppressPackageStartupMessages({
  library(data.table)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

#' Colapsa un data.table de anotaciones a un registro por ID.
#'
#' @param dt        data.table cruda. Debe tener chr, start, end, ID.
#' @param keep_meta si TRUE, conserva metadatos concatenando valores únicos por
#'                  ID (más lento). Si FALSE (def.), esquema mínimo chr/start/end/ID.
#' @param id_col    nombre de la columna identificador (def. "ID").
#' @return data.table con un registro por ID (chr, start, end, ID [, metadatos]).
collapse_by_id <- function(dt, keep_meta = FALSE, id_col = "ID") {
  dt <- as.data.table(dt)
  stopifnot(all(c("chr", "start", "end", id_col) %in% names(dt)))
  if (id_col != "ID") setnames(dt, id_col, "ID")

  dt[, start := as.integer(start)]
  dt[, end   := as.integer(end)]

  # Saneamiento: descartar coordenadas inválidas
  n0 <- nrow(dt)
  dt <- dt[!is.na(start) & !is.na(end) & end >= start]
  if (nrow(dt) < n0)
    .msg("  saneamiento: descartadas ", n0 - nrow(dt), " filas inválidas.")

  if (!keep_meta) {
    collapsed <- dt[, .(chr = chr[1L], start = min(start), end = max(end)), by = ID]
  } else {
    # Conserva metadatos concatenando valores únicos por ID (trazabilidad).
    # Todas las columnas de metadatos crudos son categóricas (biosample,
    # ontologías, tipo de muestra): concatenar valores únicos es lo correcto.
    meta_cols <- setdiff(names(dt), c("chr", "start", "end", "ID"))
    .collapse_u <- function(x, sep = ";")
      paste(sort(unique(as.character(x[!is.na(x) & x != ""]))), collapse = sep)
    collapsed <- dt[, c(
      .(chr = chr[1L], start = min(start), end = max(end),
        n_rows_collapsed = .N),               # nº de filas crudas que se colapsaron
      lapply(.SD, .collapse_u)
    ), by = ID, .SDcols = meta_cols]
  }

  # Longitud recalculada de las coordenadas colapsadas (útil aguas abajo)
  collapsed[, length := end - start + 1L]

  setcolorder(collapsed, c("chr", "start", "end", "ID"))
  setkey(collapsed, chr, start, end)
  collapsed[]
}

#' Carga un fichero crudo de cromosoma y lo colapsa por ID.
#'
#' @param path      ruta al fichero (tsv/tsv.gz/bed...). Se lee con fread.
#' @param keep_meta pasar a collapse_by_id.
#' @param id_col    columna identificador.
#' @param ...       argumentos extra para fread (sep, header...).
#' @return data.table colapsada por ID.
load_and_collapse <- function(path, keep_meta = FALSE, id_col = "ID", ...) {
  .msg("Cargando ", path, "...")
  raw <- fread(path, ...)
  .msg("  filas crudas: ", nrow(raw))
  out <- collapse_by_id(raw, keep_meta = keep_meta, id_col = id_col)
  .msg("  tras colapso por ID: ", nrow(out), " entidades únicas.")
  out
}

## ============================================================================
## EJEMPLO DE USO
## ----------------------------------------------------------------------------
## source("p1_load_collapse.R")
## tad_unique <- load_and_collapse("data/TADs/chr8.tsv.gz")
## enh_unique <- load_and_collapse("data/enhancers/chr8.tsv.gz")
## # Con metadatos (más lento, conserva biosample/cell_line para trazabilidad):
## # tad_unique <- load_and_collapse("data/TADs/chr8.tsv.gz", keep_meta = TRUE)
## ============================================================================
