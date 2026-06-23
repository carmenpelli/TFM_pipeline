################################################################################
##  p1_load_collapse.R — Carga y colapso por identificador
##  --------------------------------------------------------------------------
##  Este script constituye la etapa inicial del pipeline. Carga anotaciones de
##  TADs o CRMs por cromosoma y las consolida en un único registro por
##  identificador.
##
##  Para cada identificador se conserva el cromosoma observado, la coordenada
##  inicial mínima y la coordenada final máxima. Las filas con coordenadas
##  ausentes o inconsistentes se descartan.
##
##  Dependencias:
##    data.table.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

#' Colapsa un data.table de anotaciones a un registro por ID.
#'
#' @param dt        data.table de entrada. Debe tener chr, start, end, ID.
#' @param keep_meta si TRUE, conserva metadatos concatenando valores únicos por
#'                  ID (más lento). Si FALSE (por defecto), esquema mínimo chr/start/end/ID.
#' @param id_col    nombre de la columna identificador (por defecto "ID").
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
    # Las columnas de metadatos suelen ser categóricas (biosample,
    # ontologías, tipo de muestra); se concatenan valores únicos.
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

#' Carga un fichero de entrada de cromosoma y lo colapsa por ID.
#'
#' @param path      ruta al fichero (tsv/tsv.gz/bed...). Se lee con fread.
#' @param keep_meta pasar a collapse_by_id.
#' @param id_col    columna identificador.
#' @param ...       argumentos extra para fread (sep, header...).
#' @return data.table colapsada por ID.
load_and_collapse <- function(path, keep_meta = FALSE, id_col = "ID", ...) {
  .msg("Cargando ", path, "...")
  raw <- fread(path, ...)
  .msg("  filas de entrada: ", nrow(raw))
  out <- collapse_by_id(raw, keep_meta = keep_meta, id_col = id_col)
  .msg("  tras colapso por ID: ", nrow(out), " entidades únicas.")
  out
}

## ============================================================================
## Ejemplo de uso
## ----------------------------------------------------------------------------
## source("p1_load_collapse.R")
## tad_unique <- load_and_collapse("data/TADs/chr8.tsv.gz")
## enh_unique <- load_and_collapse("data/enhancers/chr8.tsv.gz")
## # Con metadatos, conservando biosample/cell_line para trazabilidad:
## # tad_unique <- load_and_collapse("data/TADs/chr8.tsv.gz", keep_meta = TRUE)
## ============================================================================
