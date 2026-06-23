################################################################################
##  REDUCCIÓN DEFINITIVA DE REDUNDANCIA DE TADs
##  --------------------------------------------------------------------------
##  Criterio fijado tras la fase exploratoria (validado en chr8):
##    - Aristas: reciprocal_overlap >= 0.80  (recíproco PURO; sin Simpson,
##      porque en TADs Simpson capta jerarquía sub-TAD/TAD, no redundancia).
##    - Clustering: componentes conexas (igraph::components), igual que el
##      pipeline de CRMs. Validado: con_nucleo=1, sin chaining espurio.
##
##  Filosofía: NO crea TADs biológicos nuevos. Colapsa anotaciones
##  técnicamente redundantes de múltiples bases de datos / biosamples,
##  preservando TRAZABILIDAD COMPLETA.
##
##  Salida por cada TAD consenso:
##    - AMBAS coordenadas: consensus (max start/min end) y union (min start/max end)
##    - has_core: TRUE si existe núcleo común (consensus válido)
##    - recommended_coord: "consensus" si has_core, si no "union"
##      (pensado para escalar a todos los cromosomas: si un cluster no tiene
##       núcleo común, se usa la unión)
##    - metadatos de trazabilidad CONCATENADOS (valores únicos, ;-separados)
##    - lista de IDs originales del cluster
##  Y ADEMÁS una tabla de mapeo ID_original -> cluster_id con los metadatos
##  originales intactos (trazabilidad pura y completa).
##
##  Dependencias: data.table, GenomicRanges, igraph
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(igraph)
})

.msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

## ----------------------------------------------------------------------------
## Métricas de pares (recíproco). Copia mínima autocontenida.
## ----------------------------------------------------------------------------
.compute_pairs_recip <- function(dt) {
  dt <- as.data.table(dt)
  gr <- GRanges(seqnames = dt$chr,
                ranges   = IRanges(start = dt$start, end = dt$end),
                ID       = dt$ID)
  hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
  qh <- queryHits(hits); sh <- subjectHits(hits)
  keep <- qh < sh
  qh <- qh[keep]; sh <- sh[keep]
  if (length(qh) == 0L)
    return(data.table(id_i = character(), id_j = character(),
                      reciprocal_overlap = numeric()))
  st <- start(gr); en <- end(gr); len <- width(gr)
  ov <- pmin(en[qh], en[sh]) - pmax(st[qh], st[sh]) + 1L
  ov[ov < 0L] <- 0L
  data.table(
    id_i = mcols(gr)$ID[qh],
    id_j = mcols(gr)$ID[sh],
    reciprocal_overlap = pmin(ov / len[qh], ov / len[sh])
  )
}

## ----------------------------------------------------------------------------
## Helper: concatena valores únicos de una columna, separando por ";" y
## tratando los multivaluados que ya vienen con ";" (los desdobla antes).
## Ignora el placeholder "-" salvo que sea el único valor.
## ----------------------------------------------------------------------------
.collapse_unique <- function(x, sep = ";") {
  # Desdobla valores que ya vienen concatenados con ";"
  parts <- unlist(strsplit(as.character(x), sep, fixed = TRUE))
  parts <- trimws(parts)
  parts <- parts[parts != "" & !is.na(parts)]
  # Quita el placeholder "-" si hay valores reales
  real <- parts[parts != "-"]
  if (length(real) > 0L) parts <- real
  uq <- unique(parts)
  if (length(uq) == 0L) return(NA_character_)
  paste(uq, collapse = sep)
}

## ============================================================================
## FUNCIÓN PRINCIPAL: reduce_redundancy_tads()
## ============================================================================

#' Reduce la redundancia de TADs por componentes conexas a recíproco >= thr.
#'
#' @param tad_unique data.table de TADs colapsados por ID. Debe contener al
#'        menos chr,start,end,ID. Las columnas de metadatos presentes en
#'        'meta_cols' se concatenan como trazabilidad.
#' @param thr_recip  umbral de reciprocal_overlap para crear arista (def. 0.80).
#' @param meta_cols  columnas de trazabilidad a concatenar (valores únicos).
#'        Por defecto las de BioGateway/TADHS. Las ausentes se ignoran.
#' @param sum_cols   columnas numéricas a SUMAR por cluster (no concatenar).
#'        Por defecto n_rows_collapsed: el total refleja cuántas filas crudas
#'        representó el cluster. En la salida aparecen con prefijo "sum_".
#' @param id_sep     separador para la lista de IDs originales y metadatos.
#' @return lista con:
#'   tad_reduced : data.table, 1 fila por TAD consenso, con AMBAS coordenadas,
#'                 has_core, recommended_coord, metadatos concatenados,
#'                 columnas sum_* agregadas, n_entities, original_ids, longitudes.
#'   mapping     : data.table ID_original -> cluster_id + metadatos originales
#'                 (trazabilidad completa).
#'   summary     : data.table resumen de la reducción.
reduce_redundancy_tads <- function(tad_unique,
                                   thr_recip = 0.80,
                                   meta_cols = c("biosample_name",
                                                 "cell_line_CLO",
                                                 "cell_type_CL",
                                                 "anatomical_structures_UBERON",
                                                 "BTO",
                                                 "biological_sample_type"),
                                   sum_cols  = c("n_rows_collapsed"),
                                   id_sep = ";") {

  dt <- as.data.table(copy(tad_unique))
  stopifnot(all(c("chr", "start", "end", "ID") %in% names(dt)))
  dt[, start := as.integer(start)]
  dt[, end   := as.integer(end)]
  dt[, length := end - start + 1L]

  # Sólo metadatos que existan realmente en la tabla
  meta_present <- intersect(meta_cols, names(dt))
  meta_missing <- setdiff(meta_cols, names(dt))
  if (length(meta_missing) > 0L)
    .msg("Aviso: columnas de metadatos ausentes (se ignoran): ",
         paste(meta_missing, collapse = ", "))

  # Columnas numéricas a SUMAR por cluster (p.ej. n_rows_collapsed)
  sum_present <- intersect(sum_cols, names(dt))
  sum_missing <- setdiff(sum_cols, names(dt))
  if (length(sum_missing) > 0L)
    .msg("Aviso: columnas de suma ausentes (se ignoran): ",
         paste(sum_missing, collapse = ", "))
  # Aseguramos tipo numérico en las columnas a sumar
  for (cc in sum_present) dt[, (cc) := as.numeric(get(cc))]

  n_input <- nrow(dt)
  .msg("Reducción TADs: ", n_input, " TADs de entrada. Umbral recíproco = ",
       thr_recip, ".")

  ## --- 1) Pares + aristas ---------------------------------------------------
  pairs <- .compute_pairs_recip(dt)
  edges <- if (nrow(pairs) > 0L) pairs[reciprocal_overlap >= thr_recip] else pairs
  .msg("Pares solapados: ", nrow(pairs), " | aristas (>=", thr_recip, "): ",
       nrow(edges), ".")

  ## --- 2) Grafo + componentes conexas --------------------------------------
  edf <- if (nrow(edges) > 0L)
    as.data.frame(edges[, .(id_i, id_j)], stringsAsFactors = FALSE) else
    data.frame(id_i = character(), id_j = character(), stringsAsFactors = FALSE)
  g <- graph_from_data_frame(
    d = edf, directed = FALSE,
    vertices = data.frame(name = dt$ID, stringsAsFactors = FALSE)
  )
  comp <- igraph::components(g)
  membership <- data.table(ID = names(comp$membership),
                           cluster_id = as.integer(comp$membership))

  n_clusters <- comp$no
  .msg("Clusters (entidades no redundantes): ", n_clusters,
       " (reducción ", round(100 * (1 - n_clusters / n_input), 1), "%).")

  ## --- 3) Anotar y construir consenso/union por cluster --------------------
  d <- merge(dt, membership, by = "ID")

  # Coordenadas: consensus (max start/min end) y union (min start/max end)
  coords <- d[, .(
    chr             = chr[1L],
    consensus_start = max(start),
    consensus_end   = min(end),
    union_start     = min(start),
    union_end       = max(end),
    n_entities      = .N
  ), by = cluster_id]

  # has_core y coordenada recomendada
  coords[, has_core := consensus_end >= consensus_start]
  coords[, consensus_length := pmax(consensus_end - consensus_start + 1L, 0L)]
  coords[, union_length     := union_end - union_start + 1L]
  coords[, boundary_variability := union_length - consensus_length]
  coords[, recommended_coord := ifelse(has_core, "consensus", "union")]
  # Coordenada efectiva lista para usar (consensus si hay núcleo; si no, union)
  coords[, repr_start := ifelse(has_core, consensus_start, union_start)]
  coords[, repr_end   := ifelse(has_core, consensus_end,   union_end)]

  ## --- 4) representative_id: longitud más cercana a la mediana del cluster --
  rep_dt <- d[, {
    med <- median(length)
    .(representative_id = ID[which.min(abs(length - med))])
  }, by = cluster_id]

  ## --- 5) Lista de IDs originales por cluster -------------------------------
  ids_dt <- d[, .(original_ids = paste(ID, collapse = id_sep),
                  n_original_ids = .N), by = cluster_id]

  ## --- 6) Metadatos de trazabilidad concatenados (valores únicos) ----------
  if (length(meta_present) > 0L) {
    meta_dt <- d[, lapply(.SD, .collapse_unique),
                 by = cluster_id, .SDcols = meta_present]
  } else {
    meta_dt <- unique(d[, .(cluster_id)])
  }

  ## --- 6b) Columnas numéricas SUMADAS por cluster (p.ej. n_rows_collapsed) --
  ## Se suma sobre las filas crudas que cada cluster representa, de modo que
  ## el total refleja cuántas anotaciones originales había antes del colapso.
  if (length(sum_present) > 0L) {
    sum_dt <- d[, lapply(.SD, sum, na.rm = TRUE),
                by = cluster_id, .SDcols = sum_present]
    # Prefijo para dejar claro que es una suma agregada
    setnames(sum_dt, sum_present, paste0("sum_", sum_present))
  } else {
    sum_dt <- NULL
  }

  ## --- 7) Ensamblar tad_reduced --------------------------------------------
  parts <- list(coords, rep_dt, ids_dt, meta_dt)
  if (!is.null(sum_dt)) parts <- c(parts, list(sum_dt))
  tad_reduced <- Reduce(function(a, b) merge(a, b, by = "cluster_id"), parts)

  # Orden de columnas: identidad, coordenadas, flags, métricas, trazas, metadatos
  front <- c("cluster_id", "representative_id", "chr",
             "repr_start", "repr_end", "recommended_coord", "has_core",
             "consensus_start", "consensus_end",
             "union_start", "union_end",
             "consensus_length", "union_length", "boundary_variability",
             "n_entities", "n_original_ids",
             paste0("sum_", sum_present), "original_ids")
  front <- intersect(front, names(tad_reduced))
  setcolorder(tad_reduced, c(front, setdiff(names(tad_reduced), front)))
  setkey(tad_reduced, cluster_id)

  ## --- 8) Tabla de mapeo (trazabilidad pura: original -> cluster) ----------
  map_cols <- c("ID", "cluster_id", "chr", "start", "end", "length",
                meta_present, sum_present)
  map_cols <- intersect(map_cols, names(d))
  mapping <- d[, ..map_cols]
  setnames(mapping, "ID", "original_id")
  setkey(mapping, cluster_id, original_id)

  ## --- 9) Resumen -----------------------------------------------------------
  summary_dt <- data.table(
    chr              = dt$chr[1L],
    thr_recip        = thr_recip,
    n_input          = n_input,
    n_clusters       = n_clusters,
    n_singletons     = coords[n_entities == 1L, .N],
    n_merged         = coords[n_entities > 1L, .N],
    max_cluster_size = max(coords$n_entities),
    n_without_core   = coords[has_core == FALSE, .N],   # clusters que usan union
    reduction_ratio  = 1 - n_clusters / n_input
  )
  .msg("Clusters sin núcleo común (usan union): ", summary_dt$n_without_core, ".")

  list(
    tad_reduced = tad_reduced[],
    mapping     = mapping[],
    summary     = summary_dt[]
  )
}

## ----------------------------------------------------------------------------
## Guardado de resultados (TSV.GZ, coherente con tu estilo fwrite).
## ----------------------------------------------------------------------------

#' Guarda los resultados de la reducción de TADs en disco.
#'
#' @param reduction salida de reduce_redundancy_tads().
#' @param chr       etiqueta de cromosoma para el nombre de archivo.
#' @param output_dir carpeta de salida.
save_tad_reduction <- function(reduction, chr,
                               output_dir = "results/intermediate") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  f1 <- file.path(output_dir, paste0("tad_reduced_", chr, ".tsv.gz"))
  f2 <- file.path(output_dir, paste0("tad_mapping_", chr, ".tsv.gz"))
  f3 <- file.path(output_dir, paste0("tad_reduction_summary_", chr, ".tsv"))
  fwrite(reduction$tad_reduced, f1, sep = "\t")
  fwrite(reduction$mapping,     f2, sep = "\t")
  fwrite(reduction$summary,     f3, sep = "\t")
  .msg("Guardado: ", f1)
  .msg("Guardado: ", f2)
  .msg("Guardado: ", f3)
  invisible(c(f1, f2, f3))
}

## ============================================================================
## EJEMPLO DE USO
## ----------------------------------------------------------------------------
## source("tad_reduce_final.R")
##
## red <- reduce_redundancy_tads(tad_unique, thr_recip = 0.80)
##
## print(red$summary)
## head(red$tad_reduced)     # 1 fila por TAD consenso, ambas coordenadas
## head(red$mapping)         # trazabilidad: original_id -> cluster_id
##
## # Guardar:
## save_tad_reduction(red, chr = "chr8")
##
## # --- Para escalar a todos los cromosomas ---------------------------------
## # Como cada archivo es de un cromosoma, basta iterar:
## # for (chr in chroms) {
## #   tad_unique <- collapse_redundant_annotations(load_tad_chr(chr), "ID")
## #   red <- reduce_redundancy_tads(tad_unique, thr_recip = 0.80)
## #   save_tad_reduction(red, chr = chr)
## # }
## # Revisa summary$n_without_core por cromosoma: si es alto en alguno,
## # indica clusters sin núcleo común (ahí se usa la coordenada union).
## ============================================================================
