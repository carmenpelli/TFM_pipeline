################################################################################
##  COMPARACIÓN DE CRITERIOS PARA REGIONES REGULATORIAS DENSAS (DRRs)
##  --------------------------------------------------------------------------
##  Enfrenta dos formas de DEFINIR las regiones, ambas terminando en las MISMAS
##  4 clases (Simple / Compact / Dense_complex / Extended_complex) para una
##  comparación limpia:
##
##    Criterio A — PROXIMIDAD (el antiguo, 06_detect_dense_regulatory_regions.R):
##      Agrupa CRMs consecutivos con gap <= max_gap (12.5 kb).
##      [Aviso: con datos densos (99.98% solapando) tiende a producir regiones
##       gigantes; esta comparación sirve para CONFIRMARLO con números.]
##
##    Criterio B — APILAMIENTO (el nuevo, Opción A):
##      Define regiones como BLOQUES de solapamiento contiguo (GenomicRanges
##      reduce). La "densidad" surge de cuántos CRMs se apilan, no de proximidad.
##
##  Input: crm_reduced (de reduce_redundancy_crms).
##    support = n_entities (CRMs reducidos del cluster)
##    coordenada = repr_start / repr_end
##
##  IMPORTANTE (biología, del documento): son CANDIDATOS estructurales /
##  regiones reguladoras densas. NUNCA super-enhancers confirmados.
##
##  Dependencias: data.table, GenomicRanges
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

if (!exists(".msg")) .msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

## ----------------------------------------------------------------------------
## Clasificación en 4 clases (idéntica lógica al script antiguo).
## Se factoriza para aplicarla igual a ambos criterios.
## ----------------------------------------------------------------------------
.classify_drr <- function(drr) {
  med_density <- median(drr$density_support_crms_per_kb, na.rm = TRUE)
  drr[, candidate_class := fifelse(
    n_consensus_crms <= 2 & support_crms <= 5,
    "Simple_DRR",
    fifelse(
      drr_length <= 12500 & support_crms > 5,
      "Compact_DRR",
      fifelse(
        density_support_crms_per_kb >= med_density,
        "Dense_complex_DRR",
        "Extended_complex_DRR"
      )
    )
  )]
  drr[]
}

## ----------------------------------------------------------------------------
## Estadísticas comunes de una tabla de DRRs (para comparar criterios).
## ----------------------------------------------------------------------------
.drr_stats <- function(drr, label) {
  cls <- drr[, .(n = .N,
                 median_len = median(drr_length),
                 median_support = median(support_crms),
                 median_density = median(density_support_crms_per_kb)),
             by = candidate_class][order(candidate_class)]
  overall <- data.table(
    criterion        = label,
    n_drr            = nrow(drr),
    max_drr_length   = max(drr$drr_length),
    median_drr_length= median(drr$drr_length),
    n_simple         = drr[candidate_class == "Simple_DRR", .N],
    n_compact        = drr[candidate_class == "Compact_DRR", .N],
    n_dense_complex  = drr[candidate_class == "Dense_complex_DRR", .N],
    n_extended       = drr[candidate_class == "Extended_complex_DRR", .N],
    n_extensive_ovl  = drr[candidate_class == "Extensive_overlap", .N]
  )
  list(by_class = cls[], overall = overall[])
}

## ============================================================================
## CRITERIO A — PROXIMIDAD (stitching por gap <= max_gap)
## ============================================================================
detect_drr_proximity <- function(crm_reduced,
                                 max_gap = 12500L,
                                 min_consensus_crms = 2L) {
  dt <- as.data.table(copy(crm_reduced))
  dt[, s := repr_start]
  dt[, e := repr_end]
  dt[, support := n_entities]
  dt <- dt[!is.na(s) & !is.na(e)]
  setorder(dt, chr, s, e)

  # Encadenamiento por gap (por cromosoma)
  dt[, previous_end := data.table::shift(e), by = chr]
  dt[, gap_to_previous := as.numeric(s - previous_end - 1L)]
  dt[is.na(gap_to_previous), gap_to_previous := Inf]
  dt[, new_drr := gap_to_previous > max_gap]
  dt[, drr_index := cumsum(new_drr), by = chr]
  dt[, drr_internal_id := paste(chr, drr_index, sep = "__")]

  drr <- dt[, .(
    chr = chr[1L],
    drr_start = min(s),
    drr_end   = max(e),
    n_consensus_crms = .N,
    support_crms = sum(support, na.rm = TRUE)
  ), by = drr_internal_id]

  drr[, drr_length := drr_end - drr_start + 1L]
  drr[, density_support_crms_per_kb := support_crms / (drr_length / 1000)]
  drr <- drr[n_consensus_crms >= min_consensus_crms]
  drr <- .classify_drr(drr)
  drr[, drr_id := paste0("DRRp_", sprintf("%08d", .I))]
  drr[]
}

## ============================================================================
## CRITERIO B — APILAMIENTO (bloques de solapamiento contiguo)
## ============================================================================
## Define cada región como un bloque contiguo donde los CRMs se solapan/tocan.
## La densidad (support/kb) discrimina dentro de esos bloques. Como TODO se
## solapa, habrá pocos bloques grandes; lo que varía es CUÁNTOS CRMs se apilan
## en cada tramo, capturado por support_crms y la densidad.
##
## Filtro de tamaño (literatura de super-enhancers, mediana SE 8.7-19 kb,
## dominios SE hasta ~90 kb): las regiones > max_size_bp (def. 100 kb) se
## marcan como "extensive_overlap" y se EXCLUYEN de la clasificación/validación
## (posibles artefactos de la alta densidad de solapamiento), pero se CONSERVAN.
detect_drr_stacking <- function(crm_reduced,
                                min_consensus_crms = 2L,
                                max_block_gap = 0L,
                                max_size_bp = 100000L) {
  dt <- as.data.table(copy(crm_reduced))
  dt[, s := repr_start]
  dt[, e := repr_end]
  dt[, support := n_entities]
  dt <- dt[!is.na(s) & !is.na(e)]

  gr <- GRanges(dt$chr, IRanges(dt$s, dt$e))
  # reduce() fusiona rangos que se solapan o están a <= max_block_gap.
  # max_block_gap=0 => sólo fusiona los que se solapan o tocan exactamente.
  blocks <- reduce(gr, min.gapwidth = max_block_gap + 1L)

  # Asignar cada CRM a su bloque
  ov <- findOverlaps(gr, blocks)
  dt[, block_id := NA_integer_]
  dt[queryHits(ov), block_id := subjectHits(ov)]

  drr <- dt[!is.na(block_id), .(
    chr = chr[1L],
    drr_start = min(s),
    drr_end   = max(e),
    n_consensus_crms = .N,
    support_crms = sum(support, na.rm = TRUE)
  ), by = block_id]

  drr[, drr_length := drr_end - drr_start + 1L]
  drr[, density_support_crms_per_kb := support_crms / (drr_length / 1000)]
  drr <- drr[n_consensus_crms >= min_consensus_crms]

  # Marca de tamaño: candidatas (<= max_size_bp) vs extensas (> max_size_bp)
  drr[, size_class := fifelse(drr_length <= max_size_bp,
                              "candidate", "extensive_overlap")]
  n_ext <- drr[size_class == "extensive_overlap", .N]
  .msg("Apilamiento: ", nrow(drr), " regiones | extensas (>",
       max_size_bp, " pb) marcadas: ", n_ext,
       " (", round(100 * n_ext / nrow(drr), 2), "%).")

  # Sólo las candidatas (<= max_size_bp) se clasifican en las 4 clases.
  # Las extensas conservan candidate_class = "Extensive_overlap".
  drr_cand <- .classify_drr(drr[size_class == "candidate"])
  drr_ext  <- drr[size_class == "extensive_overlap"]
  if (nrow(drr_ext) > 0L) drr_ext[, candidate_class := "Extensive_overlap"]

  drr <- rbind(drr_cand, drr_ext, fill = TRUE)
  drr[, drr_id := paste0("DRRs_", sprintf("%08d", .I))]
  drr[]
}

## ============================================================================
## COMPARACIÓN DE LOS DOS CRITERIOS
## ============================================================================
compare_drr_criteria <- function(crm_reduced, max_gap = 12500L,
                                 min_consensus_crms = 2L) {
  .msg("Criterio A (proximidad, max_gap=", max_gap, ")...")
  drr_prox <- detect_drr_proximity(crm_reduced, max_gap = max_gap,
                                   min_consensus_crms = min_consensus_crms)
  sp <- .drr_stats(drr_prox, "A_proximity")

  .msg("Criterio B (apilamiento, bloques de solapamiento)...")
  drr_stack <- detect_drr_stacking(crm_reduced,
                                   min_consensus_crms = min_consensus_crms)
  ss <- .drr_stats(drr_stack, "B_stacking")

  comparison <- rbind(sp$overall, ss$overall)

  list(
    proximity        = drr_prox,
    stacking         = drr_stack,
    comparison       = comparison[],
    by_class_prox    = sp$by_class,
    by_class_stack   = ss$by_class
  )
}

## ============================================================================
## EJEMPLO DE USO
## ----------------------------------------------------------------------------
## source("drr_compare_criteria.R")
##
## cmp <- compare_drr_criteria(red_crm$crm_reduced, max_gap = 12500)
##
## print(cmp$comparison)        # tamaños y nº de regiones por criterio
## print(cmp$by_class_prox)     # clases del criterio proximidad
## print(cmp$by_class_stack)    # clases del criterio apilamiento
##
## # Clave a mirar:
## #   max_drr_length -> ¿alguna región es de megabases? (señal de chaining)
## #   median_drr_length -> ¿tamaños biológicamente plausibles? (SEs: decenas kb)
## #   distribución de clases -> ¿discrimina o todo cae en una clase?
## #
## # Luego, el criterio elegido se valida contra SEdb (drr_start/drr_end como
## # candidate_start/candidate_end en tu script de validación).
## ============================================================================
