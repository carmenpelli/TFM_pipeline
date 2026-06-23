################################################################################
##  main.R — Orquestador principal del pipeline
##  ==========================================================================
##  Este script ejecuta el pipeline completo para un cromosoma, incluyendo:
##
##    - Reducción de redundancia de TADs y CRMs.
##    - Asignación de CRMs a un marco topológico reducido.
##    - Reducción intra-TAD de CRMs.
##    - Detección de regiones reguladoras densas (DRRs).
##    - Comparación de las DRRs con anotaciones de super-enhancers de SEdb.
##    - Propagación de la información topológica a las DRRs.
##
##  Las DRRs se interpretan como regiones estructuralmente densas definidas por
##  el pipeline, no como super-enhancers funcionalmente validados. La comparación
##  con SEdb se utiliza como análisis de concordancia posicional y génica.
##
##  Orden general del pipeline:
##    P1   carga y colapso por identificador de TADs y CRMs
##    P2   reducción de redundancia de TADs
##    P6   asignación de CRMs al TAD reducido más específico
##    P7   reducción intra-TAD de CRMs
##    DRR  detección y clasificación de regiones reguladoras densas
##    SEdb comparación posicional y génica con SEdb
##    TOP  anotación topológica de DRRs
##
##  Uso en terminal:
##    Rscript main.R chr8
##
##  Uso desde R:
##    .args <- "chr8"
##    source("main.R")
################################################################################

## ---------------------------------------------------------------------------
## 0. Configuración general
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

# Cromosoma a procesar.
# - Desde terminal: Rscript main.R chr8 [force]
# - Desde R: definir .args antes de source("main.R"), por ejemplo:
#       .args <- "chr8"
#       .args <- c("chr8", "force")
if (!exists(".args")) .args <- commandArgs(trailingOnly = TRUE)
CHR <- if (length(.args) >= 1) .args[1] else "chr8"

# Recálculo forzado: si TRUE, ignora checkpoints y recalcula todas las etapas.
FORCE_RECOMPUTE <- length(.args) >= 2 && tolower(.args[2]) %in% c("force","true","1")

# Rutas relativas a la raíz del proyecto.
PATHS <- list(
  data_dir   = "data",
  tad_dir    = "data/tad_per_chr",
  enh_dir    = "data/enh_per_chr",
  results    = "results",
  R_dir      = "R",
  sedb_bed   = "data/SEdb_Human_SE.bed",
  enh2gene   = "data/enh2gene.tsv.gz",
  human_genes= "data/human_genes.tsv"
)

# Directorio de checkpoints para permitir la reanudación del análisis.
options(tfm.checkpoint_dir = file.path(PATHS$results, "checkpoints"))

# Parámetros principales del pipeline.
PARAMS <- list(
  tad_thr_recip   = 0.80,
  crm_thr_recip   = 0.50,
  crm_jaccard_thr = 0.70,
  crm_simpson_thr = 0.99,
  drr_max_gap     = 12500L,
  drr_max_size_bp = 100000L,
  drr_min_members = 2L
)

dir.create(file.path(PATHS$results, "intermediate"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PATHS$results, "lineB", CHR),    recursive = TRUE, showWarnings = FALSE)

.msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
.msg("===== PIPELINE TFM | cromosoma: ", CHR, " =====")

## ---------------------------------------------------------------------------
## 1. Carga de funciones del pipeline
## ---------------------------------------------------------------------------

# Los scripts se cargan en orden para respetar sus dependencias internas.
pipeline_scripts <- c(
  "checkpoint.R",
  "p1_load_collapse.R",
  "tad_reduce_final.R",
  "crm_explore.R",
  "crm_antichaining_fast.R",
  "crm_reduce_final.R",
  "crm_reduce_intratad.R",
  "step6_assign_crm_subtad.R",
  "drr_compare_criteria.R",
  "drr_sedb_validation.R",
  "drr_gene_analysis.R",
  "drr_size_control_cases.R",
  "drr_lineB_stats_figs.R",
  "drr_propagate_subtad.R"
)

for (s in pipeline_scripts) source(file.path(PATHS$R_dir, s))
.msg("Funciones cargadas (", length(pipeline_scripts), " scripts).")

## ---------------------------------------------------------------------------
## 2. P1 — Carga y colapso por identificador
## ---------------------------------------------------------------------------

## Se parte de los ficheros de entrada por cromosoma y se obtiene un registro
## único por identificador, conservando metadatos cuando están disponibles.
.msg("P1: carga y colapso por ID de ", CHR, "...")

tad_unique <- cache_step("tad_unique", CHR, force = FORCE_RECOMPUTE,
  load_and_collapse(
    file.path(PATHS$tad_dir, paste0(CHR, ".tsv.gz")),
    keep_meta = TRUE))

enh_unique <- cache_step("enh_unique", CHR, force = FORCE_RECOMPUTE,
  load_and_collapse(
    file.path(PATHS$enh_dir, paste0(CHR, ".tsv.gz")),
    keep_meta = TRUE))

# enh2gene contiene relaciones globales CRM-gen; se filtra a los CRMs presentes
# en el cromosoma procesado.
enh2gene_path <- PATHS$enh2gene

if (file.exists(enh2gene_path)) {
  enh2gene_all <- fread(enh2gene_path)
  enh2gene <- enh2gene_all[crm_ID %in% enh_unique$ID]
  .msg("enh2gene: ", nrow(enh2gene_all), " pares totales -> ",
       nrow(enh2gene), " en ", CHR, ".")
} else {
  enh2gene <- NULL
  .msg("Advertencia: no se encontró ", enh2gene_path,
       "; la comparación génica se omitirá.")
}

.msg("TADs únicos: ", nrow(tad_unique), " | CRMs únicos: ", nrow(enh_unique))

## ---------------------------------------------------------------------------
## 3. P2 — Reducción de redundancia de TADs
## ---------------------------------------------------------------------------

.msg("P2: reducción de TADs...")

red_tad <- cache_step("red_tad", CHR, force = FORCE_RECOMPUTE,
  reduce_redundancy_tads(tad_unique, thr_recip = PARAMS$tad_thr_recip))

save_tad_reduction(red_tad, chr = CHR)

.msg("  TADs: ", nrow(tad_unique), " -> ", nrow(red_tad$tad_reduced))

## ---------------------------------------------------------------------------
## 4. P6 — Asignación de CRMs al TAD reducido más específico
## ---------------------------------------------------------------------------

.msg("P6: asignación CRM -> sub-TAD...")

a6 <- cache_step("a6_crm_subtad", CHR, force = FORCE_RECOMPUTE,
  annotate_crm_subtad(
    crm_reduced  = enh_unique,
    tad_reduced  = red_tad$tad_reduced,
    crm_id_col   = "ID",
    crm_start_col = "start", crm_end_col = "end",
    tad_start_col = "repr_start", tad_end_col = "repr_end"))

fwrite(a6$crm_tad, file.path(PATHS$results, "intermediate",
                             paste0("crm_subtad_", CHR, ".tsv.gz")), sep = "\t")

.msg("  CRMs asignados: ", a6$summary$n_assigned, " | TADs usados: ",
     a6$summary$n_tads_used)

## ---------------------------------------------------------------------------
## 5. P7 — Reducción intra-TAD de CRMs
## ---------------------------------------------------------------------------

## La comparación CRM-CRM se realiza dentro de cada TAD reducido, manteniendo
## la coherencia topológica del análisis.
.msg("P7: reducción de CRMs intra-TAD...")

red_crm <- cache_step("red_crm_intratad", CHR, force = FORCE_RECOMPUTE,
  reduce_redundancy_crms_intratad(
    crm_unique = enh_unique, crm_tad = a6$crm_tad,
    crm_id_in_unique = "ID", crm_id_in_tad = "crm_id",
    thr_recip   = PARAMS$crm_thr_recip,
    jaccard_thr = PARAMS$crm_jaccard_thr,
    simpson_thr = PARAMS$crm_simpson_thr))

fwrite(red_crm$crm_reduced, file.path(PATHS$results, "intermediate",
                                      paste0("crm_reduced_intratad_", CHR, ".tsv.gz")), sep = "\t")

.msg("  CRMs: ", a6$summary$n_assigned, " -> ", nrow(red_crm$crm_reduced))

## ---------------------------------------------------------------------------
## 6. Detección y clasificación de regiones reguladoras densas
## ---------------------------------------------------------------------------

.msg("DRR: detección de regiones reguladoras densas...")

cmp <- cache_step("drr_cmp", CHR, force = FORCE_RECOMPUTE,
  compare_drr_criteria(red_crm$crm_reduced, max_gap = PARAMS$drr_max_gap,
                       min_consensus_crms = PARAMS$drr_min_members))

drr <- cmp$stacking

fwrite(drr, file.path(PATHS$results, "lineB", CHR,
                      paste0("drr_", CHR, ".tsv.gz")), sep = "\t")

.msg("  DRRs candidatas: ",
     drr[candidate_class != "Extensive_overlap", .N],
     " | extensas: ", drr[candidate_class == "Extensive_overlap", .N])

## ---------------------------------------------------------------------------
## 7. Comparación con SEdb
## ---------------------------------------------------------------------------

.msg("SEdb: comparación posicional y génica...")

val <- cache_step("sedb_val", CHR, force = FORCE_RECOMPUTE,
  run_sedb_validation(
    drr = drr, crm_reduced = red_crm$crm_reduced,
    mapping = red_crm$mapping,
    enh2gene = enh2gene,
    bed_path = PATHS$sedb_bed,
    chr = CHR,
    crm_id_col = "crm_ID",
    gene_col   = "hgnc_symbol_target_genes"))

fwrite(val$positional$by_class, file.path(PATHS$results, "lineB", CHR,
       paste0("sedb_by_class_", CHR, ".tsv")), sep = "\t")

if (!is.null(val$genic))
  fwrite(val$genic$summary, file.path(PATHS$results, "lineB", CHR,
         paste0("sedb_gene_validation_", CHR, ".tsv")), sep = "\t")

## ---------------------------------------------------------------------------
## 8. Control por tamaño, modelo estadístico y figuras
## ---------------------------------------------------------------------------

.msg("Modelo: control por tamaño y comparación entre clases...")

sm  <- sedb_size_matched(drr, val$positional$overlap)
mod <- model_density_controlled(sm$per_drr, ref_class = "Extended_complex_DRR")

if (!is.null(mod$coef_table))
  fwrite(mod$coef_table, file.path(PATHS$results, "lineB", CHR,
         paste0("sedb_model_coef_", CHR, ".tsv")), sep = "\t")

fig_control_bin(sm$per_drr,
  file = file.path(PATHS$results, "lineB", CHR, paste0("figB_control_", CHR, ".png")))

fig_bins(sm$per_drr,
  file = file.path(PATHS$results, "lineB", CHR, paste0("figB_bins_", CHR, ".png")))

## ---------------------------------------------------------------------------
## 9. Anotación topológica de DRRs
## ---------------------------------------------------------------------------

.msg("Topología: propagación de sub-TAD a DRR...")

prop <- propagate_subtad_to_drr(drr = drr, crm_reduced = red_crm$crm_reduced,
                                crm_tad = a6$crm_tad)

fwrite(prop$by_class, file.path(PATHS$results, "lineB", CHR,
       paste0("drr_subtad_byclass_", CHR, ".tsv")), sep = "\t")

## ---------------------------------------------------------------------------
## 10. Exportación de objetos clave y resumen final
## ---------------------------------------------------------------------------

saveRDS(list(red_tad = red_tad, a6 = a6, red_crm = red_crm, cmp = cmp,
             val = val, sm = sm, mod = mod, prop = prop, PARAMS = PARAMS),
        file.path(PATHS$results, "lineB", CHR, paste0("pipeline_objects_", CHR, ".rds")))

.msg("===== PIPELINE ", CHR, " COMPLETADO =====")

if (!is.null(mod$coef_table)) print(mod$coef_table) else
  .msg("Modelo no ajustado en ", CHR, " por señal SEdb insuficiente; ",
       "el contraste se evaluará a nivel global. Diagnóstico: n_SE>0 = ",
       mod$diagnostics$n_nonzero, " de ", mod$diagnostics$n_obs, " DRRs.")

## ===========================================================================
## Referencia para ejecución secuencial por cromosomas
## ---------------------------------------------------------------------------
## chrs <- paste0("chr", c(1:22, "X"))
## for (c in chrs) system(paste("Rscript main.R", c))
## ===========================================================================
