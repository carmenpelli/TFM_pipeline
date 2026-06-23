################################################################################
##  main.R — ORQUESTADOR DEL PIPELINE TFM
##  ==========================================================================
##  Reducción de redundancia de TADs y CRMs (BioGateway) + detección y
##  validación de regiones reguladoras densas / candidatos estructurales a
##  super-enhancer (Línea B), con coherencia topológica intra-TAD.
##
##  Ejecuta el pipeline COMPLETO para un cromosoma. Estructurado para escalar
##  a todos los cromosomas (ver bucle al final, comentado).
##
##  IMPORTANTE (biología): las regiones densas son CANDIDATOS ESTRUCTURALES,
##  no super-enhancers confirmados (sin H3K27ac/MED1/BRD4). SEdb se usa como
##  validación cruzada posicional (recuperación de SE anotados), no como
##  verdad experimental independiente.
##
##  ORDEN DEL PIPELINE (fiel al diseño):
##    P1  colapso por ID            (TADs y CRMs por separado)
##    P2  reducción redundancia TAD (recíproco 0.80 + componentes conexas)
##    P6  asignación CRM -> sub-TAD más específico (criterio bibliográfico)
##    P7  reducción redundancia CRM INTRA-TAD (= comparación CRM-CRM)
##    LB  detección DRRs (apilamiento + filtro 100 kb + 4 clases)
##    LB  validación SEdb (posicional + control de tamaño + génica)
##    LB  anotación topológica de DRRs (sub-TAD)
##
##  USO:
##    Rscript main.R chr8
##    # o dentro de R:  source("main.R") tras fijar CHR abajo.
################################################################################

## ---------------------------------------------------------------------------
## 0. CONFIGURACIÓN
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

# Cromosoma a procesar (argumento de línea de comandos o valor por defecto)
# Cromosoma a procesar.
# - Desde terminal:  Rscript main.R chr8 [force]
# - Desde la consola de R:  define .args antes de source("main.R"), p.ej.
#       .args <- "chr8"          ; source("main.R")
#       .args <- c("chr8","force"); source("main.R")
if (!exists(".args")) .args <- commandArgs(trailingOnly = TRUE)
CHR <- if (length(.args) >= 1) .args[1] else "chr8"

# Recálculo forzado: si TRUE, ignora checkpoints y recalcula todo.
# Puedes pasarlo como 2º argumento: Rscript main.R chr8 force
FORCE_RECOMPUTE <- length(.args) >= 2 && tolower(.args[2]) %in% c("force","true","1")

# Rutas (relativas a la raíz del proyecto TFM_pipeline; ajusta a tus datos)
PATHS <- list(
  data_dir   = "data",                  # entradas
  tad_dir    = "data/tad_per_chr",      # TADs crudos por cromosoma: chrN.tsv.gz
  enh_dir    = "data/enh_per_chr",      # CRMs crudos por cromosoma: chrN.tsv.gz
  results    = "results",               # salidas
  R_dir      = "R",                     # scripts del pipeline
  sedb_bed   = "data/SEdb_Human_SE.bed",   # global
  enh2gene   = "data/enh2gene.tsv.gz",     # global (se filtra por chr)
  human_genes= "data/human_genes.tsv"      # global
)
# Carpeta de checkpoints (ejecución recuperable)
options(tfm.checkpoint_dir = file.path(PATHS$results, "checkpoints"))

# Parámetros del pipeline (todos justificados en la memoria)
PARAMS <- list(
  tad_thr_recip   = 0.80,    # reducción TAD: recíproco para "mismo cuerpo de TAD"
  crm_thr_recip   = 0.50,    # reducción CRM: recíproco (estándar de facto)
  crm_jaccard_thr = 0.70,    # criterio compuesto CRM
  crm_simpson_thr = 0.99,    # criterio compuesto CRM
  drr_max_gap     = 12500L,  # (solo criterio proximidad, comparación)
  drr_max_size_bp = 100000L, # filtro de tamaño DRR (dominios SE en literatura)
  drr_min_members = 2L
)

dir.create(file.path(PATHS$results, "intermediate"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PATHS$results, "lineB", CHR),    recursive = TRUE, showWarnings = FALSE)

.msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
.msg("===== PIPELINE TFM | cromosoma: ", CHR, " =====")

## ---------------------------------------------------------------------------
## 1. CARGA DE FUNCIONES DEL PIPELINE
## ---------------------------------------------------------------------------
# El ORDEN importa: las dependencias primero.
pipeline_scripts <- c(
  "checkpoint.R",              # utilidad de ejecución recuperable
  "p1_load_collapse.R",        # P1: carga + colapso por ID
  "tad_reduce_final.R",        # reducción TAD
  "crm_explore.R",             # compute_crm_edges_chunked, etc.
  "crm_antichaining_fast.R",   # anti_chaining_fast
  "crm_reduce_final.R",        # reduce_redundancy_crms
  "crm_reduce_intratad.R",     # reduce_redundancy_crms_intratad
  "step6_assign_crm_subtad.R", # annotate_crm_subtad (P6)
  "drr_compare_criteria.R",    # detect_drr_stacking, compare_drr_criteria
  "drr_sedb_validation.R",     # run_sedb_validation
  "drr_gene_analysis.R",       # human_genes, promotores, normalización
  "drr_size_control_cases.R",  # control por tamaño + casos destacados
  "drr_lineB_stats_figs.R",    # modelo estadístico + figuras
  "drr_propagate_subtad.R"     # anotación topológica DRR
)
for (s in pipeline_scripts) source(file.path(PATHS$R_dir, s))
.msg("Funciones cargadas (", length(pipeline_scripts), " scripts).")

## ---------------------------------------------------------------------------
## 2. P1 — CARGA Y COLAPSO POR ID (desde ficheros CRUDOS)
## ---------------------------------------------------------------------------
## Parte de los ficheros crudos del cromosoma (IDs repetidos en varias filas)
## y los colapsa a un registro por ID. keep_meta=TRUE conserva biosample/etc.
.msg("P1: carga y colapso por ID de ", CHR, "...")

tad_unique <- cache_step("tad_unique", CHR, force = FORCE_RECOMPUTE,
  load_and_collapse(
    file.path(PATHS$tad_dir, paste0(CHR, ".tsv.gz")),
    keep_meta = TRUE))           # conserva metadatos de biosample/ontología (trazabilidad)

enh_unique <- cache_step("enh_unique", CHR, force = FORCE_RECOMPUTE,
  load_and_collapse(
    file.path(PATHS$enh_dir, paste0(CHR, ".tsv.gz")),
    keep_meta = TRUE))           # idem para CRMs (más lento sobre ~1,3M filas)

# enh2gene: relación CRM original -> gen (para la validación génica).
# El fichero NO es específico de cromosoma; se filtra por los crm_ID de este chr
# (los CRMs colapsados de enh_unique tienen su ID y su chr).
enh2gene_path <- PATHS$enh2gene
if (file.exists(enh2gene_path)) {
  enh2gene_all <- fread(enh2gene_path)
  # Filtrar a los CRMs presentes en este cromosoma (enh_unique$ID)
  enh2gene <- enh2gene_all[crm_ID %in% enh_unique$ID]
  .msg("enh2gene: ", nrow(enh2gene_all), " pares totales -> ",
       nrow(enh2gene), " en ", CHR, ".")
} else {
  enh2gene <- NULL
  .msg("  AVISO: no se encontró ", enh2gene_path,
       " — la validación génica se omitirá.")
}

.msg("TADs únicos: ", nrow(tad_unique), " | CRMs únicos: ", nrow(enh_unique))

## ---------------------------------------------------------------------------
## 3. P2 — REDUCCIÓN DE REDUNDANCIA DE TADs
## ---------------------------------------------------------------------------
.msg("P2: reducción de TADs...")
red_tad <- cache_step("red_tad", CHR, force = FORCE_RECOMPUTE,
  reduce_redundancy_tads(tad_unique, thr_recip = PARAMS$tad_thr_recip))
save_tad_reduction(red_tad, chr = CHR)
.msg("  TADs: ", nrow(tad_unique), " -> ", nrow(red_tad$tad_reduced))

## ---------------------------------------------------------------------------
## 4. P6 — ASIGNACIÓN CRM -> SUB-TAD MÁS ESPECÍFICO
##    (sobre los CRMs colapsados por ID; criterio: TAD más pequeño que contiene)
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
## 5. P7 — REDUCCIÓN DE REDUNDANCIA DE CRMs **INTRA-TAD**
##    (la comparación CRM-CRM y la reducción son el mismo paso, dentro del TAD)
## ---------------------------------------------------------------------------
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
## 6. LÍNEA B — DETECCIÓN DE DRRs (apilamiento + filtro de tamaño + 4 clases)
## ---------------------------------------------------------------------------
.msg("LB: detección de DRRs...")
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
## 7. LÍNEA B — VALIDACIÓN SEdb (posicional; génica si hay enh2gene)
## ---------------------------------------------------------------------------
.msg("LB: validación SEdb...")
val <- cache_step("sedb_val", CHR, force = FORCE_RECOMPUTE,
  run_sedb_validation(
    drr = drr, crm_reduced = red_crm$crm_reduced,
    mapping = red_crm$mapping,    # ahora la reducción intra-TAD SÍ devuelve mapping
    enh2gene = enh2gene,          # NULL si no se encontró el fichero (génica se omite)
    bed_path = PATHS$sedb_bed,
    chr = CHR,                       # <-- cromosoma real (antes filtraba SIEMPRE a chr8)
    crm_id_col = "crm_ID",            # ajusta al nombre real en tu enh2gene
    gene_col   = "hgnc_symbol_target_genes"))
fwrite(val$positional$by_class, file.path(PATHS$results, "lineB", CHR,
       paste0("sedb_by_class_", CHR, ".tsv")), sep = "\t")
if (!is.null(val$genic))
  fwrite(val$genic$summary, file.path(PATHS$results, "lineB", CHR,
         paste0("sedb_gene_validation_", CHR, ".tsv")), sep = "\t")

## ---------------------------------------------------------------------------
## 8. LÍNEA B — CONTROL POR TAMAÑO + MODELO ESTADÍSTICO + FIGURAS
## ---------------------------------------------------------------------------
.msg("LB: control por tamaño y modelo...")
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
## 9. LÍNEA B — ANOTACIÓN TOPOLÓGICA DE DRRs (sub-TAD)
## ---------------------------------------------------------------------------
.msg("LB: anotación topológica DRR -> sub-TAD...")
prop <- propagate_subtad_to_drr(drr = drr, crm_reduced = red_crm$crm_reduced,
                                crm_tad = a6$crm_tad)
fwrite(prop$by_class, file.path(PATHS$results, "lineB", CHR,
       paste0("drr_subtad_byclass_", CHR, ".tsv")), sep = "\t")

## ---------------------------------------------------------------------------
## 10. GUARDAR OBJETOS CLAVE Y RESUMEN
## ---------------------------------------------------------------------------
saveRDS(list(red_tad = red_tad, a6 = a6, red_crm = red_crm, cmp = cmp,
             val = val, sm = sm, mod = mod, prop = prop, PARAMS = PARAMS),
        file.path(PATHS$results, "lineB", CHR, paste0("pipeline_objects_", CHR, ".rds")))

.msg("===== PIPELINE ", CHR, " COMPLETADO =====")
if (!is.null(mod$coef_table)) print(mod$coef_table) else
  .msg("Modelo degenerado en ", CHR, " (señal SEdb insuficiente); ",
       "se analizará a nivel global. Diagnóstico: n_SE>0 = ",
       mod$diagnostics$n_nonzero, " de ", mod$diagnostics$n_obs, " DRRs.")

## ===========================================================================
## ESCALADO A TODOS LOS CROMOSOMAS (descomentar para usar)
## ---------------------------------------------------------------------------
## chrs <- paste0("chr", c(1:22, "X"))
## for (c in chrs) system(paste("Rscript main.R", c))
## # o envolver los pasos 2-10 en una función run_chr(CHR) y aplicar lapply.
## ===========================================================================
