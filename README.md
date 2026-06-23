# Pipeline TFM — Reducción de redundancia (TADs/CRMs) y regiones reguladoras densas

Pipeline en R para (1) reducir la redundancia de anotaciones de TADs y CRMs
integradas de múltiples bases de datos (BioGateway) y (2) detectar y validar
regiones reguladoras densas / candidatos estructurales a super-enhancer
(Línea B), con coherencia topológica respecto a TADs.

> **Nota biológica importante:** las regiones densas detectadas son
> **candidatos estructurales**, no super-enhancers confirmados. No hay datos
> de H3K27ac/MED1/BRD4. SEdb se usa como validación cruzada **posicional**
> (recuperación de super-enhancers anotados), no como verdad experimental
> independiente — SEdb y los CRMs comparten parcialmente fuentes.

## Estructura de carpetas

```
TFM_pipeline/
├── main.R                      # orquestador (ejecuta un cromosoma de principio a fin)
├── README.md                   # este archivo
├── R/                          # funciones del pipeline (las que main.R hace source())
│   ├── p1_load_collapse.R          # P1: carga ficheros crudos + colapso por ID
│   ├── tad_reduce_final.R          # P2: reducción redundancia TAD
│   ├── crm_explore.R               # utilidades de aristas/solapes CRM
│   ├── crm_antichaining_fast.R     # anti-chaining (R puro, rápido)
│   ├── crm_reduce_final.R          # reduce_redundancy_crms (criterio compuesto)
│   ├── crm_reduce_intratad.R       # P7: reducción CRM intra-TAD
│   ├── step6_assign_crm_subtad.R   # P6: asignación CRM -> sub-TAD específico
│   ├── drr_compare_criteria.R      # LB: detección DRRs (apilamiento + clases)
│   ├── drr_sedb_validation.R       # LB: validación SEdb posicional + génica
│   ├── drr_gene_analysis.R         # LB: human_genes, promotores, normalización
│   ├── drr_size_control_cases.R    # LB: control por tamaño + casos destacados
│   ├── drr_lineB_stats_figs.R      # LB: modelo estadístico + figuras
│   └── drr_propagate_subtad.R      # LB: anotación topológica DRR -> sub-TAD
├── R_exploratory/              # análisis exploratorio (no se ejecutan en main.R)
│   ├── tad_threshold_sweep.R       # barrido de umbrales TAD
│   ├── tad_cluster_tools.R         # inspección de clusters TAD
│   ├── crm_inspect.R               # diagnóstico chaining vs densidad
│   ├── crm_compare_strategies.R    # conexas vs anti-chaining
│   ├── crm_antichaining_cohesion.R # cohesión de clusters
│   └── crm_dense_regions.R         # 1er criterio DRR por proximidad (descartado)
├── data/                       # entradas (créala con tus ficheros)
│   ├── tad_per_chr/<chr>.tsv.gz      # TADs CRUDOS por cromosoma (el P1 los colapsa)
│   ├── enh_per_chr/<chr>.tsv.gz      # CRMs CRUDOS por cromosoma (el P1 los colapsa)
│   ├── enh2gene.tsv.gz              # relación CRM->gen GLOBAL (se filtra por chr)
│   ├── SEdb_Human_SE.bed            # super-enhancers SEdb (global)
│   └── human_genes.tsv             # genes: chr,start,end,symbol,strand (global)
└── results/                    # se crea automáticamente al ejecutar
    ├── intermediate/           # reducciones y asignaciones por cromosoma
    └── lineB/<chr>/            # DRRs, validación SEdb, figuras, objetos .rds
```

## Cómo ejecutar

Un cromosoma de principio a fin:

```bash
Rscript main.R chr8
```

### Ejecución recuperable (checkpoints)

Cada paso costoso (colapso, reducción TAD, P6, reducción CRM intra-TAD,
detección DRRs, validación SEdb) guarda su resultado en
`results/checkpoints/<paso>_<chr>.rds`. En una **segunda ejecución**, esos
pasos se **cargan desde disco** en lugar de recalcularse — así el pipeline
retoma donde se quedó si falló, o reusa lo ya hecho.

Forzar recálculo completo (ignorar checkpoints):

```bash
Rscript main.R chr8 force
```

Borrar checkpoints de un cromosoma desde R:

```r
source("R/checkpoint.R")
options(tfm.checkpoint_dir = "results/checkpoints")
clear_checkpoints("chr8")                       # todos
clear_checkpoints("chr8", "red_crm_intratad")   # solo uno
```

Además, cada paso escribe sus tablas finales en `results/intermediate/` y
`results/lineB/<chr>/`, y al final un `.rds` con todos los objetos del cromosoma.

### Escalado a todos los cromosomas

Para procesar los 23 cromosomas en paralelo (2 a la vez, según CPUs):

```bash
Rscript run_all_chromosomes.R              # todos, reusando checkpoints
Rscript run_all_chromosomes.R force        # ignorando checkpoints
Rscript run_all_chromosomes.R "" chr1,chr2 # solo algunos
```

Cada cromosoma corre en su propio proceso R aislado (`Rscript main.R chrN`),
con su log en `results/logs/<chr>.log`. Los cromosomas se ordenan intercalando
grande+pequeño para equilibrar la memoria de cada tanda. Política **fail-stop**:
si un cromosoma falla, se deja terminar su tanda y no se lanzan más (para
investigar). Reanudar tras arreglar el problema reusa los checkpoints de los
cromosomas ya completados.

## Orden del pipeline (y por qué)

| Paso | Qué hace | Decisión clave |
|------|----------|----------------|
| **P1** | Carga ficheros crudos + colapso por ID (TADs y CRMs por separado) | Agrupa filas con mismo ID: min(start), max(end). Sanea coordenadas inválidas |
| **P2** | Reducción redundancia TAD | Recíproco ≥ 0.80 + componentes conexas. 0.80 = "mismo cuerpo de TAD" (bibliografía) |
| **P6** | Asignación CRM → sub-TAD | TAD **más pequeño que contiene** el CRM (sub-dominio específico, bibliografía). Evita el sesgo hacia TADs gigantes del "máximo solape" |
| **P7** | Reducción CRM **intra-TAD** | La comparación CRM-CRM y la reducción son el mismo paso, **dentro de cada TAD** (coherencia topológica). Criterio: recíproco ≥ 0.50 Y (Jaccard ≥ 0.70 O Simpson ≥ 0.99) + **anti-chaining** |
| **LB** | Detección DRRs | **Apilamiento** (bloques de solapamiento), NO proximidad. Filtro de tamaño 100 kb (dominios SE en literatura). 4 clases por densidad |
| **LB** | Validación SEdb | Posicional + **control por tamaño** (modelo binomial negativa: clase + log-longitud) + génica |
| **LB** | Anotación topológica | Cada DRR → sub-TAD(s); % en un único sub-TAD |

## Decisiones metodológicas justificadas

- **Reducción TAD a 0.80 recíproco:** umbral que la literatura asocia a "mismo
  cuerpo de TAD"; validado (sin chaining, todos los clusters con núcleo).
- **Reducción CRM con anti-chaining:** las componentes conexas producían un
  megacluster por la altísima densidad de solapamiento (99.98% de CRMs solapan).
  El anti-chaining garantiza núcleo común en cada cluster.
- **CRM intra-TAD vs global:** la reducción global y la intra-TAD difieren <0.1%,
  confirmando que la redundancia de CRMs es esencialmente intra-dominio. Se usa
  la intra-TAD por coherencia con el diseño.
- **DRRs por apilamiento, no proximidad:** con 99.98% de solapamiento, el
  criterio de stitching (12.5 kb, ROSE) colapsa el cromosoma entero en una
  región. El apilamiento + filtro de tamaño produce regiones plausibles.
- **Filtro 100 kb:** la literatura sitúa la mediana de SE en 8.7–19 kb y los
  dominios SE hasta ~90 kb; >100 kb se marca "Extensive_overlap" (se conserva,
  no se valida).
- **Control por tamaño en SEdb:** el enriquecimiento bruto de Dense_complex es
  en parte efecto del tamaño. El modelo controlado por longitud da IRR ≈ 7.8
  (p < 1e-67) para Dense_complex vs Extended, confirmando señal de densidad
  independiente del tamaño.

## Resultados de referencia (chr8)

- TADs: 8372 → 1580 (reducción 81.1%)
- CRMs: 1.343.498 → ~531.000 (reducción ~60%, intra-TAD)
- DRRs: ~4200 candidatas + ~300 extensas
- Validación SEdb: Dense_complex IRR ≈ 7.8 (controlado por tamaño)
- Topología: ~93% de las DRRs en un único sub-TAD

## Dependencias R

`data.table`, `GenomicRanges`, `IRanges`, `igraph` (TAD), `MASS` (glm.nb),
`ggplot2` (figuras). Conflicto conocido: `shift()` de data.table se enmascara
por IRanges — usar siempre `data.table::shift` explícito.
