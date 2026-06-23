# TFM_pipeline

Repositorio asociado al Trabajo Fin de Máster en Bioinformática.

Este repositorio contiene un pipeline en R para reducir la redundancia estructural de anotaciones genómicas de TADs y CRMs, generar CRMs consenso y detectar regiones reguladoras densas (DRRs) dentro de un marco topológico reducido. Las DRRs se analizan posteriormente mediante comparación posicional y génica con anotaciones de super-enhancers recopiladas en SEdb.

Las regiones densas detectadas se interpretan como **candidatos estructurales** definidos por acumulación local de CRMs consenso, no como super-enhancers funcionalmente validados. El pipeline no utiliza señales homogéneas de H3K27ac, MED1, BRD4 ni otros datos funcionales equivalentes. Por tanto, la comparación con SEdb se plantea como una evaluación de concordancia entre representaciones, no como una validación experimental independiente.

## Estructura del repositorio

```text
TFM_pipeline/
├── main.R                      # Ejecución principal para un cromosoma
├── run_all_chromosomes.R       # Ejecución por cromosomas
├── aggregate_global.R          # Agregación global de resultados cromosómicos
├── README.md
├── R/                          # Funciones principales del pipeline
│   ├── p1_load_collapse.R
│   ├── tad_reduce_final.R
│   ├── crm_explore.R
│   ├── crm_antichaining_fast.R
│   ├── crm_reduce_final.R
│   ├── crm_reduce_intratad.R
│   ├── step6_assign_crm_subtad.R
│   ├── drr_compare_criteria.R
│   ├── drr_sedb_validation.R
│   ├── drr_gene_analysis.R
│   ├── drr_size_control_cases.R
│   ├── drr_lineB_stats_figs.R
│   └── drr_propagate_subtad.R
├── R_exploratory/              # Scripts exploratorios no ejecutados por main.R
├── data/                       # Directorio esperado para datos de entrada
└── results/                    # Directorio de salida generado durante la ejecución
````

## Datos de entrada

Los datos de entrada no se incluyen en este repositorio debido a su tamaño y a las condiciones de redistribución de algunos recursos externos.

La estructura esperada es:

```text
data/
├── tad_per_chr/
│   └── <chr>.tsv.gz
├── enh_per_chr/
│   └── <chr>.tsv.gz
├── enh2gene.tsv.gz
├── SEdb_Human_SE.bed
└── human_genes.tsv
```

Los ficheros por cromosoma contienen las anotaciones de TADs y CRMs utilizadas como entrada del pipeline. Los ficheros globales incluyen relaciones CRM-gen, anotaciones de super-enhancers de SEdb y coordenadas génicas humanas.

## Ejecución

Para ejecutar el pipeline sobre un cromosoma:

```bash
Rscript main.R chr8
```

Para forzar el recálculo completo e ignorar resultados intermedios guardados:

```bash
Rscript main.R chr8 force
```

Para ejecutar varios cromosomas mediante el script general:

```bash
Rscript run_all_chromosomes.R
```

También pueden indicarse cromosomas concretos:

```bash
Rscript run_all_chromosomes.R "" chr1,chr2
```

## Resultados intermedios y checkpoints

El pipeline guarda resultados intermedios en formato RDS para evitar recalcular pasos costosos. Estos archivos se generan localmente en:

```text
results/checkpoints/
```

Los resultados finales por cromosoma se almacenan en:

```text
results/lineB/<chr>/
```

Los archivos intermedios de gran tamaño no están versionados en GitHub.

## Flujo general del pipeline

| Paso      | Descripción                                    | Criterio principal                                              |
| --------- | ---------------------------------------------- | --------------------------------------------------------------- |
| P1        | Carga y normalización inicial de TADs y CRMs   | Colapso por identificador y depuración de coordenadas           |
| P2        | Reducción de redundancia de TADs               | Solapamiento recíproco ≥ 0.80                                   |
| P6        | Asignación de CRMs a TADs reducidos            | TAD reducido más específico disponible                          |
| P7        | Reducción intra-TAD de CRMs                    | Solapamiento recíproco ≥ 0.50 y Jaccard ≥ 0.70 o Simpson ≥ 0.99 |
| DRR       | Detección de regiones reguladoras densas       | Agrupación local de CRMs consenso                               |
| SEdb      | Comparación con anotaciones de super-enhancers | Concordancia posicional y génica                                |
| Topología | Caracterización topológica de DRRs             | Asignación de DRRs a TADs reducidos                             |

## Decisiones metodológicas principales

La reducción de TADs se realiza mediante solapamiento recíproco para consolidar anotaciones topológicas altamente coincidentes sin inferir nuevos dominios tridimensionales.

La reducción de CRMs se aplica dentro de cada TAD reducido, combinando solapamiento recíproco, índice de Jaccard y coeficiente de Simpson. Esta estrategia permite agrupar regiones estructuralmente redundantes manteniendo la trazabilidad con las anotaciones originales.

Las DRRs se definen como agrupaciones locales de CRMs consenso dentro del marco topológico reducido. Su clasificación se basa en propiedades estructurales como longitud, número de CRMs consenso, soporte acumulado y densidad.

La comparación con SEdb se utiliza para evaluar la concordancia entre las DRRs y anotaciones previamente recopiladas de super-enhancers. Debido al solapamiento parcial entre fuentes de anotación, esta comparación no se interpreta como una validación funcional independiente.

## Resultados de referencia

En el análisis desarrollado para el TFM, el pipeline permitió reducir de forma sustancial la redundancia de las anotaciones de TADs y CRMs y generar un conjunto de CRMs consenso trazables. A partir de estos consensos se identificaron DRRs, que fueron clasificadas y comparadas con anotaciones de SEdb.

Los resultados concretos dependen del conjunto de datos de entrada y de los cromosomas procesados. Las tablas y figuras finales se generan localmente en el directorio `results/`.

## Dependencias

El pipeline se desarrolló en R. Las principales dependencias son:

* data.table
* GenomicRanges
* IRanges
* igraph
* MASS
* ggplot2

Se recomienda registrar la versión de R y de los paquetes utilizados mediante `sessionInfo()` para facilitar la reproducibilidad.

## Reproducibilidad

El análisis se organiza mediante scripts modulares y ejecución por cromosoma. Los parámetros principales se encuentran definidos explícitamente en el código, y los resultados intermedios se guardan para permitir la reanudación del análisis.

Los datos de entrada y los resultados intermedios pesados no se distribuyen en este repositorio. Para reproducir el análisis, deben colocarse los ficheros requeridos en la estructura indicada y ejecutar los scripts correspondientes.
