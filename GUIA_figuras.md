# Guía de figuras del TFM — `figuras_tfm.R`

## Cómo ejecutar

```r
# 1) Coloca figuras_tfm.R donde tengas la carpeta results/
# 2) Apunta a tu carpeta de resultados (una de las dos vías):
Sys.setenv(TFM_RESULTS = "ruta/a/results")   # o edita RESULTS_DIR en el script
# 3) Genera todas las figuras:
source("figuras_tfm.R")
main(igv_engine = "both")     # "ggplot", "gviz" o "both"
# Salida en figs_out/
```

Para una figura suelta, sin generar todas:

```r
source("figuras_tfm.R")
D <- load_all("ruta/a/results")
print(fig06_forest(D))                       # la figura estrella
save_fig(fig06_forest(D), "forest.png", 7.5, 7.5)
```

Para el IGV de cualquier región:

```r
reg <- load_region("chr1", 149740000, 149960000, "ruta/a/results")
igv_gviz(reg, "igv_chr1.png")        # versión Gviz
print(igv_ggplot(reg))               # versión ggplot2
```

## Catálogo de figuras (recomendación)

| Fig | Función | Contenido | Sugerencia |
|-----|---------|-----------|------------|
| F1  | `fig01_reduccion`      | % reducción TAD/CRM por cromosoma | **Resultados** |
| F1b | `fig01b_reduccion_abs` | Conteos absolutos antes/después   | Suplementario (alternativa a F1) |
| F2  | `fig02_tam_cluster`    | Tamaño de clúster de reducción    | Suplementario |
| F3  | `fig03_conteos`        | DRRs por clase (genoma)           | **Resultados** |
| F3b | `fig03b_conteos_chr`   | Composición de clases por cromosoma | Suplementario |
| F4  | `fig04_tamanos`        | Tamaños por clase + banda SE lit. | **Resultados** |
| F5  | `fig05_recuperacion`   | Recuperación SEdb posicional+génica | **Resultados** |
| F6  | `fig06_forest`         | Forest plot IRR (chr + global)    | **Resultados (principal)** |
| F7  | `fig07_topologia`      | % en sub-TAD único por clase      | **Resultados** |
| F8  | `fig08_igv_*`          | Visualización IGV (Gviz/ggplot/panel) | **Resultados** |
| **F9**  | `fig09_flujo`      | **Esquema cuantitativo del pipeline** (flujo) | **Resultados** |
| **F10** | `fig10_soporte`    | **Soporte de los CRMs consenso** (donut) | **Resultados** |
| **F11** | `fig11_longitud`   | **Longitud CRMs originales vs consenso** | **Resultados** |
| **F12a**| `fig12_crm_por_tad`| **Carga reguladora por sub-TAD** (histograma) | Resultados/Supl. |
| **F12b**| `fig12b_crm_por_tad_chr` | **CRMs por sub-TAD antes/después** por cromosoma | **Resultados** |
| **F13** | `fig13_contraste`  | **Panel IGV de 4 regiones de contraste** | **Resultados** |

Para Resultados, una selección contundente que equilibra los tres bloques:
F9 (flujo), F1, F10 (soporte), F11 (longitud), F12b (carga por TAD), F3, F4,
F5, F6, F7, y F8/F13 (IGV + contraste). El resto a suplementario.

### Notas sobre las figuras nuevas
- **F9** lee solo los `summary` de los `.rds`; números globales del genoma.
- **F10** usa `n_entities` de los consenso; donut con bins 1 / 2–5 / 6–20 / >20.
  Variante de barra apilada: `fig10_soporte(D, "bar")`.
- **F11** necesita los **CRMs originales** en `data/enh_per_chr/chrN.tsv.gz`. El
  runner los busca en `<results>/../data/enh_per_chr`. Si no los encuentra, cae
  automáticamente a una variante alternativa (longitud del consenso por nivel de
  soporte), que también es válida. Comprueba que tus .tsv.gz tienen columnas
  `chr,start,end,ID`; si el ID está en otra columna, ajusta `idc` en la función.
- **F12** en dos variantes (`a` histograma global, `b` boxplot por cromosoma
  antes/después). Elige la que prefieras.
- **F13** usa `CONTRAST_REGIONS` (4 regiones predefinidas que cubren el rango de
  comportamiento, incluido un caso denso con recuperación SEdb nula). Edítalas
  libremente.

## Decisiones de diseño codificadas
- Paleta por clase fija (gradiente frío→cálido con la densidad): editable en
  `CLASS_COLOR`.
- Orden de clases siempre Simple → Compact → Extended → Dense complex.
- Banda gris en F4 = rango de mediana de SE en literatura (8,7–19 kb).
- F6 usa el modelo global de `results/global/` si existe; si no, solo por cromosoma.

## Reproducibilidad
Todo se deriva de `results/lineB/<chr>/pipeline_objects_<chr>.rds` y, para el
modelo global, de `results/global/model_coef_global.tsv`. No hay tablas
intermedias externas: basta el directorio `results/`.

## Aviso técnico
- `drr_id` se repite entre cromosomas (cada uno numera desde DRRs_00000001).
  Para cruces entre cromosomas usa la clave compuesta `chr + drr_id`.
- karyoploteR es una alternativa válida al IGV; no se incluye porque Gviz cubre
  el mismo propósito y es más fácil de instalar.