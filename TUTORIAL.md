# Tutorial: ensamble conformacional de LBPCheY con BioEmu

Este tutorial reproduce paso a paso lo que se hizo para obtener los resultados del
repositorio: desde el modelo de AlphaFold3 hasta la película en PyMOL. En cada paso
se indica qué script se ejecuta, qué necesita de entrada, qué produce y cuánto tardó
en el equipo original.

**Equipo original:** Linux, GPU NVIDIA RTX 6000 Ada (48 GB), driver con CUDA 13.0,
conda (miniforge).

## Resumen del pipeline

```
AlphaFold Server ──► fold_lbpchey_model_0.cif
                            │
                            ▼  (secuencia)
              scripts/run_bioemu.sh
              ├─ 1. MSA (ColabFold) + embeddings AF2
              ├─ 2. muestreo BioEmu ─────────► topology.pdb + samples.xtc        (backbone)
              └─ 3. sidechain_relax ─────────► samples_md_equil.pdb/.xtc         (átomos pesados)
                            │
              ┌─────────────┴──────────────┐
              ▼                            ▼
   scripts/analiza_ensamble.py    scripts/prepara_visualizacion.py
   RMSD, apertura, figuras,       ensamble alineado y ordenado
   representativas                para PyMOL
              └─────────────┬──────────────┘
                            ▼
              visualizacion/ver_ensamble.pml  (PyMOL)
```

| Paso | Script | Tiempo aproximado |
|---|---|---|
| 0 | Predicción en AlphaFold Server | minutos (web) |
| 1 | Entorno conda | ~10 min, una sola vez |
| 2–4 | `run_bioemu.sh` | ~4 h 15 min en total |
| 5 | `analiza_ensamble.py` | < 1 min |
| 6 | `prepara_visualizacion.py` | < 1 min |
| 7 | `ver_ensamble.pml` | inmediato |

## Rutas usadas

Los scripts tienen rutas absolutas fijas. Si las tuyas son otras, edita las
variables al inicio de cada script.

| Variable | Ruta original | Qué hay ahí |
|---|---|---|
| `FOLD` | `~/fold_lbpchey` | Salida descomprimida de AlphaFold Server |
| `BASE` | `~/cursos/proteinas/quimerasLBPChey` | Carpeta de trabajo |
| `OUT` | `$BASE/bioemu` | Salidas de BioEmu |

---

## Paso 0. Predicción con AlphaFold3

**Entrada:** la secuencia de LBPCheY (339 aa). La solicitud enviada está en
`alphafold3/fold_lbpchey_job_request.json`: una cadena, con plantillas y una semilla.

**Cómo:** se envía el trabajo en https://alphafoldserver.com, se descarga el `.zip` y
se descomprime en `~/fold_lbpchey/`.

**Salida:** 5 modelos (`fold_lbpchey_model_{0..4}.cif`), sus confianzas (`*_summary_confidences_*.json`,
`*_full_data_*.json`), MSAs y plantillas. Solo se usa el **modelo 0**, que es el de mejor
`ranking_score` (pTM = 0.81). Está en `alphafold3/`.

> Los archivos de AlphaFold Server están sujetos a sus
> [términos de uso](https://alphafoldserver.com/output-terms)
> (`alphafold3/TERMS_OF_USE_AlphaFold_Server.md`).

---

## Paso 1. Entorno

```bash
conda create -y -n bioemu python=3.11
conda activate bioemu
export PYTHONNOUSERSITE=1
python -m pip install "bioemu[md]" mdtraj matplotlib
```

`run_bioemu.sh` hace esto automáticamente si el entorno no existe o si le falta `bioemu`.

> ⚠️ **Usa `python -m pip`, no `pip`.** En el equipo original, el `pip` que aparecía en
> el `PATH` era `~/.local/bin/pip`, que corresponde al Python 3.12 del usuario y no al
> del entorno. Con ese `pip` los paquetes se instalan fuera del entorno y luego
> aparece `ModuleNotFoundError: No module named 'bioemu'`.
> `PYTHONNOUSERSITE=1` evita que el entorno lea paquetes de `~/.local`.

Para comprobar que la GPU funciona:

```bash
python -c "import torch, bioemu; print(torch.__version__, torch.cuda.is_available())"
# 2.14.1+cu130 True
```

Versiones usadas: `bioemu` 1.4.1, `torch` 2.14.1, `mdtraj` 1.11, `openmm` 8.6.1.

---

## Pasos 2–4. `scripts/run_bioemu.sh`

```bash
bash scripts/run_bioemu.sh                 # 300 muestras
NSAMPLES=1000 BATCH=50 bash scripts/run_bioemu.sh  # opcional: más muestras, lotes más grandes
```

Son más de 4 horas, así que conviene correrlo dentro de `tmux` o `screen`, o con
`nohup ... > run_bioemu.log 2>&1 &`.

**Entrada:** `~/fold_lbpchey/fold_lbpchey_model_0.cif`.

**Parámetros (variables de entorno):**

| Variable | Default | Significado |
|---|---|---|
| `NSAMPLES` | 300 | Estructuras a generar, antes del filtro |
| `BATCH` | 10 | Tamaño de lote **para una proteína de 100 aa**. BioEmu lo escala por (100/L)², así que con 339 aa queda en 1. Con 48 GB de VRAM puedes subirlo a ~50. |

### 2. Secuencia, MSA y embeddings

El script extrae la secuencia del modelo 0 con Biopython, verifica que tenga 339
residuos y la escribe en `$OUT/lbpchey_model0.fasta`.

BioEmu **no usa las coordenadas** del modelo, solo la secuencia. El modelo de AF3
sirve aquí para obtener la secuencia y en el paso 5 como referencia.

En la primera ejecución, BioEmu:
- envía la secuencia al servidor MSA de ColabFold (`api.colabfold.com`), para lo que
  necesita internet;
- descarga los pesos de AlphaFold2 (~3.5 GB) y calcula los embeddings.

Tardó ~20 min porque JAX corrió en CPU: el `jaxlib` que se instaló no tiene CUDA.
Instalar `jax[cuda12]` lo acelera. Los embeddings quedan en `~/.bioemu_embeds_cache`,
así que las siguientes ejecuciones con la misma secuencia se saltan este paso.

### 3. Muestreo

```bash
python -m bioemu.sample --sequence $OUT/lbpchey_model0.fasta \
    --num_samples 300 --batch_size_100 10 --output_dir $OUT/lbpchey
```

El modelo de difusión `bioemu-v1.1` genera estructuras **solo de backbone**,
~3.5 s por estructura con lote de 1 (unos 17 min las 300). Al final BioEmu descarta
las estructuras con choques o con la cadena rota: aquí **quedaron 246 de 300**.
El filtro se puede desactivar con `--filter_samples=False`.

**Salida** en `$OUT/lbpchey/`:

| Archivo | Contenido |
|---|---|
| `topology.pdb` | Topología (backbone) |
| `samples.xtc` | Las 246 estructuras |
| `batch_*.npz` | Lotes crudos; sirven para reanudar si se interrumpe |

Si vuelves a correr el script con la misma `--output_dir`, BioEmu solo genera las
muestras que falten.

### 4. Cadenas laterales y minimización

```bash
python -m bioemu.sidechain_relax --pdb-path $OUT/lbpchey/topology.pdb \
    --xtc-path $OUT/lbpchey/samples.xtc --outpath $OUT/lbpchey
```

1. **HPacker** reconstruye las cadenas laterales. La primera vez se instala en
   `~/.hpacker_venv`. Tardó ~21 s por estructura (**1 h 30 min**).
   Salida: `samples_sidechain_rec.pdb` / `.xtc`.
2. **OpenMM** solvata cada estructura (amber99sb + TIP3P), restringe el backbone y
   minimiza. Tardó ~25–45 s por estructura (**2 h 07 min**). La mayor parte del tiempo
   se va en preparar el sistema en un solo hilo de CPU, por eso la GPU marca 0% casi
   siempre. Salida: `samples_md_equil.pdb` / `.xtc`, con las **246 estructuras con
   átomos pesados**.

Si este paso falla, el script avisa y continúa: el análisis del paso 5 solo usa el
backbone.

Al terminar el log dice:
```
 Listo. Resultados en: .../bioemu/lbpchey
```

---

## Paso 5. Análisis: `scripts/analiza_ensamble.py`

```bash
python scripts/analiza_ensamble.py
```

**Entrada:**
- `$OUT/lbpchey/samples.xtc` + `topology.pdb` (ensamble de backbone);
- `~/fold_lbpchey/fold_lbpchey_model_0.cif` como referencia. Si no lo puede leer, usa el
  primer frame del ensamble.

**Qué calcula**, solo con Cα (las fronteras están al inicio del script):

| Medida | Definición |
|---|---|
| RMSD lóbulo LBP | Residuos 1–140 + 244–339, alineando sobre esos mismos residuos |
| RMSD lóbulo CheY | Residuos 141–243 |
| RMSD global | Toda la proteína |
| Apertura | Distancia entre los centros de masa de los Cα de LBP y de CheY |

Si los dos lóbulos tienen RMSD bajo y el global es alto, cada lóbulo se pliega bien
pero su orientación relativa cambia.

**Salida en pantalla:**
```
ensamble: 246 estructuras, 339 residuos
referencia: modelo de AlphaFold

--- RMSD respecto a la referencia ---
lobulo LBP  (1-140, 244-339)         2.30 +/-  1.87 A   [0.76 - 18.71]
lobulo CheY (141-243)                0.81 +/-  0.30 A   [0.40 - 2.35]
estructura completa                  3.56 +/-  2.32 A   [0.79 - 21.37]

--- Orientacion relativa ---
distancia entre centros             27.86 +/-  2.04 A   [22.95 - 36.17]
```

**Archivos** (en `$BASE`):

| Archivo | Contenido |
|---|---|
| `rmsd_lobulos.png` | Histogramas de RMSD por lóbulo y global |
| `apertura_lobulos.png` | Histograma de la distancia entre lóbulos |
| `repr_cerrada.pdb` | Estructura con la menor distancia entre lóbulos (22.9 Å) |
| `repr_mediana.pdb` | Estructura con la distancia mediana (27.7 Å); no es un promedio de coordenadas |
| `repr_abierta.pdb` | Estructura con la mayor distancia (36.2 Å) |

Las representativas tienen solo backbone porque salen del ensamble del paso 3.
En este repositorio, las figuras están en `figuras/` y las representativas en `bioemu/`.

> Para comparar con la referencia, el script toma los índices de Cα **por separado**
> del ensamble y del modelo de AF3, porque el modelo tiene todos los átomos y el
> ensamble solo backbone. Si se usaran los mismos índices para ambos, el RMSD sería
> incorrecto.

---

## Paso 6. Preparar el ensamble para PyMOL: `scripts/prepara_visualizacion.py`

```bash
python scripts/prepara_visualizacion.py
```

**Entrada:** `samples_md_equil.pdb` / `.xtc` (paso 4) y el modelo 0 de AF3.

**Qué hace:**
1. Alinea las 246 estructuras sobre el lóbulo **LBP** del modelo de AF3, para que en
   la película LBP quede fijo y se vea moverse a CheY.
2. Las ordena de la **más cerrada a la más abierta** según la distancia entre los
   centros de los lóbulos.
3. Escribe un PDB multimodelo comprimido.

**Salida** (en `$BASE/visualizacion/`; en este repositorio están copiados en `bioemu/`):

| Archivo | Contenido |
|---|---|
| `ensamble_lbpchey.pdb.gz` | 246 modelos ordenados (~10 MB; 51 MB sin comprimir) |
| `orden_apertura.txt` | `estado_pelicula  muestra_original  distancia_centros_A` |

```
# estado_pelicula muestra_original distancia_centros_A
1 200 23.08
2 59 24.22
...
246 101 36.22
```

Estas distancias se calculan sobre el ensamble minimizado, así que difieren unas
décimas de las del paso 5, que usa el backbone sin minimizar.

---

## Paso 7. Ver en PyMOL: `visualizacion/ver_ensamble.pml`

Desde la raíz del repositorio:

```bash
pymol visualizacion/ver_ensamble.pml
```

**Entrada:** `bioemu/ensamble_lbpchey.pdb.gz`, `alphafold3/fold_lbpchey_model_0.cif` y
`bioemu/repr_*.pdb`, con rutas relativas a la raíz del repositorio.

**Qué hace:** carga todo y alinea el modelo de AF3 y las representativas sobre LBP.
Luego pone el fondo negro, asigna estructura secundaria (`dss`), muestra todo en
cartoon y arma una película de 246 frames.

| Objeto | Color | Al abrir |
|---|---|---|
| `ensamble`: LBP | gris | visible |
| `ensamble`: CheY | naranja | visible |
| `af3` | azul | visible |
| `cerrada` | verde claro | apagado |
| `mediana` | beige (wheat) | apagado |
| `abierta` | salmón | apagado |

Comandos útiles:

```
mplay                      # reproducir la película (cerrada → abierta)
mstop                      # detener
set all_states, on         # ver los 246 estados encimados
enable mediana             # encender una representativa
disable ensamble           # apagar el ensamble para comparar representativas
```

> Los estados son **muestras independientes** del equilibrio, no una trayectoria
> en el tiempo. Ordenarlos por apertura muestra la tendencia, pero puede haber saltos
> entre estados con apertura parecida y CheY girado hacia otro lado.

> **Detalle de PyMOL:** en un `.pml` no pongas comentarios `#` al final de una línea
> `load`, porque PyMOL los toma como parte del nombre del objeto. Tampoco pongas `;`
> en los comentarios, porque PyMOL lo interpreta como separador de comandos.

---

## Paso 8. Poblaciones y energías libres: `scripts/poblaciones_energias.py`

```bash
python scripts/poblaciones_energias.py
```

BioEmu no calcula energías. Las muestras valen todas lo mismo, así que la población
de un estado es la fracción de muestras que caen en él, y su energía libre relativa es
ΔG = −kT·ln(pᵢ/p₁), con T = 300 K.

**Entrada:** `samples.xtc` + `topology.pdb` (ensamble de backbone) y el modelo 0 de AF3.
Requiere `scikit-learn` (`python -m pip install scikit-learn`).

**Qué hace:**
1. **Plegamiento por lóbulo:** calcula la fracción de contactos nativos Cα–Cα respecto
   a AF3, igual que `bioemu/training/foldedness.py`, y considera plegado ≥ 0.65. Si no hay
   muestras desplegadas, reporta solo la cota ΔG > kT·ln(N).
2. **Conformaciones:**
   - alinea todo sobre LBP;
   - hace PCA de las coordenadas Cα de CheY;
   - agrupa con k-means sobre 3 componentes, con `N_CONF = 3`.
   - Elige k mirando la silueta: 0.42 con k = 2 y 0.41 con k = 3.
3. **Por conformación:** población, ΔG con IC 95% (2000 remuestreos de bootstrap),
   distancia entre lóbulos, desplazamiento de CheY, RMSD frente a AF3 y la estructura
   más cercana al centroide.

**Salida en pantalla (resumida):**
```
LBP   plegadas 246/246 ...   dG_despl = > 3.3 (cota: 0 de 246 desplegadas) kcal/mol
CheY  plegadas 246/246 ...   dG_despl = > 3.3 (cota: 0 de 246 desplegadas) kcal/mol
conf    n           poblacion %             dG kcal/mol        dist A  desp CheY A  RMSD glob
1     172     69.9 [64.2-75.6]       0.00 [-0.00-0.00]    27.4 ± 1.7        5.1        2.8
2      54     22.0 [17.1-27.2]       0.69 [0.52-0.89]    29.2 ± 2.2       10.5        4.7
3      20      8.1 [ 4.9-11.8]       1.28 [1.04-1.59]    28.3 ± 2.2       15.6        7.3
el modelo de AlphaFold3 cae en la conformacion 1
```

**Archivos** (en `$BASE/energias/`; en este repositorio, en `energias/`):

| Archivo | Contenido |
|---|---|
| `conformaciones.csv` | Una fila por conformación: población, ΔG con IC, descriptores |
| `asignacion_muestras.txt` | Conformación, distancia, desplazamiento y contactos nativos de cada muestra |
| `conformacion_{1,2,3}.pdb` | Estructura representativa de cada conformación (backbone, alineada sobre LBP) |
| `energia_libre.png` | Mapa de energía libre en PC1–PC2 y ΔG por conformación |

> **Límites:** con N muestras, la población mínima es 1/N y el ΔG máximo medible es
> kT·ln(N): ~3.3 kcal/mol con 246 muestras. Para resolver estados más raros, genera más
> muestras.

---

## Problemas encontrados y soluciones

| Síntoma | Causa | Solución |
|---|---|---|
| `No module named 'bioemu'` aunque "se instaló" | `pip` del `PATH` era el de `~/.local` (Python 3.12) | `python -m pip` + `PYTHONNOUSERSITE=1` |
| `conda activate` falla con `unbound variable` | `set -u` en el script | Usar `set -eo pipefail`, sin `-u` |
| Embeddings muy lentos (~20 min) | `jaxlib` sin CUDA, corre en CPU | Opcional: `python -m pip install "jax[cuda12]"` |
| `Using batch size 1` | `batch_size_100` se escala por (100/L)² | Subir `BATCH` si sobra VRAM |
| RMSD incorrecto contra AF3 | Índices de átomos del ensamble aplicados al `.cif` | Seleccionar los Cα por separado en cada estructura |
| La película de PyMOL no se movía | Comentario `#` al final de `load` cambiaba el nombre del objeto | Comentarios en líneas aparte |
