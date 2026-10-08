# BioEmu

Cómo funciona BioEmu por dentro: qué representa, qué red usa, cómo se entrenó, cómo
genera muestras y qué limitaciones tiene. Este documento trata solo del método. Para
ver cómo se aplicó a LBPCheY, consulta [TUTORIAL.md](TUTORIAL.md).

Cada afirmación lleva su fuente:
- **[P]** el artículo.
- **[R]** el README.
- **[MC]** la ficha del modelo.
- **[CFG]** la configuración del checkpoint.
- **[CODE]** el código. Revisamos `bioemu` 1.4.1, la versión instalada.

La lista completa está en [Referencias](#referencias). Si el artículo y el código no
coinciden, se indica.

---

## 1. Qué problema resuelve

Las proteínas no tienen una sola estructura. Pasan por un **ensamble de
conformaciones** con poblaciones dadas por la termodinámica (la distribución de
Boltzmann). Hay dos formas habituales de estudiarlo, y ninguna basta:

- **AlphaFold y similares** predicen una estructura, normalmente la más probable.
  No dicen cuánto se mueve la proteína ni con qué frecuencia visita otros estados.
- **La dinámica molecular (MD)** sí muestrea el ensamble, pero los cambios
  funcionales ocurren en escalas de µs a ms. Muestrearlos bien cuesta semanas o meses
  de GPU por proteína.

BioEmu (*Biomolecular Emulator*) es un **modelo generativo** que aprende a producir
muestras **estadísticamente independientes** de ese ensamble de equilibrio a partir
de la secuencia. Reporta una velocidad entre cuatro y cinco órdenes de magnitud mayor
que la MD [P]: 10 000 estructuras independientes en minutos a pocas horas en una GPU
[P], o "miles de estructuras por hora" [HF].

---

## 2. Visión general

```
 secuencia ──► MSA ──► AlphaFold2 (Evoformer) ──► embeddings single (L×384) y pair (L×L×128)
                                                         │  (se calculan una vez y se guardan)
                                                         ▼
 ruido ──► [ red de puntuación (IPA, 8 capas) ] ◄── tiempo de difusión t
   ▲              │ predice la "score" de posiciones y orientaciones
   │              ▼
   └──── solver DPM de 2.º orden, t: 0.99 → 0.001 ────► marcos de residuo (posición + rotación)
                                                         │
                                                         ▼
                              backbone N, CA, C, O ──► filtro de calidad ──► ensamble
```

Cada muestra parte de ruido distinto y se condiciona con los mismos embeddings. Por
eso las muestras son **independientes entre sí**: no forman una trayectoria.

---

## 3. Representación: la proteína como marcos rígidos

BioEmu no trabaja con átomos. Representa cada residuo como un **marco rígido**,
igual que AlphaFold2 [P]:

- **Una posición** en ℝ³, el Cα. El modelo trabaja en nanómetros [CODE `chemgraph.py`].
- **Una orientación**, una rotación en SO(3) dada como matriz 3×3, que define cómo
  están colocados N, Cα y C.

Una proteína de L residuos es entonces un conjunto de L pares (posición, rotación). La
difusión ocurre en ese espacio.

**Al final, de marcos a átomos.** Con cada marco se colocan N, Cα y C con la geometría
ideal de AlphaFold2. El oxígeno del carbonilo se pone en el plano del enlace peptídico,
a 1.23 Å del C [CODE `convert_chemgraph.py`, `compute_backbone`, `_adjust_oxygen_pos`].

El resultado tiene **solo átomos de backbone**. "No se modelan explícitamente ni las
cadenas laterales ni los hidrógenos" [P].

---

## 4. Condicionamiento: embeddings de AlphaFold2

BioEmu no lee la secuencia directamente. Usa las **representaciones internas del
Evoformer de AlphaFold2** [P]:

| Representación | Dimensión | Qué contiene |
|---|---|---|
| *single* | L × 384 | Un vector por residuo |
| *pair* | L × L × 128 | Un vector por par de residuos: información coevolutiva y de contactos |

Detalles [P] [CODE `colabfold_inline/model_runner.py`]:
- Se usa el **modelo 3 de AlphaFold2**, sin plantillas y **sin reciclado**
  (`num_recycle=0`). El código lo indica: *"Hardcoded to monomer alphafold2, model_3,
  num_recycle=0"*.
- Solo se usan los embeddings; la estructura que predice AF2 se descarta.
- Se calculan **una sola vez** por secuencia y se guardan en caché
  (`~/.bioemu_embeds_cache/<sha256>_single.npy` y `_pair.npy`) [CODE `get_embeds.py`].
- Los pesos de AF2 (~3.5 GB) se descargan la primera vez [R].

**El MSA, según la fuente:**
- **Artículo:** hhblits contra Uniclust30 (versión de agosto de 2018) [P].
- **Código publicado:** el servidor MMseqs2 de ColabFold, o un `.a3m` del usuario
  [R] [CODE `get_embeds.py`]. Si das tu propio MSA, el código advierte que
  *"might result in suboptimal performance"* si no se generó como los de ColabFold.

Como el MSA alimenta los embeddings, su calidad afecta todo lo que viene después.

---

## 5. La red de puntuación (*score model*)

La red recibe los marcos con ruido, el tiempo de difusión *t* y los embeddings.
Predice hacia dónde mover cada posición y cada orientación para quitar ruido. Su
arquitectura **"se parece a los módulos de estructura de AlphaFold2 y del
Distributional Graphormer (DiG)"** [P] [MC].

**Entrada** [CODE `models.py`, `DistributionalGraphormer.forward`]:
1. El embedding *single* se proyecta a 512 dimensiones y se le suma un embedding
   sinusoidal del tiempo *t*. Así la red sabe cuánto ruido hay.
2. El embedding *pair* se proyecta a 256 dimensiones y se le suma un sesgo según la
   distancia en secuencia entre residuos (64 cubetas, hasta 128 posiciones).

**Núcleo: atención de puntos invariante (IPA).** La red apila **8 capas** de IPA
alternadas con normalización y MLP [CFG] [CODE `structure_module.py`]. La IPA es el
mecanismo del módulo de estructura de AF2 (Alg. 22 del suplemento de AF2): cada residuo
atiende a los demás usando puntos expresados en su propio marco local.

| Hiperparámetro (v1.1) | Valor [CFG] |
|---|---|
| Capas | 8 |
| Cabezas de atención | 32 |
| Dimensión del modelo | 512 |
| Dimensión *pair* | 256 |
| Capa oculta del MLP | 1024 |
| Parámetros entrenables | 31.4 M (v1.0/v1.1), 35.7 M (v1.2) [MC] |

**Salida: dos puntuaciones** por residuo:
- **Traslación:** un vector 3D. Es **equivariante**: si giras la proteína, el vector
  gira con ella.
- **Rotación:** un vector en representación eje-ángulo. Es **invariante** a la pose
  global.

Gracias a esto, el modelo no depende de la orientación con que esté colocada la
proteína en el espacio [CODE `models.py`, docstring de `DistributionalGraphormer`].

---

## 6. El proceso de difusión

Un modelo de difusión aprende a **revertir un proceso que agrega ruido**:

- **Hacia adelante** (entrenamiento): se toma una estructura real y se le agrega ruido
  poco a poco hasta que queda aleatoria.
- **Hacia atrás** (muestreo): se parte de ruido y se quita paso a paso con la red,
  hasta obtener una estructura nueva.

BioEmu lo hace en paralelo sobre **posiciones** y **orientaciones**. Cada residuo se
perturba de manera independiente [P].

### 6.1 Posiciones: SDE que preserva la varianza (VP) con calendario coseno

La posición de cada residuo sigue una ecuación diferencial estocástica de tipo
*variance-preserving*, con calendario de ruido coseno [P] [CFG `CosineVPSDE`, s = 0.008]:

```
x_t = α(t) · x_0 + σ(t) · ε,      ε ~ N(0, I)
α(t) = cos( (t + s)/(1 + s) · π/2 ) / cos( s/(1 + s) · π/2 ),     σ(t)² = 1 − α(t)²
```

- En **t = 0** las posiciones son las de la estructura real.
- En **t = 1** son ruido gaussiano puro.
- El calendario coseno agrega ruido despacio al principio, lo que conserva la
  estructura fina durante más pasos [CODE `sde_lib.py`].

### 6.2 Orientaciones: difusión en SO(3)

Una rotación no se puede sumar con ruido gaussiano, porque SO(3) es una variedad
curva. BioEmu usa un **movimiento browniano en SO(3)**, cuya distribución marginal es
la **IGSO(3)**, la gaussiana isotrópica sobre rotaciones [CODE `so3_sde.py`, `DiGSO3SDE`]:
- **Calendario:** geométrico [P], con σ entre 0.02 y 2.33 [CFG].
- **Tipo de proceso:** *variance exploding* [CODE].
- **Densidades:** la IGSO(3) y su gradiente no tienen forma cerrada simple, así que
  se precalculan en tablas por expansión en series. Se guardan en `~/sampling_so3_cache`.
  Eso es lo que hace el mensaje `Computing igso3_expansion` al arrancar.

---

## 7. Muestreo: de ruido a estructura

Al generar una muestra [CODE `denoiser.py`, `dpm_solver`]:

1. Se crea un "grafo de contexto" con la secuencia y los embeddings. Las posiciones y
   orientaciones empiezan vacías [CODE `sample.py`, `get_context_chemgraph`].
2. Se sortean posiciones gaussianas y rotaciones uniformes, que es el estado en t ≈ 1.
3. Se integra la ecuación inversa desde **t = 0.99** hasta **t = 0.001** con un
   **solver de segundo orden**. En cada paso la red se evalúa dos veces: una en el
   punto actual y otra en un punto intermedio (DPM-Solver-2, Lu et al. 2022).
4. Las orientaciones se actualizan con pasos propios en SO(3), usando la puntuación de
   rotación.

**Número de pasos según la fuente:**

| Fuente | Solver | Pasos |
|---|---|---|
| Artículo [P] | "sampler de segundo orden" | 100 |
| Código, opción por defecto (`denoiser_type="dpm"`) [CODE `config/denoiser/dpm.yaml`] | DPM-Solver, determinista (`noise: 0`) | 50 |
| Código, alternativa (`heun.yaml`) | Heun, estocástico (`noise: 0.5`) | 100 |

**Lotes.** El costo de memoria crece con L². Por eso el tamaño de lote se escala así
[CODE `sample.py`]:

```
batch_size = max(1, int(batch_size_100 × (100 / L)²))
```

`batch_size_100` es el lote que usarías con una proteína de 100 residuos.

**Semillas.** Cada lote usa la semilla `base_seed + número de muestras ya generadas`.
Así se puede reanudar una corrida o repetirla de forma reproducible [CODE `sample.py`].

**Tiempos de referencia** con 1000 muestras en una A100 80 GB y `batch_size_100=20` [R]:

| Longitud | Tiempo |
|---|---|
| 100 aa | ~4 min |
| 300 aa | ~40 min |
| 600 aa | ~150 min |

---

## 8. Filtrado de muestras no físicas

Algunas muestras salen con la cadena rota o con choques. Por defecto
(`filter_samples=True`) se descartan, así que **obtienes menos muestras de las que
pediste** [R]. Una muestra se rechaza si cumple cualquiera de estos criterios
[CODE `convert_chemgraph.py`, `filter_unphysical_traj`]:

| Criterio | Umbral |
|---|---|
| Distancia Cα–Cα entre residuos consecutivos | > 4.5 Å (cadena rota) |
| Distancia C–N entre residuos consecutivos (enlace peptídico) | > 2.0 Å |
| Dos átomos de residuos distintos (sin contar vecinos i+1 e i+2) | < 1.0 Å (choque) |

### 8.1 *Steering* físico (opcional)

En lugar de solo descartar, se puede **guiar el muestreo** para que salgan menos
muestras defectuosas [R] [CODE `config/steering/physical_steering.yaml`]. Durante el
tramo final (t de 0.1 a 0), se agregan potenciales que penalizan:
- distancias Cα–Cα lejos de 3.8 Å;
- pares de átomos a menos de 4.1 Å, sin contar vecinos cercanos en secuencia.

Se mantienen varias "partículas" por muestra (5 por defecto) y se remuestrean con
**Monte Carlo secuencial (SMC)**, según sus pesos. Con 3–10 partículas por muestra
"se reducen mucho" las muestras no físicas [R]. Se activa con:

```bash
python -m bioemu.sample ... --denoiser_config src/bioemu/config/steering/physical_steering.yaml
```

---

## 9. Entrenamiento en tres etapas

Lo que hace que BioEmu reproduzca **ensambles de equilibrio**, y no solo estructuras
plausibles, es cómo se entrenó.

### Etapa 1: preentrenamiento con la AlphaFold Database
- **Agrupamiento:** ~93 M secuencias de la AFDB se agruparon al 80% de identidad
  (~1.4 M clusters) y se filtraron a ~50 000 clusters estructuralmente diversos [P].
  La ficha del modelo cita 161 000 estructuras de la AFDB [MC].
- **Aumento de datos:** en cada ejemplo se elige un cluster al azar y luego una
  estructura dentro de él. Así una misma secuencia queda asociada a **varias
  conformaciones**, y el modelo aprende que una secuencia corresponde a una
  distribución de estructuras, no a una sola [P].
- **Objetivo:** el entrenamiento de difusión estándar (*score matching*). La red
  aprende a predecir el ruido que se agregó.

### Etapa 2: ajuste fino con dinámica molecular
Más de **200 ms de MD de todos los átomos**; 216 ms para v1.0/v1.1 [P] [MC]:

| Conjunto | Tamaño |
|---|---|
| Dominios CATH | ~1100 dominios de 50–200 aa, 1–5 µs cada uno; 46 ms en total |
| Octapéptidos | 1100 × 5 µs; 8 ms |
| MEGAsim (simulaciones propias) | 271 proteínas silvestres y 21 458 mutantes, ~1–1.5 µs cada una |
| Otros sistemas públicos | SARS-CoV-2, DDR1, SETD8, entre otros |

Campos de fuerza: AMBER ff99SB-ILDN, ff14SB y ff99SB-disp (este último para regiones
desordenadas) [P].

Las simulaciones no llegan al equilibrio por sí solas. Cuando fue posible, los datos de
MD se **repesaron hacia el equilibrio** con **modelos de Markov** o con pesos derivados
de datos experimentales [P]. Así el modelo no hereda el sesgo de dónde empezó cada
simulación.

### Etapa 3: ajuste con datos experimentales (PPFT)
*Property-Prediction Fine-Tuning* ajusta el modelo para que la **estabilidad de
plegamiento** de sus ensambles coincida con la medida en el laboratorio [P]
[CODE `training/loss.py`, `training/foldedness.py`].

- **Datos:** el conjunto MEGAscale (Tsuboyama et al.), con cientos de miles de ΔG
  experimentales. Son 502 000 en v1.1 y 1.3 M en v1.2 [MC].
- **Qué cuenta como plegado:** la **fracción de contactos nativos** (FNC) respecto a
  una estructura de referencia. El código convierte la FNC en una "plegabilidad" entre
  0 y 1 con una sigmoide, `foldedness = σ(2·k·(FNC − umbral))` [CODE]. El artículo fija
  el estado plegado en FNC ≥ 0.65 [P].
- **Cómo se entrena:**
  1. Para cada proteína se generan varias muestras con un muestreo rápido de pocos
     pasos (8 pasos según [P]).
  2. Se calcula su plegabilidad media y se compara con la que se espera del ΔG medido.
  3. El error se propaga hacia atrás a través del propio muestreo.

  **No se necesitan estructuras** de los estados desplegados, solo el número
  experimental.

Esta etapa es la que permite que BioEmu prediga **energías libres relativas**.

### Versiones del modelo [R] [MC]

| Checkpoint | Datos | Nota |
|---|---|---|
| `bioemu-v1.0` | AFDB + MD (216 ms) + 19 k ΔG | El del preprint |
| `bioemu-v1.1` | AFDB + MD (216 ms) + 502 k ΔG | El del artículo en *Science*; es el **predeterminado** |
| `bioemu-v1.2` | AFDB + MD (145.4 ms) + 1.3 M ΔG | Agrega embeddings extra de tipo de residuo y de pares; 35.7 M parámetros |

---

## 10. Qué tan bien funciona

El artículo define la **cobertura** como el porcentaje de estructuras de referencia
para las que BioEmu genera muestras cercanas, dentro de un umbral de RMSD [P]. Resultados:

| Prueba | Resultado [P] | Ficha del modelo v1.1 [MC] |
|---|---|---|
| Movimientos de dominio (≤ 3 Å) | 85% | 83% |
| Desplegamiento local: estado plegado / desplegado | 72% / 74% | 70% / 82% |
| Bolsillos crípticos: holo / apo (1.5 Å) | 85% / 49% | 88% / 55% |
| Energías libres vs. MD de CATH (MAE) | 0.91 kcal/mol | — |
| Energías libres vs. *fast folders* de DESRES (MAE) | 0.74 kcal/mol | — |
| Estabilidad experimental | MAE < 0.8 kcal/mol, Spearman > 0.65 | ~0.9 kcal/mol, Spearman 0.6 |

Los valores de la ficha del modelo difieren un poco de los del artículo. Reflejan la
versión publicada del checkpoint. Las comparaciones con otros métodos (AlphaFlow,
AFCluster, MD) están en el repositorio
[bioemu-benchmarks](https://github.com/microsoft/bioemu-benchmarks).

---

## 11. Después del muestreo: cadenas laterales

BioEmu entrega solo backbone. El paquete incluye un paso opcional
(`pip install bioemu[md]`) [R] [CODE `sidechain_relax.py`]:

1. **HPacker** coloca las cadenas laterales sobre cada backbone. Se instala en un
   entorno aparte.
2. **OpenMM** prepara cada estructura:
   - **Sistema:** AMBER99SB, agua TIP3P y una caja con 1 nm de margen.
   - **Minimización:** se mantiene el backbone restringido y se minimiza para quitar
     choques.
   - **Equilibrado (opcional):** NVT/NPT corto.

Esto **no es una simulación de equilibrio**. Solo hace químicamente razonables las
muestras de BioEmu.

---

## 12. Limitaciones

Según los autores [P] [R] [MC]:

- **Una sola cadena:** no modela complejos, multímeros, ligandos, membranas ni
  cofactores. Unir cadenas con un *linker* "no funcionó bien".
- **Condición fija:** ensambles a ~300 K, la temperatura de la MD de entrenamiento [P]. No se puede cambiar la temperatura, el pH ni
  la concentración de sales.
- **Solo backbone:** las cadenas laterales se agregan después y no salen del modelo
  generativo.
- **Sin cinética:** las muestras son independientes. No dan tiempos, velocidades ni
  rutas de transición, así que **ordenarlas no produce una trayectoria**.
- **Hereda aproximaciones:** de los embeddings de AF2 (y por lo tanto del MSA) y de los
  campos de fuerza usados en la MD de entrenamiento.
- **Longitud:** el artículo no fija un máximo. El entrenamiento con MD usó dominios de
  50–200 aa, el README reporta tiempos hasta 600 aa, y el costo crece con L².

---

## 13. Resumen

| Aspecto | BioEmu |
|---|---|
| Entrada | Secuencia de una cadena (con MSA de ColabFold o un `.a3m` propio) |
| Condicionamiento | Embeddings *single* y *pair* del Evoformer de AF2 (modelo 3, sin reciclado) |
| Representación | Un marco rígido por residuo: posición en ℝ³ y rotación en SO(3) |
| Red | Transformer de IPA estilo DiG/AF2: 8 capas, ~31 M parámetros |
| Difusión | VP-SDE coseno para las posiciones; IGSO(3) con calendario geométrico para las rotaciones |
| Muestreo | DPM-Solver de 2.º orden, 50 pasos (código) o 100 (artículo) |
| Entrenamiento | AFDB → MD repesada (>200 ms) → PPFT con ΔG experimentales |
| Salida | Muestras independientes de backbone (N, Cα, C, O), filtradas |
| Velocidad | Miles de estructuras por hora en una GPU |

---

## Referencias

- **[P]** Lewis S., Hempel T., Jiménez-Luna J., et al. *Scalable emulation of protein
  equilibrium ensembles with generative deep learning*. **Science** 389 (2025).
  doi:[10.1126/science.adv9817](https://doi.org/10.1126/science.adv9817).
  Preprint: [bioRxiv 10.1101/2024.12.05.626885](https://www.biorxiv.org/content/10.1101/2024.12.05.626885v2.full).
- **[R]** README del repositorio: <https://github.com/microsoft/bioemu>
- **[MC]** Ficha del modelo: <https://github.com/microsoft/bioemu/blob/main/MODEL_CARD.md>
- **[HF]** Pesos y ficha en Hugging Face: <https://huggingface.co/microsoft/bioemu>
- **[CFG]** Configuración del checkpoint v1.1:
  <https://huggingface.co/microsoft/bioemu/resolve/main/checkpoints/bioemu-v1.1/config.yaml>
- **[CODE]** Código fuente de `bioemu` 1.4.1, en `src/bioemu/` del repositorio:
  - `models.py`, `structure_module.py`: la red de puntuación.
  - `sde_lib.py`, `so3_sde.py`: los procesos de difusión.
  - `denoiser.py`, `config/denoiser/`: los solvers de muestreo.
  - `sample.py`, `get_embeds.py`, `colabfold_inline/`: el muestreo y los embeddings.
  - `convert_chemgraph.py`: el backbone y el filtro.
  - `training/`: la pérdida de PPFT.
  - `steering/`, `config/steering/`: el *steering*.
  - `sidechain_relax.py`: cadenas laterales y relajación.
- Arquitectura base, Distributional Graphormer (DiG): Zheng S. et al. *Predicting
  equilibrium distributions for molecular systems with deep learning*. Nature Machine
  Intelligence (2024). [arXiv:2306.05445](https://arxiv.org/abs/2306.05445).
- AlphaFold2: Jumper J. et al. *Highly accurate protein structure prediction with
  AlphaFold*. Nature 596, 583–589 (2021).
- DPM-Solver: Lu C. et al. *DPM-Solver: A Fast ODE Solver for Diffusion Probabilistic
  Model Sampling in Around 10 Steps*. NeurIPS (2022).
  [arXiv:2206.00927](https://arxiv.org/abs/2206.00927).
- *Steering*: [arXiv:2501.06848](https://arxiv.org/abs/2501.06848), citado en [R].
- MEGAscale: Tsuboyama K. et al. *Mega-scale experimental analysis of protein folding
  stability in biology and design*. Nature 620, 434–444 (2023).
