# Ensamble conformacional de la quimera LBPCheY con BioEmu

Ensamble de equilibrio de la proteína quimérica **LBPCheY** (339 residuos), en la que
el dominio CheY está insertado en la proteína de unión a leucina (LBP), generado con
[BioEmu](https://github.com/microsoft/bioemu) a partir de la secuencia del modelo 0
de AlphaFold3 (AlphaFold Server).

La pregunta: los dos lóbulos se pliegan bien por separado, pero ¿se mueven uno
respecto al otro?

> **Aviso de términos de uso.** Los archivos de `alphafold3/` son salida de AlphaFold
> Server, y los de `bioemu/` y `figuras/` derivan de ella (usan su secuencia y su
> modelo como referencia). Se distribuyen bajo y sujetos a los
> [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms)
> (copia en `alphafold3/TERMS_OF_USE_AlphaFold_Server.md`): solo uso no comercial,
> sin usarlos para entrenar modelos de predicción de estructura ni en sistemas
> automatizados de predicción de unión a ligandos. Modificaciones respecto a la
> salida original: ninguna en el `.cif`; el ensamble es un derivado generado con BioEmu.

## Dominios (numeración 1-indexada)

| Región | Residuos |
|---|---|
| LBP N-terminal | 1–140 |
| CheY | 141–243 |
| LBP C-terminal | 244–339 |

## Pipeline

1. **AlphaFold3** (AlphaFold Server): modelo 0, pTM = 0.81 (`alphafold3/`).
2. **BioEmu v1.4.1**, checkpoint `bioemu-v1.1`, 300 muestras con MSA de ColabFold.
   El filtro de calidad de BioEmu deja **246**.
3. **Cadenas laterales** con HPacker y **minimización local** con OpenMM
   (amber99sb + TIP3P, backbone restringido), con `bioemu.sidechain_relax`.
4. **Análisis** (`scripts/analiza_ensamble.py`): RMSD de Cα por lóbulo contra el modelo
   de AF3 y distancia entre los centros de masa de los dos lóbulos.

`scripts/run_bioemu.sh` corre los pasos 2 y 3. Las rutas del script son las
del servidor original; ajústalas antes de usarlo. El paso a paso completo, con
entradas, salidas y tiempos, está en **[TUTORIAL.md](TUTORIAL.md)**. Cómo funciona
BioEmu por dentro está explicado en **[BIOEMU.md](BIOEMU.md)**.

## Resultados

| Medida (Å) | Media ± sd | Rango |
|---|---|---|
| RMSD lóbulo CheY | 0.81 ± 0.30 | 0.40 – 2.35 |
| RMSD lóbulo LBP | 2.30 ± 1.87 | 0.76 – 18.71 |
| RMSD global | 3.56 ± 2.32 | 0.79 – 21.37 |
| Distancia entre centros de los lóbulos | 27.9 ± 2.0 | 23.0 – 36.2 |

- **CheY** se mantiene prácticamente rígido y muy parecido al modelo de AF3.
- **LBP** también se conserva en general, pero tiene unas cuantas muestras muy alejadas
  del modelo.
- **Entre los lóbulos** hay movimiento: el RMSD global es mayor que el de cada lóbulo y
  la distancia entre centros varía unos 13 Å.

![RMSD por lóbulo](figuras/rmsd_lobulos.png)
![Apertura entre lóbulos](figuras/apertura_lobulos.png)

Las estadísticas se calcularon con el ensamble de backbone de BioEmu (antes de la
minimización) y los centros de masa con los Cα.

## Poblaciones y energías libres por conformación

BioEmu solo entrega estructuras: muestras independientes, todas con el mismo peso. Las
poblaciones se obtienen contando cuántas muestras caen en cada conformación, y las
energías libres con ΔG = −kT·ln(pᵢ/p₁), a 300 K (kT = 0.596 kcal/mol). El cálculo está en
`scripts/poblaciones_energias.py`.

**Plegamiento.** Se usó la fracción de contactos nativos respecto al modelo de AF3,
con la misma definición que BioEmu usa en su entrenamiento (≥ 0.65 = plegado). LBP y
CheY están plegados en las **246/246** muestras. Por eso solo hay una cota: ΔG de
desplegamiento > 3.3 kcal/mol.

**Conformaciones.** Se agruparon según la posición de CheY tras alinear sobre LBP (PCA +
k-means, k = 3). Los intervalos son IC 95% por bootstrap:

| Conf. | Muestras | Población | ΔG (kcal/mol) | Distancia entre lóbulos | Desplazamiento de CheY vs AF3 | RMSD global vs AF3 |
|---|---|---|---|---|---|---|
| 1 | 172 | 69.9% (64–76) | 0 (referencia) | 27.4 ± 1.7 Å | 5.1 Å | 2.8 Å |
| 2 | 54 | 22.0% (17–27) | +0.69 (0.52–0.89) | 29.2 ± 2.2 Å | 10.5 Å | 4.7 Å |
| 3 | 20 | 8.1% (5–12) | +1.28 (1.04–1.59) | 28.3 ± 2.2 Å | 15.6 Å | 7.3 Å |

- **Conformación 1:** el modelo de AlphaFold3 cae en ella.
- **Conformaciones 2 y 3:** están a menos de ~2 kT y son accesibles a temperatura ambiente.
- **Tipo de movimiento:** la distancia entre lóbulos casi no cambia; lo que varía es la
  posición de CheY, que se desplaza o gira hasta 16 Å respecto a AF3.

![Energía libre por conformación](energias/energia_libre.png)

**Precauciones:**
- **Número de conformaciones:** la silueta es ~0.4 tanto con k = 2 como con k = 3, así
  que las posiciones de CheY forman más un continuo que estados separados. Con k = 2 queda
  77% / 23%, con ΔG = 0.71 kcal/mol.
- **Resolución:** con 246 muestras, una muestra equivale a 0.4%, y el ΔG máximo medible
  es ~3.3 kcal/mol.
- **Filtro:** BioEmu descartó 54 de 300 muestras, lo que puede sesgar las poblaciones.
- **Error del modelo:** BioEmu reporta ~1 kcal/mol de error en energías libres, así que
  la diferencia de 0.7 kcal/mol está dentro de ese margen.

## Archivos

| Ruta | Contenido |
|---|---|
| `alphafold3/` | Modelo 0 (`.cif`), resumen de confianza, solicitud del trabajo y términos de uso |
| `bioemu/ensamble_lbpchey.pdb.gz` | 246 estructuras minimizadas (átomos pesados) en un PDB multimodelo, alineadas sobre LBP y **ordenadas de cerrada a abierta** |
| `bioemu/orden_apertura.txt` | Para cada estado: muestra original y distancia entre lóbulos |
| `bioemu/repr_{cerrada,mediana,abierta}.pdb` | Estructuras representativas (solo backbone) |
| `bioemu/lbpchey_model0.fasta` | Secuencia usada |
| `visualizacion/ver_ensamble.pml` | Script de PyMOL |
| `energias/` | Poblaciones y ΔG por conformación (`conformaciones.csv`), asignación de cada muestra, estructura representativa de cada conformación y figura |
| `scripts/` | `run_bioemu.sh`, `analiza_ensamble.py`, `prepara_visualizacion.py` y `poblaciones_energias.py` |
| `BIOEMU.md` | Cómo funciona BioEmu por dentro, con referencias |
| `TUTORIAL.md` | Cómo se ejecutó todo el pipeline, paso a paso |
| `presentacion/` | Presentación sobre BioEmu y este caso (`.pptx` y `.pdf`) |

## Visualizar en PyMOL

Desde la raíz del repositorio:

```bash
pymol visualizacion/ver_ensamble.pml
```

Fondo negro, todo en cartoon. LBP en gris, CheY en naranja, el modelo de AF3 en azul.
Las representativas (`cerrada`, `mediana`, `abierta`) empiezan apagadas.

- `mplay` / `mstop` reproducen la película de cerrada a abierta.
- `set all_states, on` muestra los 246 estados encimados.

Los estados son muestras independientes del equilibrio, no una trayectoria
temporal, así que la película puede tener saltos.
