# =====================================================================
#  Visualizacion del ensamble BioEmu de LBPCheY
#  Uso (desde la raiz del repositorio):  pymol visualizacion/ver_ensamble.pml
#  (las rutas son relativas a la raiz del repositorio)
# =====================================================================

# --- estructuras --------------------------------------------------------
# ensamble: 246 estados alineados sobre LBP y ordenados de cerrada a abierta (ver orden_apertura.txt)
# af3: modelo 0 de AlphaFold3
load bioemu/ensamble_lbpchey.pdb.gz, ensamble
load alphafold3/fold_lbpchey_model_0.cif, af3
load bioemu/repr_cerrada.pdb, cerrada
load bioemu/repr_mediana.pdb, mediana
load bioemu/repr_abierta.pdb, abierta

# alinear AF3 y las representativas sobre el lobulo LBP del ensamble
select lbp,  resi 1-140 or resi 244-339
select chey, resi 141-243
super af3     and lbp and name CA, ensamble and lbp and name CA
super cerrada and lbp and name CA, ensamble and lbp and name CA
super mediana and lbp and name CA, ensamble and lbp and name CA
super abierta and lbp and name CA, ensamble and lbp and name CA

# --- estilo -------------------------------------------------------------
hide everything
dss
show cartoon
set cartoon_transparency, 0
bg_color black

color grey80, ensamble and lbp
color orange, ensamble and chey
color marine, af3
color palegreen, cerrada
color wheat,     mediana
color salmon,    abierta

# solo el ensamble y AF3 visibles al inicio
disable cerrada
disable mediana
disable abierta

deselect
orient ensamble

# --- pelicula -------------------------------------------------------------
# mplay / mstop para reproducir, o los botones abajo a la derecha.
# "set all_states, on" muestra los 246 estados encimados.
mset 1 -246
set movie_fps, 15
