#!/usr/bin/env python3
"""
Prepara el ensamble minimizado de BioEmu para verlo en PyMOL.

- alinea todas las estructuras sobre el lobulo LBP del modelo de AlphaFold
- las ordena de la mas cerrada a la mas abierta (distancia entre los
  centros de los lobulos LBP y CheY)
- escribe un PDB multimodelo comprimido y una tabla con el orden

Uso:  python prepara_visualizacion.py
"""
import gzip, os, shutil
import numpy as np, mdtraj as md

BASE = os.path.expanduser("~/cursos/proteinas/quimerasLBPChey")
ENS  = f"{BASE}/bioemu/lbpchey"
AF   = os.path.expanduser("~/fold_lbpchey/fold_lbpchey_model_0.cif")
OUT  = f"{BASE}/visualizacion"

# mismas fronteras que analiza_ensamble.py, en resid de mdtraj (0-indexado)
LBP  = "name CA and (resid 0 to 139 or resid 243 to 338)"
CHEY = "name CA and resid 140 to 242"

os.makedirs(OUT, exist_ok=True)
traj = md.load(f"{ENS}/samples_md_equil.xtc", top=f"{ENS}/samples_md_equil.pdb")
ref  = md.load(AF)

# alinear sobre LBP para que lo que se mueva en la pelicula sea CheY
traj.superpose(ref, atom_indices=traj.topology.select(LBP),
               ref_atom_indices=ref.topology.select(f"protein and {LBP}"))

lbp, chey = traj.topology.select(LBP), traj.topology.select(CHEY)
dist = np.linalg.norm(traj.xyz[:, lbp].mean(1) - traj.xyz[:, chey].mean(1), axis=1) * 10
orden = np.argsort(dist)

pdb = f"{OUT}/ensamble_lbpchey.pdb"
traj[orden].save_pdb(pdb)
with open(pdb, "rb") as fi, gzip.open(pdb + ".gz", "wb", compresslevel=9) as fo:
    shutil.copyfileobj(fi, fo)
os.remove(pdb)

np.savetxt(f"{OUT}/orden_apertura.txt",
           np.c_[np.arange(1, len(orden) + 1), orden + 1, dist[orden]],
           fmt=["%d", "%d", "%.2f"],
           header="estado_pelicula muestra_original distancia_centros_A")

print(f"{traj.n_frames} estructuras, de {dist.min():.1f} a {dist.max():.1f} A")
print(f"escrito {pdb}.gz y {OUT}/orden_apertura.txt")
