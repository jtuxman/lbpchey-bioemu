#!/usr/bin/env python3
"""
Analisis del ensamble de BioEmu para la quimera LBPCheY.

Pregunta central: los dos lobulos plegan bien por separado, pero
se mueven uno respecto al otro?  Eso es lo que el articulo observo
en el cristal (orientacion distinta a la del diseno) y lo que el
DSC sugiere (desplegamiento en dos etapas).
"""
import os, numpy as np, mdtraj as md
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.expanduser("~/cursos/proteinas/quimerasLBPChey")
ENS  = f"{BASE}/bioemu/lbpchey"
AF   = os.path.expanduser("~/fold_lbpchey/fold_lbpchey_model_0.cif")

# fronteras de la quimera (1-indexado, como el FASTA)
LBP_N  = (1, 140)
CHEY   = (141, 243)
LBP_C  = (244, 339)

# ---------------------------------------------------------------
traj = md.load(f"{ENS}/samples.xtc", top=f"{ENS}/topology.pdb")
print(f"ensamble: {traj.n_frames} estructuras, {traj.n_residues} residuos")

def ca(t, lo, hi):
    return t.topology.select(f"protein and name CA and resid {lo-1} to {hi-1}")

sel_lbp  = np.concatenate([ca(traj, *LBP_N), ca(traj, *LBP_C)])
sel_chey = ca(traj, *CHEY)

# referencia = el modelo de AlphaFold, si mdtraj puede leer el cif;
# si no, el primer frame del ensamble
try:
    ref = md.load(AF)
    print("referencia: modelo de AlphaFold")
except Exception:
    ref = traj[0]
    print("referencia: primer frame del ensamble (no pude leer el .cif)")

# la referencia tiene todos los atomos: sus indices de CA son otros
ref_lbp  = np.concatenate([ca(ref, *LBP_N), ca(ref, *LBP_C)])
ref_chey = ca(ref, *CHEY)
ref_all  = ref.topology.select("protein and name CA")

# ---------------------------------------------------------------
# 1. Cada lobulo por separado: plegaron bien?
# ---------------------------------------------------------------
rmsd_lbp  = md.rmsd(traj, ref, atom_indices=sel_lbp,  ref_atom_indices=ref_lbp)  * 10
rmsd_chey = md.rmsd(traj, ref, atom_indices=sel_chey, ref_atom_indices=ref_chey) * 10

# 2. Todo junto: el ensamble global se mantiene?
sel_all   = traj.topology.select("protein and name CA")
rmsd_glob = md.rmsd(traj, ref, atom_indices=sel_all, ref_atom_indices=ref_all) * 10

# 3. Distancia entre centros de masa de los dos lobulos
#    (si varia mucho, los lobulos se abren y cierran)
com_lbp  = md.compute_center_of_mass(traj.atom_slice(sel_lbp))
com_chey = md.compute_center_of_mass(traj.atom_slice(sel_chey))
dist = np.linalg.norm(com_lbp - com_chey, axis=1) * 10

# ---------------------------------------------------------------
def resumen(nombre, v, unidad="A"):
    print(f"{nombre:<34} {v.mean():6.2f} +/- {v.std():5.2f} {unidad}"
          f"   [{v.min():.2f} - {v.max():.2f}]")

print("\n--- RMSD respecto a la referencia ---")
resumen("lobulo LBP  (1-140, 244-339)", rmsd_lbp)
resumen("lobulo CheY (141-243)",        rmsd_chey)
resumen("estructura completa",          rmsd_glob)
print("\n--- Orientacion relativa ---")
resumen("distancia entre centros",      dist)

print("\nLectura: si los dos primeros son bajos y el tercero alto,")
print("cada fragmento pliega bien pero el ensamble global es flexible.")

# ---------------------------------------------------------------
# graficas
# ---------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for a, (v, t) in zip(ax, [(rmsd_lbp,  "RMSD lobulo LBP"),
                          (rmsd_chey, "RMSD lobulo CheY"),
                          (rmsd_glob, "RMSD global")]):
    a.hist(v, bins=40, color="#4a7c59", edgecolor="white")
    a.set_xlabel("RMSD (A)"); a.set_title(t); a.set_ylabel("estructuras")
plt.tight_layout(); plt.savefig(f"{BASE}/rmsd_lobulos.png", dpi=200)

plt.figure(figsize=(6, 4))
plt.hist(dist, bins=40, color="#7b4397", edgecolor="white")
plt.xlabel("distancia entre centros de los lobulos (A)")
plt.ylabel("estructuras"); plt.tight_layout()
plt.savefig(f"{BASE}/apertura_lobulos.png", dpi=200)

# ---------------------------------------------------------------
# estructuras representativas: la mas cerrada, la mediana, la mas abierta
# ---------------------------------------------------------------
orden = np.argsort(dist)
for etiqueta, idx in [("cerrada",  orden[0]),
                      ("mediana",  orden[len(orden)//2]),
                      ("abierta",  orden[-1])]:
    traj[idx].save_pdb(f"{BASE}/repr_{etiqueta}.pdb")
    print(f"guardado repr_{etiqueta}.pdb  (distancia {dist[idx]:.1f} A)")

print(f"\ngraficas en {BASE}/rmsd_lobulos.png y apertura_lobulos.png")
