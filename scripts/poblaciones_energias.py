#!/usr/bin/env python3
"""
Poblaciones y energias libres por conformacion del ensamble de BioEmu.

BioEmu solo entrega estructuras: muestras independientes del ensamble de
equilibrio, todas con el mismo peso. Las poblaciones salen de contar
cuantas muestras caen en cada conformacion, y las energias libres de

    dG_i = -kT ln(p_i / p_ref)        (T = 300 K, la de entrenamiento de BioEmu)

Pasos:
  1. Plegamiento de cada lobulo con la fraccion de contactos nativos (FNC),
     con la misma definicion que BioEmu usa en su entrenamiento
     (bioemu/training/foldedness.py), respecto al modelo 0 de AF3.
  2. Conformaciones = posicion de CheY respecto a LBP: se alinea todo sobre
     LBP, se hace PCA de las coordenadas Ca de CheY y k-means sobre los
     3 primeros componentes. k se fija con N_CONF.
  3. Por conformacion: poblacion, dG con IC 95% (bootstrap), descriptores
     y una estructura representativa (la mas cercana al centroide).

Uso:  python poblaciones_energias.py
"""
import os
from itertools import combinations
import numpy as np, mdtraj as md
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.expanduser("~/cursos/proteinas/quimerasLBPChey")
ENS  = f"{BASE}/bioemu/lbpchey"
AF   = os.path.expanduser("~/fold_lbpchey/fold_lbpchey_model_0.cif")
OUT  = f"{BASE}/energias"
N_CONF = 3
N_BOOT = 2000
KT = 0.0019872041 * 300.0          # kcal/mol a 300 K
FNC_PLEGADO = 0.65                 # umbral de plegado del articulo de BioEmu

LBP  = "name CA and (resid 0 to 139 or resid 243 to 338)"
CHEY = "name CA and resid 140 to 242"
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(0)

traj = md.load(f"{ENS}/samples.xtc", top=f"{ENS}/topology.pdb")
ref  = md.load(AF)
iL, iC = traj.topology.select(LBP), traj.topology.select(CHEY)
rL, rC = ref.topology.select(f"protein and {LBP}"), ref.topology.select(f"protein and {CHEY}")
rA, iA = ref.topology.select("protein and name CA"), traj.topology.select("name CA")
n = traj.n_frames
print(f"ensamble: {n} estructuras")

# ------------------------------------------------------------------ 1. plegamiento
def fnc(idx_t, idx_r):
    """Fraccion de contactos nativos Ca-Ca, como bioemu/training/foldedness.py."""
    r = ref.xyz[0, idx_r] * 10
    pares = np.array([(i, j) for i, j in combinations(range(len(r)), 2) if abs(i - j) > 3])
    d0 = np.linalg.norm(r[pares[:, 0]] - r[pares[:, 1]], axis=1)
    pares, d0 = pares[d0 <= 10.0], d0[d0 <= 10.0]
    x = traj.xyz[:, idx_t] * 10
    d = np.linalg.norm(x[:, pares[:, 0]] - x[:, pares[:, 1]], axis=2)
    return (0.5 * (1 + np.tanh(-2.5 * (d - 1.2 * d0)))).mean(1)   # = expit(-5(d - 1.2 d0))

fnc_lbp, fnc_chey = fnc(iL, rL), fnc(iC, rC)

def dg_plegamiento(f):
    plegado = f >= FNC_PLEGADO
    nf, nu = plegado.sum(), (~plegado).sum()
    if nu == 0:   # ninguna muestra desplegada: solo hay cota inferior
        return nf, nu, f"> {KT * np.log(n):.1f} (cota: 0 de {n} desplegadas)"
    return nf, nu, f"{-KT * np.log(nu / nf):.2f}"

# ------------------------------------------------------------------ 2. conformaciones
traj.superpose(ref, atom_indices=iL, ref_atom_indices=rL)
X = (traj.xyz[:, iC] * 10).reshape(n, -1)
pca = PCA(3).fit(X)
Y = pca.transform(X)
km = KMeans(N_CONF, n_init=50, random_state=0).fit(Y)
sil = silhouette_score(Y, km.labels_)

# ordenar conformaciones de mayor a menor poblacion
orden = np.argsort(-np.bincount(km.labels_, minlength=N_CONF))
lab = np.empty(n, int)
for nuevo, viejo in enumerate(orden):
    lab[km.labels_ == viejo] = nuevo
cent = km.cluster_centers_[orden]

# descriptores por estructura
com = lambda x, idx: x.xyz[:, idx].mean(1) * 10
dist = np.linalg.norm(com(traj, iL) - com(traj, iC), axis=1)
desp_chey = np.sqrt(((traj.xyz[:, iC] - ref.xyz[0, rC]) ** 2).sum(-1).mean(-1)) * 10  # CheY tras alinear LBP
rmsd_glob = md.rmsd(traj, ref, atom_indices=iA, ref_atom_indices=rA) * 10
rmsd_lbp  = md.rmsd(traj, ref, atom_indices=iL, ref_atom_indices=rL) * 10
rmsd_chey = md.rmsd(traj, ref, atom_indices=iC, ref_atom_indices=rC) * 10
traj.superpose(ref, atom_indices=iL, ref_atom_indices=rL)   # md.rmsd no mueve traj, por si acaso

# donde cae el modelo de AF3
y_af = pca.transform((ref.xyz[0, rC] * 10).reshape(1, -1))
conf_af = int(np.argmin(np.linalg.norm(cent - y_af, axis=1)))

# ------------------------------------------------------------------ 3. poblaciones y dG
cuentas = np.bincount(lab, minlength=N_CONF)
p = cuentas / n
ref_i = 0                                    # referencia: la conformacion mas poblada
dG = -KT * np.log(p / p[ref_i]) + 0.0   # + 0.0 evita "-0.00"
dG[ref_i] = 0.0
boot_p = np.array([np.bincount(rng.choice(lab, n), minlength=N_CONF) / n for _ in range(N_BOOT)])
with np.errstate(divide="ignore"):
    boot_dG = -KT * np.log(boot_p / boot_p[:, [ref_i]])
p_lo, p_hi = np.percentile(boot_p, [2.5, 97.5], axis=0)
g_lo, g_hi = np.nanpercentile(np.where(np.isfinite(boot_dG), boot_dG, np.nan), [2.5, 97.5], axis=0)

filas = []
for k in range(N_CONF):
    m = lab == k
    medoide = np.where(m)[0][np.argmin(np.linalg.norm(Y[m] - cent[k], axis=1))]
    traj[medoide].save_pdb(f"{OUT}/conformacion_{k + 1}.pdb")
    filas.append(dict(
        conf=k + 1, n=int(cuentas[k]), p=100 * p[k], p_lo=100 * p_lo[k], p_hi=100 * p_hi[k],
        dG=dG[k], g_lo=g_lo[k], g_hi=g_hi[k],
        dist=dist[m].mean(), dist_sd=dist[m].std(),
        desp=desp_chey[m].mean(), glob=rmsd_glob[m].mean(),
        lbp=rmsd_lbp[m].mean(), chey=rmsd_chey[m].mean(),
        fl=fnc_lbp[m].mean(), fc=fnc_chey[m].mean(), medoide=int(medoide) + 1))

# ------------------------------------------------------------------ salida
cols = ["conf", "n", "p", "p_lo", "p_hi", "dG", "g_lo", "g_hi", "dist", "dist_sd",
        "desp", "glob", "lbp", "chey", "fl", "fc", "medoide"]
with open(f"{OUT}/conformaciones.csv", "w") as fh:
    fh.write("conformacion,n_muestras,poblacion_pct,pob_IC95_bajo,pob_IC95_alto,dG_kcal_mol,"
             "dG_IC95_bajo,dG_IC95_alto,dist_lobulos_A,dist_lobulos_sd_A,desplaz_CheY_vs_AF3_A,"
             "RMSD_global_A,RMSD_LBP_A,RMSD_CheY_A,FNC_LBP,FNC_CheY,muestra_representativa\n")
    for f in filas:
        fh.write(",".join(f"{f[c]:.3f}" if isinstance(f[c], float) else str(f[c]) for c in cols) + "\n")
np.savetxt(f"{OUT}/asignacion_muestras.txt",
           np.c_[np.arange(1, n + 1), lab + 1, dist, desp_chey, fnc_lbp, fnc_chey],
           fmt=["%d", "%d", "%.2f", "%.2f", "%.3f", "%.3f"],
           header="muestra conformacion dist_lobulos_A desplaz_CheY_A FNC_LBP FNC_CheY")

print(f"\nT = 300 K, kT = {KT:.3f} kcal/mol, k = {N_CONF} (silueta {sil:.2f}), bootstrap {N_BOOT}")
print("\n--- Plegamiento de cada lobulo (FNC >= 0.65 respecto a AF3) ---")
for nom, f in [("LBP", fnc_lbp), ("CheY", fnc_chey)]:
    nf, nu, g = dg_plegamiento(f)
    print(f"{nom:<5} plegadas {nf}/{n}, FNC media {f.mean():.3f} (min {f.min():.2f})   dG_despl = {g} kcal/mol")

print("\n--- Conformaciones (posicion de CheY respecto a LBP) ---")
print(f"{'conf':<5}{'n':>4}{'poblacion %':>22}{'dG kcal/mol':>24}{'dist A':>14}{'desp CheY A':>13}{'RMSD glob':>11}")
for f in filas:
    print(f"{f['conf']:<5}{f['n']:>4}{f['p']:>9.1f} [{f['p_lo']:4.1f}-{f['p_hi']:4.1f}]"
          f"{f['dG']:>11.2f} [{f['g_lo']:4.2f}-{f['g_hi']:4.2f}]"
          f"{f['dist']:>8.1f} ± {f['dist_sd']:.1f}{f['desp']:>11.1f}{f['glob']:>11.1f}")
print(f"\nel modelo de AlphaFold3 cae en la conformacion {conf_af + 1}")
print(f"resolucion: 1 muestra = {100 / n:.2f}% ; dG maximo medible = {KT * np.log(n):.1f} kcal/mol")

# ------------------------------------------------------------------ figura
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
H, xe, ye = np.histogram2d(Y[:, 0], Y[:, 1], bins=25)
with np.errstate(divide="ignore"):
    F = -KT * np.log(H.T / H.max())
im = ax[0].imshow(np.ma.masked_invalid(np.where(np.isfinite(F), F, np.nan)), origin="lower",
                  extent=[xe[0], xe[-1], ye[0], ye[-1]], aspect="auto", cmap="viridis")
colores = ["#E8792B", "#2A9D8F", "#5B6B78", "#9B5DE5", "#F15BB5"]
for k in range(N_CONF):
    ax[0].scatter(Y[lab == k, 0], Y[lab == k, 1], s=8, c=colores[k], edgecolor="none", alpha=0.8,
                  label=f"conf {k + 1} ({100 * p[k]:.0f}%)")
ax[0].scatter(*y_af[0, :2], marker="*", s=250, c="white", edgecolor="black", label="AF3 modelo 0")
ax[0].set_xlabel(f"PC1 ({100 * pca.explained_variance_ratio_[0]:.0f}%)")
ax[0].set_ylabel(f"PC2 ({100 * pca.explained_variance_ratio_[1]:.0f}%)")
ax[0].set_title("Energia libre de la posicion de CheY (kcal/mol)")
plt.colorbar(im, ax=ax[0]); ax[0].legend(fontsize=8, loc="best")

x = np.arange(N_CONF)
ax[1].bar(x, [f["dG"] for f in filas], color=colores[:N_CONF],
          yerr=[[f["dG"] - f["g_lo"] for f in filas], [f["g_hi"] - f["dG"] for f in filas]], capsize=6)
for i, f in enumerate(filas):
    ax[1].text(i, f["g_hi"] + 0.05, f"{f['p']:.0f}%", ha="center")
ax[1].set_ylim(0, max(f["g_hi"] for f in filas) * 1.18)
ax[1].set_xticks(x, [f"conf {k + 1}" for k in range(N_CONF)])
ax[1].set_ylabel("dG respecto a conf 1 (kcal/mol)")
ax[1].set_title("Energia libre relativa por conformacion (IC 95%)")
plt.tight_layout(); plt.savefig(f"{OUT}/energia_libre.png", dpi=200)
print(f"\nresultados en {OUT}/")
