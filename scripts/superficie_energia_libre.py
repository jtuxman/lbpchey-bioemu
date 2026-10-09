#!/usr/bin/env python3
"""
Superficies de energia libre del ensamble de BioEmu, al estilo de la Fig. 3
del articulo de BioEmu (Lewis et al., Science 2025): F = -kT ln p, en kcal/mol,
con curvas de nivel.

En el articulo los ejes son las dos componentes TICA mas lentas, construidas
con simulaciones de MD. Aqui no hay MD (y TICA requiere series de tiempo), asi
que se usan dos alternativas:

  a) PC1-PC2 de las coordenadas Ca de CheY tras alinear sobre LBP
     (la misma proyeccion que poblaciones_energias.py)
  b) variables fisicas: distancia entre los centros de los lobulos y
     angulo de rotacion de CheY respecto al modelo de AF3

La densidad se estima con un kernel gaussiano (KDE) porque 246 muestras son
pocas para un histograma 2D. Las regiones con F > F_MAX quedan en blanco.

Uso:  python superficie_energia_libre.py
"""
import os
import numpy as np, mdtraj as md
from scipy.stats import gaussian_kde
from scipy.spatial.transform import Rotation
from sklearn.decomposition import PCA
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.expanduser("~/cursos/proteinas/quimerasLBPChey")
ENS  = f"{BASE}/bioemu/lbpchey"
AF   = os.path.expanduser("~/fold_lbpchey/fold_lbpchey_model_0.cif")
OUT  = f"{BASE}/energias"
KT    = 0.0019872041 * 300.0       # kcal/mol a 300 K
F_MAX = 4.0                        # kcal/mol; por encima no hay muestras suficientes
NIVELES = np.arange(0, F_MAX + 0.01, 0.5)

LBP  = "name CA and (resid 0 to 139 or resid 243 to 338)"
CHEY = "name CA and resid 140 to 242"
os.makedirs(OUT, exist_ok=True)

traj = md.load(f"{ENS}/samples.xtc", top=f"{ENS}/topology.pdb")
ref  = md.load(AF)
iL, iC = traj.topology.select(LBP), traj.topology.select(CHEY)
rL, rC = ref.topology.select(f"protein and {LBP}"), ref.topology.select(f"protein and {CHEY}")
traj.superpose(ref, atom_indices=iL, ref_atom_indices=rL)

# conformacion de cada muestra (de poblaciones_energias.py), si existe
f_asig = f"{OUT}/asignacion_muestras.txt"
conf = np.loadtxt(f_asig, dtype=int, usecols=1) if os.path.exists(f_asig) else None

# ---------------------------------------------------------------- variables colectivas
X = (traj.xyz[:, iC] * 10).reshape(traj.n_frames, -1)
x_af = (ref.xyz[0, rC] * 10).reshape(1, -1)
pca = PCA(2).fit(X)
pc, pc_af = pca.transform(X), pca.transform(x_af)[0]

def centro(x, idx):
    return x.xyz[:, idx].mean(1) * 10

dist = np.linalg.norm(centro(traj, iL) - centro(traj, iC), axis=1)
dist_af = np.linalg.norm(centro(ref, rL) - centro(ref, rC), axis=1)[0]

# giro de CheY respecto a AF3 (con LBP ya alineado): Kabsch sobre los Ca de CheY.
# La magnitud de una rotacion es siempre >= 0 y su distribucion tiene un factor
# geometrico (hay mas formas de girar 20 grados que 0), asi que se usa un angulo
# CON SIGNO: la proyeccion del vector de rotacion sobre su eje principal (PCA).
A = ref.xyz[0, rC] * 10
A0 = A - A.mean(0)
rotvec = np.empty((traj.n_frames, 3))
for k in range(traj.n_frames):
    B = traj.xyz[k, iC] * 10
    rot, _ = Rotation.align_vectors(B - B.mean(0), A0)      # AF3 -> muestra
    rotvec[k] = rot.as_rotvec()
eje = PCA(1).fit(rotvec)
ang = np.degrees(rotvec @ eje.components_[0])               # 0 = orientacion de AF3
if np.median(ang) < 0:
    ang = -ang
var_eje = 100 * eje.explained_variance_ratio_[0]
magn = np.degrees(np.linalg.norm(rotvec, axis=1))

# ---------------------------------------------------------------- superficie
def superficie(ax, x, y, x_af, y_af, xlabel, ylabel, titulo, leyenda="lower left"):
    kde = gaussian_kde(np.vstack([x, y]))
    px, py = 0.15 * np.ptp(x), 0.15 * np.ptp(y)
    gx, gy = np.meshgrid(np.linspace(x.min() - px, x.max() + px, 200),
                         np.linspace(y.min() - py, y.max() + py, 200))
    dens = kde(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)
    F = -KT * np.log(dens / dens.max())
    F = np.ma.masked_greater(F, F_MAX)
    cf = ax.contourf(gx, gy, F, levels=NIVELES, cmap="viridis")
    ax.contour(gx, gy, F, levels=NIVELES, colors="k", linewidths=0.4, alpha=0.6)
    if conf is not None:
        colores = {1: "#E8792B", 2: "#FFFFFF", 3: "#E63946"}
        for c in sorted(set(conf)):
            m = conf == c
            ax.scatter(x[m], y[m], s=7, c=colores.get(c, "grey"), edgecolor="k", linewidths=0.2,
                       label=f"conformación {c}", zorder=3)
    else:
        ax.scatter(x, y, s=5, c="white", edgecolor="k", linewidths=0.2, zorder=3)
    ax.scatter(x_af, y_af, marker="*", s=320, c="gold", edgecolor="k", zorder=4,
               label="AlphaFold3 (modelo 0)")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(titulo)
    ax.legend(fontsize=8, loc=leyenda, framealpha=0.85)
    return cf

fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))
v = pca.explained_variance_ratio_ * 100
superficie(ax[0], pc[:, 0], pc[:, 1], pc_af[0], pc_af[1],
           f"PC1 ({v[0]:.0f}%)", f"PC2 ({v[1]:.0f}%)",
           "a) Posición de CheY: componentes principales", leyenda="upper right")
cf = superficie(ax[1], dist, ang, dist_af, 0.0,
                "Distancia entre centros de los lóbulos (Å)", f"Giro de CheY respecto a AF3 (°, eje principal, {var_eje:.0f}%)",
                "b) Variables físicas")
cb = fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.02)
cb.set_label("Energía libre (kcal/mol)")
fig.suptitle(f"Superficie de energía libre de LBPCheY (BioEmu, {traj.n_frames} muestras, 300 K)",
             fontsize=13)
plt.savefig(f"{OUT}/superficie_energia_libre.png", dpi=200, bbox_inches="tight")

# perfiles 1D
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
for a, x, x_af, lab in [(ax[0], dist, dist_af, "Distancia entre centros de los lóbulos (Å)"),
                        (ax[1], ang, 0.0, "Giro de CheY respecto a AF3 (°, eje principal)")]:
    kde = gaussian_kde(x)
    g = np.linspace(x.min() - 0.1 * np.ptp(x), x.max() + 0.1 * np.ptp(x), 400)
    F = -KT * np.log(kde(g) / kde(g).max())
    m = F <= F_MAX
    a.plot(g[m], F[m], color="#2A9D8F", lw=2.5)
    a.fill_between(g[m], F[m], F_MAX, color="#2A9D8F", alpha=0.12)
    a.axvline(x_af, color="goldenrod", ls="--", lw=1.5, label="AlphaFold3")
    a.set_xlabel(lab); a.set_ylabel("Energía libre (kcal/mol)"); a.set_ylim(0, F_MAX)
    a.legend(fontsize=9)
fig.suptitle("Perfiles de energía libre 1D", fontsize=12)
plt.tight_layout(); plt.savefig(f"{OUT}/perfiles_energia_libre.png", dpi=200)

print(f"distancia AF3 {dist_af:.1f} A | ensamble {dist.mean():.1f} ± {dist.std():.1f} A")
print(f"giro de CheY sobre el eje principal ({var_eje:.0f}% de la varianza): "
      f"mediana {np.median(ang):.1f}°, rango {ang.min():.0f} a {ang.max():.0f}°")
print(f"magnitud total de la rotacion: mediana {np.median(magn):.1f}°, max {magn.max():.0f}°")
if conf is not None:
    for c in sorted(set(conf)):
        m = conf == c
        print(f"conformacion {c}: giro {ang[m].mean():.1f} ± {ang[m].std():.1f}°, "
              f"magnitud {np.median(magn[m]):.1f}° (mediana)")
print(f"figuras en {OUT}/superficie_energia_libre.png y perfiles_energia_libre.png")
