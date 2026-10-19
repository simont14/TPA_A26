# -*- coding: utf-8 -*-
"""
Sismique refraction sur le Grand Axe: temps de premiere arrivee aux trois
geophones pour chaque position de la plaque, ajustement d'une droite brisee
(onde directe + onde refractee) et calcul de l'epaisseur de la couche de surface.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from signaux_geo import BASE, DISTANCES, GEOPHONES, temps_arrivee, trace_moyenne

OUT_DIR = os.environ.get("ANALYSE_OUT_DIR", os.path.join(BASE, "figures"))
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    "font.size": 11,
    "figure.figsize": (8, 5),
    "axes.grid": True,
    "grid.alpha": 0.3,
})

COULEURS = {3: "C0", 1: "C1", 2: "C2"}


def fr(x, fmt="g"):
    """Virgule decimale (convention francaise) pour les etiquettes."""
    return format(x, fmt).replace(".", ",")


def use_fr_ticks(ax):
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos=None: fr(x)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos=None: fr(x)))


def droite_brisee(x, s1, s2, xc):
    """Temps de premiere arrivee: droite de pente s1 partant de l'origine jusqu'a
    xc, puis pente s2 (ms, m)."""
    return np.where(x < xc, s1 * x, s1 * xc + s2 * (x - xc))


# ---------------------------------------------------------------- reperage
mesures = []      # (position, canal, x, t ou None)
traces = {}
for p, d in enumerate(DISTANCES):
    for canal, (nom, x0) in GEOPHONES.items():
        t, z = trace_moyenne(p, canal)
        traces[(p, canal)] = (t, z)
        mesures.append((p, canal, d + x0, temps_arrivee(t, z)))

X = np.array([m[2] for m in mesures if m[3] is not None])
T = np.array([m[3] for m in mesures if m[3] is not None])
C = np.array([m[1] for m in mesures if m[3] is not None])

# ---------------------------------------------------------------- ajustement
popt, pcov = curve_fit(droite_brisee, X, T, p0=[2.0, 1.0, 2.0])
perr = np.sqrt(np.diag(pcov))
s1, s2, xc = popt
ds1, ds2, dxc = perr
residus = T - droite_brisee(X, *popt)
rms = np.sqrt(np.mean(residus**2))
r2 = 1 - np.sum(residus**2) / np.sum((T - T.mean())**2)

V1, V2 = 1e3 / s1, 1e3 / s2               # m/s (pentes en ms/m)
dV1, dV2 = 1e3 * ds1 / s1**2, 1e3 * ds2 / s2**2

e = xc / 2 * np.sqrt((V2 - V1) / (V2 + V1))
de_dxc = e / xc
de_dV1 = -e * V2 / (V2**2 - V1**2)
de_dV2 = e * V1 / (V2**2 - V1**2)
de = np.sqrt((de_dxc * dxc)**2 + (de_dV1 * dV1)**2 + (de_dV2 * dV2)**2)

# Profondeur minimale du roc: une onde refractee sur un roc de vitesse V3 a une
# profondeur h1 + h2 arriverait au temps
#   t = x/V3 + 2 h1 sqrt(1/V1^2 - 1/V3^2) + 2 h2 sqrt(1/V2^2 - 1/V3^2).
# Elle n'est pas observee avant le dernier geophone, donc ce temps est plus grand
# que le temps mesure au point le plus loin.
x_max = X.max()
t_max = droite_brisee(x_max, *popt)
prof_min = {}
for V3 in (2000.0, 4000.0):
    a1 = 2 * np.sqrt(1 / V1**2 - 1 / V3**2) * 1e3   # ms/m
    a2 = 2 * np.sqrt(1 / V2**2 - 1 / V3**2) * 1e3
    h2 = (t_max - x_max / V3 * 1e3 - a1 * e) / a2
    prof_min[V3] = e + h2

print("Temps d'arrivee (ms):")
for p, canal, x, t in mesures:
    print(f"  d={DISTANCES[p]:.1f} m  {GEOPHONES[canal][0]}  x={x:.2f} m  t="
          + ("rejete" if t is None else f"{t:.2f}"))
print(f"N = {len(X)} points, RMS = {rms:.2f} ms, R2 = {r2:.4f}")
print(f"V1 = {V1:.0f} +- {dV1:.0f} m/s, V2 = {V2:.0f} +- {dV2:.0f} m/s")
print(f"Dc = {xc:.2f} +- {dxc:.2f} m")
print(f"e = {e:.3f} +- {de:.3f} m (relatif {de/e:.3f}, Dc {dxc/xc:.3f})")
print(f"x_max = {x_max:.2f} m, t(x_max) = {t_max:.2f} ms")
for V3, h in prof_min.items():
    print(f"profondeur minimale du roc (V3 = {V3:.0f} m/s): {h:.3f} m")

# ---------------------------------------------------------------- figure: section
fig, ax = plt.subplots(figsize=(8, 6))
echelle = 0.25    # demi-amplitude d'une trace, en m de distance
for p, d in enumerate(DISTANCES):
    t, z = traces[(p, 2)]
    x = d + GEOPHONES[2][1]
    fen = (t > 0) & (t < 25)
    z_aff = z / np.abs(z[fen]).max() * echelle
    ax.plot(t, x + z_aff, color="0.25", lw=0.8)
    t_p = dict(((m[0], m[1]), m[3]) for m in mesures)[(p, 2)]
    if t_p is not None:
        ax.plot(t_p, x, "o", color="C3", ms=5, zorder=5)
x_mod = np.linspace(0, X.max() + 0.3, 200)
ax.plot(droite_brisee(x_mod, *popt), x_mod, "--", color="C3", lw=1.2,
        label="Ajustement des temps d'arrivée")
ax.set_xlim(-1, 25)
ax.set_ylim(-0.3, 6.0)
ax.set_xlabel("Temps après l'impact (ms)")
ax.set_ylabel("Distance plaque-géophone (m)")
ax.legend(loc="lower right")
use_fr_ticks(ax)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "01_section_sismique.png"), dpi=200)
plt.close(fig)

# ---------------------------------------------------------------- figure: dromochronique
fig, ax = plt.subplots(figsize=(8, 5))
for canal, (nom, _) in GEOPHONES.items():
    sel = C == canal
    ax.errorbar(X[sel], T[sel], yerr=rms, fmt="o", ms=5, capsize=3,
                color=COULEURS[canal], label=f"Géophone {nom[1]}")
x_mod = np.linspace(0, X.max() + 0.2, 300)
ax.plot(x_mod, droite_brisee(x_mod, *popt), "k-", lw=1.3, label="Ajustement")
ax.axvline(xc, color="0.4", ls="--", lw=1)
ax.text(xc + 0.08, 1.0, f"$D_c$ = {fr(xc, '.2f')} m", color="0.3")
ax.text(0.12, 3.7, f"$V_1$ = {fr(V1, '.0f')} m/s")
ax.text(3.55, 5.9, f"$V_2$ = {fr(V2, '.0f')} m/s")
ax.set_xlim(0, 5.6)
ax.set_ylim(0, 11.5)
ax.set_xlabel("Distance plaque-géophone (m)")
ax.set_ylabel("Temps de première arrivée (ms)")
ax.legend(loc="upper left")
use_fr_ticks(ax)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "02_dromochronique.png"), dpi=200)
plt.close(fig)
