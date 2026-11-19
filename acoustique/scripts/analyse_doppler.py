# -*- coding: utf-8 -*-
"""
Effet Doppler d'une voiture qui passe devant le microphone en jouant un sinus
de 440 Hz. Produit le tableau de resultats et les figures du rapport.

Pour chaque passage, v et d sont mesures a partir des trois regions choisies
dans fit_manuel.py (regions_manuelles.json): les queues donnent v, le centre
donne d (methode decrite dans doppler_commun.py). Un passage absent du fichier
prend les regions de depart automatiques.
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from doppler_commun import (BASE, C, DC, DT_C, F0, T_C, analyse_regions, charger_passages,
                            lire_regions, modele, regions_depart, vitesses_queues)

OUT_DIR = os.environ.get("ANALYSE_OUT_DIR", os.path.join(BASE, "figures"))
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
VIRG = FuncFormatter(lambda x, _: f"{x:g}".replace(".", ",").replace("-", "−"))

# ---------------------------------------------------------------- analyse
manuelles = lire_regions()
passages = charger_passages()
for p in passages:
    R, tp0 = manuelles.get(p["nom"]), None
    p["source"] = "manuel" if R else "auto"
    if R is None:
        R, tp0 = regions_depart(p)
    r = analyse_regions(p, R, tp0=tp0)
    if r is None:
        raise SystemExit(f"{p['nom']}: pas assez de points dans une des regions")
    p.update(R=R, r=r, pA=np.array([r["v"], r["d"], r["tp"]]), f0=r["f0"], evA=r["ev"], edA=r["ed"],
             selA=r["sel"], vq=vitesses_queues(p, R, r))

print(f"c = {C:.1f} +- {DC:.1f} m/s  (T = {T_C} +- {DT_C} degC)")
print(f"{'source':>6s} {'v visee':>7s} {'d galon':>7s} | {'v mesuree':>14s} {'d mesuree':>14s} "
      f"{'f0':>6s} {'rms q/c':>10s} {'n a/c/e':>14s} conv | v approche / eloignement (f0 = 440 Hz)")
for p in passages:
    r = p["r"]
    print(f"{p['source']:>6s} {p['vn']:7d} {p['d']:7.2f} | {r['v']*3.6:6.1f} +- {r['ev']*3.6:5.1f} "
          f"{r['d']:6.2f} +- {r['ed']:5.2f} {r['f0']:6.1f} {r['rms_q']:4.2f}/{r['rms_c']:4.2f} "
          f"{r['na']:4d}/{r['nc']:4d}/{r['ne']:4d} {r['conv']} | {p['vq'][0]*3.6:5.1f} / {p['vq'][1]*3.6:5.1f}")

# ---------------------------------------------------------------- figure 1: spectrogrammes
# Une fenetre de ~0,37 s donne une raie fine a 50 km/h sans trop etaler la chute.
def spectrogramme(ax, p, demi=3.5, n=16384, hop=512):
    x = p["x"].astype(float)
    fs, tp = p["fs"], p["pA"][2]
    i0 = max(0, int((tp - demi) * fs) - n // 2)
    i1 = min(len(x), int((tp + demi) * fs) + n // 2)
    seg = x[i0:i1] - x[i0:i1].mean()
    deb = np.arange(0, len(seg) - n, hop)
    tr = np.lib.stride_tricks.sliding_window_view(seg, n)[deb] * np.hanning(n)
    S = 20 * np.log10(np.abs(np.fft.rfft(tr, 4 * n, axis=1)) + 1e-9)
    fr = np.fft.rfftfreq(4 * n, 1 / fs)
    m = (fr > 420) & (fr < 465)  # centre sur la raie a 50 km/h
    tt = (i0 + deb + n / 2) / fs - tp
    S = S[:, m] - S[:, m].max()
    im = ax.pcolormesh(tt, fr[m], S.T, shading="gouraud", cmap="magma", vmin=-45, vmax=0, rasterized=True)
    ax.set_xlim(-demi, demi)
    ax.set_xlabel("Temps par rapport au passage (s)")
    ax.set_ylabel("Fréquence (Hz)")
    ax.grid(False)
    ax.xaxis.set_major_formatter(VIRG)
    return im


def trouve(vn, d):
    return next(p for p in passages if p["vn"] == vn and abs(p["d"] - d) < 0.01)


fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
for ax, d in zip(axs, (3.53, 8.97)):
    im = spectrogramme(ax, trouve(50, d))
    ax.text(0.97, 0.94, f"50 km/h, $d$ = {d:.2f} m".replace(".", ","), transform=ax.transAxes,
            ha="right", va="top", color="w", fontsize=10)
axs[1].set_ylabel("")
fig.colorbar(im, ax=axs, label="Niveau relatif (dB)", pad=0.02)
fig.savefig(os.path.join(OUT_DIR, "01_spectrogrammes.png"), dpi=200, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- figure 2: courbes f(t) et fits
# Une case par passage (ligne = distance, colonne = vitesse visee): lisible en noir et blanc.
vitesses = sorted({p["vn"] for p in passages})
distances = (3.53, 8.97)
fig, axs = plt.subplots(len(distances), len(vitesses), figsize=(10, 4.6), sharex=True, sharey=True)
for i, d in enumerate(distances):
    for j, vn in enumerate(vitesses):
        ax, p = axs[i, j], trouve(vn, d)
        v, dd, tp = p["pA"]
        s = p["selA"]
        ax.plot(p["t"][s] - tp, p["f"][s], ".", ms=1.5, color="0.6")
        tt = np.linspace(p["t"][s].min(), p["t"][s].max(), 400)
        ax.plot(tt - tp, modele(tt, v, dd, tp, p["f0"]), color="k", lw=1.4)
        ax.axhline(F0, color="k", ls=":", lw=0.8)
        ax.text(0.95, 0.93, f"$v$ = {vn} km/h", transform=ax.transAxes, ha="right", va="top", fontsize=9)
        ax.text(0.05, 0.05, f"$d$ = {d:.2f} m".replace(".", ","), transform=ax.transAxes, fontsize=9)
        ax.tick_params(labelsize=9)
axs[0, 0].set_xlim(-3.5, 3.5)
axs[0, 0].set_ylim(414, 474)
axs[0, 0].xaxis.set_major_formatter(VIRG)
fig.supxlabel("Temps par rapport au passage (s)", fontsize=11)
fig.supylabel("Fréquence reçue (Hz)", fontsize=11)
fig.tight_layout(w_pad=0.4, h_pad=0.4)
fig.savefig(os.path.join(OUT_DIR, "02_courbes_doppler.png"), dpi=200)
plt.close(fig)

# ---------------------------------------------------------------- figure 3: vitesse mesuree vs visee
fig, ax = plt.subplots(figsize=(5.6, 4.6))
# Marqueur plein contre marqueur vide: les deux distances restent distinctes en noir et blanc.
for d, mk, c, fc in ((3.53, "o", "C0", "C0"), (8.97, "s", "C1", "w")):
    ps = [p for p in passages if abs(p["d"] - d) < 0.01]
    ax.errorbar([p["vn"] for p in ps], [p["pA"][0] * 3.6 for p in ps], yerr=[p["evA"] * 3.6 for p in ps],
                fmt=mk, color=c, mfc=fc, mew=1.5, capsize=3, ms=6, label=f"$d$ = {d:.2f} m".replace(".", ","))
vv = np.array([0, 75])
ax.plot(vv, vv, color="0.4", ls="--", lw=1, label="$v$ mesurée = $v$ visée")
ax.set_xlim(0, 75)
ax.set_ylim(0, 75)
ax.set_aspect("equal")
ax.set_xlabel("Vitesse visée au compteur (km/h)")
ax.set_ylabel("Vitesse mesurée par effet Doppler (km/h)")
ax.xaxis.set_major_formatter(VIRG)
ax.yaxis.set_major_formatter(VIRG)
ax.legend(fontsize=9, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "03_vitesses.png"), dpi=200)
plt.close(fig)
