# -*- coding: utf-8 -*-
"""
Effet Doppler d'une voiture qui passe devant le sonometre en jouant un sinus
de 440 Hz. Pour chaque passage: suivi de la frequence du pic dans le
spectrogramme, puis curve fit du modele Doppler (source en ligne droite a
vitesse constante, observateur fixe, delai de propagation inclus).

Deux fits par passage, sur la zone ou la vitesse est constante:
  A) d fixee a la distance mesuree au ruban  -> vitesse v
  B) d libre                                  -> v et d a partir de l'audio seul

Les noms de fichiers donnent la vitesse visee (km/h) et la distance (cm).
"""
import glob
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.io import wavfile
from scipy.ndimage import median_filter
from scipy.optimize import least_squares

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
OUT_DIR = os.environ.get("ANALYSE_OUT_DIR", os.path.join(BASE, "figures"))
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
VIRG = FuncFormatter(lambda x, _: f"{x:g}".replace(".", ",").replace("-", "−"))

# ---------------------------------------------------------------- constantes
T_C, DT_C = 1.0, 1.0                     # temperature de l'air [degC]
C = 331.3 * np.sqrt(1 + T_C / 273.15)    # vitesse du son [m/s]
DC = 0.606 * DT_C                        # incertitude sur c [m/s]
F0 = 440.0                               # frequence jouee [Hz]
L_ZONE = 10.0                            # zone a vitesse constante: |x| <= 10 m
NFFT, HOP, PAD = 4096, 256, 4            # trame de 93 ms, pas de 5,8 ms
BANDE = (400.0, 485.0)                   # recherche du pic [Hz] (evite la raie fixe a 385 Hz)
SNR_MIN = 15.0                           # pic / mediane de la bande [dB]
CORR = np.sqrt(NFFT / HOP)               # 16 trames successives partagent les memes donnees


# ---------------------------------------------------------------- suivi du pic
def suivi(x, fs):
    x = x.astype(float) - x.mean()
    w = np.hanning(NFFT)
    n = NFFT * PAD
    fr = np.fft.rfftfreq(n, 1 / fs)
    m = np.where((fr > BANDE[0]) & (fr < BANDE[1]))[0]
    debuts = np.arange(0, len(x) - NFFT, HOP)
    tr = np.lib.stride_tricks.sliding_window_view(x, NFFT)[debuts]
    X = np.abs(np.fft.rfft(tr * w, n, axis=1))[:, m[0] - 1:m[-1] + 2]
    sub = X[:, 1:-1]
    k = np.argmax(sub, 1)
    i, r = k + 1, np.arange(len(k))
    a, b, c = (np.log(X[r, i + j] + 1e-12) for j in (-1, 0, 1))
    dp = 0.5 * (a - c) / (a - 2 * b + c)            # interpolation parabolique
    f = fr[m[0] - 1 + i] + dp * fs / n
    snr = 20 * np.log10(sub[r, k] / np.median(sub, 1))
    t = (debuts + NFFT / 2) / fs
    g = snr > SNR_MIN
    return t[g], f[g]


# ---------------------------------------------------------------- modele
def u_emission(t, v, d, tp):
    """Temps d'emission (relatif au passage tp) du son recu au temps t."""
    tau = t - tp
    a = C * C - v * v
    return (C * C * tau - np.sqrt(C**4 * tau**2 - a * (C * C * tau**2 - d * d))) / a


def modele(t, v, d, tp, f0=F0):
    u = u_emission(t, v, d, tp)
    x = v * u
    return f0 / (1 + v * x / (C * np.sqrt(x * x + d * d)))


def zone_monotone(t, f, tp):
    """A vitesse constante, f(t) ne fait que diminuer: on garde la partie entre
    le maximum avant le passage et le minimum apres (tolerance 1,5 Hz)."""
    lis = median_filter(f, size=int(0.3 / np.median(np.diff(t))) | 1, mode="nearest")
    ia = np.argmax(np.where(t < tp, lis, -np.inf))
    i0 = ia
    while i0 > 0 and lis[i0 - 1] >= lis[ia] - 1.5:
        i0 -= 1
    ie = np.argmin(np.where(t > tp, lis, np.inf))
    i1 = ie
    while i1 < len(t) - 1 and lis[i1 + 1] <= lis[ie] + 1.5:
        i1 += 1
    return (t >= t[i0]) & (t <= t[i1])


def ajuste(t, f, mono, p, d_fixe=None):
    """Curve fit sur la zone |x| <= L_ZONE (recalculee avec les parametres)."""
    for _ in range(40):
        sel = mono & (np.abs(p[0] * u_emission(t, *p)) <= L_ZONE)
        if d_fixe is None:
            res = lambda q: modele(t[sel], *q) - f[sel]
            q0, lb, ub = p, [0.3, 0.3, t.min()], [40, 30, t.max()]
        else:
            res = lambda q: modele(t[sel], q[0], d_fixe, q[1]) - f[sel]
            q0, lb, ub = p[[0, 2]], [0.3, t.min()], [40, t.max()]
        r = least_squares(res, q0, loss="soft_l1", f_scale=1.0, bounds=(lb, ub))
        pn = r.x if d_fixe is None else np.array([r.x[0], d_fixe, r.x[1]])
        fini = np.allclose(pn, p, rtol=1e-5, atol=1e-5)
        p = pn
        if fini:
            break
    cov = np.linalg.inv(r.jac.T @ r.jac) * np.sum(r.fun**2) / (sel.sum() - len(q0))
    e = np.sqrt(np.diag(cov)) * CORR
    if d_fixe is not None:
        e = np.array([e[0], 0.0, e[1]])
    rms = np.sqrt(np.mean(r.fun**2))
    return p, e, sel, rms


# ---------------------------------------------------------------- analyse
passages = []
for chemin in sorted(glob.glob(os.path.join(DATA, "*.wav"))):
    nom = os.path.basename(chemin)
    vn, dcm = (int(s) for s in re.findall(r"_(\d+)(?:kmh|cm)", nom))
    fs, x = wavfile.read(chemin)
    t, f = suivi(x, fs)
    # passage approximatif: pente la plus negative de f(t) lissee sur 0,5 s, loin des bords
    lis = median_filter(f, size=int(0.5 / np.median(np.diff(t))) | 1, mode="nearest")
    pente = np.gradient(lis, t)
    pente[(t < t.min() + 1.5) | (t > t.max() - 1.5)] = 0
    tp0 = t[np.argmin(pente)]
    mono = zone_monotone(t, f, tp0)
    d = dcm / 100
    pA, eA, selA, rmsA = ajuste(t, f, mono, np.array([vn / 3.6, d, tp0]), d_fixe=d)
    pB, eB, selB, rmsB = ajuste(t, f, mono, pA.copy())
    pA2, _, _, _ = ajuste(t, f, mono, pA.copy(), d_fixe=d + 0.9)   # source au centre de la voiture
    evA = np.hypot(eA[0], pA[0] * DC / C)        # propagation: fit et vitesse du son
    passages.append(dict(nom=nom, vn=vn, d=d, fs=fs, x=x, t=t, f=f, pA=pA, evA=evA,
                         selA=selA, rmsA=rmsA, pB=pB, eB=eB, rmsB=rmsB, vA2=pA2[0], nA=selA.sum()))

print(f"c = {C:.1f} +- {DC:.1f} m/s  (T = {T_C} +- {DT_C} degC)")
print(f"{'v visee':>7s} {'d':>5s} | {'v (d fixe)':>13s} {'rms':>4s} {'n':>4s} | {'v (d+0,9)':>9s} | {'v (d libre)':>14s} {'d (libre)':>12s}")
for p in sorted(passages, key=lambda p: (p["d"], p["vn"])):
    vB, dB = p["pB"][0] * 3.6, p["pB"][1]
    print(f"{p['vn']:7d} {p['d']:5.2f} | {p['pA'][0]*3.6:6.1f} +- {p['evA']*3.6:4.1f} {p['rmsA']:4.2f} {p['nA']:4d} | "
          f"{p['vA2']*3.6:9.1f} | {vB:6.1f} +- {p['eB'][0]*3.6:5.1f} {dB:5.2f} +- {p['eB'][1]:5.2f}")

# ---------------------------------------------------------------- figure 1: spectrogrammes
def spectrogramme(ax, p, demi=3.5):
    x = p["x"].astype(float)
    fs, tp = p["fs"], p["pA"][2]
    n, hop = 8192, 512
    i0 = max(0, int((tp - demi) * fs) - n // 2)
    i1 = min(len(x), int((tp + demi) * fs) + n // 2)
    seg = x[i0:i1] - x[i0:i1].mean()
    deb = np.arange(0, len(seg) - n, hop)
    tr = np.lib.stride_tricks.sliding_window_view(seg, n)[deb] * np.hanning(n)
    S = 20 * np.log10(np.abs(np.fft.rfft(tr, 4 * n, axis=1)) + 1e-9)
    fr = np.fft.rfftfreq(4 * n, 1 / fs)
    m = (fr > 380) & (fr < 500)
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
cmap = plt.get_cmap("viridis")
vitesses = sorted({p["vn"] for p in passages})
coul = {v: cmap(i / (len(vitesses) - 1) * 0.9) for i, v in enumerate(vitesses)}
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
for ax, d in zip(axs, (3.53, 8.97)):
    for p in sorted((p for p in passages if abs(p["d"] - d) < 0.01), key=lambda p: p["vn"]):
        v, dd, tp = p["pA"]
        s = p["selA"]
        ax.plot(p["t"][s] - tp, p["f"][s], ".", ms=1.5, color=coul[p["vn"]], alpha=0.35)
        tt = np.linspace(p["t"][s].min(), p["t"][s].max(), 400)
        ax.plot(tt - tp, modele(tt, v, dd, tp), color=coul[p["vn"]], lw=1.8, label=f"{p['vn']} km/h")
    ax.axhline(F0, color="0.4", ls="--", lw=0.8)
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(415, 470)
    ax.xaxis.set_major_formatter(VIRG)
    ax.set_xlabel("Temps par rapport au passage (s)")
    ax.text(0.03, 0.04, f"$d$ = {d:.2f} m".replace(".", ","), transform=ax.transAxes)
    ax.legend(fontsize=9, loc="upper right", title="Vitesse visée", title_fontsize=9)
axs[0].set_ylabel("Fréquence reçue (Hz)")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "02_courbes_doppler.png"), dpi=200)
plt.close(fig)

# ---------------------------------------------------------------- figure 3: vitesse mesuree vs visee
fig, ax = plt.subplots(figsize=(5.6, 4.6))
for d, mk, c in ((3.53, "o", "C0"), (8.97, "s", "C1")):
    ps = [p for p in passages if abs(p["d"] - d) < 0.01]
    ax.errorbar([p["vn"] for p in ps], [p["pA"][0] * 3.6 for p in ps], yerr=[p["evA"] * 3.6 for p in ps],
                fmt=mk, color=c, capsize=3, ms=6, label=f"$d$ = {d:.2f} m".replace(".", ","))
vv = np.array([0, 75])
ax.fill_between(vv, np.maximum((vv - 4) / 1.1, 0), vv, color="0.85", lw=0, label="Plage du règlement n° 39")
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
