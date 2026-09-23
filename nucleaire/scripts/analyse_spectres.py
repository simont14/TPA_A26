# -*- coding: utf-8 -*-
"""
Calibration en energie (canal -> keV) d'un scintillateur a partir de spectres
gamma de sources connues, puis trace des spectres (coups/s vs energie).
"""
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPEC_DIR = os.path.join(BASE, "data")
OUT_DIR = os.environ.get("ANALYSE_OUT_DIR", os.path.join(BASE, "figures"))
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    "font.size": 11,
    "figure.figsize": (8, 5),
    "axes.grid": True,
    "grid.alpha": 0.3,
})


def fr(x, fmt="g"):
    """Formate un nombre avec une virgule decimale (convention francaise),
    pour rester coherent avec le texte du rapport LaTeX."""
    return format(x, fmt).replace(".", ",")


def fr_sci(x, sig=3):
    """Formate un nombre en notation scientifique a-la-LaTeX (a * 10^n),
    avec virgule decimale, pour les legendes matplotlib (mathtext)."""
    if x == 0:
        return "0"
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / 10**exp
    return f"{fr(mant, f'.{sig - 1}f')} \\times 10^{{{exp}}}"


def _fr_tick(x, pos=None):
    return fr(x)


def signed_term(x, fmt="0.2f"):
    """Terme '+ 1,23' ou '- 1,23' pour construire une somme sans double
    signe disgracieux (evite '+ -3,01')."""
    s = fr(abs(x), fmt)
    return f"+ {s}" if x >= 0 else f"- {s}"


def use_fr_ticks(ax):
    """Applique le formatage a virgule decimale sur les axes lineaires
    (les axes logarithmiques gardent leur formatage standard)."""
    if ax.get_xscale() == "linear":
        ax.xaxis.set_major_formatter(FuncFormatter(_fr_tick))
    if ax.get_yscale() == "linear":
        ax.yaxis.set_major_formatter(FuncFormatter(_fr_tick))


def parse_spe(path):
    with open(path, "r", errors="ignore") as f:
        lines = [l.rstrip("\n") for l in f]

    idx = {l.strip(): i for i, l in enumerate(lines)}

    # Live time / real time
    i_t = idx["$MEAS_TIM:"]
    live_time, real_time = (float(x) for x in lines[i_t + 1].split())

    # Data block
    i_d = idx["$DATA:"]
    first_ch, last_ch = (int(x) for x in lines[i_d + 1].split())
    n = last_ch - first_ch + 1
    counts = np.array([float(lines[i_d + 2 + k]) for k in range(n)])

    # Existing MCA calibration (only used as an initial guess to locate peaks)
    a0, a1 = 0.0, 1.0
    if "$MCA_CAL:" in idx:
        i_c = idx["$MCA_CAL:"]
        vals = lines[i_c + 2].split()
        a0, a1 = float(vals[0]), float(vals[1])

    return counts, live_time, real_time, (a0, a1)


def gauss_lin(x, amp, mu, sigma, m, b):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + m * x + b


def locate_peak(channels, counts, guess_channel, search_half=160, smooth=7,
                 min_prom_frac=0.03):
    """Coarse localisation: smooth the region around `guess_channel` and pick
    the most significant local maximum near the expected position."""
    lo = max(0, int(guess_channel - search_half))
    hi = min(len(counts), int(guess_channel + search_half))
    y = counts[lo:hi]
    if len(y) < 10 or y.max() <= 0:
        return None
    kernel = np.ones(smooth) / smooth
    ys = np.convolve(y, kernel, mode="same")
    peaks, props = find_peaks(ys, prominence=max(ys.max() * min_prom_frac, 1))
    if len(peaks) == 0:
        return lo + int(np.argmax(ys))
    center = guess_channel - lo
    max_prom = props["prominences"].max()
    candidates = [p for p, pr in zip(peaks, props["prominences"])
                  if pr >= 0.15 * max_prom]
    best = min(candidates, key=lambda p: abs(p - center))
    return lo + int(best)


def fit_peak(channels, counts, guess_channel, half_win=35):
    mu_coarse = locate_peak(channels, counts, guess_channel)
    if mu_coarse is None:
        return None

    lo = max(0, int(mu_coarse - half_win))
    hi = min(len(counts), int(mu_coarse + half_win))
    x = channels[lo:hi]
    y = counts[lo:hi]
    if len(x) < 8 or y.max() <= 0:
        return None

    amp0 = max(y.max() - np.median(y), 1)
    b0 = np.median(y)
    sigma0 = half_win / 5.0

    bounds_lo = [0, mu_coarse - half_win * 0.5, 1.0, -np.inf, 0]
    bounds_hi = [np.inf, mu_coarse + half_win * 0.5, half_win, np.inf, np.inf]

    try:
        popt, pcov = curve_fit(
            gauss_lin, x, y,
            p0=[amp0, mu_coarse, sigma0, 0.0, b0],
            bounds=(bounds_lo, bounds_hi),
            maxfev=40000,
        )
    except RuntimeError:
        return None

    amp, mu, sigma, m, b = popt
    if amp <= 0 or sigma <= 0:
        return None
    mu_err = np.sqrt(pcov[1, 1]) if np.all(np.isfinite(pcov)) else np.nan
    fwhm = 2.3548 * sigma
    return {"mu": mu, "mu_err": mu_err, "sigma": sigma, "fwhm": fwhm,
            "amp": amp, "x": x, "y": y, "fit": gauss_lin(x, *popt)}


SOURCES = {
    "BA_133_1": {"label": "Ba-133", "lines": [30.85, 81.0, 302.85, 356.02], "xlim": (0, 450)},
    "CO_157_1": {"label": "Co-57", "lines": [122.06], "xlim": (0, 150)},
    "CO_60_1":  {"label": "Co-60", "lines": [1173.2, 1332.5], "xlim": (0, 1500)},
    "CS_137_1": {"label": "Cs-137", "lines": [32.0, 661.7], "xlim": (0, 800)},
    "NA_22_1":  {"label": "Na-22", "lines": [511.0, 1274.5], "xlim": (0, 1500)},
}
BACKGROUND_FILE = "BG"

# ---------------------------------------------------------------------------
# 1) Charger tous les spectres
# ---------------------------------------------------------------------------
data = {}
for key in list(SOURCES) + [BACKGROUND_FILE]:
    path = os.path.join(SPEC_DIR, key + ".Spe")
    counts, live, real, cal0 = parse_spe(path)
    channels = np.arange(len(counts))
    data[key] = {"counts": counts, "channels": channels,
                 "live": live, "real": real, "cal0": cal0}
    print(f"{key}: {len(counts)} canaux, temps vif={live:.0f}s, temps reel={real:.0f}s")

# ---------------------------------------------------------------------------
# 2) Identifier les photopics et ajuster une gaussienne pour chaque source
# ---------------------------------------------------------------------------
cal_points = []  # (channel, energy, channel_err, source_label)
peak_results = {}  # key -> list of dict(energy, fit result)

for key, info in SOURCES.items():
    d = data[key]
    a0, a1 = d["cal0"]
    peaks_here = []
    for E in info["lines"]:
        guess_ch = (E - a0) / a1
        half_win = int(np.clip(0.035 * guess_ch + 25, 25, 90))
        res = fit_peak(d["channels"], d["counts"], guess_ch, half_win=half_win)
        if res is None:
            print(f"  !! pic non trouve pour {info['label']} @ {E} keV (guess canal {guess_ch:.0f})")
            continue
        cal_points.append((res["mu"], E, res["mu_err"], info["label"]))
        peaks_here.append({"energy": E, **res})
        print(f"  {info['label']:8s} {E:8.2f} keV -> canal {res['mu']:8.2f} +/- {res['mu_err']:.2f} "
              f"(FWHM={res['fwhm']:.1f} canaux)")
    peak_results[key] = peaks_here

# ---------------------------------------------------------------------------
# 3) Regression E(canal): lineaire ET quadratique
#
# Les incertitudes statistiques des centroides (ajustement gaussien) sont
# beaucoup plus petites que la dispersion reelle des points autour de la
# droite (non-linearites, asymetrie des pics, etc.). Une regression ponderee
# par ces erreurs statistiques accorderait donc un poids demesure aux pics
# les plus etroits: on utilise une regression non ponderee (moindres carres
# ordinaires), plus robuste ici.
#
# On compare un modele lineaire et un modele quadratique: aux basses energies
# (points ~30 keV), la reponse d'un scintillateur/PM s'ecarte souvent d'une
# droite parfaite. Le modele quadratique n'est retenu que s'il ameliore
# nettement le RMSE.
# ---------------------------------------------------------------------------
chans = np.array([p[0] for p in cal_points])
energies = np.array([p[1] for p in cal_points])
chan_errs = np.array([p[2] if np.isfinite(p[2]) and p[2] > 0 else 1.0 for p in cal_points])
labels = [p[3] for p in cal_points]

n = len(chans)


def fit_stats(coef, x, y):
    pred = np.polyval(coef, x)
    resid = y - pred
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1 - ss_res / ss_tot
    rmse = np.sqrt(ss_res / len(x))
    return pred, resid, r2, rmse


lin_coef, lin_cov = np.polyfit(chans, energies, 1, cov=True)
lin_pred, lin_resid, lin_r2, lin_rmse = fit_stats(lin_coef, chans, energies)

quad_coef, quad_cov = np.polyfit(chans, energies, 2, cov=True)
quad_pred, quad_resid, quad_r2, quad_rmse = fit_stats(quad_coef, chans, energies)

print("\nCalibration lineaire:    E = {:.5e}*canal + {:.3f}  (R^2={:.6f}, RMSE={:.2f} keV)".format(
    lin_coef[0], lin_coef[1], lin_r2, lin_rmse))
print("Calibration quadratique: E = {:.5e}*canal^2 + {:.5e}*canal + {:.3f}  (R^2={:.6f}, RMSE={:.2f} keV)".format(
    quad_coef[0], quad_coef[1], quad_coef[2], quad_r2, quad_rmse))

USE_QUAD = quad_rmse < 0.7 * lin_rmse
if USE_QUAD:
    model_coef, pred, resid, r2, rmse = quad_coef, quad_pred, quad_resid, quad_r2, quad_rmse
    model_label = (f"$E = {fr_sci(quad_coef[0])} \\cdot \\mathrm{{canal}}^2 "
                   f"{signed_term(quad_coef[1], '.4f')} \\cdot \\mathrm{{canal}} "
                   f"{signed_term(quad_coef[2], '.2f')}$\n"
                   f"$R^2$ = {fr(r2, '.5f')}, RMSE = {fr(rmse, '.1f')} keV")
    print(f"\n-> Modele retenu: QUADRATIQUE (amelioration nette du RMSE, non-linearite detectee)")
else:
    model_coef, pred, resid, r2, rmse = lin_coef, lin_pred, lin_resid, lin_r2, lin_rmse
    model_label = (f"$E = {fr(lin_coef[0], '.4f')} \\cdot \\mathrm{{canal}} "
                   f"{signed_term(lin_coef[1], '.2f')}$\n"
                   f"$R^2$ = {fr(r2, '.5f')}, RMSE = {fr(rmse, '.1f')} keV")
    print(f"\n-> Modele retenu: LINEAIRE (le quadratique n'ameliore pas assez le fit)")


def chan_to_E(ch):
    return np.polyval(model_coef, ch)


# ---------------------------------------------------------------------------
# 4) Graphique de calibration
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(6.6, 5.4), sharex=True,
    gridspec_kw={"height_ratios": [3, 1]}
)

unique_labels = sorted(set(labels))
cmap = plt.get_cmap("tab10")
color_of = {lab: cmap(i) for i, lab in enumerate(unique_labels)}

for lab in unique_labels:
    idxs = [i for i, l in enumerate(labels) if l == lab]
    ax1.errorbar(chans[idxs], energies[idxs], xerr=chan_errs[idxs],
                 fmt="o", ms=7, capsize=3, color=color_of[lab], label=lab)

xfit = np.linspace(0, max(chans) * 1.08, 200)
if USE_QUAD:
    ax1.plot(xfit, np.polyval(lin_coef, xfit), color="0.6", ls=":", lw=1.2,
              label=f"Ajust. linéaire (RMSE={fr(lin_rmse, '.1f')} keV)")
ax1.plot(xfit, chan_to_E(xfit), "k--", lw=1.3, label=model_label)
ax1.set_ylabel("Énergie (keV)")
ax1.set_title("Étalonnage en énergie du scintillateur")
ax1.legend(loc="upper left", fontsize=9)
use_fr_ticks(ax1)

local_slope = np.polyval(np.polyder(model_coef), chans)
ax2.axhline(0, color="k", lw=1)
ax2.axhspan(-rmse, rmse, color="gray", alpha=0.15, label=f"±RMSE ({fr(rmse, '.1f')} keV)")
for lab in unique_labels:
    idxs = [i for i, l in enumerate(labels) if l == lab]
    ax2.errorbar(chans[idxs], resid[idxs], yerr=np.abs(chan_errs[idxs] * local_slope[idxs]),
                 fmt="o", ms=6, capsize=3, color=color_of[lab])
ax2.set_ylabel("Résidu (keV)")
ax2.set_xlabel("Numéro de canal")
ax2.legend(loc="upper right", fontsize=8)
use_fr_ticks(ax2)

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "01_calibration.png"), dpi=200)
plt.close(fig)

# ---------------------------------------------------------------------------
# 5) Spectres individuels (coups/s vs energie) avec pics identifies
# ---------------------------------------------------------------------------
def plot_spectrum(key, title, filename, xlim, mark_peaks=None, bg_key=None, yscale="linear"):
    d = data[key]
    E = chan_to_E(d["channels"])
    cps = d["counts"] / d["live"]

    if bg_key is not None:
        # Meme MCA, meme gain, meme decoupage en canaux pour toutes les
        # acquisitions (le gain n'a pas ete retouche pendant la seance):
        # la soustraction se fait donc canal a canal, sans interpolation,
        # directement sur les taux de comptage (coups/s - coups/s = taux
        # de comptage net, valide car ce sont deux processus de Poisson
        # independants).
        cps_bg = data[bg_key]["counts"] / data[bg_key]["live"]
        cps = cps - cps_bg

    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.step(E, cps, where="mid", color="tab:blue", lw=1.0, label=title)

    mask = (E >= xlim[0]) & (E <= xlim[1])
    ymax = cps[mask].max() if np.any(mask) else cps.max()

    ax.set_xlim(*xlim)
    ax.set_yscale(yscale)
    if yscale == "log":
        floor = cps[(cps > 0) & mask].min() if np.any((cps > 0) & mask) else 1e-3
        ax.set_ylim(floor * 0.5, ymax * 5)
    else:
        ax.set_ylim(0, ymax * 1.45)

    if mark_peaks:
        counts = d["counts"]
        smooth_kernel = np.ones(5) / 5
        counts_smooth = np.convolve(counts, smooth_kernel, mode="same")
        for pk in mark_peaks:
            Epk = pk["energy"]
            mu = pk.get("mu")
            sigma = pk.get("sigma")
            if mu is not None and sigma is not None:
                # Snap the marker to the apex of a lightly smoothed version of
                # the histogram near the fitted centroid, so the dotted line
                # lands visually centered on the peak tip rather than on a
                # single noisy bin or on the literature energy (which can
                # differ slightly from the fit due to calibration residuals).
                lo_c = max(0, int(round(mu - 1.5 * sigma)))
                hi_c = min(len(counts), int(round(mu + 1.5 * sigma)) + 1)
                idx = lo_c + int(np.argmax(counts_smooth[lo_c:hi_c]))
                E_line, local_max = E[idx], cps[idx]
            else:
                window = (E > Epk - 20) & (E < Epk + 20)
                if np.any(window):
                    idx = np.flatnonzero(window)[np.argmax(cps[window])]
                    E_line, local_max = E[idx], cps[idx]
                else:
                    E_line, local_max = Epk, ymax

            label_y = local_max * 2.2 if yscale == "log" else local_max + ymax * 0.08
            ax.axvline(E_line, color="red", ls=":", lw=1, alpha=0.7)
            ax.annotate(f"{fr(Epk, '.2f')} keV",
                        xy=(E_line, local_max), xytext=(E_line, label_y),
                        ha="center", va="bottom",
                        fontsize=8, color="red",
                        bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1))

    ax.set_xlabel("Énergie (keV)")
    ax.set_ylabel("Taux de comptage net (coups/s)" if bg_key is not None
                  else "Taux de comptage (coups/s)")
    ax.set_title(title, pad=14)
    use_fr_ticks(ax)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, filename), dpi=200)
    plt.close(fig)


plot_num = 2
for key, info in SOURCES.items():
    fname = f"{plot_num:02d}_spectre_{info['label'].replace('-', '')}.png"
    plot_spectrum(key, f"Spectre {info['label']}", fname, info["xlim"],
                  mark_peaks=peak_results[key], bg_key=BACKGROUND_FILE)
    plot_num += 1

bg_xlim = (0, chan_to_E(data[BACKGROUND_FILE]["channels"].max()))
plot_spectrum(BACKGROUND_FILE, "Bruit de fond (background)",
              f"{plot_num:02d}_spectre_background.png", bg_xlim, yscale="log")
plot_num += 1

# ---------------------------------------------------------------------------
# 6) Identification de la raie dominante du bruit de fond (K-40, 1460.8 keV)
# ---------------------------------------------------------------------------
K40_ENERGY = 1460.8  # keV, litterature (IAEA)

d_bg = data[BACKGROUND_FILE]
E_bg_grid = chan_to_E(d_bg["channels"])
guess_ch_k40 = float(np.interp(K40_ENERGY, E_bg_grid, d_bg["channels"]))
k40_fit = fit_peak(d_bg["channels"], d_bg["counts"], guess_ch_k40, half_win=220)

if k40_fit is not None:
    slope_k40 = np.polyval(np.polyder(model_coef), k40_fit["mu"])
    k40_E_mesuree = chan_to_E(k40_fit["mu"])
    k40_E_err = k40_fit["mu_err"] * slope_k40
    print(f"\nBruit de fond: pic dominant a {k40_E_mesuree:.1f} +/- {k40_E_err:.1f} keV "
          f"(K-40 attendu: {K40_ENERGY:.1f} keV)")
else:
    print("\nBruit de fond: pic du K-40 non trouve automatiquement.")

fig, ax = plt.subplots(figsize=(8.5, 5))
E_bg = E_bg_grid
cps_bg = d_bg["counts"] / d_bg["live"]
ax.step(E_bg, cps_bg, where="mid", color="0.35", lw=0.8)
ax.set_yscale("log")
mask_bg = (E_bg >= bg_xlim[0]) & (E_bg <= bg_xlim[1])
floor_bg = cps_bg[(cps_bg > 0) & mask_bg].min()
ax.set_ylim(floor_bg * 0.5, cps_bg[mask_bg].max() * 5)
ax.set_xlim(*bg_xlim)
if k40_fit is not None:
    idx_k40 = np.argmin(np.abs(E_bg - k40_E_mesuree))
    y_k40 = cps_bg[max(0, idx_k40 - 3):idx_k40 + 4].max()
    ax.axvline(k40_E_mesuree, color="red", ls=":", lw=1.2)
    ax.annotate(
        f"$^{{40}}$K, {fr(k40_E_mesuree, '.0f')} keV",
        xy=(k40_E_mesuree, y_k40),
        xytext=(k40_E_mesuree - 430, y_k40 * 12),
        fontsize=9, color="red", ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color="red", lw=1),
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5),
    )
ax.set_xlabel("Énergie (keV)")
ax.set_ylabel("Taux de comptage (coups/s)")
ax.set_title("Bruit de fond ambiant: identification du $^{40}$K", pad=14)
use_fr_ticks(ax)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, f"{plot_num:02d}_background_K40.png"), dpi=200)
plt.close(fig)
plot_num += 1

# ---------------------------------------------------------------------------
# 7) Anatomie d'un spectre gamma: front Compton et pic de retrodiffusion
#    (illustre sur le spectre du Cs-137, qui n'a qu'une seule raie a 661.7 keV)
# ---------------------------------------------------------------------------
MC2 = 511.0  # keV, energie de masse au repos de l'electron
E_CS137 = 661.7
E_compton_edge = 2 * E_CS137**2 / (MC2 + 2 * E_CS137)
E_backscatter = E_CS137 - E_compton_edge

print(f"\nCs-137: front Compton attendu a {E_compton_edge:.1f} keV, "
      f"pic de retrodiffusion attendu a {E_backscatter:.1f} keV")

# Positions mesurees sur le spectre (point mi-hauteur du front, centroide de
# la bosse de retrodiffusion), telles que rapportees dans le texte du
# rapport: valeurs fixes, non re-derivees ici, pour eviter toute divergence
# avec la discussion.
E_compton_edge_mesure = 462.0
E_backscatter_mesure = 195.0

d_cs = data["CS_137_1"]
E_cs = chan_to_E(d_cs["channels"])
cps_cs = d_cs["counts"] / d_cs["live"] - data[BACKGROUND_FILE]["counts"] / data[BACKGROUND_FILE]["live"]

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.step(E_cs, cps_cs, where="mid", color="tab:blue", lw=1.1)
ax.set_xlim(0, 750)
mask_cs = (E_cs >= 0) & (E_cs <= 750)
ymax_cs = cps_cs[mask_cs].max()
ax.set_ylim(0, ymax_cs * 1.5)

photopeak_cs = [p for p in peak_results["CS_137_1"] if abs(p["energy"] - E_CS137) < 1][0]


def _mark(E_line, label, color, xytext, ls=":"):
    idx = np.argmin(np.abs(E_cs - E_line))
    y_local = cps_cs[max(0, idx - 3):idx + 4].max()
    ax.axvline(E_line, color=color, ls=ls, lw=1.3, alpha=0.8)
    ax.annotate(
        label, xy=(E_line, y_local), xytext=xytext,
        fontsize=8, color=color, ha="center", va="bottom",
        arrowprops=dict(arrowstyle="->", color=color, lw=1),
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1),
    )


# Photopic: theorie et mesure coincident (c'est la calibration elle-meme).
_mark(E_CS137, "Photopic (effet\nphotoélectrique)\n661,7 keV", "tab:red",
      (E_CS137, cps_cs[np.argmin(np.abs(E_cs - E_CS137))] + 0.16 * ymax_cs))

# Front Compton et retrodiffusion: ligne pointillee = position theorique
# (eq. 1/2), fleche separee = position mesuree sur le spectre.
_mark(E_compton_edge, f"Front Compton\nthéorique (éq. 1)\n{fr(E_compton_edge, '.0f')} keV",
      "tab:purple", (E_compton_edge + 90, ymax_cs * 0.62), ls=":")
_mark(E_compton_edge_mesure, f"mesuré\n$\\approx${fr(E_compton_edge_mesure, '.0f')} keV",
      "tab:purple", (E_compton_edge_mesure + 90, ymax_cs * 0.30), ls="-.")

_mark(E_backscatter, f"Pic de rétrodiffusion\nthéorique (éq. 2)\n{fr(E_backscatter, '.0f')} keV",
      "tab:green", (E_backscatter - 60, ymax_cs * 0.95), ls=":")
_mark(E_backscatter_mesure, f"mesuré\n$\\approx${fr(E_backscatter_mesure, '.0f')} keV",
      "tab:green", (E_backscatter_mesure - 60, ymax_cs * 0.65), ls="-.")

ax.axvspan(0, E_compton_edge, color="tab:orange", alpha=0.08)
ax.text(320, ymax_cs * 1.42, "continuum Compton",
        fontsize=9, color="tab:orange", ha="center", style="italic")

ax.set_xlabel("Énergie (keV)")
ax.set_ylabel("Taux de comptage net (coups/s)")
ax.set_title("Anatomie du spectre du Cs-137: diffusion Compton et rétrodiffusion", pad=14)
use_fr_ticks(ax)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, f"{plot_num:02d}_compton_cs137.png"), dpi=200)
plt.close(fig)
plot_num += 1

# ---------------------------------------------------------------------------
# 8) Resolution en energie (FWHM/E) du scintillateur en fonction de l'energie
# ---------------------------------------------------------------------------
res_E, res_pct = [], []
for key, peaks in peak_results.items():
    for p in peaks:
        slope_here = np.polyval(np.polyder(model_coef), p["mu"])
        fwhm_keV = p["fwhm"] * slope_here
        res_E.append(p["energy"])
        res_pct.append(100 * fwhm_keV / p["energy"])
res_E = np.array(res_E)
res_pct = np.array(res_pct)

# Regression en loi de puissance R(%) = A * E^b via les logarithmes.
b, log_a = np.polyfit(np.log(res_E), np.log(res_pct), 1)
A = np.exp(log_a)

fig, ax = plt.subplots(figsize=(7.5, 5.2))
ax.scatter(res_E, res_pct, s=45, color="tab:blue", zorder=3, label="Pics ajustés")
E_fit_grid = np.linspace(res_E.min() * 0.8, res_E.max() * 1.1, 200)
ax.plot(E_fit_grid, A * E_fit_grid**b, "k--", lw=1.3,
        label=f"$R = {fr(A, '.0f')} \\cdot E^{{{fr(b, '.2f')}}}$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Énergie (keV)")
ax.set_ylabel("Résolution $R = \\mathrm{FWHM}/E$ (%)")
ax.set_title("Résolution en énergie du scintillateur NaI(Tl)", pad=14)
ax.legend(fontsize=9)
use_fr_ticks(ax)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, f"{plot_num:02d}_resolution.png"), dpi=200)
plt.close(fig)
plot_num += 1

print(f"\nResolution: R(%) = {A:.0f} * E^{b:.3f}  (E en keV)")

print("\nGraphiques sauvegardes dans:", OUT_DIR)
