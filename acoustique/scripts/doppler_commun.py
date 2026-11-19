# -*- coding: utf-8 -*-
"""
Code commun a analyse_doppler.py et fit_manuel.py: lecture des passages,
suivi de la frequence, modele Doppler et mesure de v et d par regions.

Methode (une region choisie a la main pour chaque information):
  - vitesse: les deux queues de la courbe (approche et eloignement), ou la
    frequence tend vers f0/(1 -+ v/c). Le rapport des deux niveaux donne v,
    et leur moyenne harmonique donne la frequence de la source f0.
  - distance: la region centrale, ou la pente de la chute (~ v^2/d) donne d,
    avec v et f0 fixes.
Les queues n'atteignent jamais tout a fait leur plateau (cos(theta) < 1), et
cet ecart depend de d. Les deux etapes sont donc alternees jusqu'a ce que v et
d ne changent plus: queues -> (v, f0) avec d fixe, centre -> (d, tp) avec v fixe.
"""
import glob
import json
import os
import re

import numpy as np
from scipy.io import wavfile
from scipy.ndimage import median_filter
from scipy.optimize import least_squares

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
REGIONS_JSON = os.path.join(BASE, "scripts", "regions_manuelles.json")

# ---------------------------------------------------------------- constantes
T_C, DT_C = 1.0, 1.0                     # temperature de l'air [degC]
C = 331.3 * np.sqrt(1 + T_C / 273.15)    # vitesse du son [m/s]
DC = 0.606 * DT_C                        # incertitude sur c [m/s]
F0 = 440.0                               # frequence jouee [Hz]
L_ZONE = 10.0                            # zone de depart du fit automatique: |x| <= 10 m
NFFT, HOP, PAD = 4096, 256, 4            # trame de 93 ms, pas de 5,8 ms
BANDE = (400.0, 485.0)                   # recherche du pic [Hz] (evite la raie fixe a 385 Hz)
SNR_MIN = 15.0                           # pic / mediane de la bande [dB]
CORR = np.sqrt(NFFT / HOP)               # 16 trames successives partagent les memes donnees
EXCLUS = {"440Hz_55kmh_353cm.wav"}       # premier passage, test seulement (pas fait aux deux distances)
D_MAX = 30.0                             # limite permise pour d [m]
CLES = ("a0", "a1", "c0", "c1", "e0", "e1")   # bornes des regions approche, centre, eloignement


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


def position(t, v, d, tp):
    """Position x de la source [m] au moment ou elle emet le son recu au temps t."""
    return v * u_emission(t, v, d, tp)


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


def ajuste_auto(t, f, mono, p, d_fixe):
    """Curve fit automatique (d fixee) sur la zone |x| <= L_ZONE, recalculee avec les
    parametres. Sert seulement a placer les regions de depart."""
    for _ in range(40):
        sel = mono & (np.abs(position(t, *p)) <= L_ZONE)
        res = lambda q: modele(t[sel], q[0], d_fixe, q[1]) - f[sel]
        r = least_squares(res, p[[0, 2]], loss="soft_l1", f_scale=1.0,
                          bounds=([0.3, t.min()], [40, t.max()]))
        pn = np.array([r.x[0], d_fixe, r.x[1]])
        fini = np.allclose(pn, p, rtol=1e-5, atol=1e-5)
        p = pn
        if fini:
            break
    return p


# ---------------------------------------------------------------- mesure par regions
def _ls(res, q0, lb, ub):
    q0 = np.clip(q0, np.array(lb) + 1e-9, np.array(ub) - 1e-9)
    return least_squares(res, q0, loss="soft_l1", f_scale=1.0, bounds=(lb, ub)).x


def _fit_centre(tc, fc, v, f0, d, tp):
    """(d, tp) sur la region centrale, v et f0 fixes."""
    return _ls(lambda q: modele(tc, v, q[0], q[1], f0) - fc, [d, tp], [0.3, tp - 2], [D_MAX, tp + 2])


def _fit_queues(tq, fq, d, tp, v, f0):
    """(v, f0) sur les deux queues, d et tp fixes."""
    return _ls(lambda q: modele(tq, q[0], d, tp, q[1]) - fq, [v, f0], [0.3, 400.0], [60.0, 480.0])


def _incertitudes(res, q):
    """Incertitudes des parametres (matrice de covariance, jacobien numerique),
    multipliees par CORR, et ecart quadratique moyen des residus."""
    q = np.asarray(q, float)
    r = res(q)
    n, k = len(r), len(q)
    J = np.empty((n, k))
    for j in range(k):
        h = 1e-6 * max(abs(q[j]), 1.0)
        dq = np.zeros(k)
        dq[j] = h
        J[:, j] = (res(q + dq) - res(q - dq)) / (2 * h)
    e = np.full(k, np.nan)
    if n > k:
        try:
            cov = np.linalg.inv(J.T @ J) * np.sum(r**2) / (n - k)
            e = np.sqrt(np.abs(np.diag(cov))) * CORR
        except np.linalg.LinAlgError:
            pass
    return e, float(np.sqrt(np.mean(r**2)))


def selections(p, R):
    t = p["t"]
    return ((t >= R["a0"]) & (t <= R["a1"]), (t >= R["c0"]) & (t <= R["c1"]),
            (t >= R["e0"]) & (t <= R["e1"]))


def analyse_regions(p, R, tp0=None):
    """Mesure v et d d'un passage a partir des trois regions R (temps absolus).
    Retourne None s'il manque des points dans une region."""
    t, f = p["t"], p["f"]
    sa, sc, se = selections(p, R)
    if min(sa.sum(), sc.sum(), se.sum()) < 5:
        return None
    tc, fc = t[sc], f[sc]
    tq, fq = np.r_[t[sa], t[se]], np.r_[f[sa], f[se]]
    fa, fe = float(np.median(f[sa])), float(np.median(f[se]))
    # depart: plateaux pris tels quels
    v = max(C * (fa - fe) / (fa + fe), 0.5)
    f0 = 2 * fa * fe / (fa + fe)
    d = p["d"]
    tp = 0.5 * (R["c0"] + R["c1"]) if tp0 is None else tp0
    conv = False
    for _ in range(200):
        d_n, tp = _fit_centre(tc, fc, v, f0, d, tp)
        v_n, f0 = _fit_queues(tq, fq, d_n, tp, v, f0)
        conv = abs(v_n - v) < 1e-5 and abs(d_n - d) < 1e-5
        v, d = v_n, d_n
        if conv:
            break
    # incertitudes: fit de chaque region, plus le couplage entre v et d
    (ev_fit, _), rms_q = _incertitudes(lambda q: modele(tq, q[0], d, tp, q[1]) - fq, [v, f0])
    (ed_fit, _), rms_c = _incertitudes(lambda q: modele(tc, v, q[0], q[1], f0) - fc, [d, tp])
    hv, hd = 0.01 * v, 0.01 * d
    dd_dv = (_fit_centre(tc, fc, v + hv, f0, d, tp)[0] - _fit_centre(tc, fc, v - hv, f0, d, tp)[0]) / (2 * hv)
    dv_dd = (_fit_queues(tq, fq, d + hd, tp, v, f0)[0] - _fit_queues(tq, fq, d - hd, tp, v, f0)[0]) / (2 * hd)
    ev = float(np.sqrt(ev_fit**2 + (dv_dd * ed_fit)**2 + (v * DC / C)**2))
    ed = float(np.sqrt(ed_fit**2 + (dd_dv * ev_fit)**2))
    return dict(v=float(v), d=float(d), tp=float(tp), f0=float(f0), ev=ev, ed=ed,
                fa=fa, fe=fe, rms_q=rms_q, rms_c=rms_c, na=int(sa.sum()), nc=int(sc.sum()),
                ne=int(se.sum()), conv=conv, sel=sa | sc | se)


def vitesses_queues(p, R, r):
    """Vitesse donnee par chaque queue separement, avec f0 = 440 Hz (sinus numerique)
    et d, tp mesures. Un ecart entre les deux montre que la vitesse a change."""
    sa, _, se = selections(p, R)
    t, f = p["t"], p["f"]
    out = []
    for s in (sa, se):
        q = _ls(lambda q: modele(t[s], q[0], r["d"], r["tp"]) - f[s], [r["v"]], [0.1], [60.0])
        out.append(float(q[0]))
    return out


# ---------------------------------------------------------------- passages
def charger_passages():
    """Lit les .wav, suit la frequence et estime le moment du passage.
    Les noms de fichiers donnent la vitesse visee (km/h) et la distance (cm)."""
    passages = []
    for chemin in sorted(glob.glob(os.path.join(DATA, "*.wav"))):
        nom = os.path.basename(chemin)
        if nom in EXCLUS:
            continue
        vn, dcm = (int(s) for s in re.findall(r"_(\d+)(?:kmh|cm)", nom))
        fs, x = wavfile.read(chemin)
        t, f = suivi(x, fs)
        # passage approximatif: pente la plus negative de f(t) lissee sur 0,5 s, loin des bords
        lis = median_filter(f, size=int(0.5 / np.median(np.diff(t))) | 1, mode="nearest")
        pente = np.gradient(lis, t)
        pente[(t < t.min() + 1.5) | (t > t.max() - 1.5)] = 0
        tp0 = t[np.argmin(pente)]
        passages.append(dict(nom=nom, vn=vn, d=dcm / 100, fs=fs, x=x, t=t, f=f, tp0=tp0,
                             mono=zone_monotone(t, f, tp0)))
    passages.sort(key=lambda p: (p["d"], p["vn"]))
    return passages


K_CENTRE = 0.77    # centre: |t - tp| < 0,77 tau (pente > 50 % de son maximum)
K_QUEUE = (2.0, 5.0)   # queues: 2 tau < |t - tp| < 5 tau (pente < 15 %, cos(theta) de 0,89 a 0,98)
T_MIN = 0.3        # duree minimale d'une queue [s]
N_MIN = 10         # nombre minimal de points suivis dans une queue


def regions_depart(p):
    """Regions placees avec des criteres objectifs, sans utiliser la vitesse visee
    ni la distance au galon (temps absolus), et moment du passage estime.

    1. Partie ou la frequence diminue (zone_monotone): avant, la voiture
       accelerait encore, apres, elle freinait deja.
    2. Une marche douce f = m - A tanh((t - tp)/w) y est ajustee. Elle ne sert
       qu'a trouver le moment du passage tp et la duree de la chute. Pour une
       courbe Doppler, tau = d/v vaut environ 0,86 w.
    3. Centre: |t - tp| < 0,77 tau, la ou la pente depasse la moitie de son
       maximum (|x| < 0,77 d). Queues: de 2 tau a 5 tau de chaque cote
       (|x| de 2 d a 5 d), sans sortir de la partie ou la frequence diminue."""
    t, f, mono = p["t"], p["f"], p["mono"]
    tm, fm = t[mono], f[mono]
    q0 = [np.median(fm), 0.5 * (fm.max() - fm.min()), p["tp0"], 0.3]
    lb, ub = [380.0, 0.0, tm.min(), 0.02], [500.0, 60.0, tm.max(), 20.0]
    q = least_squares(lambda q: q[0] - q[1] * np.tanh((tm - q[2]) / q[3]) - fm,
                      np.clip(q0, np.array(lb) + 1e-9, np.array(ub) - 1e-9),
                      loss="soft_l1", f_scale=1.0, bounds=(lb, ub)).x
    tp, tau = q[2], 0.86 * q[3]
    m0, m1 = tm.min(), tm.max()
    a1, e0 = tp - K_QUEUE[0] * tau, tp + K_QUEUE[0] * tau
    a0, e1 = max(tp - K_QUEUE[1] * tau, m0), min(tp + K_QUEUE[1] * tau, m1)
    R = dict(a0=min(a0, a1 - T_MIN), a1=a1, c0=tp - K_CENTRE * tau, c1=tp + K_CENTRE * tau,
             e0=e0, e1=max(e1, e0 + T_MIN))
    # une queue avec trop peu de points suivis (signal faible) est allongee vers l'exterieur
    avant, apres = t[t <= R["a1"]], t[t >= R["e0"]]
    if np.sum(avant >= R["a0"]) < N_MIN and len(avant) >= N_MIN:
        R["a0"] = avant[-N_MIN]
    if np.sum(apres <= R["e1"]) < N_MIN and len(apres) >= N_MIN:
        R["e1"] = apres[N_MIN - 1]
    R = {k: float(np.clip(x, t.min(), t.max())) for k, x in R.items()}
    return R, float(tp)


def lire_regions():
    """Regions choisies a la main {nom: {a0, a1, c0, c1, e0, e1}} (temps absolus [s])."""
    if not os.path.exists(REGIONS_JSON):
        return {}
    with open(REGIONS_JSON, encoding="utf-8") as fh:
        return {nom: R for nom, R in json.load(fh).items() if all(k in R for k in CLES)}


def ecrire_regions(regions):
    regions = {nom: {k: round(float(R[k]), 4) for k in CLES} for nom, R in regions.items()}
    with open(REGIONS_JSON, "w", encoding="utf-8") as fh:
        json.dump(regions, fh, indent=2, ensure_ascii=False)
