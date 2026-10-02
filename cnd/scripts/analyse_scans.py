#!/usr/bin/env python3
"""
analyse_scans.py : comparaison de deux acquisitions TOPAZ du bloc LASERAX avec le CAD officiel.

  Scan_Bloc_Laser_3.UVData       balayage a la main, position donnee par l'horloge (20 Hz, 0,25 mm/position)
  Encodeur_Bloc_Laser_6.UVData   balayage a la main, position donnee par l'encodeur (0,25 mm/position)
  Laserax_CAD.step               piece usinee : plaque de 9,5 mm + lettres en relief de 5 mm

La sonde est posee sur la face plane. Sans lettre, l'echo de fond sort a l'epaisseur de la plaque (9,5 mm).
Sous une lettre, le son continue dans la lettre et l'echo sort au dessus de la lettre (14,5 mm).

Lancement :  python analyse_scans.py      (ecrit les figures dans ../figures et imprime les valeurs du rapport)
"""
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.interpolate import BSpline
from scipy.signal import correlate
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import least_squares

import uvdata_3d_gui as U

ICI = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ICI, '..', 'data', 'Data_Session_2')
FIG = os.path.join(ICI, '..', 'figures')

SCANS = [('Scan_Bloc_Laser_3.UVData', 'Horloge'), ('Encodeur_Bloc_Laser_6.UVData', 'Encodeur')]
LETTRES = ['L', 'A', 'S', 'E', 'R', 'A', 'X']

Y_FACE, Y_BASE, Y_HAUT = 9.5, 19.0, 24.0       # face sondee, dessus de la plaque, dessus des lettres (CAD, mm)
V_NOMINALE = 6320.0                             # vitesse longitudinale de l'aluminium dans UltraVision (m/s)
FENETRE_FOND = (7.5, 12.0)                      # fenetres de recherche des echos (profondeur, mm)
FENETRE_LETTRE = (12.5, 17.0)
SEUIL_SIGNAL = 10.0                             # % ecran : en dessous, A-scan sans couplage
DX_CAD = 0.05                                   # pas du masque CAD (mm)
# faisceaux retenus pour les ajustements : la sonde (38,4 mm) deborde du bloc (24 mm) et l'echo de fond
# baisse deja a environ 1 mm du bord du bloc
BANDE = lambda S: (S['z_cad'] > -10.6) & (S['z_cad'] < 10.9)

plt.rcParams.update({'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 9,
                     'savefig.dpi': 200, 'savefig.bbox': 'tight'})
C_HORLOGE, C_ENCODEUR = '#D9822B', '#2F5D8C'           # orange / bleu : paire lisible en daltonisme


# =====================================================================  CAD (STEP)
def lire_step(path):
    data = open(path, encoding='utf-8', errors='ignore').read().split('DATA;')[1]
    ents = {}
    for m in re.finditer(r'#(\d+)\s*=\s*([A-Z_0-9]*)\s*\((.*?)\);\s*(?=#|ENDSEC|\()', data, re.S):
        ents[int(m.group(1))] = (m.group(2), m.group(3).replace('\n', ''))
    return ents


def _refs(s):
    return [int(x) for x in re.findall(r'#(\d+)', s)]


def _pt(ents, i):
    return np.array([float(x) for x in re.search(r'\(([^()]*)\)\s*$', ents[i][1]).group(1).split(',')])


def aretes_step(ents):
    """Toutes les aretes du solide, echantillonnees en polylignes 3D."""
    vtx = lambda i: _pt(ents, _refs(ents[i][1])[0])
    segs = []
    for t, a in ents.values():
        if t != 'EDGE_CURVE':
            continue
        r = _refs(a)
        A, B = vtx(r[0]), vtx(r[1])
        kind, arg = ents[r[2]]
        if kind == 'CIRCLE':
            ax = _refs(ents[_refs(arg)[0]][1])
            o, nz, nx = _pt(ents, ax[0]), _pt(ents, ax[1]), _pt(ents, ax[2])
            R, ny = float(arg.split(',')[-1]), np.cross(nz, nx)
            ang = lambda p: np.arctan2((p - o) @ ny, (p - o) @ nx)
            meme_sens = a.strip().endswith('.T.')             # .F. : l'arete parcourt le cercle a l'envers
            t0, t1 = (ang(A), ang(B)) if meme_sens else (ang(B), ang(A))
            if t1 <= t0 + 1e-9:
                t1 += 2 * np.pi
            ts = np.linspace(t0, t1, 24)
            arc = o + R * (np.outer(np.cos(ts), nx) + np.outer(np.sin(ts), ny))
            segs.append(arc if meme_sens else arc[::-1])
        elif kind.startswith('B_SPLINE'):
            deg = int(arg.split(',')[1])
            cps = np.array([_pt(ents, k) for k in _refs(arg) if ents[k][0] == 'CARTESIAN_POINT'])
            nums = re.findall(r'\(([-0-9.,E+\s]*)\)', arg)
            T = np.repeat([float(x) for x in nums[-1].split(',')], [int(x) for x in nums[-2].split(',')])
            P = BSpline(T, cps, deg)(np.linspace(T[deg], T[-deg - 1], 400))
            ia, ib = np.argmin(np.linalg.norm(P - A, axis=1)), np.argmin(np.linalg.norm(P - B, axis=1))
            seg = P[ia:ib + 1] if ia <= ib else P[ib:ia + 1][::-1]
            segs.append(np.vstack([A, seg[::8], B]))
        else:
            segs.append(np.array([A, B]))
    return segs


def contours_lettres(segs, y=Y_BASE):
    """Contours fermes (x, z) des lettres a la base des lettres (y = 19 mm)."""
    pool = [s[:, [0, 2]] for s in segs
            if np.all(abs(s[:, 1] - y) < 1e-3) and np.all(abs(s[:, 2]) < 10.5) and np.all(abs(s[:, 0]) < 69)]
    proche = lambda a, b: np.hypot(*(a - b)) < 1e-3
    used, loops = [False] * len(pool), []
    for i in range(len(pool)):
        if used[i]:
            continue
        used[i], loop, change = True, list(pool[i]), True
        while change and not proche(loop[0], loop[-1]):
            change = False
            for j, s in enumerate(pool):
                if used[j]:
                    continue
                if proche(loop[-1], s[0]):
                    loop += list(s[1:])
                elif proche(loop[-1], s[-1]):
                    loop += list(s[::-1][1:])
                else:
                    continue
                used[j] = change = True
                break
        loops.append(np.array(loop))
    return loops


def masque(loops, xs, zs):
    X, Z = np.meshgrid(xs, zs, indexing='ij')
    pts = np.c_[X.ravel(), Z.ravel()]
    M = np.zeros(len(pts), bool)
    for l in loops:
        M ^= Path(l).contains_points(pts)          # pair-impair : les trous du A et du R restent vides
    return M.reshape(X.shape)


def cad():
    segs = aretes_step(lire_step(os.path.join(DATA, 'Laserax_CAD.step')))
    loops = contours_lettres(segs)
    xs, zs = np.arange(-72, 72, DX_CAD), np.arange(-13, 13, DX_CAD)
    return dict(segs=segs, loops=loops, xs=xs, zs=zs, M=masque(loops, xs, zs))


# =====================================================================  acquisitions
def _echo(vol, d, lo, hi):
    """Amplitude max et profondeur (centre de la zone > 50 % du max) de l'echo dans [lo, hi] mm."""
    k = np.nonzero((d >= lo) & (d <= hi))[0]
    w = vol[..., k]
    amp = w.max(-1)
    m = w >= 0.5 * amp[..., None]
    prof = (m * w * d[k]).sum(-1) / np.maximum((m * w).sum(-1), 1e-9)
    return amp, prof


def charger(nom):
    info, vols = U.read_acquisition(os.path.join(DATA, nom))
    vol, meta = vols[0]
    filled = meta['filled'].copy()
    if info['internal_clock']:
        keep = np.nonzero(filled)[0]
        v, f = vol[keep[0]:keep[-1] + 1], filled[keep[0]:keep[-1] + 1]
        trous = []
    else:
        pos = U.prepare_positions(vol, filled, interpolate=True)
        v, f = pos['v'], filled[pos['first']:pos['last'] + 1]
        trous = pos['gaps']
    dz = meta['sample_period_s'] * meta['velocity_m_s'] / 2 * 1000
    d = np.arange(v.shape[2]) * dz
    a_f, p_f = _echo(v, d, *FENETRE_FOND)
    a_l, p_l = _echo(v, d, *FENETRE_LETTRE)
    signal = (a_f + a_l) > SEUIL_SIGNAL
    r = np.where(signal, a_l / np.maximum(a_f + a_l, 1e-9), np.nan)     # 0 : plaque, 1 : lettre
    # chute de l'echo de fond : 0 sur la plaque, 1 sous une lettre. La reference est prise faisceau par
    # faisceau, parce que la calibration de sensibilite n'est pas appliquee aux donnees.
    ref = np.nanpercentile(np.where(signal, a_f, np.nan), 90, axis=0)
    q = np.where(signal, np.clip(1 - a_f / np.maximum(ref, 1e-9), 0, 1), np.nan)
    return dict(nom=nom, info=info, meta=meta, vol=v, filled=f, trous=trous, dz=dz, d=d,
                ds=meta['scan_res_mm'], di=meta['index_res_mm'], a_f=a_f, p_f=p_f, a_l=a_l, p_l=p_l, r=r, q=q,
                v=meta['velocity_m_s'], rate=info['rate_hz'])


# =====================================================================  recalage sur le CAD
def recaler(S, C):
    """Retournements et decalages entiers qui superposent le mieux la carte r au masque CAD
    echantillonne sur la meme grille (0,25 mm x 0,6 mm). Le pas de balayage reste le pas nominal."""
    ds, di = S['ds'], S['di']
    xs, zs = np.arange(-75, 75, ds), np.arange(-15, 15, di)
    ix = np.clip(np.round((xs - C['xs'][0]) / DX_CAD).astype(int), 0, len(C['xs']) - 1)
    iz = np.clip(np.round((zs - C['zs'][0]) / DX_CAD).astype(int), 0, len(C['zs']) - 1)
    ref = C['M'][np.ix_(ix, iz)].astype(float)
    ref[(xs < C['xs'][0]) | (xs > C['xs'][-1])] = 0
    ref -= ref.mean()
    r = np.nan_to_num(S['q'] - np.nanmean(S['q']))
    best = None
    for fx in (1, -1):
        for fz in (1, -1):
            img = r[::fx, ::fz]
            c = correlate(ref, img, mode='full', method='fft')
            k = np.unravel_index(np.argmax(c), c.shape)
            if best is None or c[k] > best[0]:
                best = (c[k], fx, fz, k[0] - (img.shape[0] - 1), k[1] - (img.shape[1] - 1))
    _, fx, fz, ox, oz = best
    nx, ny = S['q'].shape
    # coordonnees CAD de chaque position du scan et de chaque faisceau
    i = np.arange(nx)
    j = np.arange(ny)
    xi = xs[0] + ((nx - 1 - i) if fx < 0 else i) * ds + ox * ds
    zj = zs[0] + ((ny - 1 - j) if fz < 0 else j) * di + oz * di
    S.update(fx=fx, fz=fz, x_cad=xi, z_cad=zj)
    return S


def iou(S, C, seuil=0.5):
    """Recouvrement (intersection / union) entre les lettres mesurees et le CAD, sur la zone balayee."""
    ix = np.clip(np.round((S['x_cad'] - C['xs'][0]) / DX_CAD).astype(int), 0, len(C['xs']) - 1)
    iz = np.clip(np.round((S['z_cad'] - C['zs'][0]) / DX_CAD).astype(int), 0, len(C['zs']) - 1)
    ref = C['M'][np.ix_(ix, iz)]
    ok = ~np.isnan(S['q']) & (S['z_cad'] > -11.5)[None, :] & (S['z_cad'] < 11.5)[None, :]
    mes = (np.nan_to_num(S['q']) >= seuil) & ok
    ref = ref & ok
    return (mes & ref).sum() / max((mes | ref).sum(), 1)


# =====================================================================  mesures
def segmenter(S):
    """Repere grossierement les sept lettres le long du balayage (point de depart de l'ajustement)."""
    bande = (S['z_cad'] > -9.7) & (S['z_cad'] < 10.4)       # faisceaux qui passent sous les lettres
    prof = np.nan_to_num(S['r'])[:, bande].max(1)
    on = prof >= 0.5
    runs, i = [], 0
    while i < len(on):
        if on[i]:
            j = i
            while j + 1 < len(on) and on[j + 1]:
                j += 1
            runs.append((i, j))
            i = j + 1
        else:
            i += 1
    runs = [ru for ru in runs if (ru[1] - ru[0] + 1) * S['ds'] > 2.0]
    if len(runs) != 7:
        raise RuntimeError('%s : %d lettres trouvees au lieu de 7' % (S['nom'], len(runs)))
    if S['fx'] < 0:
        runs = runs[::-1]                                    # ordre L A S E R A X
    return runs


def gabarit(C, sx, sz, k=None, LC=None):
    """Masque CAD (d'une lettre, ou de toutes) floute par un faisceau gaussien d'ecarts-types sx, sz (mm)."""
    M = C['M'].astype(float)
    if k is not None:
        keep = (C['xs'] >= LC[k]['s0'] - DX_CAD / 2) & (C['xs'] <= LC[k]['s1'] + DX_CAD / 2)
        M[~keep] = 0
    T = gaussian_filter(M, (sx / DX_CAD, sz / DX_CAD), mode='constant')
    return RegularGridInterpolator((C['xs'], C['zs']), T, bounds_error=False, fill_value=0.0)


def _modele(p, T, X, Z, cx, cz):
    ax, bx, az, bz = p
    return T(np.c_[(cx + (X - cx - bx) / ax).ravel(), (cz + (Z - cz - bz) / az).ravel()])


def _ajuster(q, T, X, Z, cx, cz, p0, libres=(True, True, True, True)):
    """Moindres carres sur (etirement, decalage) en x et en z ; amplitude et fond sont lineaires."""
    p0 = np.array(p0, float)
    lib = np.array(libres)
    ok = ~np.isnan(q).ravel()
    y = np.nan_to_num(q).ravel()[ok]

    def res(pl):
        p = p0.copy()
        p[lib] = pl
        m = _modele(p, T, X, Z, cx, cz)[ok]
        A = np.c_[m, np.ones_like(m)]
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        return A @ coef - y

    sol = least_squares(res, p0[lib], x_scale=np.array([0.1, 1.0, 0.1, 1.0])[lib], diff_step=1e-3)
    p = p0.copy()
    p[lib] = sol.x
    return p, np.sqrt(np.mean(sol.fun ** 2))


def flou(S, C):
    """Largeur du faisceau : ecarts-types du flou gaussien qui rend le mieux la carte mesuree a partir du CAD."""
    band = BANDE(S)
    X, Z = np.meshgrid(S['x_cad'], S['z_cad'][band], indexing='ij')
    q = S['q'][:, band]
    best = None
    for sx in np.arange(0.5, 4.01, 0.25):
        for sz in np.arange(0.25, 2.51, 0.25):
            T = gabarit(C, sx, sz)
            m = _modele((1, 0, 1, 0), T, X, Z, 0, 0)
            ok = ~np.isnan(q).ravel()
            A = np.c_[m[ok], np.ones(ok.sum())]
            coef, *_ = np.linalg.lstsq(A, q.ravel()[ok], rcond=None)
            e = np.sqrt(np.mean((A @ coef - q.ravel()[ok]) ** 2))
            if best is None or e < best[0]:
                best = (e, sx, sz)
    return best[1], best[2]


def lettres_mesurees(S, C, LC, sx, sz):
    """Pour chaque lettre, le contour CAD floute est etire et deplace jusqu'a superposer la carte mesuree.
    L'etirement donne la longueur (le long du balayage) et la hauteur (le long de la sonde) de la lettre."""
    runs = segmenter(S)
    band = BANDE(S)
    xi = S['x_cad']
    out = []
    for k, (a, b) in enumerate(runs):
        L = LC[k]
        cx, cz, W = L['centre'], L['zc'], L['largeur']
        tronq = (a == 0) or (b == len(xi) - 1)
        x0, x1 = sorted((xi[a], xi[b]))
        ax0 = max((x1 - x0 + S['ds']) / (W - 2.0), 0.2)       # le flou raccourcit la zone > 50 %
        bx0 = (x0 + x1) / 2 - cx
        if tronq:
            out.append(dict(tronquee=True))
            continue
        lo, hi = cx + bx0 - ax0 * (W / 2 + 2.2), cx + bx0 + ax0 * (W / 2 + 2.2)
        sel = (xi >= lo) & (xi <= hi) & (abs(xi) < 64)       # pas les bouts du bloc
        X, Z = np.meshgrid(xi[sel], S['z_cad'][band], indexing='ij')
        T = gabarit(C, sx, sz, k, LC)
        p, e = _ajuster(S['q'][np.ix_(sel, band)], T, X, Z, cx, cz, (ax0, bx0, 1.0, 0.0))
        out.append(dict(tronquee=False, ax=p[0], az=p[2], centre=cx + p[1], zc=cz + p[3],
                        largeur=p[0] * W, hauteur=p[2] * L['hauteur'], rms=e))
    return out


def lettres_cad(C):
    xs, zs, M = C['xs'], C['zs'], C['M']
    col = M.any(1)
    out, i = [], 0
    while i < len(col):
        if col[i]:
            j = i
            while j + 1 < len(col) and col[j + 1]:
                j += 1
            rows = np.nonzero(M[i:j + 1].any(0))[0]
            out.append(dict(s0=xs[i], s1=xs[j], largeur=xs[j] - xs[i] + DX_CAD, centre=(xs[i] + xs[j]) / 2,
                            hauteur=zs[rows[-1]] - zs[rows[0]] + DX_CAD, zc=(zs[rows[-1]] + zs[rows[0]]) / 2))
            i = j + 1
        else:
            i += 1
    return out


def _stat(x):
    """Mediane et ecart-type robuste (1,4826 x ecart absolu median) : insensible aux A-scans de bord."""
    m = np.median(x)
    return m, 1.4826 * np.median(abs(x - m)), x.size


def profondeurs(S):
    """Profondeur de l'echo de fond (plaque) et du dessus des lettres, sur les A-scans sans ambiguite."""
    q = S['q']
    plaque = (q < 0.05) & ~np.isnan(q)
    lettre = (q > 0.8) & (S['a_l'] > SEUIL_SIGNAL)
    plaque &= S['filled'][:, None]                            # pas de A-scan interpole
    lettre &= S['filled'][:, None]
    # deuxieme echo de fond (aller-retour de plus dans la plaque) : epaisseur sans le zero de profondeur
    a2, p2 = _echo(S['vol'], S['d'], 17.5, S['d'][-1])
    ok2 = plaque & (a2 > SEUIL_SIGNAL)
    fond, lettre_ = _stat(S['p_f'][plaque]), _stat(S['p_l'][lettre])
    return dict(fond=fond, lettre=lettre_, relief=lettre_[0] - fond[0], double=_stat(p2[ok2] - S['p_f'][ok2]),
                a_fond=np.median(S['a_f'][plaque]), a_lettre=np.median(S['a_l'][lettre]))


# =====================================================================  figures
def _carte(ax, S, C, lettre):
    """Carte de la chute de l'echo de fond dans les coordonnees du CAD, contour CAD superpose."""
    ox, oz = np.argsort(S['x_cad']), np.argsort(S['z_cad'])
    img = S['q'][np.ix_(ox, oz)]
    x, z = S['x_cad'][ox], S['z_cad'][oz]
    ext = [x[0] - S['ds'] / 2, x[-1] + S['ds'] / 2, z[0] - S['di'] / 2, z[-1] + S['di'] / 2]
    im = ax.imshow(img.T, origin='lower', extent=ext, aspect='equal', cmap='viridis', vmin=0, vmax=1,
                   interpolation='nearest')
    for l in C['loops']:
        ax.plot(l[:, 0], l[:, 1], color='white', lw=0.6)
    ax.text(0.006, 0.93, lettre, transform=ax.transAxes, color='white', va='top', fontsize=9)
    return im


def _axes_bloc(ax):
    ax.set_xlim(-72, 72)
    ax.set_ylim(-13, 13)
    ax.invert_yaxis()                                         # vue de dessus : les lettres se lisent
    ax.set_ylabel('$z$ (mm)')


def fig_cscans(SH, SE, C, path):
    fig, axs = plt.subplots(3, 1, figsize=(7.0, 4.1), sharex=True, gridspec_kw=dict(hspace=0.08))
    ext = [C['xs'][0], C['xs'][-1], C['zs'][0], C['zs'][-1]]
    axs[0].imshow(C['M'].T, origin='lower', extent=ext, aspect='equal', cmap='Greys', vmin=0, vmax=1.6)
    axs[0].add_patch(plt.Rectangle((-70, -12), 140, 24, fill=False, color='0.4', lw=0.8))
    axs[0].text(0.006, 0.93, '(a)', transform=axs[0].transAxes, va='top', fontsize=9)
    _carte(axs[1], SH, C, '(b)')
    im = _carte(axs[2], SE, C, '(c)')
    vides = SE['x_cad'][~SE['filled']]
    axs[2].plot(vides, np.full(vides.size, -12.4), '|', color='#F28E2B', ms=5, mew=1.0)
    for a in axs:
        _axes_bloc(a)
    axs[2].set_xlabel('$x$ (mm)')
    cb = fig.colorbar(im, ax=axs, shrink=0.9, pad=0.015)
    cb.set_label("Chute de l'écho de fond $q$")
    fig.savefig(path)
    plt.close(fig)


def fig_ascan(SH, SE, path):
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    for S, c, lab in ((SH, C_HORLOGE, 'horloge'), (SE, C_ENCODEUR, 'encodeur')):
        q = S['q']
        plaque = (q < 0.05) & S['filled'][:, None]
        lettre = (q > 0.8) & (S['a_l'] > SEUIL_SIGNAL) & S['filled'][:, None]
        ax.plot(S['d'], S['vol'][plaque].mean(0), '-', color=c, lw=1.2, label='Plaque, %s' % lab)
        ax.plot(S['d'], S['vol'][lettre].mean(0), '--', color=c, lw=1.2, label='Lettre, %s' % lab)
    for y, t in ((Y_BASE - Y_FACE, 'plaque'), (Y_HAUT - Y_FACE, 'lettre'),
                 (2 * (Y_BASE - Y_FACE), r'$2\times$ plaque')):
        ax.axvline(y, color='0.45', ls=':', lw=0.9)
        ax.text(y - 0.25, 97, t, fontsize=8, color='0.3', rotation=90, va='top', ha='right')
    ax.set_xlabel('Profondeur (mm)')
    ax.set_ylabel('Amplitude moyenne (% écran)')
    ax.set_xlim(0, SH['d'][-1])
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(0.62, 0.86))
    fig.savefig(path)
    plt.close(fig)


def fig_positions(SH, SE, LC, path):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.9),
                                 gridspec_kw=dict(width_ratios=(1.15, 1), wspace=0.32))
    xc = np.array([l['centre'] for l in LC])
    for S, c, m, lab in ((SH, C_HORLOGE, 'o', 'Horloge'), (SE, C_ENCODEUR, 's', 'Encodeur')):
        k = [i for i, l in enumerate(S['L']) if not l['tronquee']]
        s = np.array([S['L'][i]['centre'] for i in k])
        a1.plot(xc[k], s - xc[k], m + '-', color=c, ms=5, lw=0.8, label=lab)
        a2.plot(np.array(k) + (0.1 if m == 's' else -0.1), [S['L'][i]['ax'] for i in k], m, color=c, ms=5,
                label=lab)
    a1.axhline(0, color='0.5', lw=0.8, ls=':')
    a1.set_xlabel('Centre de la lettre dans le CAD (mm)')
    a1.set_ylabel('Centre mesuré $-$ CAD (mm)')
    a1.legend(frameon=False, loc='lower left')
    a2.axhline(1, color='0.5', lw=0.8, ls=':')
    a2.set_xticks(range(7))
    a2.set_xticklabels(LETTRES)
    a2.set_ylabel('Longueur mesurée / CAD')
    a2.set_ylim(0, 1.6)
    a1.text(-0.2, 1.03, '(a)', transform=a1.transAxes, fontsize=9)
    a2.text(-0.2, 1.03, '(b)', transform=a2.transAxes, fontsize=9)
    fig.savefig(path)
    plt.close(fig)


def fig_3d(S, path, elev=50, azim=-80):
    """Surface reconstruite : profondeur de l'echo le plus profond (dessus de la plaque ou d'une lettre)."""
    q = np.nan_to_num(S['q'])
    h = np.where(q >= 0.5, S['p_l'], S['p_f'])
    ix = np.nonzero(abs(S['x_cad']) < 66)[0]                  # sans les bouts du bloc
    iz = np.nonzero(abs(S['z_cad'] - 0.35) < 10.3)[0]         # sans les bords du bloc
    ox, oz = ix[np.argsort(S['x_cad'][ix])], iz[np.argsort(S['z_cad'][iz])]
    X, Z = np.meshgrid(S['x_cad'][ox], S['z_cad'][oz], indexing='ij')
    H = gaussian_filter(h[np.ix_(ox, oz)], (0.8, 0.5))
    fig = plt.figure(figsize=(9.0, 4.0))
    ax = fig.add_subplot(111, projection='3d')
    H = np.clip(H, 9, 15)
    ax.plot_surface(X, Z, H, cmap='viridis', vmin=9, vmax=15.5, rstride=1, cstride=1, linewidth=0,
                    antialiased=False, shade=True)
    ax.set_box_aspect((np.ptp(X), np.ptp(Z), 3.0 * 6))    # relief exagere 3 fois (un zoom couperait la surface)
    ax.set_xlim(X.min(), X.max())
    ax.set_zlim(9, 15)                                        # profondeur vue de la sonde : lettres en relief
    ax.set_ylim(10.5, -10.5)
    ax.set_xlabel('$x$ (mm)', labelpad=14)
    ax.set_ylabel('')
    ax.set_zlabel('')
    ax.set_zticks([9.5, 14.5])
    ax.set_yticks([])
    ax.tick_params(axis='z', pad=1)
    ax.tick_params(axis='y', pad=0)
    ax.view_init(elev=elev, azim=azim)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    fig.savefig(path)
    plt.close(fig)
    from PIL import Image                                     # une vue 3D laisse beaucoup de blanc autour
    im = Image.open(path).convert('RGB')
    a = np.asarray(im).min(2) < 245
    r, c = np.nonzero(a.any(1))[0], np.nonzero(a.any(0))[0]
    im.crop((max(c[0] - 15, 0), max(r[0] - 15, 0), c[-1] + 15, r[-1] + 15)).save(path)


# =====================================================================  main
def ajuster(L, LC):
    """Centre mesure = a * centre CAD + b (moindres carres), sur les lettres completes."""
    k = [i for i, l in enumerate(L) if not l['tronquee']]
    xc = np.array([LC[i]['centre'] for i in k])
    s = np.array([L[i]['centre'] for i in k])
    a, b = np.polyfit(xc, s, 1)
    res = s - (a * xc + b)
    return a, b, np.sqrt(np.mean(res ** 2)), res


def main():
    os.makedirs(FIG, exist_ok=True)
    C = cad()
    LC = lettres_cad(C)
    print('CAD : %d contours, %d lettres' % (len(C['loops']), len(LC)))
    for k, l in enumerate(LC):
        print('  %s  x = %7.2f a %7.2f  centre %7.2f  largeur %5.2f  hauteur %5.2f'
              % (LETTRES[k], l['s0'], l['s1'], l['centre'], l['largeur'], l['hauteur']))
    res = {}
    for nom, lab in SCANS[::-1]:                               # encodeur d'abord : il fixe la largeur du faisceau
        S = recaler(charger(nom), C)
        if lab == 'Encodeur':
            sx, sz = flou(S, C)
            print('\nflou du faisceau (encodeur) : sigma_x = %.2f mm, sigma_z = %.2f mm (FWHM %.1f et %.1f mm)'
                  % (sx, sz, 2.355 * sx, 2.355 * sz))
        S['L'] = lettres_mesurees(S, C, LC, sx, sz)
        S['prof'] = profondeurs(S)
        S['iou'] = iou(S, C)
        res[lab] = S
        I, P = S['info'], S['prof']
        print('\n=====', nom, '(%s)' % lab)
        print('  positions %d (pas %.3f mm, %.1f mm), faisceaux %d (pas %.2f mm), echantillons %d (%.3f mm)'
              % (S['q'].shape[0], S['ds'], (S['q'].shape[0] - 1) * S['ds'], S['q'].shape[1], S['di'],
                 S['vol'].shape[2], S['dz']))
        print('  vitesse %.0f m/s, frequence %.1f Hz (%s), retournements fx=%d fz=%d, x de %.1f a %.1f'
              % (S['v'], S['rate'], I['rate_mode'], S['fx'], S['fz'], S['x_cad'][0], S['x_cad'][-1]))
        if S['trous']:
            n_t = sum(S['trous'])
            print('  trous : %d series %s, %d positions vides sur %d (%.1f %%), plus long %d positions (%.2f mm)'
                  % (len(S['trous']), S['trous'], n_t, len(S['filled']), 100 * n_t / len(S['filled']),
                     max(S['trous']), max(S['trous']) * S['ds']))
            print('  vitesse max sans trou : %.1f mm/s' % (S['ds'] * S['rate']))
            vides = np.nonzero(~S['filled'])[0]
            print('  positions vides (x CAD) : %s' % np.round(S['x_cad'][vides], 1))
        print('  fond    %.3f +- %.3f mm (N=%d), amplitude mediane %.0f %%' % (P['fond'] + (P['a_fond'],)))
        print('  lettre  %.3f +- %.3f mm (N=%d), amplitude mediane %.0f %%' % (P['lettre'] + (P['a_lettre'],)))
        print('  relief  %.3f mm' % P['relief'])
        print('  2e echo - 1er echo : %.3f +- %.3f mm (N=%d)' % P['double'])
        print('  avec %.0f m/s : fond %.3f, lettre %.3f, relief %.3f, double %.3f' % (
            V_NOMINALE, P['fond'][0] * V_NOMINALE / S['v'], P['lettre'][0] * V_NOMINALE / S['v'],
            P['relief'] * V_NOMINALE / S['v'], P['double'][0] * V_NOMINALE / S['v']))
        print('  IoU lettres/CAD (recalage global, pas nominal) : %.3f' % S['iou'])
        a, b, rms, rr = ajuster(S['L'], LC)
        S['fit'] = (a, b, rms, rr)
        print('  centres : pente %.4f, rms residus %.2f mm, residus %s' % (a, rms, np.round(rr, 2)))
        for k, l in enumerate(S['L']):
            if l['tronquee']:
                print('  %s  coupee au debut ou a la fin du balayage' % LETTRES[k])
                continue
            print('  %s  centre %7.2f  longueur %5.2f (CAD %5.2f, ecart %+5.2f)  hauteur %5.2f (CAD %5.2f, '
                  'ecart %+5.2f)  etirements %.3f %.3f  rms %.3f'
                  % (LETTRES[k], l['centre'], l['largeur'], LC[k]['largeur'], l['largeur'] - LC[k]['largeur'],
                     l['hauteur'], LC[k]['hauteur'], l['hauteur'] - LC[k]['hauteur'], l['ax'], l['az'], l['rms']))
        if lab == 'Horloge':
            print('  vitesse de la main par lettre (pas nominal / etirement) : %s mm/s'
                  % [round(S['ds'] * S['rate'] / l['ax'], 1) for l in S['L'] if not l['tronquee']])
            k = [i for i, l in enumerate(S['L']) if not l['tronquee']]
            xc = np.array([LC[i]['centre'] for i in k])
            sc = np.array([S['L'][i]['centre'] for i in k])
            dt = np.diff(sc) / S['ds'] / S['rate']
            print('  vitesse de la main entre centres : %s mm/s (temps %s s)'
                  % (np.round(np.diff(xc) / dt, 1), np.round(dt, 2)))
            print('  duree du balayage : %.1f s' % (S['q'].shape[0] / S['rate']))

    SH, SE = res['Horloge'], res['Encodeur']
    fig_cscans(SH, SE, C, os.path.join(FIG, '01_cscans_cad.png'))
    fig_ascan(SH, SE, os.path.join(FIG, '02_ascans_profondeur.png'))
    fig_positions(SH, SE, LC, os.path.join(FIG, '03_positions_lettres.png'))
    fig_3d(SE, os.path.join(FIG, '04_surface_3d_encodeur.png'))
    fig_3d(SH, os.path.join(FIG, '05_surface_3d_horloge.png'))


if __name__ == '__main__':
    main()
