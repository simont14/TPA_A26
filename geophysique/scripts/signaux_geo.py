# -*- coding: utf-8 -*-
"""
Lecture des enregistrements de l'oscilloscope, moyenne des trois tirs de chaque
position et reperage automatique de la premiere arrivee. Partage entre
analyse_refraction.py et pointage_manuel.py.
"""
import json
import os
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCOPE_DIR = os.path.join(BASE, "data", "scope")
POINTAGE_JSON = os.environ.get("POINTAGE_JSON", os.path.join(BASE, "data", "pointage_manuel.json"))

# Distance bord de plaque -> geophone 1 pour chaque groupe de trois tirs
# (scope_0-2: 0 m, scope_3-5: 0,3 m, ...), voir data/Lab_Geo_Notes.txt
DISTANCES = [0.0, 0.3, 0.6, 0.9, 1.5, 2.1, 2.7, 3.3, 3.9, 5.0]
# Canal de l'oscilloscope -> (nom du geophone, distance centre de plaque -> geophone
# quand le bord de la plaque touche le geophone 1). Le centre de la plaque est a
# 13 cm de son bord et les geophones sont espaces de 10 cm.
GEOPHONES = {3: ("G1", 0.13), 1: ("G2", 0.23), 2: ("G3", 0.33)}

SMOOTH_MS = 0.5   # largeur de la moyenne glissante
SEUIL = 2.0       # seuil de detection du creux, en ecarts-types du bruit
DUREE_MS = 0.5    # duree minimale sous le seuil
SEUIL_POSITIF = 6.0  # bosse positive qui, si elle precede le creux, fait rejeter la mesure
BRUIT = 1.0       # niveau de retour pour remonter au debut du creux
DT_COMMUN = 0.05  # grille de temps commune pour moyenner les tirs (ms)


def load_scope(i):
    """Retourne t (ms) et le tableau des 4 canaux (V)."""
    d = np.loadtxt(os.path.join(SCOPE_DIR, f"scope_{i}.csv"), delimiter=",", skiprows=2)
    return d[:, 0] * 1e3, d[:, 1:]


def lisser(z):
    w = int(round(SMOOTH_MS / DT_COMMUN))
    return np.convolve(z, np.ones(w) / w, "same")


def tirs_position(position, canal):
    """Les trois tirs d'une position, chacun normalise par l'ecart-type de son
    bruit avant l'impact, sur une grille de temps commune (ms)."""
    t_c = np.arange(-24, 48, DT_COMMUN)
    tirs = []
    for j in range(3):
        t, v = load_scope(3 * position + j)
        y = v[:, canal - 1]
        pre = t < -1
        y = (y - y[pre].mean()) / y[pre].std()
        tirs.append(np.interp(t_c, t, y, left=np.nan, right=np.nan))
    tirs = np.array(tirs)
    ok = ~np.isnan(tirs).any(axis=0)
    return t_c[ok], tirs[:, ok]


def trace_moyenne(position, canal):
    """Moyenne des trois tirs, lissee puis renormalisee par le bruit."""
    t, tirs = tirs_position(position, canal)
    z = lisser(tirs.mean(axis=0))
    z /= z[(t < -1) & (t > -20)].std()
    return t, z


def premiere_sequence(masque, n):
    """Premier indice k tel que masque[k:k+n] est vrai partout, sinon None."""
    run = np.convolve(masque.astype(int), np.ones(n, int), "valid")
    idx = np.where(run == n)[0]
    return idx[0] if len(idx) else None


def temps_arrivee(t, z):
    """Debut du premier creux negatif. Retourne None si une bosse positive
    apparait avant le creux (premiere arrivee noyee dans le bruit)."""
    n = int(round(DUREE_MS / DT_COMMUN))
    k_neg = premiere_sequence((z < -SEUIL) & (t > 0), n)
    k_pos = premiere_sequence((z > SEUIL_POSITIF) & (t > 0), n)
    if k_neg is None or (k_pos is not None and k_pos < k_neg):
        return None
    k = k_neg
    while k > 0 and z[k] < -BRUIT:
        k -= 1
    return t[k]


def lire_pointage():
    """Pointage manuel sauvegarde par pointage_manuel.py, ou None s'il n'existe pas."""
    if not os.path.exists(POINTAGE_JSON):
        return None
    with open(POINTAGE_JSON, encoding="utf-8") as f:
        return json.load(f)
