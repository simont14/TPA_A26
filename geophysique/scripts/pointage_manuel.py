# -*- coding: utf-8 -*-
"""
Pointage manuel des premieres arrivees.

Affiche, pour chaque position de la plaque et chaque geophone (30 signaux), la
moyenne des trois tirs. Un clic gauche place le debut de l'onde. Le resultat est
sauvegarde dans data/pointage_manuel.json a chaque modification, et relu au
prochain lancement pour reprendre ou on etait.

Commandes:
  clic gauche         placer le debut de l'onde (dans un des deux graphiques)
  fleche droite / n   signal suivant (marque le signal courant comme verifie)
  fleche gauche / p   signal precedent
  r                   rejeter / accepter (trop de bruit)
  a                   revenir au pointage automatique
  i                   afficher / cacher les trois tirs individuels
  q ou fermer         sauvegarder et quitter
La barre d'outils (loupe, deplacement) reste utilisable pour zoomer. Quand la
loupe ou le deplacement est actif, les clics ne placent pas de point.
"""
import json
import os
import sys
import numpy as np
import matplotlib

# Les fleches et 'r' servent a la navigation dans cet outil, pas a matplotlib
for cle, touches in [("keymap.back", ["left", "c"]), ("keymap.forward", ["right", "v"]),
                     ("keymap.home", ["r"]), ("keymap.quit", ["q"])]:
    matplotlib.rcParams[cle] = [k for k in matplotlib.rcParams[cle] if k not in touches]
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from signaux_geo import (DISTANCES, GEOPHONES, POINTAGE_JSON, SMOOTH_MS, lisser,
                         lire_pointage, temps_arrivee, tirs_position)

ORDRE_GEOPHONES = sorted(GEOPHONES, key=lambda c: GEOPHONES[c][0])   # G1, G2, G3

# ---------------------------------------------------------------- donnees
print("Chargement des 30 signaux...")
signaux = []
for p, d in enumerate(DISTANCES):
    for canal in ORDRE_GEOPHONES:
        nom, x0 = GEOPHONES[canal]
        t, tirs = tirs_position(p, canal)
        bruit = (t < -1) & (t > -20)
        # toutes les courbes sont divisees par le meme facteur (bruit de la moyenne
        # lissee) pour garder la meme echelle d'amplitude
        brut = tirs.mean(axis=0)
        echelle = lisser(brut)[bruit].std()
        z = lisser(brut) / echelle
        brut = brut / echelle
        indiv = np.array([lisser(y) for y in tirs]) / echelle
        auto = temps_arrivee(t, z)
        signaux.append({
            "position": p, "d_m": d, "geophone": int(nom[1]), "canal": canal,
            "x_m": round(d + x0, 3), "t": t, "z": z, "brut": brut, "indiv": indiv,
            "auto_t_ms": None if auto is None else round(float(auto), 3),
        })

# Etat du pointage: repris du JSON s'il existe, sinon pointage automatique
etat = [{"t_ms": s["auto_t_ms"], "rejete": False, "verifie": False} for s in signaux]
ancien = lire_pointage()
if ancien is not None:
    par_cle = {(e["position"], e["geophone"]): e for e in ancien["pointes"]}
    for s, e in zip(signaux, etat):
        a = par_cle.get((s["position"], s["geophone"]))
        if a is not None:
            e.update(t_ms=a["t_ms"], rejete=a["rejete"], verifie=a["verifie"])
    print(f"Pointage existant repris: {POINTAGE_JSON}")


def sauvegarder():
    pointes = []
    for s, e in zip(signaux, etat):
        pointes.append({
            "position": s["position"], "d_m": s["d_m"], "geophone": s["geophone"],
            "canal": s["canal"], "x_m": s["x_m"],
            "t_ms": None if e["t_ms"] is None else round(float(e["t_ms"]), 3),
            "rejete": e["rejete"], "verifie": e["verifie"], "auto_t_ms": s["auto_t_ms"],
        })
    contenu = {
        "description": "Temps de premiere arrivee pointes a la main sur la moyenne des "
                       "trois tirs de chaque position (t en ms apres l'impact, x en m "
                       "entre le centre de la plaque et le geophone).",
        "lissage_ms": SMOOTH_MS,
        "pointes": pointes,
    }
    with open(POINTAGE_JSON, "w", encoding="utf-8") as f:
        json.dump(contenu, f, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------- interface
fig = plt.figure(figsize=(13, 8.5))
ax_vue = fig.add_axes([0.06, 0.55, 0.90, 0.33])
ax_zoom = fig.add_axes([0.06, 0.14, 0.90, 0.33])
texte = fig.text(0.06, 0.95, "", fontsize=13, weight="bold")
aide = fig.text(0.06, 0.925, "Clic gauche: début de l'onde | ←/→: précédent/suivant | "
                "r: rejeter | a: auto | i: tirs individuels | q: quitter",
                fontsize=9, color="0.35")
courant = {"i": 0, "indiv": False}


def bouton(x, label, action):
    b = Button(fig.add_axes([x, 0.03, 0.13, 0.055]), label)
    b.on_clicked(lambda evt: action())
    return b


def dessiner(garder_zoom=False):
    i = courant["i"]
    s, e = signaux[i], etat[i]
    t, z = s["t"], s["z"]
    lim_vue = ax_vue.get_xlim() if garder_zoom else (-2, 25)
    lim_zoom = ax_zoom.get_xlim() if garder_zoom else None

    for ax in (ax_vue, ax_zoom):
        ax.clear()
        ax.grid(alpha=0.3)
        ax.axhline(0, color="0.6", lw=0.6)
        ax.axvline(0, color="0.6", lw=0.6, ls=":")
        if courant["indiv"]:
            for y in s["indiv"]:
                ax.plot(t, y, color="0.75", lw=0.6)
        ax.plot(t, s["brut"], color="0.55", lw=0.5)
        ax.plot(t, z, color="C0", lw=1.4)
        if s["auto_t_ms"] is not None:
            ax.axvline(s["auto_t_ms"], color="C2", ls="--", lw=1)
        if e["t_ms"] is not None and not e["rejete"]:
            ax.axvline(e["t_ms"], color="C3", lw=1.8)

    # vue d'ensemble
    fen = (t > lim_vue[0]) & (t < lim_vue[1])
    amp = max(6.0, 1.1 * np.abs(z[fen]).max())
    ax_vue.set_xlim(*lim_vue)
    ax_vue.set_ylim(-amp, amp)
    ax_vue.set_ylabel("Signal / bruit")

    # zoom autour du point (ou du pointage auto, ou du debut de la fenetre)
    centre = e["t_ms"] if e["t_ms"] is not None else (s["auto_t_ms"] or 3.0)
    if lim_zoom is None:
        lim_zoom = (centre - 4, centre + 3)
    fen = (t > lim_zoom[0]) & (t < centre + 0.5)
    amp = max(6.0, 1.4 * np.abs(z[fen]).max()) if fen.any() else 6.0
    for niv in (-3, -1, 1, 3):
        ax_zoom.axhline(niv, color="0.7", lw=0.6, ls="--")
    ax_zoom.set_xlim(*lim_zoom)
    ax_zoom.set_ylim(-amp, amp)
    ax_zoom.set_xlabel("Temps après l'impact (ms)")
    ax_zoom.set_ylabel("Signal / bruit (zoom)")

    if e["rejete"]:
        statut, couleur = "REJETÉ (trop de bruit)", "C3"
    elif e["t_ms"] is None:
        statut, couleur = "aucun point: cliquer ou rejeter", "C1"
    else:
        statut, couleur = f"début = {e['t_ms']:.2f} ms", "k"
    n_verif = sum(x["verifie"] for x in etat)
    texte.set_text(f"[{i + 1}/30]  plaque à d = {s['d_m']:.1f} m, géophone {s['geophone']} "
                   f"(x = {s['x_m']:.2f} m)   |   {statut}   |   vérifiés: {n_verif}/30"
                   + ("   ✔" if e["verifie"] else ""))
    texte.set_color(couleur)
    fig.canvas.draw_idle()


def aller(pas):
    etat[courant["i"]]["verifie"] = True
    courant["i"] = (courant["i"] + pas) % len(signaux)
    sauvegarder()
    dessiner()


def rejeter():
    e = etat[courant["i"]]
    e["rejete"] = not e["rejete"]
    e["verifie"] = True
    sauvegarder()
    dessiner(garder_zoom=True)


def auto():
    e = etat[courant["i"]]
    e["t_ms"], e["rejete"] = signaux[courant["i"]]["auto_t_ms"], False
    sauvegarder()
    dessiner()


def basculer_indiv():
    courant["indiv"] = not courant["indiv"]
    dessiner(garder_zoom=True)


def quitter():
    sauvegarder()
    plt.close(fig)


def sur_clic(evt):
    if evt.button != 1 or evt.inaxes not in (ax_vue, ax_zoom) or evt.xdata is None:
        return
    barre = fig.canvas.toolbar
    if barre is not None and barre.mode:
        return
    e = etat[courant["i"]]
    e["t_ms"], e["rejete"], e["verifie"] = float(evt.xdata), False, True
    sauvegarder()
    dessiner(garder_zoom=True)


def sur_touche(evt):
    actions = {"right": lambda: aller(1), "n": lambda: aller(1),
               "left": lambda: aller(-1), "p": lambda: aller(-1),
               "r": rejeter, "a": auto, "i": basculer_indiv, "q": quitter}
    if evt.key in actions:
        actions[evt.key]()


boutons = [bouton(0.06 + 0.15 * k, label, action) for k, (label, action) in enumerate([
    ("◀ Précédent", lambda: aller(-1)), ("Suivant ▶", lambda: aller(1)),
    ("Rejeter / accepter", rejeter), ("Pointage auto", auto),
    ("Tirs individuels", basculer_indiv), ("Terminer", quitter)])]
fig.canvas.mpl_connect("button_press_event", sur_clic)
fig.canvas.mpl_connect("key_press_event", sur_touche)
fig.canvas.mpl_connect("close_event", lambda evt: sauvegarder())
fig.legend(handles=[plt.Line2D([], [], color="C0", lw=1.4, label="Moyenne des 3 tirs (lissée)"),
                    plt.Line2D([], [], color="0.55", lw=0.5, label="Moyenne non lissée"),
                    plt.Line2D([], [], color="C2", ls="--", label="Pointage automatique"),
                    plt.Line2D([], [], color="C3", lw=1.8, label="Ton point")],
           loc="lower right", bbox_to_anchor=(0.96, 0.882), ncol=4, fontsize=9, frameon=False)

dessiner()
plt.show()

restants = [f"d={s['d_m']:.1f} m G{s['geophone']}" for s, e in zip(signaux, etat) if not e["verifie"]]
print(f"Pointage sauvegarde: {POINTAGE_JSON}")
if restants:
    print(f"{len(restants)} signal(aux) pas encore verifie(s): " + ", ".join(restants))
else:
    print("Les 30 signaux ont ete verifies.")
