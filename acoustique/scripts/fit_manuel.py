# -*- coding: utf-8 -*-
"""
Outil pour choisir, sur le spectrogramme, les trois regions qui servent a
mesurer v et d pour chaque passage (methode decrite dans doppler_commun.py).

    python fit_manuel.py

  - region bleue (approche) et region orange (eloignement): les queues de la
    courbe, la ou la frequence se stabilise. Elles donnent la vitesse.
  - region blanche (centre): la chute. Elle donne la distance.
v et d sont recalcules des qu'une region change (au relachement du curseur).
Les regions sauvegardees sont celles du rapport: analyse_doppler.py refait le
meme calcul a partir de regions_manuelles.json.

Courbes: cyan = modele avec v, d et f0 mesures, vert pointille = theorie avec
la vitesse visee, le d du galon et 440 Hz. Lignes horizontales cyan: plateaux
f0/(1 -+ v/c) que la courbe atteindrait loin du micro.

Curseurs: clic ou glisser, ou taper la valeur dans la boite a droite puis Entree.
Clavier (cliquer dans le graphique d'abord):
  haut/bas      choisir le curseur actif (etiquette en rouge)
  gauche/droite ajuster le curseur actif d'un pas (Maj: 10 pas)
  c             recentrer la vue sur le passage
  PgUp/PgDn     passage precedent/suivant
"""
import os
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, CheckButtons, Slider, TextBox

from doppler_commun import (C, CLES, F0, REGIONS_JSON, analyse_regions, charger_passages,
                            ecrire_regions, lire_regions, modele, position, regions_depart, selections)

for cle, touches in (("keymap.back", ("left", "c")), ("keymap.forward", ("right",))):
    plt.rcParams[cle] = [k for k in plt.rcParams[cle] if k not in touches]

print("Lecture des passages et regions de depart...")
PASSAGES = charger_passages()
DEPART = {}
for p in PASSAGES:
    DEPART[p["nom"]], p["tref"] = regions_depart(p)     # tref: zero de l'axe du temps
SAUVEES = lire_regions()
REGIONS = {p["nom"]: dict(SAUVEES.get(p["nom"], DEPART[p["nom"]])) for p in PASSAGES}

# curseur: (etiquette, pas du clavier, format de la boite). Temps relatifs a tref.
CURSEURS = {
    "a0": ("approche: début", 0.01, "%.2f"),
    "a1": ("approche: fin", 0.01, "%.2f"),
    "c0": ("centre: début", 0.01, "%.2f"),
    "c1": ("centre: fin", 0.01, "%.2f"),
    "e0": ("éloignement: début", 0.01, "%.2f"),
    "e1": ("éloignement: fin", 0.01, "%.2f"),
    "zoom": ("largeur vue ± (s)", 0.1, "%.1f"),
}
ORDRE = list(CURSEURS)
COUL = {"a": "deepskyblue", "c": "white", "e": "orange"}


def spectrogramme(p):
    """Spectrogramme de tout l'enregistrement, normalise instant par instant
    (0 dB au maximum de chaque trame) pour que la raie du sinus reste visible
    meme loin du micro. La trame est plus longue pour les passages lents, ou la
    chute de frequence est petite et lente."""
    x, fs, tref = p["x"].astype(float), p["fs"], p["tref"]
    duree = 2 * p["d"] / (max(p["vn"], 5) / 3.6)          # duree approximative de la chute [s]
    n = int(np.clip(2 ** np.round(np.log2(duree / 6 * fs)), 8192, 32768))
    hop = n // 16
    seg = x - x.mean()
    deb = np.arange(0, len(seg) - n, hop)
    fr = np.fft.rfftfreq(4 * n, 1 / fs)
    m = (fr >= 380) & (fr <= 500)
    w = np.hanning(n)
    vues = np.lib.stride_tricks.sliding_window_view(seg, n)
    S = np.empty((len(deb), m.sum()))
    for k in range(0, len(deb), 100):
        S[k:k + 100] = np.abs(np.fft.rfft(vues[deb[k:k + 100]] * w, 4 * n, axis=1))[:, m]
    S = 20 * np.log10(S + 1e-9)
    S -= S.max(axis=1, keepdims=True)
    tt = (deb + n / 2) / fs - tref
    return (tt[0], tt[-1], fr[m][0], fr[m][-1]), S.T


class Outil:
    def __init__(self):
        self.i = 0
        self.actif = "a0"
        self.bloque = False
        self.modifie = False
        self.a_calculer = False
        self.voir_points = False
        self.spectros = {}
        self.res = None
        self.spans = []
        self.demi = 3.0

        self.fig = plt.figure(figsize=(14, 9))
        self.fig.canvas.manager.set_window_title("Doppler: régions pour v et d")
        self.ax = self.fig.add_axes([0.06, 0.38, 0.60, 0.57])

        self.sl, self.tb = {}, {}
        for k, nom in enumerate(ORDRE):
            lab, pas, fmt = CURSEURS[nom]
            y = 0.30 - 0.038 * k
            s = Slider(self.fig.add_axes([0.15, y, 0.41, 0.024]), lab, 0.0, 1.0, valinit=0.5)
            s.valtext.set_visible(False)
            if nom[0] in COUL:
                s.poly.set_color({"a": "deepskyblue", "c": "0.55", "e": "orange"}[nom[0]])
            s.on_changed(lambda val, nom=nom: self.curseur(nom, val))
            self.sl[nom] = s
            tb = TextBox(self.fig.add_axes([0.58, y - 0.003, 0.08, 0.03]), "", initial="")
            tb.on_submit(lambda txt, nom=nom: self.texte(nom, txt))
            self.tb[nom] = tb

        self.info = self.fig.text(0.70, 0.95, "", va="top", family="monospace", fontsize=9)
        self.check = CheckButtons(self.fig.add_axes([0.70, 0.40, 0.13, 0.045]), ["points suivis"], [False])
        self.check.on_clicked(self.bascule_points)

        def bouton(rect, texte, fn):
            b = Button(self.fig.add_axes(rect), texte)
            b.on_clicked(lambda _: fn())
            return b

        self.btn = [
            bouton([0.70, 0.33, 0.13, 0.045], "< Précédent", lambda: self.va(-1)),
            bouton([0.84, 0.33, 0.13, 0.045], "Suivant >", lambda: self.va(+1)),
            bouton([0.70, 0.27, 0.13, 0.045], "Recentrer (c)", self.recentre),
            bouton([0.84, 0.27, 0.13, 0.045], "Régions de départ", self.reinit),
            bouton([0.70, 0.21, 0.13, 0.045], "Sauvegarder", self.sauve),
            bouton([0.84, 0.21, 0.13, 0.045], "Sauver + figures", self.genere),
        ]
        self.fig.text(0.70, 0.17,
                      "Bleu: approche, blanc: centre, orange: éloignement\n"
                      "Boîtes de texte: taper la valeur puis Entrée\n"
                      "Clavier (cliquer dans le graphique d'abord):\n"
                      "  haut/bas: choisir le curseur\n"
                      "  gauche/droite: ajuster (Maj: x10)\n"
                      "  c: recentrer   PgUp/PgDn: passage",
                      va="top", fontsize=9)
        self.fig.canvas.mpl_connect("key_press_event", self.touche)
        self.fig.canvas.mpl_connect("button_release_event", self.relache)
        self.fig.canvas.mpl_connect("close_event", self.fermeture)
        self.dessine_passage()

    # ------------------------------------------------------------ etat
    @property
    def p(self):
        return PASSAGES[self.i]

    @property
    def R(self):
        return REGIONS[self.p["nom"]]

    def plages(self):
        p = self.p
        a, b = p["t"].min() - p["tref"], p["t"].max() - p["tref"]
        pl = {k: (a, b) for k in CLES}
        pl["zoom"] = (0.5, max(-a, b))
        return pl

    def calcule(self):
        self.res = analyse_regions(self.p, self.R, tp0=self.res["tp"] if self.res else None)
        self.a_calculer = False

    # ------------------------------------------------------------ dessin
    def dessine_passage(self):
        p = self.p
        ax = self.ax
        ax.clear()
        self.spans = []
        if p["nom"] not in self.spectros:
            self.spectros[p["nom"]] = spectrogramme(p)
        ext, S = self.spectros[p["nom"]]
        ax.imshow(S, origin="lower", aspect="auto", extent=ext, cmap="magma", vmin=-10, vmax=0,
                  interpolation="bilinear")
        self.pts_hors, = ax.plot([], [], ".", ms=2, color="0.6")
        self.pts_in, = ax.plot([], [], ".", ms=2.5, color="lime")
        self.theorie, = ax.plot([], [], "--", color="lime", lw=1.3,
                                label=f"théorie: {p['vn']} km/h, d = {p['d']:.2f} m, 440 Hz")
        self.courbe, = ax.plot([], [], "-", color="cyan", lw=1.5, label="modèle avec v, d, f0 mesurés")
        self.pa = ax.axhline(np.nan, color="cyan", lw=0.8, ls="--", alpha=0.8)
        self.pe = ax.axhline(np.nan, color="cyan", lw=0.8, ls="--", alpha=0.8)
        ax.axhline(F0, color="w", ls=":", lw=0.6, alpha=0.6)
        demi_f = max(1.6 * F0 * (p["vn"] / 3.6) / C, 10.0)
        ax.set_ylim(F0 - demi_f, F0 + demi_f)
        ax.set_ylabel("Fréquence (Hz)")
        ax.set_xlabel("Temps (s)")
        ax.legend(loc="upper right", fontsize=8)
        ax.set_title(f"Passage {self.i + 1}/{len(PASSAGES)}: {p['vn']} km/h visés, "
                     f"d = {p['d']:.2f} m   ({p['nom']})", fontsize=10)
        R = self.R
        self.demi = max(2.0, 1.2 * max(p["tref"] - R["a0"], R["e1"] - p["tref"]))
        for nom, (a, b) in self.plages().items():
            sl = self.sl[nom]
            sl.valmin, sl.valmax = a, b
            sl.ax.set_xlim(a, b)
        self.res = None
        self.calcule()
        self.synchro()
        self.maj()
        self.recentre()

    def recentre(self):
        c = (self.res["tp"] if self.res else self.p["tref"]) - self.p["tref"]
        self.ax.set_xlim(c - self.demi, c + self.demi)
        self.fig.canvas.draw_idle()

    def synchro(self):
        tref, R = self.p["tref"], self.R
        valeurs = {k: R[k] - tref for k in CLES}
        valeurs["zoom"] = self.demi
        self.bloque = True
        for nom, val in valeurs.items():
            sl = self.sl[nom]
            sl.set_val(float(np.clip(val, sl.valmin, sl.valmax)))
        self.bloque = False
        for nom in ORDRE:
            self.maj_texte(nom)
        self.colore()

    def maj_texte(self, nom):
        tb = self.tb[nom]
        if tb.capturekeystrokes:        # ne pas ecraser ce que l'utilisateur tape
            return
        bloque, self.bloque = self.bloque, True
        tb.set_val(CURSEURS[nom][2] % self.sl[nom].val)
        self.bloque = bloque

    def texte(self, nom, txt):
        if self.bloque:
            return
        sl = self.sl[nom]
        try:
            val = float(txt.strip().replace(",", ".").replace("±", ""))
        except ValueError:
            self.maj_texte(nom)
            return
        val = float(np.clip(val, sl.valmin, sl.valmax))
        if val == sl.val:
            self.maj_texte(nom)
        else:
            sl.set_val(val)
        self.relache(None)

    def colore(self):
        for nom, s in self.sl.items():
            s.label.set_color("C3" if nom == self.actif else "k")
            s.label.set_fontweight("bold" if nom == self.actif else "normal")

    def maj(self):
        p, R, r = self.p, self.R, self.res
        t, f, tref = p["t"], p["f"], p["tref"]
        for sp in self.spans:
            sp.remove()
        self.spans = [self.ax.axvspan(R[z + "0"] - tref, R[z + "1"] - tref, color=COUL[z],
                                      alpha=0.13 if z != "c" else 0.10, lw=0) for z in "ace"]
        for z in "ace":
            for k in "01":
                self.spans.append(self.ax.axvline(R[z + k] - tref, color=COUL[z], lw=1))
        sa, sc, se = selections(p, R)
        w = sa | sc | se
        vis = np.full(len(t), self.voir_points)
        self.pts_hors.set_data(t[vis & ~w] - tref, f[vis & ~w])
        self.pts_in.set_data(t[vis & w] - tref, f[vis & w])
        tt = np.linspace(t.min(), t.max(), 4000)
        if r is None:
            self.courbe.set_data([], [])
            self.theorie.set_data([], [])
            self.pa.set_ydata([np.nan] * 2)
            self.pe.set_ydata([np.nan] * 2)
            lignes = ["Pas assez de points dans une région", "(au moins 5 par région)."]
        else:
            v, d, tp, f0 = r["v"], r["d"], r["tp"], r["f0"]
            self.courbe.set_data(tt - tref, modele(tt, v, d, tp, f0))
            self.theorie.set_data(tt - tref, modele(tt, p["vn"] / 3.6, p["d"], tp))
            self.pa.set_ydata([f0 / (1 - v / C)] * 2)
            self.pe.set_ydata([f0 / (1 + v / C)] * 2)
            xs = {k: position(np.array(R[k]), v, d, tp) for k in CLES}
            lignes = [
                "Vitesse (queues bleue et orange):",
                f"  v  = {v * 3.6:6.1f} ± {r['ev'] * 3.6:.1f} km/h  (visée {p['vn']})",
                f"  plateaux mesurés: {r['fa']:.1f} / {r['fe']:.1f} Hz",
                f"  f0 de la source:  {f0:.1f} Hz",
                "Distance (centre blanc):",
                f"  d  = {d:6.2f} ± {r['ed']:.2f} m    (galon {p['d']:.2f})",
                "",
                f"Écart quadratique: queues {r['rms_q']:.2f} Hz,",
                f"                   centre {r['rms_c']:.2f} Hz",
                f"Points: {r['na']} / {r['nc']} / {r['ne']}",
                "Position de la voiture aux bords (m):",
                f"  approche     {xs['a0']:+6.1f} à {xs['a1']:+6.1f}",
                f"  centre       {xs['c0']:+6.1f} à {xs['c1']:+6.1f}",
                f"  éloignement  {xs['e0']:+6.1f} à {xs['e1']:+6.1f}",
            ]
            if not r["conv"]:
                lignes.append("ATTENTION: v et d n'ont pas convergé")
        lignes += ["", "Modifié, pas sauvegardé" if self.modifie else "Sauvegardé"]
        self.info.set_text("\n".join(lignes))
        self.fig.canvas.draw_idle()

    # ------------------------------------------------------------ actions
    def curseur(self, nom, val):
        if self.bloque:
            return
        self.actif = nom
        self.colore()
        self.maj_texte(nom)
        if nom == "zoom":
            self.demi = float(val)
            self.recentre()
            return
        self.R[nom] = float(self.p["tref"] + val)
        self.modifie = True
        self.a_calculer = True
        self.maj()              # deplace les regions tout de suite, le calcul vient au relachement

    def relache(self, _):
        if self.a_calculer:
            self.calcule()
            self.maj()

    def bascule_points(self, _):
        self.voir_points = not self.voir_points
        self.maj()

    def va(self, pas):
        self.relache(None)
        self.i = (self.i + pas) % len(PASSAGES)
        self.dessine_passage()

    def reinit(self):
        REGIONS[self.p["nom"]] = dict(DEPART[self.p["nom"]])
        self.modifie = True
        self.calcule()
        self.synchro()
        self.maj()

    def sauve(self):
        ecrire_regions(REGIONS)
        self.modifie = False
        print(f"Sauvegardé dans {REGIONS_JSON}")
        self.maj()

    def genere(self):
        self.sauve()
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "analyse_doppler.py")
        print("Génération des figures et du tableau...")
        r = subprocess.run([sys.executable, script], capture_output=True, text=True)
        print(r.stdout, r.stderr)

    def touche(self, ev):
        k = ev.key or ""
        if any(tb.capturekeystrokes for tb in self.tb.values()):   # en train de taper dans une boite
            return
        if k in ("up", "down"):
            j = ORDRE.index(self.actif)
            self.actif = ORDRE[(j + (1 if k == "down" else -1)) % len(ORDRE)]
            self.colore()
            self.fig.canvas.draw_idle()
        elif k in ("left", "right", "shift+left", "shift+right"):
            s = self.sl[self.actif]
            pas = CURSEURS[self.actif][1] * (10 if k.startswith("shift") else 1)
            s.set_val(float(np.clip(s.val + (pas if k.endswith("right") else -pas), s.valmin, s.valmax)))
            self.relache(None)
        elif k == "c":
            self.recentre()
        elif k == "pagedown":
            self.va(+1)
        elif k == "pageup":
            self.va(-1)

    def fermeture(self, _):
        if self.modifie:
            ecrire_regions(REGIONS)
            print(f"Sauvegardé à la fermeture dans {REGIONS_JSON}")


if __name__ == "__main__":
    outil = Outil()
    plt.show()
