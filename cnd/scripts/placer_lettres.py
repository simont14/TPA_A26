#!/usr/bin/env python3
"""
placer_lettres.py : outil (a usage unique) pour placer a la main le contour CAD de chaque lettre sur les deux
acquisitions du bloc LASERAX (horloge et encodeur).

Chaque lettre a quatre parametres, independants d'une lettre a l'autre, ce qui permet une distorsion qui change
le long du balayage :
    x_mesure = cx + bx + ax * (x_cad - cx)        (ax : etirement le long du balayage, bx : decalage)
    z_mesure = cz + bz + az * (z_cad - cz)        (az : etirement le long de la sonde,  bz : decalage)
ou (cx, cz) est le centre de la lettre dans le CAD. Au depart, les parametres viennent de l'ajustement automatique
d'analyse_scans.py (ou du fichier deja enregistre, si on reprend le travail).

Souris
    glisser (bouton gauche)        deplace la lettre choisie
    clic sur une lettre            la choisit
    molette                        etire / comprime le long du balayage (x)
    Maj + molette                  etire / comprime le long de la sonde (z)
Clavier
    fleches                        deplacent de 0,1 mm (Maj : 1 mm)
    a / d                          etirement x -1 % / +1 %   (Maj : 0,2 %)
    s / w                          etirement z -1 % / +1 %   (Maj : 0,2 %)
    n / p                          lettre suivante / precedente
    z                              zoom sur la lettre / vue complete
    Ctrl + s                       enregistrer

Le resultat est enregistre dans ../data/Data_Session_2/placement_lettres.json (aussi a la fermeture de la fenetre).
analyse_scans.py le relit pour tracer les contours et mesurer les lettres.

Lancement :  python placer_lettres.py
"""
import json
import os
import datetime

import numpy as np
import analyse_scans as A                # impose le mode Agg (figures sans fenetre) a son import
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons, Slider, Button, CheckButtons

plt.switch_backend('TkAgg')              # cet outil a besoin d'une vraie fenetre

SORTIE = os.path.join(A.DATA, 'placement_lettres.json')
SIGMA_X, SIGMA_Z = 1.5, 1.0          # flou du faisceau trouve par analyse_scans.flou() sur l'acquisition encodeur
NOMS = [lab for _, lab in A.SCANS]   # 'Horloge', 'Encodeur'

for k in [k for k in plt.rcParams if k.startswith('keymap.')]:
    plt.rcParams[k] = []             # les raccourcis de matplotlib entrent en conflit avec ceux de l'outil


def reg(S):
    return dict(x_cad0=float(S['x_cad'][0]), dx=float(S['x_cad'][1] - S['x_cad'][0]),
                z_cad0=float(S['z_cad'][0]), dz=float(S['z_cad'][1] - S['z_cad'][0]))


def depart(S, C, LC, lab, sauve):
    """Parametres de depart : fichier enregistre si la grille correspond, sinon ajustement automatique."""
    if sauve and lab in sauve.get('scans', {}):
        d = sauve['scans'][lab]
        if all(abs(d['recalage'][k] - v) < 1e-6 for k, v in reg(S).items()):
            return [dict(ax=l['ax'], bx=l['bx'], az=l['az'], bz=l['bz'], tronquee=l['tronquee']) for l in d['lettres']]
        print('  %s : le recalage a change depuis l\'enregistrement, on repart de l\'ajustement automatique.' % lab)
    auto = A.lettres_mesurees(S, C, LC, SIGMA_X, SIGMA_Z)
    out = []
    for L, a in zip(LC, auto):
        if a['tronquee']:
            out.append(dict(ax=1.0, bx=0.0, az=1.0, bz=0.0, tronquee=True))
        else:
            out.append(dict(ax=a['ax'], bx=a['centre'] - L['centre'], az=a['az'], bz=a['zc'] - L['zc'], tronquee=False))
    return out


class Outil:
    def __init__(self):
        print('Lecture du CAD et des deux acquisitions...')
        self.C = A.cad()
        self.LC = A.lettres_cad(self.C)
        self.grp = A.contours_par_lettre(self.C, self.LC)
        sauve = json.load(open(SORTIE, encoding='utf-8')) if os.path.exists(SORTIE) else None
        self.S, self.P, self.P0 = {}, {}, {}
        for nom, lab in A.SCANS:
            S = A.recaler(A.charger(nom), self.C)
            self.S[lab] = S
            self.P0[lab] = depart(S, self.C, self.LC, lab, None)
            self.P[lab] = depart(S, self.C, self.LC, lab, sauve)
            print('  %s pret' % lab)
        self.lab, self.k, self.zoom = 'Encodeur', 0, False
        self.drag = None
        self.construire()

    # ------------------------------------------------------------ fenetre
    def construire(self):
        self.fig = plt.figure(figsize=(14, 7.6))
        self.fig.canvas.manager.set_window_title('Placement des lettres sur les acquisitions')
        self.ax = self.fig.add_axes([0.05, 0.36, 0.93, 0.6])
        self.im = None
        self.lignes = [[self.ax.plot([], [], color='white', lw=0.9, ls=(0, (3, 2)))[0] for _ in g] for g in self.grp]
        self.info = self.fig.text(0.05, 0.315, '', family='monospace', fontsize=9.5)
        self.fig.text(0.05, 0.015, 'Glisser : déplacer   molette : étirer en x   Maj + molette : étirer en z   '
                      'flèches : 0,1 mm (Maj : 1 mm)   a/d, s/w : étirements   n/p : lettre   z : zoom   Ctrl+s : enregistrer',
                      fontsize=8.5, color='0.35')

        axs = self.fig.add_axes([0.05, 0.08, 0.09, 0.2], title='Acquisition')
        self.r_scan = RadioButtons(axs, NOMS, active=NOMS.index(self.lab))
        self.r_scan.on_clicked(self.choisir_scan)
        axl = self.fig.add_axes([0.155, 0.08, 0.07, 0.2], title='Lettre')
        self.r_let = RadioButtons(axl, ['%d %s' % (i + 1, c) for i, c in enumerate(A.LETTRES)], active=0)
        self.r_let.on_clicked(lambda t: self.choisir_lettre(int(t.split()[0]) - 1))

        def sl(y, nom, a, b):
            s = Slider(self.fig.add_axes([0.33, y, 0.36, 0.03]), nom, a, b, valinit=0)
            s.on_changed(self.slider)
            return s
        self.s_bx = sl(0.25, 'Décalage x (mm)', -30, 30)
        self.s_ax = sl(0.20, 'Étirement x', 0.3, 1.8)
        self.s_bz = sl(0.15, 'Décalage z (mm)', -5, 5)
        self.s_az = sl(0.10, 'Étirement z', 0.7, 1.3)

        axc = self.fig.add_axes([0.73, 0.19, 0.12, 0.08])
        self.c_tr = CheckButtons(axc, ['Lettre coupée'], [False])
        self.c_tr.on_clicked(self.coupee)
        self.b_zoom = Button(self.fig.add_axes([0.73, 0.10, 0.12, 0.05]), 'Zoom lettre')
        self.b_zoom.on_clicked(lambda e: self.basculer_zoom())
        self.b_reset = Button(self.fig.add_axes([0.865, 0.19, 0.115, 0.05]), 'Ajustement auto')
        self.b_reset.on_clicked(lambda e: self.reinit())
        self.b_save = Button(self.fig.add_axes([0.865, 0.10, 0.115, 0.05]), 'Enregistrer')
        self.b_save.on_clicked(lambda e: self.enregistrer())

        c = self.fig.canvas
        c.mpl_connect('button_press_event', self.press)
        c.mpl_connect('button_release_event', lambda e: setattr(self, 'drag', None))
        c.mpl_connect('motion_notify_event', self.move)
        c.mpl_connect('scroll_event', self.scroll)
        c.mpl_connect('key_press_event', self.key)
        c.mpl_connect('close_event', lambda e: self.enregistrer())
        self.dessiner_carte()
        self.sync()

    def dessiner_carte(self):
        S = self.S[self.lab]
        ox, oz = np.argsort(S['x_cad']), np.argsort(S['z_cad'])
        img = S['q'][np.ix_(ox, oz)]
        x, z = S['x_cad'][ox], S['z_cad'][oz]
        ext = [x[0] - S['ds'] / 2, x[-1] + S['ds'] / 2, z[0] - S['di'] / 2, z[-1] + S['di'] / 2]
        if self.im is not None:
            self.im.remove()
        self.im = self.ax.imshow(img.T, origin='lower', extent=ext, aspect='equal', cmap='viridis', vmin=0, vmax=1,
                                 interpolation='nearest', zorder=0)
        self.ax.set_xlabel('$x$ (mm)')
        self.ax.set_ylabel('$z$ (mm)')
        self.cadrer()

    def cadrer(self):
        if self.zoom:
            L, p = self.LC[self.k], self.P[self.lab][self.k]
            c = L['centre'] + p['bx']
            demi = 0.5 * L['largeur'] * p['ax'] + 6
            self.ax.set_xlim(c - demi, c + demi)
        else:
            self.ax.set_xlim(-73, 73)
        self.ax.set_ylim(13, -13)                    # vue de dessus : les lettres se lisent

    # ------------------------------------------------------------ etat
    def sync(self):
        """Contours, curseurs et texte d'etat a partir des parametres."""
        P = self.P[self.lab]
        for k, (g, lines) in enumerate(zip(self.grp, self.lignes)):
            L, p = self.LC[k], P[k]
            for loop, ln in zip(g, lines):
                ln.set_data(L['centre'] + p['bx'] + p['ax'] * (loop[:, 0] - L['centre']),
                            L['zc'] + p['bz'] + p['az'] * (loop[:, 1] - L['zc']))
                ln.set_linewidth(1.8 if k == self.k else 0.9)
                ln.set_alpha(0.45 if p['tronquee'] and k != self.k else 1.0)
        p = P[self.k]
        self._muet = True
        for s, v in ((self.s_bx, p['bx']), (self.s_ax, p['ax']), (self.s_bz, p['bz']), (self.s_az, p['az'])):
            s.set_val(float(np.clip(v, s.valmin, s.valmax)))
        if self.c_tr.get_status()[0] != p['tronquee']:
            self.c_tr.set_active(0)
        self._muet = False
        L = self.LC[self.k]
        self.info.set_text('%s, lettre %d (%s) : centre %.2f mm, longueur %.2f mm (CAD %.2f), hauteur %.2f mm (CAD %.2f)%s'
                           % (self.lab, self.k + 1, A.LETTRES[self.k], L['centre'] + p['bx'], p['ax'] * L['largeur'],
                              L['largeur'], p['az'] * L['hauteur'], L['hauteur'],
                              '   [coupée : exclue des mesures]' if p['tronquee'] else ''))
        if self.zoom:
            self.cadrer()
        self.fig.canvas.draw_idle()

    def slider(self, _):
        if getattr(self, '_muet', False):
            return
        p = self.P[self.lab][self.k]
        p.update(bx=self.s_bx.val, ax=self.s_ax.val, bz=self.s_bz.val, az=self.s_az.val)
        self.sync()

    def coupee(self, _):
        if getattr(self, '_muet', False):
            return
        self.P[self.lab][self.k]['tronquee'] = self.c_tr.get_status()[0]
        self.sync()

    def choisir_scan(self, lab):
        self.lab = lab
        self.dessiner_carte()
        self.sync()

    def choisir_lettre(self, k):
        self.k = k % len(self.LC)
        if self.r_let.value_selected != self.r_let.labels[self.k].get_text():
            self.r_let.set_active(self.k)
        self.sync()

    def basculer_zoom(self):
        self.zoom = not self.zoom
        self.b_zoom.label.set_text('Vue complète' if self.zoom else 'Zoom lettre')
        self.cadrer()
        self.fig.canvas.draw_idle()

    def reinit(self):
        self.P[self.lab][self.k] = dict(self.P0[self.lab][self.k])
        self.sync()

    def modifier(self, **d):
        p = self.P[self.lab][self.k]
        for key, v in d.items():
            p[key] += v
        p['ax'] = float(np.clip(p['ax'], 0.2, 2.5))
        p['az'] = float(np.clip(p['az'], 0.5, 1.5))
        self.sync()

    # ------------------------------------------------------------ souris et clavier
    def lettre_sous(self, x):
        P = self.P[self.lab]
        for k, L in enumerate(self.LC):
            c, demi = L['centre'] + P[k]['bx'], 0.5 * L['largeur'] * P[k]['ax']
            if abs(x - c) <= demi:
                return k
        return None

    def press(self, e):
        if e.inaxes is not self.ax or e.button != 1:
            return
        k = self.lettre_sous(e.xdata)
        if k is not None and k != self.k:
            self.choisir_lettre(k)
        p = self.P[self.lab][self.k]
        self.drag = (e.xdata, e.ydata, p['bx'], p['bz'])

    def move(self, e):
        if not self.drag or e.inaxes is not self.ax or e.xdata is None:
            return
        x0, z0, bx0, bz0 = self.drag
        p = self.P[self.lab][self.k]
        p['bx'], p['bz'] = bx0 + (e.xdata - x0), bz0 + (e.ydata - z0)
        self.sync()

    def scroll(self, e):
        f = 0.01 * (1 if e.button == 'up' else -1)
        if e.key == 'shift':
            self.modifier(az=f)
        else:
            self.modifier(ax=f)

    def key(self, e):
        k = e.key or ''
        fin = k.startswith('shift+') or k.isupper()
        b = k.replace('shift+', '').lower()
        pas, et = (1.0, 0.002) if fin else (0.1, 0.01)
        actions = {'left': dict(bx=-pas), 'right': dict(bx=pas), 'up': dict(bz=-pas), 'down': dict(bz=pas),
                   'a': dict(ax=-et), 'd': dict(ax=et), 's': dict(az=-et), 'w': dict(az=et)}
        if k == 'ctrl+s':
            self.enregistrer()
        elif b in actions:
            self.modifier(**actions[b])
        elif b == 'n':
            self.choisir_lettre(self.k + 1)
        elif b == 'p':
            self.choisir_lettre(self.k - 1)
        elif b == 'z':
            self.basculer_zoom()

    # ------------------------------------------------------------ sortie
    def enregistrer(self):
        out = dict(description='Placement manuel du contour CAD de chaque lettre. x_mesure = cx + bx + ax * (x_cad - cx), '
                               'z_mesure = cz + bz + az * (z_cad - cz), en mm, dans les coordonnees du recalage global.',
                   date=datetime.datetime.now().isoformat(timespec='seconds'), scans={})
        for nom, lab in A.SCANS:
            S = self.S[lab]
            out['scans'][lab] = dict(fichier=nom, recalage=reg(S), lettres=[
                dict(lettre=A.LETTRES[k], cx=L['centre'], cz=L['zc'], largeur_cad=L['largeur'], hauteur_cad=L['hauteur'],
                     **{key: float(self.P[lab][k][key]) for key in ('ax', 'bx', 'az', 'bz')},
                     tronquee=bool(self.P[lab][k]['tronquee']))
                for k, L in enumerate(self.LC)])
        with open(SORTIE, 'w', encoding='utf-8') as f:
            json.dump(out, f, indent=2, ensure_ascii=False)
        print('Enregistré :', SORTIE)
        self.fig.suptitle('Enregistré à %s' % datetime.datetime.now().strftime('%H:%M:%S'), fontsize=9, x=0.9, y=0.995,
                          color='0.4')
        self.fig.canvas.draw_idle()


if __name__ == '__main__':
    Outil()
    plt.show()
