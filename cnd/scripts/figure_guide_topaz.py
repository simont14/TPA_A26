#!/usr/bin/env python3
"""
figure_guide_topaz.py : image de la page couverture du guide Topaz (protocole_topaz).

Surface 3D du bloc LASERAX reconstruite a partir du balayage avec encodeur, sans axes. Reprend le calcul et
le rendu de analyse_scans.py (figure 04, panneau b), sans modifier les figures du rapport.

Lancement :  python figure_guide_topaz.py      (ecrit ../figures/guide_topaz_scan3d.png)
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyse_scans as A

SORTIE = os.path.join(A.FIG, 'guide_topaz_scan3d.png')


def rendu(X, Z, H, path, elev=40, azim=-75):
    fig = plt.figure(figsize=(10.0, 4.6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Z, np.clip(H, 9, 15), cmap='viridis', vmin=9, vmax=15.5, rstride=1, cstride=1,
                    linewidth=0, antialiased=False, shade=True)
    ax.set_box_aspect((A.X3D[1] - A.X3D[0], A.Z3D[1] - A.Z3D[0], 3.0 * 6))
    ax.set_xlim(*A.X3D)
    ax.set_ylim(A.Z3D[1], A.Z3D[0])
    ax.set_zlim(9, 15)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    fig.savefig(path, dpi=300, facecolor='white')
    plt.close(fig)
    im = A._rogner([Image.open(path).convert('RGB')], marge=10)[0]
    im.save(path)


def main():
    C = A.cad()
    S = A.recaler(A.charger('Encodeur_Bloc_Laser_6.UVData'), C)
    rendu(*A.surface_mesuree(S), SORTIE)
    print('ecrit', SORTIE)


if __name__ == '__main__':
    main()
