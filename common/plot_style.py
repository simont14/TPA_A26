"""
Style matplotlib personnalisé pour des graphiques de qualité scientifique.

Usage :

    from common.plot_style import set_style, error_plot

    set_style()  # à appeler une fois, en haut du script/notebook

    fig, ax = plt.subplots()
    error_plot(ax, x, y, yerr=sigma_y, xerr=sigma_x,
               xlabel="Temps", xunit="s", ylabel="Tension", yunit="V")
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

# Style de base : police lisible, grille discrète, traits nets.
STYLE = {
    "figure.figsize": (7, 5),
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.size": 12,
    "font.family": "serif",
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
    "legend.frameon": False,
    "legend.fontsize": 11,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "lines.linewidth": 1.5,
    "lines.markersize": 5,
    "errorbar.capsize": 3,
    "axes.spines.top": False,
    "axes.spines.right": False,
}


def set_style() -> None:
    """Applique le style matplotlib personnalisé à toute la session."""
    plt.rcParams.update(STYLE)


def _label_avec_unite(nom: str | None, unite: str | None) -> str | None:
    if nom is None:
        return None
    return f"{nom} [{unite}]" if unite else nom


def error_plot(
    ax,
    x,
    y,
    xerr=None,
    yerr=None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    xunit: str | None = None,
    yunit: str | None = None,
    label: str | None = None,
    fmt: str = "o",
    **kwargs,
):
    """Trace des points expérimentaux avec barres d'erreur et labels
    contenant les unités.

    Wrapper autour de `ax.errorbar` qui uniformise l'habillage des
    graphiques (labels avec unités, barres d'erreur, style de points).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt=fmt, label=label, **kwargs)

    xl = _label_avec_unite(xlabel, xunit)
    yl = _label_avec_unite(ylabel, yunit)
    if xl:
        ax.set_xlabel(xl)
    if yl:
        ax.set_ylabel(yl)
    if label:
        ax.legend()

    return ax
