"""
Fonctions de fit génériques basées sur `scipy.optimize.curve_fit`, avec
extraction propre des incertitudes sur les paramètres et quelques
statistiques utiles (chi carré réduit).

Usage :

    from common.fitting import fit_curve

    def lineaire(x, m, b):
        return m * x + b

    resultat = fit_curve(lineaire, x, y, sigma=yerr, p0=[1, 0])
    print(resultat.parametres)      # {'m': ufloat(...), 'b': ufloat(...)}
    print(resultat.chi2_reduit)
"""

from __future__ import annotations

import inspect
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import curve_fit
from uncertainties import correlated_values


@dataclass
class ResultatFit:
    """Résultat d'un ajustement de courbe."""

    fonction: Callable
    noms_parametres: list[str]
    valeurs: np.ndarray
    covariance: np.ndarray
    parametres: dict = field(init=False)
    chi2: float = field(init=False)
    chi2_reduit: float = field(init=False)
    ndl: int = field(init=False)

    def __post_init__(self):
        # Paramètres corrélés avec incertitude (objets `ufloat`), pratiques
        # pour des calculs subséquents qui propagent correctement les
        # covariances entre paramètres.
        params_correles = correlated_values(self.valeurs, self.covariance)
        self.parametres = dict(zip(self.noms_parametres, params_correles))

    def evaluer(self, x):
        """Évalue la fonction ajustée (valeurs nominales) sur `x`."""
        return self.fonction(x, *self.valeurs)


def fit_curve(
    func: Callable,
    xdata,
    ydata,
    sigma=None,
    p0=None,
    absolute_sigma: bool = True,
    **kwargs,
) -> ResultatFit:
    """Ajuste `func` aux données par moindres carrés non linéaires.

    Parameters
    ----------
    func : Callable
        Fonction du modèle, signature `func(x, param1, param2, ...)`.
    xdata, ydata : array-like
        Données à ajuster.
    sigma : array-like, optional
        Incertitudes-types sur `ydata`. Fortement recommandé pour obtenir
        un chi carré et des incertitudes de paramètres significatifs.
    p0 : array-like, optional
        Estimé initial des paramètres.
    absolute_sigma : bool
        Si True (défaut), `sigma` est interprété comme une incertitude
        absolue réelle (recommandé en physique expérimentale). Si False,
        seule l'incertitude relative entre points compte et la covariance
        est mise à l'échelle par le chi2 réduit.

    Returns
    -------
    ResultatFit
        Contient les paramètres ajustés (avec incertitudes, sous forme de
        `ufloat` corrélés), la matrice de covariance et le chi carré réduit.
    """
    xdata = np.asarray(xdata, dtype=float)
    ydata = np.asarray(ydata, dtype=float)

    noms_parametres = list(inspect.signature(func).parameters.keys())[1:]

    valeurs, covariance = curve_fit(
        func, xdata, ydata, p0=p0, sigma=sigma, absolute_sigma=absolute_sigma, **kwargs
    )

    resultat = ResultatFit(
        fonction=func,
        noms_parametres=noms_parametres,
        valeurs=valeurs,
        covariance=covariance,
    )

    residus = ydata - func(xdata, *valeurs)
    ndl = len(xdata) - len(valeurs)
    if sigma is not None:
        chi2 = float(np.sum((residus / np.asarray(sigma)) ** 2))
    else:
        chi2 = float(np.sum(residus**2))

    resultat.chi2 = chi2
    resultat.ndl = ndl
    resultat.chi2_reduit = chi2 / ndl if ndl > 0 else float("nan")

    return resultat


# Modèles courants réutilisables tels quels avec `fit_curve`.


def lineaire(x, m, b):
    return m * x + b


def exponentielle(x, a, tau, c):
    return a * np.exp(-x / tau) + c


def gaussienne(x, amplitude, mu, sigma, offset=0.0):
    return amplitude * np.exp(-((x - mu) ** 2) / (2 * sigma**2)) + offset


def sinusoide(x, amplitude, frequence, phase, offset=0.0):
    return amplitude * np.sin(2 * np.pi * frequence * x + phase) + offset
