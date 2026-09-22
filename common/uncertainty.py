"""
Propagation d'incertitude.

Deux approches sont offertes :

1. `propagate_numeric` : propagation par dérivées partielles calculées
   numériquement (différences finies centrées). Utile pour une fonction
   Python quelconque (pas besoin qu'elle soit "traçable" symboliquement).

2. Les fonctions `ufloat` / `to_ufloat` : wrapper mince autour du package
   `uncertainties`, pratique quand les variables d'entrée sont déjà
   indépendantes et qu'on veut laisser `uncertainties` gérer la propagation
   (incluant les corrélations créées par les opérations elles-mêmes).

En général : utiliser `uncertainties` pour des expressions algébriques
simples, et `propagate_numeric` quand le calcul est trop complexe pour être
réécrit avec des `ufloat` (ex. appel à `scipy`, boucles, etc.).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
from uncertainties import ufloat, unumpy
from uncertainties.core import AffineScalarFunc


def to_ufloat(valeur: float, incertitude: float):
    """Crée un nombre avec incertitude (ufloat) à partir d'une valeur et
    de son incertitude-type."""
    return ufloat(valeur, incertitude)


def valeurs_et_incertitudes(ufloats: Sequence[AffineScalarFunc]):
    """Sépare une séquence de `ufloat` en deux tableaux numpy :
    (valeurs nominales, incertitudes-types)."""
    valeurs = unumpy.nominal_values(ufloats)
    incertitudes = unumpy.std_devs(ufloats)
    return valeurs, incertitudes


def propagate_numeric(
    func: Callable[..., float],
    valeurs: Sequence[float],
    incertitudes: Sequence[float],
    correlations: np.ndarray | None = None,
    h_relatif: float = 1e-6,
) -> tuple[float, float]:
    """Propage une incertitude par dérivées partielles numériques.

    Calcule f(valeurs) et son incertitude-type via la formule générale
    de propagation :

        sigma_f^2 = sum_i (df/dx_i)^2 * sigma_i^2
                    + 2 * sum_{i<j} (df/dx_i)(df/dx_j) * cov(x_i, x_j)

    Les dérivées partielles sont estimées par différences finies centrées.

    Parameters
    ----------
    func : Callable
        Fonction scalaire prenant les variables comme arguments positionnels,
        ex. `func(x, y, z)`.
    valeurs : Sequence[float]
        Valeurs nominales des variables, dans l'ordre attendu par `func`.
    incertitudes : Sequence[float]
        Incertitudes-types associées à chaque variable.
    correlations : np.ndarray, optional
        Matrice de covariance complète (n x n). Si fournie, remplace le
        produit `incertitudes[i] * incertitudes[j]` pour les termes croisés.
        Si None, les variables sont supposées indépendantes.
    h_relatif : float
        Pas relatif utilisé pour la différence finie centrée.

    Returns
    -------
    (valeur, incertitude) : tuple[float, float]
    """
    valeurs = np.asarray(valeurs, dtype=float)
    incertitudes = np.asarray(incertitudes, dtype=float)
    n = len(valeurs)

    if correlations is None:
        cov = np.diag(incertitudes**2)
    else:
        cov = np.asarray(correlations, dtype=float)

    valeur_centrale = func(*valeurs)

    # Dérivées partielles par différences finies centrées.
    derivees = np.zeros(n)
    for i in range(n):
        h = h_relatif * valeurs[i] if valeurs[i] != 0 else h_relatif
        valeurs_plus = valeurs.copy()
        valeurs_moins = valeurs.copy()
        valeurs_plus[i] += h
        valeurs_moins[i] -= h
        derivees[i] = (func(*valeurs_plus) - func(*valeurs_moins)) / (2 * h)

    variance = derivees @ cov @ derivees
    return float(valeur_centrale), float(np.sqrt(variance))
