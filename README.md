# TPA_A26 — GPH-3000, Travaux pratiques avancés (automne 2026)

Dépôt centralisant le code Python utilisé pour analyser les données de chaque
labo : calculs, propagation d'incertitude, graphiques et ajustements de
courbes (fits).

## Structure

```
TPA_A26/
├── common/              # Code partagé entre tous les labos
│   ├── uncertainty.py    # Propagation d'incertitude (numérique + wrapper `uncertainties`)
│   ├── plot_style.py     # Style matplotlib personnalisé + helper de graphique avec barres d'erreur
│   └── fitting.py        # Wrapper autour de scipy.optimize.curve_fit (incertitudes, chi2)
│
├── nucleaire/            # Radiations nucléaires (N1-N2)
├── cnd/                  # CND ultrasons et courants induits (M1-M2)
├── geophysique/          # Géophysique
├── acoustique/            # Acoustique (A1-A2)
├── rayons-x/              # Rayons-X (RX1-RX2)
├── laser-co2/             # Laser CO2
│
├── requirements.txt
└── .gitignore
```

Chaque dossier de labo suit la même structure :

- `data/` : données brutes du labo (non modifiées après acquisition).
- `scripts/` : scripts Python d'analyse (un script par manipulation ou par
  section du labo, au besoin).
- `figures/` : graphiques générés par les scripts (idéalement reproductibles
  en relançant le script correspondant).

## Installation

```bash
python -m venv venv
source venv/bin/activate        # Windows : venv\Scripts\activate
pip install -r requirements.txt
```

## Utilisation du code partagé (`common/`)

Un script dans un dossier de labo (ex. `nucleaire/scripts/analyse_n1.py`)
importe `common` en ajoutant la racine du dépôt au path, ou en lançant les
scripts depuis la racine avec `python -m nucleaire.scripts.analyse_n1` :

```python
from common.plot_style import set_style, error_plot
from common.fitting import fit_curve, lineaire
from common.uncertainty import to_ufloat, propagate_numeric

set_style()

resultat = fit_curve(lineaire, x, y, sigma=yerr, p0=[1, 0])
print(resultat.parametres)       # {'m': ufloat(...), 'b': ufloat(...)}
print(resultat.chi2_reduit)
```

Voir les docstrings dans `common/uncertainty.py`, `common/plot_style.py` et
`common/fitting.py` pour le détail de chaque fonction.

## Ajouter un nouveau labo

1. Créer un dossier à la racine (ex. `mon-labo/`) avec les sous-dossiers
   `data/`, `scripts/`, `figures/`.
2. Ajouter un `README.md` minimal dans ce dossier décrivant le labo.
3. Réutiliser les fonctions de `common/` plutôt que de dupliquer du code de
   propagation d'incertitude, de style de graphique ou de fit.
4. Si une nouvelle dépendance Python est nécessaire, l'ajouter à
   `requirements.txt`.
