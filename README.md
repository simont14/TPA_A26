# TPA_A26 — GPH-3000, Travaux pratiques avancés (automne 2026)

Dépôt centralisant le code Python utilisé pour analyser les données de chaque
labo : calculs, propagation d'incertitude, graphiques et ajustements de
courbes (fits).

## Structure

```
TPA_A26/
├── consignes/            # Consignes générales de rédaction des rapports
│                         # (valables pour tous les labos)
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

`consignes/` contient les documents généraux du cours qui s'appliquent à tous
les rapports (consignes de rédaction, grille d'évaluation, etc.), par
opposition à `protocole/` dans chaque labo qui contient le protocole
spécifique à ce labo.

Chaque dossier de labo suit la même structure :

- `data/` : données brutes du labo (non modifiées après acquisition).
- `scripts/` : scripts Python d'analyse (un script par manipulation ou par
  section du labo, au besoin).
- `figures/` : graphiques générés par les scripts (idéalement reproductibles
  en relançant le script correspondant).
- `protocole/` : protocole(s) du labo fourni(s) par le cours (PDF, etc.).
- `rapport/` : rapport LaTeX (ou autre) du labo, incluant les fichiers générés
  à la compilation. Les figures utilisées dans le rapport sont référencées
  directement depuis `figures/` (chemin relatif, ex. `../figures/nom.png`),
  sans être dupliquées dans `rapport/`.

## Installation

```bash
python -m venv venv
source venv/bin/activate        # Windows : venv\Scripts\activate
pip install -r requirements.txt
```

## Ajouter un nouveau labo

1. Créer un dossier à la racine (ex. `mon-labo/`) avec les sous-dossiers
   `data/`, `scripts/`, `figures/`, `protocole/`, `rapport/`.
2. Ajouter un `README.md` minimal dans ce dossier décrivant le labo.
3. Si une nouvelle dépendance Python est nécessaire, l'ajouter à
   `requirements.txt`.
