# Instructions pour Claude

## Contexte du dépôt

Rapports GPH-3000 (groupe 85213, A26). Un dossier par expérience (acoustique,
cnd, geophysique, laser-co2, nucleaire, rayons-x), structuré ainsi :

- `data/` — données brutes
- `figures/` — figures générées (scripts Python)
- `protocole/` — protocole spécifique à l'expérience
- `Rapport/` — le rapport LaTeX (`.tex` + `.pdf` compilé), avec un
  sous-dossier `Rapport/scripts/` pour le code d'analyse

Compilation LaTeX via MikTeX.

## Rédaction des rapports LaTeX

Avant de rédiger ou de générer le contenu d'un rapport (`<labo>/Rapport/*.tex`),
lire dans l'ordre :

1. `consignes/Consignes_Rapports_1.pdf` et `consignes/Grille évaluation rapports.pdf`
   — structure attendue, sections obligatoires, grille de notation. Ce sont
   les consignes standard du cours.
2. Le protocole spécifique du labo concerné (`<labo>/protocole/`).

**Ajustements A26, confirmés en session par le professeur (Simon Rainville),
qui priment sur les consignes standard ci-dessus :**
- Longueur cible : **5-7 pages** (les consignes standard visent ~10 pages)
- Remise : vendredi 23h55 de la semaine suivant l'expérience, dépôt PDF sur
  Gradescope (+1 semaine si l'échéance tombe pendant la semaine de lecture)

## Conventions de rédaction

- Langue : français, registre scientifique standard (pas de québécois oral)
- Nombres en mode math : virgule décimale (`icomma` ou `siunitx` avec
  `locale=FR`), jamais de notation calculatrice type `5E-3` — utiliser
  `5\times10^{-3}`
- Unités : jamais en italique — `\text{}`, `\mathrm{}` ou `\textrm{}`, avec
  une espace insécable entre le nombre et l'unité
- Parenthèses/crochets de grande taille : `\left(`/`\right)`,
  `\left[`/`\right]`
- Toute incertitude doit être justifiée et propagée par dérivées partielles
- Éviter le ton trop lisse/générique dans les sections rédigées : varier la
  structure des paragraphes, éviter les tournures répétitives d'un rapport
  à l'autre

