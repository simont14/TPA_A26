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

## Éviter un style trop visiblement généré par IA

- Jamais de point-virgule (`;`) dans le texte courant : scinder en deux
  phrases, ou utiliser une virgule.
- Jamais de tiret cadratin/demi-cadratin (`—`, `--`, `---`) comme ponctuation
  de pause : utiliser une virgule, des parenthèses, ou une nouvelle phrase.
  Un trait d'union simple pour un mot composé (`canal-énergie`) reste correct.
- Jamais d'espace avant `:`, toujours un espace après (`titre: texte`). C'est
  l'inverse de la convention typographique française standard, donc à forcer
  explicitement : avec `babel[french]`, ajouter `\shorthandoff{;:!?}`
  **après** `\begin{document}` (pas dans le préambule, babel réactive les
  raccourcis à `\begin{document}`).
- Jamais de mot en gras au milieu d'un paragraphe ou d'un item de liste
  (pas de `\textbf{}` en pseudo-titre). Utiliser `\emph{}` pour l'emphase, et
  de vrais sous-titres (`\subsubsection*{}`, `\paragraph{}`) pour structurer.
- Pas de vocabulaire dramatique ou marketing (« un véritable outil »,
  « transformer X en Y », « à part entière », « porteur de sens/d'information »)
  et pas de motivation/narratif inventé plus grandiose que la réalité (ex. ne
  pas écrire que « le but de l'expérience était de transformer un cristal en
  véritable outil » quand l'objectif réel est simplement d'étalonner un
  détecteur). Rester factuel et sobre, comme un·e étudiant·e écrit vraiment
  un rapport.
- Éviter la tournure contrastive « ce n'était pas seulement X, mais Y » : un
  tic d'écriture IA reconnaissable. Énoncer le fait directement.

### Écrire comme un·e étudiant·e l'expliquerait, pas comme un manuel

Le `\subsection*{Résumé}` de `nucleaire/rapport/template_devoirBD.tex` (après
sa réécriture du 23 septembre 2026) est l'exemple de référence à suivre. Ce
qui en ressort :

- Ne pas empiler les détails techniques précis dans une phrase de synthèse.
  Le niveau de détail du genre « onze raies de 30,85 keV à 1332,50 keV » va
  dans les résultats, pas dans le résumé ou l'introduction : une phrase de
  synthèse dit *quoi* et *pourquoi*, pas *toutes les valeurs*.
- Préférer une phrase explicite et déroulée à une formulation compressée ou
  nominalisée. Dire « Le modèle utilisé pour faire la correspondance entre
  les canaux et l'énergie a été un modèle quadratique » plutôt que « Un
  modèle quadratique canal-énergie a été retenu ».
- Préférer un mot courant à un terme technique précis quand le mot courant
  reste correct. « Les mathématiques de la diffusion Compton » plutôt que
  « la cinématique de diffusion Compton » ; « qui correspond à » plutôt que
  « compatible avec ».
- Préférer une expression parlée à une expression soutenue pour exprimer un
  accord approximatif. « Plus ou moins la résolution du détecteur » plutôt
  que « à la résolution du détecteur près ».
- Ne pas ajouter de fioriture « poétique » même factuellement vraie (ex.
  éviter un clin d'œil du genre « ... et même de nous-mêmes » ; rester dans
  l'énumération factuelle plate : « le béton, la brique et le sol »).
- Ne jamais remplacer une anecdote concrète et vérifiable par une liste
  générique d'applications (« utilisé en sécurité, en médecine, en
  environnement... »). Chercher un exemple réel, spécifique, si possible
  historique ou surprenant, et le citer. Ex. : la centrale de Forsmark qui a
  détecté la catastrophe de Tchernobyl deux jours après qu'elle ait eu lieu,
  grâce au même principe de spectroscopie combiné à une analyse du vent.
- Ne pas présenter un résultat de façon plus définitive ou plus certaine que
  ce que le texte a réellement démontré à ce point du rapport.
- Corriger l'orthographe et la grammaire de ce qui est écrit, mais toujours
  en gardant le sens et le registre voulus : le but est un français correct
  qui sonne quand même comme écrit par une personne, pas par un logiciel.

