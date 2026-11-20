# Instructions pour Claude

## Contexte du dépôt

Rapports GPH-3000 (groupe 85213, A26). Un dossier par expérience (acoustique,
cnd, geophysique, laser-co2, nucleaire, rayons-x), structuré ainsi :

- `data/` — données brutes
- `figures/` — figures générées (scripts Python)
- `protocole/` — protocole spécifique à l'expérience
- `preparation/` — fiche de préparation du labo (`.docx`), voir plus bas
- `cahier_de_lab/` — notes et manipulations faites pendant la session de labo
- `rapport/` — le rapport LaTeX (`.tex` + `.pdf` compilé)
- `scripts/` — code d'analyse Python

Les labs déjà complétés (nucleaire, cnd, geophysique) n'ont ni
`preparation/` ni `cahier_de_lab/`. Ne pas les y ajouter.

Compilation LaTeX via MikTeX.

## Préparation d'un laboratoire

Quand l'auteur demande de « faire la préparation » d'un labo, il veut une
fiche à coller dans OneNote (page « Préparation » de la section du labo dans
le cahier TPA_Cahier_de_Lab). Référence approuvée:
`acoustique/preparation/Preparation_Acoustique.docx` (3 novembre 2026).

- Sources: tout ce qui se trouve dans `<labo>/protocole/` (protocole(s),
  page du site de cours, articles).
- Livrable: un seul fichier `<labo>/preparation/Preparation_<Labo>.docx`.
  Pas de HTML ni de texte brut à côté: le collage navigateur → OneNote rend
  mal, seul le collage depuis Word garde les puces.
- Génération: écrire le contenu en HTML simple (titres de section en `<p>`,
  listes `<ul>` imbriquées, Calibri 11, aucune couleur) dans le scratchpad,
  puis le convertir avec Word par COM (`Documents.Open(html)` puis
  `SaveAs2(docx, 16)`). Ne pas faire de script qui copie dans le
  presse-papier avec un Word caché (il reste bloqué à la fermeture).
- Sections, dans l'ordre, chacune en liste à puces:
  1. Explications: la physique et le but de chaque partie du labo
  2. Formules: en texte simple (`Lp = 20·log₁₀(p/p_réf)`), chaque paramètre
     en sous-puce avec sa signification et son unité entre crochets
  3. Équipement: avec les réglages des logiciels s'ils sont donnés
  4. Manipulations: étapes par partie, en sous-puces
  5. Définitions importantes: termes moins communs
  6. Points importants: à ne pas oublier pendant le labo
  7. Questions: celles du protocole, chacune avec une piste de réponse
  8. Pistes d'exemples intéressants: applications concrètes et historiques
     de la matière (pas une liste générique), incluant les articles fournis
- Les règles de style des rapports (pas de `;`, etc.) ne s'appliquent pas
  strictement à cette fiche, qui est un document de travail personnel.

## Rédaction des rapports LaTeX

Avant de rédiger ou de générer le contenu d'un rapport (`<labo>/rapport/*.tex`),
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
- Toute incertitude doit être justifiée et propagée par dérivées partielles,
  mais rester simple : une seule équation de propagation pour les résultats
  principaux suffit. Pas d'incertitude sur des numéros de canal (ce sont des
  entiers), pas de propagation pour chaque paramètre secondaire (ex. le
  préfacteur d'une loi de puissance), pas de sous-section « Incertitudes »
  générique dans la discussion.
- Éviter le ton trop lisse/générique dans les sections rédigées : varier la
  structure des paragraphes, éviter les tournures répétitives d'un rapport
  à l'autre
- Bibliographie : une référence qui n'est citée nulle part dans le texte ne
  doit pas apparaître. Utiliser le mécanisme de
  `geophysique/rapport/TPA_Lab3.tex` (redéfinition de `\citation` dans le
  préambule et `\refsicitee{cle}{texte}` à la place de `\bibitem{cle} texte`),
  qui masque automatiquement les entrées non citées. Dans le texte d'une
  entrée, écrire `\%` dans les URL (le texte est un argument de macro).

### Principes tirés de la révision du rapport nucléaire (24 septembre 2026)

Ces principes viennent de la révision ligne par ligne faite par l'auteur sur
`nucleaire/rapport/TPA_Lab1.tex`. Ils s'appliquent à tous les rapports.

- Rapporter ce qui a été utilisé, pas ce qui a été testé. Pas de modèle
  rejeté (ex. linéaire contre quadratique), pas de critère de sélection de
  modèle, pas de figure ou de phrase qui compare au modèle abandonné. Pour
  justifier le modèle retenu, une phrase suffit (« Le modèle quadratique est
  un bon choix, avec $R^2 = 0{,}99997$... »), sans argumentaire défensif du
  type « une amélioration suffisamment nette pour être retenue plutôt que ».
- Pas de résidus : ni colonne de résidus dans les tableaux, ni panneau de
  résidus sous les figures, sauf si le protocole le demande.
- Pas de section « Montage » qui recopie la chaîne d'acquisition du
  protocole (le lecteur a le protocole). Citer le protocole une fois dans la
  méthodologie.
- Méthodologie courte : décrire l'idée (« la région autour du canal attendu
  est analysée pour trouver un maximum local ») et non l'implémentation
  (pas de NumPy/SciPy, pas de nom de fonction ni de script, pas de formule de
  la gaussienne d'ajustement). « Curve fit » suffit pour nommer l'ajustement.
- Pas de détails entre parenthèses qui n'apportent rien (« canal par canal,
  même détecteur et même gain », « échelle logarithmique, 36 minutes
  d'acquisition » dans une légende).
- Pas de guillemets d'emphase ou d'ironie (« volontaire »).
- Couper ce qui est hors sujet ou négligeable, même si c'est vrai (ex. un
  paragraphe sur la production de paires qui « reste négligeable », un
  paragraphe sur ce qu'un détecteur HPGe aurait permis, un rappel de
  « l'indice donné dans le protocole »).
- Ne pas énumérer en texte des valeurs déjà dans un tableau (« d'environ 26
  à 31 % pour les raies sous 35 keV, jusqu'à 4,5 à 5,2 %... »). Dire la
  tendance (« diminue avec l'énergie ») et renvoyer au tableau.
- Pour chaque pic identifié dans un spectre, donner la valeur théorique ET
  la valeur mesurée, dans le texte et sur la figure.
- Légendes de figures courtes et descriptives : ce que la figure montre, pas
  l'interprétation ni la méthode (« Énergie en fonction du numéro de canal
  pour les onze raies. »). Nommer la couleur réelle (« la zone orange ») plutôt
  qu'un terme vague (« la zone ombragée »).
- Pas de phrases méta sur l'annexe (« fournies à titre complémentaire », « le
  corps du rapport est compréhensible sans s'y référer »).
- Conclusion et perspectives courtes et concrètes : ne pas redonner les
  chiffres de détail, une amélioration par phrase (« Aussi, une source
  d'étalonnage au-dessus de 1332,5 keV permettrait de réaliser un meilleur
  étalonnage. »).
- Vocabulaire simple, exemples de remplacements faits par l'auteur :
  « moindre » → « plus petite », « abruptement » → « rapidement », « en deçà »
  → « en dessous », « lequel » → « qui », « acquises individuellement » →
  « mesurées », « particulièrement instructif » → « intéressant », « On y
  distingue nettement » → « On y voit », « centroïde » (dans le texte courant)
  → « centre », « à comparer à » → « similaire à », « a occasionné une
  fraction de temps mort » → « a totalisé un temps mort », « une prochaine
  itération de l'expérience gagnerait à consigner » → « une prochaine
  expérience pourrait prendre en compte », « un continuum qui décroît
  globalement » → « une tendance qui décroît ».

## Figures générées (scripts Python)

- Jamais de titre (`ax.set_title()` / `fig.suptitle()`) sur une figure
  destinée au rapport : dans un rapport scientifique, c'est la légende
  (`\caption{}`) qui joue ce rôle. Un titre dans l'image et une légende en
  dessous font doublon. Le `label=` passé à `ax.step()`/`ax.plot()` pour la
  légende interne (ex. nom de la source) peut rester, seul `set_title` est à
  proscrire.
- Toute figure doit rester lisible imprimée en noir et blanc. Ne jamais
  distinguer des séries seulement par la couleur: marqueurs différents
  (plein contre vide), styles de ligne, ou une case par série (grille de
  petites figures) quand il y a plus de deux ou trois séries. Pour une
  carte de couleur, prendre une palette à luminosité monotone (`magma`,
  `viridis`). Vérifier en convertissant l'image en gris.
- Avec `subcaption` (`\begin{subfigure}...\caption{}...\label{fig:xxx}`), le
  `\ref{fig:xxx}` d'une sous-figure retourne déjà le numéro complet avec sa
  lettre (ex. `2a`). Ne jamais ajouter la lettre à la main après le `\ref`
  (`Figure~\ref{fig:xxx}a` devient alors « Figure 2aa »). Écrire
  `Figure~\ref{fig:xxx}` tout court, ou référencer directement le label de
  l'autre sous-figure (`fig:yyy`) pour pointer vers `2b`.
- Quand deux sous-figures sont côte à côte (`subfigure` + `\hfill`), vérifier
  visuellement (via `pdftoppm` ou équivalent) que leurs largeurs relatives
  donnent des hauteurs comparables une fois affichées : une sous-figure à
  aspect ratio différent peut paraître beaucoup plus petite que l'autre même
  avec une largeur `\textwidth` proche. Ajuster les largeurs (ou régénérer la
  figure avec un aspect ratio différent) plutôt que de laisser un
  déséquilibre visuel.
- Avant de figer une figure avec des flèches d'annotation (`ax.annotate`)
  proches les unes des autres, régénérer et regarder l'image : deux flèches
  peuvent se croiser visuellement même si le code semble correct. Inverser
  la position (`xytext`) des deux étiquettes concernées règle généralement
  le problème.

### Schémas (TikZ ou autre)

Le schéma de la Figure 1b de `geophysique/rapport/TPA_Lab3.tex` (bloc
`tikzpicture` « couches », approuvé par l'auteur le 19 octobre 2026) est
l'exemple de référence à suivre pour tout nouveau schéma. Partir de ce code.
Ce qui en ressort :

- Fond en aplats de couleur douce (`brown!20`, `brown!45`) pour les milieux,
  traits noirs fins pour les interfaces, symbole de vitesse ($V_1$, $V_2$)
  dans le coin de chaque milieu.
- Trajets en couleur (bleu plein pour l'onde directe, rouge pointillé pour
  l'onde réfractée), avec une étiquette de la même couleur, en `\small`, juste
  à côté du trajet. Les trajets se terminent tous au même capteur.
- Objets simples : rectangle gris pour la plaque, triangles noirs pour les
  capteurs avec leur nom au-dessus, grosse flèche pour la masse avec son nom
  au-dessus.
- Cotes et angles selon les règles ci-dessous.

- Toutes les cotes (distances) d'une même figure utilisent le même style :
  une ligne avec une flèche aux deux bouts (`\draw[<->]`). Ne jamais mélanger
  accolades et flèches dans une figure, c'est incohérent pour le lecteur.
- Les cotes alignées sont à la même hauteur, au-dessus des objets, avec de
  courtes lignes de rappel. Aucun bout de cote ni ligne de rappel ne doit être
  caché par un autre élément (plaque, capteur, flèche, étiquette).
- L'étiquette d'une cote est centrée au-dessus de sa ligne et ne doit pas être
  plus large que l'intervalle coté. Le schéma n'étant pas à l'échelle, élargir
  l'intervalle plutôt que de décaler l'étiquette.
- L'étiquette d'un angle est centrée sur sa bissectrice, à l'extérieur de
  l'arc, sans toucher l'arc, la normale ni le trajet. Si le même angle apparaît
  plusieurs fois (ex. angle critique à l'entrée et à la sortie), le dessiner et
  l'étiqueter à chaque endroit, de la même façon.
- Toujours regarder le rendu (`pdftoppm`) après chaque modification d'un
  schéma, en zoomant sur les zones serrées.

## Ne jamais référencer le chemin du dépôt local dans le texte du rapport

Le rapport est lu par une personne qui n'a pas accès à ce dépôt Git (ni à ses
dossiers `data/`, `figures/`, `scripts/`). Ne jamais écrire dans le texte du
rapport une phrase qui renvoie le lecteur vers un dossier ou un chemin du
projet (ex. une « Annexe » qui dit « voir le script dans le dossier
`nucleaire/scripts/` du dépôt du cours »).

Ne jamais écrire non plus le nom d'un fichier de données (`Scan_Bloc_Laser_3`,
`Encodeur_Bloc_Laser_6.UVData`, etc.) ni le nom d'un script (`analyse_scans.py`,
`uvdata_3d_gui.py`, etc.) dans le rapport, que ce soit dans le texte, un
tableau ou une légende. Un vrai rapport scientifique nomme une acquisition par
ce qu'elle est (« horloge » contre « encodeur », « source de ¹³⁷Cs ») et décrit
l'analyse par son idée (« analysé en Python »), parce que ce sont les résultats
qui comptent, pas le nom des fichiers.

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
  (pas de `\textbf{}` en pseudo-titre). Utiliser de vrais sous-titres
  (`\subsubsection*{}`, `\paragraph{}`) pour structurer.
- Jamais de mot en italique (`\emph{}`, `\textit{}`) au milieu du texte pour
  mettre en valeur un terme (ex. « le *photopic* », « une *bosse de
  rétrodiffusion* », « cède *toute* son énergie ») : un·e étudiant·e n'écrit
  pas comme ça, le terme se suffit à lui-même. L'italique reste correct
  seulement là où c'est une convention : titres d'ouvrages dans la
  bibliographie, nom du cours dans l'en-tête.
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

Le `\subsection*{Résumé}` de `nucleaire/rapport/TPA_Lab1.tex` (après
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
  détecté la catastrophe de Tchernobyl deux jours après qu'elle a eu lieu,
  grâce au même principe de spectroscopie combiné à une analyse du vent.
- Ne pas présenter un résultat de façon plus définitive ou plus certaine que
  ce que le texte a réellement démontré à ce point du rapport.
- Corriger l'orthographe et la grammaire de ce qui est écrit, mais toujours
  en gardant le sens et le registre voulus : le but est un français correct
  qui sonne quand même comme écrit par une personne, pas par un logiciel.

### Autres principes (tirés de la réécriture de l'introduction, 23 septembre 2026)

- Ne pas commencer une section par une phrase de mise en contexte générique
  et vide de contenu (« X est important dans notre société », « on sait
  tous que... »). Aller directement au premier fait concret.
- Préférer deux verbes coordonnés par « et » à une construction en gérondif
  qui compresse deux actions en une seule. « Se désintègrent et émettent »
  plutôt que « se désintègrent en émettant ».
- Couper les propositions ou phrases qui ne font que reformuler ce qui vient
  d'être dit. Ex. : ne pas ajouter « tous les événements produisent le même
  signal » juste après avoir dit qu'il n'y a « aucune information sur
  l'énergie » — c'est la même idée redite deux fois.
- Couper les adjectifs et adverbes qui ne portent pas d'information réelle.
  « Une question concrète » → « une question » ; « cet instrument
  fraîchement calibré » → « ce système calibré » ; « cinq sources
  d'énergies connues » → « cinq sources connues ».
- Préférer une formulation directe à une définition entre parenthèses quand
  c'est possible. « Convertir les numéros de canal en une énergie en keV »
  plutôt que « convertir la sortie brute de l'électronique (un numéro de
  canal) en une grandeur physique utile (une énergie en keV) ».
- Énoncer un compte plutôt que de le qualifier abstraitement : « Il y avait
  deux objectifs » plutôt que « L'objectif était double », puis énumérer
  avec « Premièrement... Deuxièmement... » plutôt que « D'abord... Ensuite... ».
- Un paragraphe qui annonce le plan du rapport doit être une suite de
  phrases courtes et déclaratives, chacune avec sa propre référence de
  section, plutôt qu'une ou deux phrases longues qui empilent plusieurs
  sections avec des clauses complexes.
- **Toute référence à une section, figure, tableau ou équation doit être une
  vraie référence LaTeX (`\label{}`/`\ref{}`), jamais un numéro écrit en
  dur.** Ajouter un `\label{sec:xxx}` juste après chaque `\section{}` ou
  `\subsection{}` visée et écrire `Section~\ref{sec:xxx}` (pas
  `section 2` ni `Section~4.2` en dur). La référence reste juste si les
  sections changent de numéro, et devient cliquable dans le PDF
  (`hyperref` est déjà chargé dans le préambule).


### Style étudiant: synthèse (rapport acoustique, 19 novembre 2026)

Tirée d'une longue série de corrections de l'auteur sur
`acoustique/rapport/TPA_Lab4.tex`. À appliquer dès le premier jet de
chaque rapport, pour que l'auteur n'ait pas à repasser chaque phrase.

Le ton visé: un·e étudiant·e qui explique son labo à un·e autre étudiant·e.
Le but d'un rapport est de partager l'information clairement, pas d'avoir
l'air savant.

- Mots simples et banals, même si c'est un peu « simplet ». « Un micro est
  positionné », « On dit que x est... », « θ est l'angle entre... »,
  « On suppose ici... », « car la théorie dit que oui », « est assez pour ».
  Pas de mot savant quand un mot courant marche (« habitacle » → « dans la
  voiture », « trame » → « morceau », « reçue » → « entendue » dans le texte
  courant).
- Ne pas avoir peur de la répétition. Redire « le point P est le point
  de... », répéter « la vitesse visée » plutôt que d'alterner avec
  « affichée », « au compteur », etc. Un même objet garde le même nom
  partout (texte, légendes, tableau, figures).
- Nommer un concept une fois, clairement, puis réutiliser ce nom: « Ce
  profil est ce qu'on appelle une courbe Doppler », « On l'appellera le
  périastre pour le reste du rapport », « les queues ». Ensuite, utiliser
  ce nom plutôt qu'une périphrase (« la courbe de la fréquence en fonction
  du temps » → « la courbe Doppler »).
- Expliquer le pourquoi en mots simples, comme à l'oral: « Pour séparer v
  et d, il faut... », « Pour vérifier ça, on peut... », « ça veut dire
  que... ». Un exemple chiffré simple vaut mieux qu'une phrase abstraite
  (« une erreur de 5 % sur v donne environ 10 % sur d »).
- Être exact sur la physique avec des mots simples: « la fréquence de son
  bruit baisse » et pas « son bruit baisse ». « Des plateaux constants en
  fréquence » plutôt que « des plateaux » tout court.
- Ne pas chercher la phrase parfaite. Une phrase un peu lourde, une
  répétition ou une construction parlée sont correctes. Corriger
  l'orthographe et les accords, mais ne jamais transformer une tournure
  simple de l'auteur en tournure plus élégante. Ne pas ajouter d'erreurs
  exprès non plus.
- Plusieurs phrases très courtes à la suite font haché: les joindre avec
  « et » ou une virgule (« Les résultats sont présentés à la Section 4,
  ils sont discutés à la Section 5 et la Section 6 conclut »).
- Anglicismes: corriger en le disant à l'auteur (« assumer » au sens de
  supposer → « supposer »).

Par section:

- Résumé: quoi et pourquoi, jamais comment (la méthode est plus loin). Un
  seul chiffre par série de mesures (écart moyen par distance), en
  commençant par la première série. Pas de sous-groupes, pas de cause.
  Ne jamais omettre en silence un cas qui ne suit pas la tendance.
- Introduction: dire les objectifs simplement, avec la raison (« déterminer
  si, en pratique, ... est assez pour ..., car la théorie dit que oui »).
- Théorie: poser les hypothèses et les symboles en phrases directes,
  définir chaque terme utilisé plus loin (queues, périastre).
- Légendes: minimales. « Schéma de l'expérience. », « Spectrogrammes des
  passages à 50 km/h, à d = 3,53 m et 8,97 m. », « Vitesse et distance pour
  les 10 passages. ». Nommer la grandeur par ce qu'elle est (« fréquence
  dominante du spectre ») et pas par l'étape d'analyse. Pas de
  gauche/droite/en haut quand chaque panneau a son étiquette. Pas de valeur
  lue sur la figure.
- Tableaux: seulement les colonnes qui servent au propos principal (la
  colonne f₀ a été retirée). Les en-têtes nomment la grandeur (« d
  réelle »), pas l'instrument.
- Discussion et conclusion: une idée par phrase, la conclusion répond aux
  objectifs de l'introduction.

Figures:

- Lisibles en noir et blanc (voir « Figures générées »).
- Étiquette des panneaux avec le symbole: « v = 30 km/h », « d = 3,53 m ».
- Pour une figure d'exemple (spectrogramme), choisir le cas le plus
  propre et centrer la courbe dans le cadre. Ajuster la fenêtre d'analyse
  et les limites d'axes au cas choisi, et montrer quelques essais à
  l'auteur plutôt que d'en imposer un.

Sources et normes:

- Vérifier qu'une norme ou un règlement cité s'applique vraiment au
  Québec/Canada. Sinon, ne pas l'utiliser comme argument.
- Ne jamais nommer le modèle de la voiture utilisée dans le rapport.

## Journal des corrections de l'auteur

Chaque correction que l'auteur demande (« Remplace X par Y ») est appliquée
telle quelle dans le `.tex`, avec l'orthographe corrigée, puis notée ici avec
ce qu'elle apprend sur son style. Relire ce journal avant de rédiger.

- Guide Topaz (12 novembre 2026): « plusieurs éléments piézoélectriques
  indépendants, alignés en réseau » → « plusieurs éléments piézoélectriques
  alignés en réseau ». La phrase suivante dit déjà « Chaque élément peut être
  excité séparément ». Leçon: ne pas annoncer dans un adjectif une idée que la
  phrase suivante explique. Une phrase de présentation décrit l'objet
  physique (combien, comment disposé), et la propriété est expliquée une seule
  fois, à l'endroit où elle sert.
- Guide Topaz (12 novembre 2026), demande générale: « réduire les termes
  trop complexes, inutiles, trop longs, trop détaillés, les informations
  triviales ». Première passe jugée encore « beaucoup trop formelle ».
  Exemple donné par l'auteur: « Déverrouiller le loquet et retirer le capot
  du connecteur. Aligner le connecteur de la sonde avec celui du panneau
  latéral, puis verrouiller le loquet » → « Aligner le connecteur et le
  verrouiller (Figure 5) ». Leçons (18 → 13 pages):
  - Ne jamais écrire une étape évidente que tout le monde fait sans y
    penser (retirer un capot pour brancher, appuyer sur Accept après chaque
    fenêtre, allumer avec le bouton marche/arrêt, remettre le capot).
  - Pas de mots techniques que personne ne dit (« loquet », « panneau
    latéral », « séquence de démarrage s'amorce »). Utiliser les mots de
    tous les jours: « le côté droit du Topaz », « allumer », « marcher »,
    « mettre X à Y », « bouts de ruban adhésif ».
  - Le moins de mots possible: un pronom plutôt que de renommer l'objet
    (« le verrouiller »), le chemin de menu en raccourci (« Dans Channel,
    Calculator, entrer... ») plutôt que « Choisir le menu Channel, puis
    appuyer sur le bouton Calculator ».
  - Regrouper en une phrase les étapes d'une liste qui ne sont que des clics
    à la suite. Une liste numérotée seulement quand l'ordre compte vraiment.
  - Une information dite deux fois est gardée à un seul endroit. Couper ce
    qui n'a pas servi (« mode non testé », « pas utilisé ici ») et les
    remarques de tableau qui redisent la valeur.
  - Ton d'un collègue qui explique à un autre, pas d'un manuel: « En jouant
    sur ces délais, on peut orienter le faisceau », « il suffit de ne pas
    aller trop vite ».
- Rapport acoustique (18 novembre 2026): passage à 55 km/h retiré partout
  (script, figures, tableau, texte). C'était un premier passage de test, fait
  à une seule distance. Leçon: un essai de test qui n'a pas de pendant dans
  le plan de mesure ne va pas dans le rapport.
- Rapport acoustique (18 novembre 2026): « mesurées au ruban » → « mesurées
  avec un galon à mesurer ». Leçon: nommer l'instrument comme l'équipe le
  nomme vraiment.
- Rapport acoustique (18 novembre 2026), Tableau 1: « $d$ au ruban » →
  « $d$ réelle », « $v$ visée » → « $v$ réelle ». Leçon: un en-tête de
  tableau dit ce que la grandeur est, pas avec quel instrument elle a été
  obtenue.
- Rapport acoustique (18 novembre 2026): « sonomètre » → « microphone »
  dans le texte et le schéma. Le sonomètre est seulement nommé une fois pour
  dire d'où vient le microphone. Leçon: nommer la pièce qui sert vraiment
  (le micro), pas l'appareil complet.
- Rapport acoustique (18 novembre 2026): étiquette « trajectoire » enlevée
  du schéma. Leçon: ne pas étiqueter ce qui est évident sur le dessin (la
  route est la trajectoire).
- Rapport acoustique (18 novembre 2026), légende du schéma: « Vue de haut du
  montage (pas à l'échelle). La voiture roulait à vitesse constante sur
  environ 10 m avant et 10 m après le point P » → « Schéma de la
  configuration de l'expérience. La voiture roule à vitesse constante v en
  ligne droite et passe par le point P, le plus proche du microphone, à une
  distance d. » La zone de 10 m avant/après a aussi été retirée du schéma, du
  résumé, de la théorie et de la méthodologie. Elle n'était pas prévue au
  protocole et sa longueur variait: on en parle seulement dans la discussion,
  comme cause des écarts avec la théorie. Leçons: une légende de schéma
  décrit le modèle idéal, sans « pas à l'échelle ». Une limite de l'expérience
  qui n'était pas planifiée ne se présente pas comme une partie du montage.
- Rapport acoustique (19 novembre 2026): spectrogrammes de 50 km/h
  (« un peu déformé »), puis de 60 km/h (« encore plus déformé »),
  remplacés par ceux de 10 km/h, puis de 30 km/h, puis de nouveau de
  50 km/h (avec une fenêtre d'analyse plus longue et l'axe centré). Leçon: pour
  une figure d'exemple, l'auteur veut une raie propre et régulière, centrée
  dans le cadre (pas une petite courbe dans le haut d'un grand cadre vide).
- Rapport acoustique (19 novembre 2026), légende de la Figure 3: « Fréquence
  reçue en fonction du temps autour du passage, à 3,53 m (en haut) et à
  8,97 m (en bas), pour chaque vitesse visée. Les points gris sont les
  fréquences suivies dans les trois régions et la ligne noire le modèle avec
  la vitesse et la distance mesurées. » → « Fréquence dominante du spectre
  en fonction du temps pour les passages à chaque vitesse et à chaque
  distance. La ligne noire est le résultat du curve fit et les lignes
  pointillées marquent 440 Hz. » Étiquettes des cases: « v = » et « d = »
  devant les valeurs. Leçons: dans une légende, nommer la grandeur par ce
  qu'elle est (fréquence dominante du spectre), pas par l'étape d'analyse
  (« suivies dans les trois régions »). Pas de positions (en haut/en bas)
  quand chaque case porte son étiquette, et l'étiquette nomme le symbole.
- Rapport acoustique (19 novembre 2026), résumé: « Quand une voiture nous
  croise, son bruit baisse d'un coup » → « ... la fréquence de son bruit
  baisse ». Leçon: être exact sur la grandeur physique. « Le bruit baisse »
  se lit comme le volume, alors que l'effet Doppler change la fréquence.
- Rapport acoustique (19 novembre 2026), résumé: « devant un microphone, à
  3,53 m ou à 8,97 m de sa trajectoire, à des vitesses visées entre 10 et
  70 km/h » → « devant un microphone à 3,53 m et 8,97 m de sa trajectoire et
  à des vitesses entre 10 km/h et 70 km/h ». Leçon: dans le résumé, garder
  la phrase simple (pas de « visées », la nuance va dans la méthodologie),
  et coordonner avec « et » plutôt qu'enchaîner des virgules.
- Rapport acoustique (19 novembre 2026), résumé: « La fréquence reçue a été
  suivie dans le spectrogramme de chacun des 10 passages. » → « Les
  spectrogrammes des 10 passages ont été enregistrés pour y analyser le
  profil des courbes Doppler. » Leçon: dans le résumé, dire le but de
  l'analyse (le profil de la courbe Doppler) plutôt que l'opération
  technique (suivre une fréquence).
- Rapport acoustique (19 novembre 2026), résumé: phrase « La vitesse a été
  mesurée sur les deux bouts de la courbe, où la fréquence se stabilise, et
  la distance sur la partie centrale, où la fréquence chute » retirée et
  fusionnée avec la précédente (« ...pour analyser le profil des courbes
  Doppler et en tirer la vitesse de la voiture et la distance entre le
  microphone et sa trajectoire »). Leçon: le résumé dit quoi on mesure, pas
  comment. La méthode est expliquée plus loin, ne pas la dire deux fois.
- Rapport acoustique (19 novembre 2026), résumé: trois phrases d'écarts de
  vitesse par groupe de passages (10 km/h à part, 30 à 70 km/h, queue
  d'approche à moins de 6 %, freinage) → « Pour les cinq passages à
  3,53 m, la vitesse mesurée s'écarte en moyenne de 30 % de la vitesse
  visée, et pour les cinq passages à 8,97 m, de 11 %. » Leçons: dans le
  résumé, un seul chiffre par série, en commençant par la première série.
  Ne pas découper en sous-groupes ni expliquer la cause, c'est le rôle de la
  discussion. Et ne jamais omettre en silence un cas qui ne suit pas la
  tendance (le passage à 10 km/h). Même chose pour la distance: « 6
  passages sur 10, incertitude de 13 à 16 % à 8,97 m, les passages où la
  vitesse changeait le plus... » → « l'écart moyen avec la distance mesurée
  avec le galon est de 29 % à 3,53 m et de 16 % à 8,97 m ».
- Rapport acoustique (19 novembre 2026), introduction: « Les résultats sont
  présentés à la Section 4. Ils sont discutés à la Section 5. La Section 6
  conclut. » → « Les résultats sont présentés à la Section 4, ils sont
  discutés à la Section 5 et la Section 6 conclut. » Leçon: la règle des
  phrases courtes pour le plan vise les phrases longues qui empilent des
  clauses. Trois phrases très courtes à la suite font haché, on peut les
  joindre en une phrase simple.
- Rapport acoustique (19 novembre 2026): règlement n° 39 de la CEE-ONU et
  zone grise retirés de la Figure 4 et de la discussion. Il ne s'applique
  pas au Canada, et l'analyse des queues montre déjà la vraie cause de
  l'écart. Leçon: une référence qui ne s'applique pas au contexte de
  l'expérience ne sert pas d'argument, même comme ordre de grandeur.
- Rapport acoustique (19 novembre 2026), légende des spectrogrammes:
  « Spectrogrammes de deux passages à 60 km/h visés, à 3,53 m (gauche) et à
  8,97 m (droite) de la trajectoire. La chute de fréquence dure environ
  0,7 s à 3,53 m et environ 1,7 s à 8,97 m. » → « Spectrogrammes des
  passages à 10 km/h, à d = 3,53 m et 8,97 m. » Leçon: légende minimale,
  pas de (gauche)/(droite) quand l'étiquette est déjà dans chaque panneau,
  pas de valeur lue sur la figure.
- Rapport acoustique (19 novembre 2026), légende du Tableau 1: « Vitesse,
  distance et fréquence de la source mesurées pour les 10 passages. » →
  « Vitesse, distance et fréquence pour les 10 passages. » Leçon: une
  légende de tableau nomme les grandeurs, les en-têtes de colonnes disent
  déjà lesquelles sont mesurées.
- Rapport acoustique (19 novembre 2026): le règlement n° 39 de la CEE-ONU
  (Europe) ne s'applique pas à une voiture vendue au Québec. Au Canada, la
  norme 101 ne fixe aucune précision pour le compteur. Leçon: vérifier
  qu'une norme citée s'applique vraiment au Québec/Canada, sinon la
  présenter comme un ordre de grandeur. Ne jamais nommer le modèle de la
  voiture dans le rapport.
