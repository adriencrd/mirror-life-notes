# La vie miroir : une bactérie que ton corps ne verrait pas

*Version accessible, sans aucun prérequis. Pour les méthodes et les chiffres,
voir [README-TECHNIQUE.md](README-TECHNIQUE.md).*

---

Imagine une bactérie qui entre dans ton corps et se multiplie tranquillement.

Pas de fièvre. Pas d'inflammation. Pas de fatigue. Ton système immunitaire ne
réagit pas, non pas parce qu'il est faible, mais parce qu'il **ne la voit
pas**. Il n'a même pas conscience qu'il y a quelque chose à voir.

Cette bactérie n'existe pas. On saurait bientôt la fabriquer.

Ce dépôt essaie de répondre à la question que personne n'a encore chiffrée :
**à quel point serions-nous vraiment aveugles ?**

---

## 1. Regarde tes mains

Ta main gauche et ta main droite sont identiques. Mêmes doigts, mêmes os, mêmes
longueurs. Et pourtant tu ne peux pas enfiler un gant droit sur ta main gauche.

Elles sont **l'image l'une de l'autre dans un miroir**. Et ça suffit à les rendre
incompatibles avec les mêmes objets.

Les molécules font exactement pareil. Beaucoup existent en version « gauche » et
version « droite ». Les chimistes appellent ça la **chiralité**, du grec
*kheir*, la main.

Maintenant, le fait vraiment étrange :

> **Tout le vivant sur Terre a choisi un seul côté.**

- Toutes tes protéines sont bâties avec des acides aminés **gauches**.
- Tous tes sucres, ton ADN, sont **droits**.
- Pareil pour les bactéries, les arbres, les champignons, les baleines, les
  moisissures du frigo.

Sans exception. Depuis quatre milliards d'années.

Et personne ne sait pourquoi. Les deux versions sont aussi stables, aussi
faciles à fabriquer. La vie a simplement pris un côté au tout début, et
absolument tout ce qui a suivi en a hérité.

---

## 2. On sait maintenant construire l'autre côté

Longtemps, c'était une curiosité de laboratoire. Plus maintenant.

Des équipes savent aujourd'hui fabriquer des protéines droites, de l'ADN gauche,
des enzymes inversées. Pièce par pièce, le jeu de construction du vivant existe
désormais **en version miroir**.

La question qui se pose donc, sérieusement : **pourrait-on assembler une
bactérie entièrement inversée ?**

Ce serait de la vraie vie. Elle mangerait, se diviserait, muterait, évoluerait.
Mais chacune de ses molécules serait le reflet des nôtres.

En décembre 2024, une quarantaine de scientifiques (biologistes de synthèse,
immunologistes, plusieurs prix Nobel) ont publié dans la revue **Science** un
avertissement accompagné d'un rapport technique de plusieurs centaines de pages.
Leur message, en substance :

> **Ne construisez pas ça. Et décidons-le maintenant, pendant que c'est encore
> impossible.**

Ce n'est pas une inquiétude vague sur « la nature qu'on ne doit pas toucher ».
C'est un argument technique, et il est glaçant de simplicité.

---

## 3. Ton système immunitaire est un trousseau de serrures

Quand une bactérie franchit ta peau, ton **système immunitaire inné** la repère
en quelques minutes. C'est la garde rapprochée, bien avant les anticorps.

Comment fait-il ? Il possède des **récepteurs** : des protéines en forme de
serrure, qui reconnaissent des morceaux caractéristiques de bactéries. Un
fragment de paroi, un bout de flagelle. Ces morceaux sont les **clés**.

```
   bactérie normale                    ton récepteur
   ┌───────────┐                       ┌──────────┐
   │  morceau  │ ─── s'emboîte ───►    │ serrure  │ ──► ALERTE ! Défense !
   │  (clé)    │                       │          │
   └───────────┘                       └──────────┘
```

Voilà le point crucial : **une serrure chirale ne reconnaît qu'une clé du bon
côté.** Exactement comme le gant.

Si la bactérie est miroir, toutes ses clés sont inversées :

```
   bactérie MIROIR                     ton récepteur (inchangé)
   ┌───────────┐                       ┌──────────┐
   │  morceau  │ ─── ne rentre ──╳     │ serrure  │ ──► ... silence ...
   │  inversé  │     pas               │          │
   └───────────┘                       └──────────┘
```

Aucun signal. Pas d'alerte, donc pas de défense. La bactérie se multiplierait
sans rencontrer de résistance.

Et ce n'est pas tout : **les antibiotiques sont eux aussi chiraux.** Ils sont
conçus pour s'emboîter dans des cibles bactériennes orientées comme les nôtres.
Sur une bactérie miroir, ils ne fonctionneraient probablement pas non plus.

Ni détection, ni traitement. C'est ce scénario qui a motivé l'alerte.

---

## 4. Le trou dans le raisonnement

Et c'est ici que ce projet commence.

Le danger a été nommé, publié dans la meilleure revue du monde, signé par des
prix Nobel. Le raisonnement est solide et tout le monde le comprend.

Mais il repose sur une **attente**, pas sur une mesure.

Personne n'a écrit : *ce récepteur-là perdrait 87 % de sa sensibilité, celui-ci
seulement 20 %.* On dit « ça ne rentrerait plus ». C'est très probablement vrai.
Mais « probablement vrai » n'est pas un résultat scientifique, et surtout : ça ne
dit pas **combien**, ni **où subsiste une faille exploitable**.

> Le danger est documenté. **Le chiffre, lui, n'existe pas.**
> À notre connaissance, aucune étude publiée ne l'a calculé, récepteur par
> récepteur.

C'est ce que ce dépôt essaie de produire. Pas un avis de plus dans le débat : un
nombre, reproductible, vérifiable, avec sa barre d'erreur.

Et, c'est peut-être le plus important, **un nombre qui a le droit de
contredire l'hypothèse de départ.**

---

## 5. Ce que fait ce projet

Tout se passe **sur ordinateur**. Aucune manipulation biologique, rien qui sorte
d'une carte graphique. C'est un point de méthode, pas une limite : on ne peut pas
étudier expérimentalement une bactérie miroir sans la fabriquer, ce qui est
précisément ce qu'il ne faut pas faire.

> **La simulation est le seul moyen de mesurer le risque sans le créer.**

La méthode tient en une phrase :

> On simule un récepteur humain face au morceau de bactérie **normal**, puis face
> au même morceau **inversé**, et on compare la force du collage.

### La mesure : le ΔΔG

On calcule une **énergie de liaison** : à quel point la clé tient dans la
serrure. Puis on fait la différence entre les deux versions.

Ce nombre s'appelle le **ΔΔG** (« delta-delta-G »). C'est le résultat central.

| Résultat | Ce que ça voudrait dire |
|---|---|
| **ΔΔG grand et positif** | Le miroir colle beaucoup moins bien → **l'inquiétude est fondée**, on serait bel et bien aveugles |
| **ΔΔG proche de zéro** | Le miroir colle presque pareil → **il reste du signal**, donc une prise pour se défendre |
| **ΔΔG négatif** | Le miroir collerait *mieux* → totalement inattendu, et passionnant |

Les trois réponses sont utiles, et c'est voulu. Une étude qui ne peut donner
qu'une seule réponse n'est pas une expérience : c'est une démonstration.

---

## 6. Comment on simule ça

**1. Trouver la serrure.** On récupère la structure 3D du récepteur humain, ici
NOD1, un détecteur de paroi bactérienne. Sa forme exacte n'a jamais été mesurée
en laboratoire : on utilise une **prédiction par intelligence artificielle**
(AlphaFold).

**2. Construire la clé, et son reflet.** On bâtit le morceau bactérien en 3D.
Pour obtenir le miroir, on applique une réflexion à **toutes** les coordonnées
d'un seul coup, comme un vrai miroir. Ça compte : en inversant les atomes un par
un, on finirait forcément par en oublier un, et une seule erreur invaliderait
tout le reste, en silence.

**3. Emboîter.** Un logiciel cherche la meilleure façon de poser la clé dans la
serrure.

**4. Laisser bouger.** Les molécules ne sont pas des pièces rigides : elles
vibrent, se tordent, l'eau les bouscule. On simule donc **50 nanosecondes** de
mouvement réel, à 37 °C, dans l'eau salée.

50 nanosecondes, c'est 50 milliardièmes de seconde. Autant dire rien. Sauf que
pour y arriver, l'ordinateur calcule les forces entre **54 738 atomes**, douze
millions et demi de fois de suite.

> Soit près de **700 milliards de positions atomiques** calculées, pour un seul
> bras. Le double pour la comparaison. En un peu moins de **6 heures**, sur une
> seule carte graphique.

**5. Mesurer.** On calcule l'énergie de collage sur des centaines d'instantanés,
pour les deux versions, et on soustrait.

---

## 7. Ce que ça permet

Le but n'est pas de répéter « on serait aveugles ». C'est de **le chiffrer**, et
un chiffre ouvre des portes qu'une intuition ne peut pas ouvrir.

**Mesurer l'ampleur, pas seulement le sens.** « Ça colle moins bien » ne se
transforme en rien d'actionnable. « Ça colle 3 kcal/mol moins bien, soit une
reconnaissance divisée par environ 150 », si.

**Classer les récepteurs.** Le projet en traite plusieurs. On peut donc les
ranger : lequel s'effondre complètement, lequel garde un peu de sensibilité.
C'est une **carte de vulnérabilité**, et elle dit par où commencer.

**Repérer ce qui résiste encore.** Chaque simulation produit la liste des acides
aminés du récepteur qui touchent encore le motif inversé, en trois catégories :
contacts *perdus*, *conservés*, et **apparus** : des contacts que le motif normal
ne faisait même pas. Ces derniers sont les plus précieux : ils signalent un point
d'accroche résiduel, donc une prise possible pour concevoir un **détecteur
artificiel** là où le corps ne verrait plus rien.

**Se tromper de façon détectable.** Si le ΔΔG ressort proche de zéro sur certains
couples, l'inquiétude devra être nuancée pour ceux-là. Le protocole est conçu
pour permettre ce démenti.

**Orienter les vraies expériences.** La paillasse est lente et coûteuse. Un
calcul qui dit « ces deux systèmes-là valent le détour, les trois autres non »
fait gagner des mois.

**Rester réutilisable.** La machine n'est pas taillée pour un cas unique :
n'importe quel couple récepteur / molécule peut y passer, avec les mêmes
garanties.

---

## 8. Où en est le projet, honnêtement

**Aucun résultat scientifique n'est encore sorti.**

La machine est construite, vérifiée de bout en bout, et rapide. Elle n'a pas
encore produit un chiffre auquel on ait le droit de croire. La différence est
essentielle, et elle est affichée ici plutôt que dissimulée.

Ce qui manque avant de pouvoir affirmer quoi que ce soit :

- **Un étalonnage.** On n'a pas encore vérifié que la méthode retrouve une valeur
  déjà connue par l'expérience. *Une balance qu'on n'a jamais tarée peut donner
  un chiffre très précis, et complètement faux.*
- **Des répétitions.** Une simulation est bruitée. Il faut la refaire plusieurs
  fois pour distinguer un vrai signal d'un coup de chance.
- **Une expérience « à blanc ».** Comparer le motif normal… avec lui-même. On
  devrait obtenir exactement zéro. L'écart au zéro donne le niveau de bruit de la
  méthode, donc le seuil en dessous duquel un résultat ne veut rien dire.
- **Une confirmation du bon endroit.** Pour NOD1, l'emplacement exact où la clé
  se pose n'a jamais été observé expérimentalement. On le déduit de la forme du
  récepteur. C'est raisonnable, ce n'est pas prouvé.

Tant que ces quatre points ne sont pas réglés, tout chiffre sorti d'ici est un
chiffre de mise au point, pas un résultat. C'est écrit noir sur blanc à chaque
étape du dépôt, y compris dans le code.

---

## 9. Ce que ce projet n'est pas

Ce dépôt ne contient **aucune information permettant de fabriquer quoi que ce
soit**. Il ne décrit pas comment construire une bactérie miroir, ni un
quelconque agent biologique.

Il simule des **récepteurs humains** face à des molécules inversées, pour mesurer
un déficit de détection. C'est un travail de **défense** : savoir précisément où
l'on est vulnérable est le préalable pour cesser de l'être, et, idéalement, pour
avoir un détecteur prêt bien avant que la question ne se pose vraiment.

---

## Pour aller plus loin

- **[README-TECHNIQUE.md](README-TECHNIQUE.md)** : la même chose, avec les
  méthodes, les paramètres et tous les chiffres.
- [docs/00-choix-methodologiques.md](docs/00-choix-methodologiques.md) : chaque
  décision, la mesure qui la justifie, le test qui la verrouille.
- [docs/01-audit-relecture-2026-07-15.md](docs/01-audit-relecture-2026-07-15.md),
  la relecture critique : ce qui casserait si un expert regardait de près.
