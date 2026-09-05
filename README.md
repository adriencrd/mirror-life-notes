# Vie miroir et immunité innée : notes de recherche

Étude ***in silico*** en cours, indépendante. Elle cherche à répondre à une
question posée publiquement mais jamais chiffrée :

> **L'immunité innée humaine reconnaîtrait-elle encore les motifs microbiens
> d'une bactérie de chiralité inversée, et de combien la reconnaissance
> chuterait-elle ?**

---

## Le point de départ

En décembre 2024, une quarantaine de chercheurs, dont plusieurs lauréats du prix
Nobel, ont publié dans *Science* un appel à ne pas construire d'organisme de
chiralité inversée, accompagné d'un rapport technique de plusieurs centaines de
pages. L'un de leurs arguments centraux est immunologique : les récepteurs de
l'immunité innée reconnaissent les motifs microbiens par complémentarité de forme
et de charge dans des sites de liaison **chiraux**. Face à des motifs inversés,
cette reconnaissance devrait s'effondrer.

L'argument est solide, et qualitatif. Il énonce une direction, pas une amplitude.
À notre connaissance, **aucune étude publiée ne donne, récepteur par récepteur,
l'ampleur de la perte** ni n'identifie les couples où subsisterait un signal
résiduel exploitable.

C'est cette quantification que ce travail vise.

---

## Par où commencer

| Document | Pour qui | Durée |
|---|---|---|
| **[README-SIMPLE.md](README-SIMPLE.md)** | Sans aucun prérequis. La chiralité par l'analogie du gant, ce qu'est la vie miroir, pourquoi le corps serait aveugle, ce que le projet permet. | ~12 min |
| **[README-TECHNIQUE.md](README-TECHNIQUE.md)** | Note technique complète : contexte, matériel, méthodes détaillées, contrôles, limites, reproductibilité, coûts mesurés. | ~30 min |
| [docs/00-choix-methodologiques.md](docs/00-choix-methodologiques.md) | Seize décisions méthodologiques, chacune avec la mesure qui la justifie et le test qui la verrouille. | référence |
| [docs/01-audit-relecture-2026-07-15.md](docs/01-audit-relecture-2026-07-15.md) | Relecture critique adversariale, statut vérifié dans le code point par point. | référence |

---

## L'observable

```
ΔΔG = ΔG_liaison(motif miroir) − ΔG_liaison(motif naturel)     [kcal/mol]
```

| Signe | Interprétation |
|---|---|
| ΔΔG ≫ 0 | perte de reconnaissance : l'immunité innée serait aveugle à ce motif |
| ΔΔG ≈ 0 | signal résiduel conservé : résultat falsifiant, et piste de contre-mesure |
| ΔΔG < 0 | liaison renforcée : inattendu, à investiguer |

Les trois issues sont exploitables. Le protocole est conçu pour pouvoir
**contredire** l'hypothèse de départ : une étude qui ne pourrait produire qu'une
seule réponse ne serait pas une expérience.

---

## La chaîne de calcul

```
Préparation      modèle AlphaFold du récepteur, découpe du domaine, protonation
   et docking    ligand construit depuis son SMILES, vérification de chiralité
                 miroir par réflexion globale des coordonnées
                 docking, paramètres strictement identiques aux deux bras
        |
Dynamique        solvatation explicite, 54 738 atomes
moléculaire      minimisation, chauffage, équilibration en pression
                 production 50 ns par bras, sur GPU
        |
MM/GBSA          retrait du solvant, réévaluation en solvant implicite
                 ΔG de chaque bras, puis ΔΔG avec son incertitude
        |
Analyse          le ligand est-il resté lié, quels contacts sont perdus,
                 conservés, ou apparus
```

Coût mesuré : environ **9 h 30 par couple** sur une RTX 4070, dont 98 % en
dynamique moléculaire.

L'apport méthodologique n'est pas l'enchaînement des outils, qui est standard,
mais le contrôle d'un invariant : garantir et **vérifier** que les deux bras de
la comparaison ne diffèrent que par la chiralité, et par rien d'autre. Cinq
mécanismes violaient cet invariant par défaut, sans lever la moindre erreur.
Chacun est documenté, mesuré et corrigé.

---

## État : aucun résultat publiable à ce jour

La chaîne est construite, vérifiée de bout en bout, couverte par 160 tests. La
toute première production a démarré le **5 septembre 2026**.

Ce qui manque avant qu'un chiffre soit interprétable, par ordre d'importance :

**Bloquants**

- **Contrôle de calibration.** La chaîne n'a jamais reproduit une énergie libre
  de liaison connue expérimentalement. Sans lui, un chiffre reste
  ininterprétable en absolu : un instrument non étalonné peut produire une valeur
  très précise et fausse.
- **Répliques et expérience nulle.** Sans ΔΔG calculé entre deux répliques du
  *même* bras naturel, le plancher de bruit de la méthode est inconnu, et aucune
  valeur ne peut être déclarée significative.

**Caveats lourds**

- **Le site de liaison de NOD1 est déduit, pas observé.** Aucun site expérimental
  n'existe pour le domaine LRR ; les résidus de contact sont extraits
  *a posteriori* de la trajectoire, ce qui est circulaire.
- **Le récepteur est protoné par heuristique**, sans calcul de pKa local, alors
  que le site compte dix résidus ionisables dans un rayon de 10 Å, dont deux
  histidines.
- **Le ΔΔG mesure la perte dans le mode de liaison du naturel**, par
  construction, et non l'impossibilité de toute liaison.

Ces limites sont énoncées ici plutôt que dissimulées : elles font partie du
travail, et aucune correction éditoriale ne les lève.

---

## Principe de travail

**Toute décision doit être justifiée par une mesure, pas par une intuition.**
Plusieurs raisonnements plausibles se sont révélés faux à la vérification :
l'origine des écarts de charges partielles, l'innocuité de la précision mixte,
l'effet de la fréquence du barostat, et jusqu'au débit de production lui-même,
longtemps mesuré dans une configuration qui n'était pas celle de la production.

**Du code qui n'a jamais été exécuté est présumé cassé.** Le module produisant
l'observable centrale du projet ne s'était jamais exécuté et portait deux erreurs
de construction, alors que 124 tests passaient au vert, tous portant sur
l'arithmétique et aucun sur la physique.

---

## Ce dépôt

Il contient **les notes et la méthodologie**. Le code source vit dans un dépôt
privé séparé et n'est pas publié à ce stade.

Ce dépôt ne contient **aucune information permettant de fabriquer quoi que ce
soit** : ni protocole de synthèse, ni méthode de construction d'un organisme
miroir ou d'un quelconque agent biologique. Il documente la simulation de
**récepteurs humains** face à des molécules inversées, pour mesurer un déficit de
détection.

C'est un travail de **défense** : savoir précisément où l'on est vulnérable est
le préalable pour cesser de l'être, et pour disposer d'un moyen de détection
avant que la question ne se pose réellement.

---

*Travail de recherche indépendant, en cours. Les documents évoluent avec le
projet.*
