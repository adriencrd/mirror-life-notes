# Vie miroir et immunité innée — notes de recherche

Notes méthodologiques d'une étude ***in silico*** en cours : **l'immunité innée
humaine reconnaîtrait-elle encore les motifs microbiens d'une bactérie
miroir — et à quel point ?**

Le danger a été posé publiquement (revue *Science*, décembre 2024, ~40 auteurs
dont plusieurs prix Nobel, appel à un moratoire). Il repose sur une attente
qualitative : les récepteurs sont chiraux, les motifs inversés ne s'y
emboîteraient plus. À notre connaissance, **aucune étude publiée ne l'a chiffré,
récepteur par récepteur.** C'est ce que ce travail essaie de produire.

---

## Par où commencer

| Document | Pour qui | Durée |
|---|---|---|
| **[README-SIMPLE.md](README-SIMPLE.md)** | Sans aucun prérequis. La chiralité par l'analogie du gant, ce qu'est la vie miroir, pourquoi le corps serait aveugle, ce que le projet permet. | ~12 min |
| **[README-TECHNIQUE.md](README-TECHNIQUE.md)** | Pipeline complète, paramètres de production, chaque décision adossée à sa mesure, limites connues. | ~25 min |
| [docs/00-choix-methodologiques.md](docs/00-choix-methodologiques.md) | Les 14 décisions méthodologiques, chacune avec la mesure qui la justifie et le test qui la verrouille. | référence |
| [docs/01-audit-relecture-2026-07-15.md](docs/01-audit-relecture-2026-07-15.md) | Relecture critique adversariale, statut vérifié dans le code point par point. | référence |

---

## L'observable

```
ΔΔG = ΔG_liaison(motif miroir) − ΔG_liaison(motif naturel)
```

| Signe | Lecture |
|---|---|
| ΔΔG ≫ 0 | perte de reconnaissance : l'immunité innée serait aveugle |
| ΔΔG ≈ 0 | signal résiduel — résultat **falsifiant**, et piste de contre-mesure |
| ΔΔG < 0 | liaison renforcée : inattendu, à investiguer |

Les trois issues sont publiables. Le protocole est conçu pour pouvoir contredire
l'hypothèse de départ.

---

## État — aucun résultat publiable à ce jour

La chaîne de calcul est construite, vérifiée de bout en bout et rapide. Elle n'a
pas encore produit de chiffre auquel on ait le droit de croire. Ce qui manque,
par ordre d'importance :

- **contrôle de calibration** — la chaîne n'a jamais reproduit un ΔG expérimental
  connu ; sans lui, un chiffre reste ininterprétable en absolu ;
- **répliques et expérience nulle** (ΔΔG naturel-vs-naturel) — le signal attendu
  est de l'ordre du bruit d'un MM/GBSA mono-trajectoire ;
- **convergence jugée sur l'énergie**, pas seulement sur le RMSD ;
- **site de liaison de NOD1** — aucun site expérimental n'est connu pour le
  domaine LRR ; il est déduit de la forme du récepteur. Raisonnable, non prouvé.

Ces limites sont énoncées ici plutôt que dissimulées : elles font partie du
travail.

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
le préalable pour cesser de l'être — et pour disposer d'un moyen de détection
avant que la question ne se pose réellement.

---

*Travail de recherche indépendant, en cours. Les documents évoluent avec le
projet.*
