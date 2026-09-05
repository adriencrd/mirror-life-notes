# Reconnaissance des motifs microbiens miroirs par l'immunité innée humaine

**Note technique complète.** Protocole, paramètres, contrôles et limites d'une
étude *in silico* quantifiant la perte de reconnaissance des motifs microbiens de
chiralité inversée par les récepteurs de l'immunité innée humaine.

Pour une introduction sans prérequis : [README-SIMPLE.md](README-SIMPLE.md).

---

## Sommaire

1. [Résumé](#1-résumé)
2. [Contexte et état de l'art](#2-contexte-et-état-de-lart)
3. [Question, hypothèse, observable](#3-question-hypothèse-observable)
4. [Matériel](#4-matériel)
5. [Méthodes](#5-méthodes)
6. [Contrôles et validations](#6-contrôles-et-validations)
7. [Limites](#7-limites)
8. [Coûts et durées mesurés](#8-coûts-et-durées-mesurés)
9. [Reproductibilité](#9-reproductibilité)
10. [Organisation du code](#10-organisation-du-code)
11. [Traçabilité](#11-traçabilité)

---

## 1. Résumé

Une bactérie de chiralité inversée, dite « miroir », présenterait à l'immunité
innée humaine l'énantiomère de chacun de ses motifs moléculaires. Les récepteurs
qui détectent ces motifs possédant des sites de liaison chiraux, la
reconnaissance devrait chuter. Cette attente est largement partagée et sert
d'argument central à un appel international au moratoire, mais elle n'a, à notre
connaissance, jamais été quantifiée récepteur par récepteur.

Ce travail construit et valide une chaîne de calcul qui produit, pour un couple
récepteur / motif donné, la quantité

```
ΔΔG = ΔG_liaison(motif miroir) − ΔG_liaison(motif naturel)
```

par docking, dynamique moléculaire tout-atome en solvant explicite, puis
estimation d'énergie libre de liaison MM/GBSA mono-trajectoire. Le couple pilote
est **NOD1 / iE-DAP**.

L'apport méthodologique principal n'est pas l'enchaînement des outils, qui est
standard, mais le **contrôle de l'invariant central** : garantir et vérifier que
les deux bras de la comparaison ne diffèrent que par la chiralité, et par rien
d'autre. Plusieurs mécanismes silencieux violent cet invariant par défaut
(section 6.2). Chacun est ici mesuré, corrigé et verrouillé par test.

**État au 5 septembre 2026 : aucun ΔΔG n'a valeur de résultat.** La chaîne est
construite, vérifiée de bout en bout, et une première production tourne. Les
contrôles qui rendraient un chiffre interprétable, en particulier une calibration
sur valeur expérimentale connue et une expérience nulle, ne sont pas faits
(section 7).

---

## 2. Contexte et état de l'art

### 2.1 Homochiralité du vivant

La vie terrestre est homochirale : la traduction ribosomique n'incorpore que des
L-acides aminés, les acides nucléiques et les sucres sont de série D. Aucune
raison physique connue n'impose ce choix, les deux séries étant énergétiquement
équivalentes.

Cette règle admet des exceptions réelles mais circonscrites, toutes non
ribosomiques : D-alanine et D-glutamate du peptidoglycane bactérien,
méso-diaminopimélate (lui-même achiral), plusieurs antibiotiques peptidiques. La
glycine, dépourvue de centre stéréogène, n'appartient à aucune série.

Ces exceptions concernent directement ce travail : le motif étudié, l'iE-DAP, en
contient deux. Une bactérie miroir les inverserait également, ce qui impose de
construire l'énantiomère **complet** du motif et non « sa version D »
(section 5.2).

### 2.2 Vie miroir et risque immunitaire

La synthèse de biomolécules de chiralité inversée est acquise pièce par pièce :
protéines de série D, acides nucléiques de série L, polymérases miroir
fonctionnelles. L'assemblage d'un organisme entièrement inversé n'est pas
réalisable aujourd'hui, mais n'est plus jugé hors d'atteinte.

En décembre 2024, un collectif d'une quarantaine de chercheurs, dont plusieurs
lauréats du prix Nobel, a publié dans *Science* un appel à ne pas poursuivre
cette voie, accompagné d'un rapport technique de plusieurs centaines de pages.
L'un des arguments centraux est immunologique : les récepteurs de l'immunité
innée reconnaissent des motifs microbiens conservés par complémentarité de forme
et de charge dans des sites de liaison chiraux. Face à des motifs inversés, cette
reconnaissance devrait s'effondrer, privant l'hôte de sa première ligne de
défense. Les antibiotiques, eux aussi majoritairement chiraux, seraient
vraisemblablement inopérants sur des cibles inversées.

### 2.3 Le point que ce travail adresse

L'argument est qualitatif. Il énonce une direction, pas une amplitude. Nous
n'avons trouvé aucune étude publiée donnant, pour un récepteur donné, l'ampleur
de la perte de reconnaissance, ni identifiant les couples où subsisterait un
signal résiduel exploitable.

C'est cette quantification que ce travail vise, avec une contrainte de méthode :
le protocole doit pouvoir **contredire** l'hypothèse de départ.

---

## 3. Question, hypothèse, observable

**Question.** De combien l'énergie libre de liaison d'un motif microbien à son
récepteur humain se dégrade-t-elle lorsque le motif est remplacé par son
énantiomère complet ?

**Hypothèse testée.** La reconnaissance chute, et le déficit est quantifiable
récepteur par récepteur.

**Observable.**

```
ΔΔG = ΔG_liaison(miroir) − ΔG_liaison(naturel)     [kcal/mol]
```

| Signe | Interprétation |
|---|---|
| ΔΔG ≫ 0 | perte de reconnaissance : l'immunité innée serait aveugle à ce motif |
| ΔΔG ≈ 0 | signal résiduel conservé : résultat falsifiant, et piste de contre-mesure |
| ΔΔG < 0 | liaison renforcée : inattendu, à investiguer |

Les trois issues sont exploitables et publiables. Un protocole qui ne pourrait
produire qu'une seule réponse ne serait pas une expérience.

**Cadre.** Travail strictement *in silico* et défensif. Aucune manipulation de
matériel biologique, aucun protocole de synthèse. La simulation est le seul moyen
d'estimer ce risque sans produire l'objet redouté.

---

## 4. Matériel

### 4.1 Couples récepteur / motif

| Couple | Récepteur | Origine de la structure | Motif | État |
|---|---|---|---|---|
| **NOD1 / iE-DAP** | domaine LRR, résidus 650 à 953 | modèle AlphaFold, UniProt Q9Y239 | iE-DAP, PubChem CID 45480617 | **pilote, production en cours** |
| TLR1/TLR2 / Pam3CSK4 | 2Z7X, humain, 2.1 Å | PDB | Pam3CSK4, CID 130704 | récepteur et ligands préparés |
| TLR4-MD2 / lipide A | 4G8A, humain, 2.4 Å | PDB | Re-LPS | structure sélectionnée |
| TLR5 / flagelline | modèle par homologie, UniProt Q9NR61 | ColabFold, non produit | FliC, 295 résidus | FASTA préparé |
| Dectine-1 / β-glucane | abandonné | 2CL8 est murin | abandonné | abandonné |

**Ordre de traitement.** Validation de la chaîne petite molécule sur le système
le plus léger, puis montée en difficulté. Cette chaîne est réutilisée sur trois
couples sur quatre. TLR5 exige un docking protéine-protéine et une chaîne
distincte.

### 4.2 Sélection des structures

Trois critères, appliqués par un audit automatisé (`inventory.py`) : origine
humaine, résolution meilleure que 2.5 Å, ligand exploitable.

**Aucune des structures initialement envisagées ne satisfaisait les trois.**
Remplacements motivés :

- **3FXI vers 4G8A** pour TLR4. 3FXI est à 3.1 Å, hors critère, et son ligand est
  éclaté en 13 résidus sur 4 chaînes, chaînes acyles mêlées à la protéine,
  liaisons acyle-sucre absentes du fichier. 4G8A est à 2.4 Å et porte le ligand
  entier dans une seule chaîne. Réserve levée par mesure : 4G8A porte les SNP
  D299G et T399I, mais ils sont à 15.2 et 17.2 Å du ligand, hors du site.
- **3V47 vers modèle ColabFold** pour TLR5. 3V47 est une chimère de TLR5b de
  poisson-zèbre, incompatible avec l'exigence d'origine humaine.
- **Dectine-1 abandonné.** 2CL8 est le domaine lectine de la protéine murine ;
  aucune structure humaine exploitable.

### 4.3 Cible NOD1

Aucune structure expérimentale du domaine LRR de NOD1 n'existe : 2DBD, 2NZ7 et
4JQW couvrent le domaine CARD. La cible est donc le modèle AlphaFold
`AF-Q9Y239-F1-model_v6`, tronqué aux résidus 650 à 953.

Qualité du modèle sur le domaine retenu, recalculée sur les 2358 atomes :

| Métrique | Valeur |
|---|---|
| pLDDT moyen | **95.0** |
| pLDDT médian | 97.0 |
| pLDDT minimum | 73.2 |
| Atomes sous 70 | **0** |

La confiance du modèle n'est pas le maillon faible de cette cible. Le choix du
site de liaison l'est (section 7.4).

### 4.4 Motif iE-DAP

γ-D-glutamyl-méso-diaminopimélate, motif minimal du peptidoglycane reconnu par
NOD1. PubChem **CID 45480617**. Le CID 194426, souvent cité, est erroné : il
porte un pont cétone au lieu de la liaison amide.

| Propriété | Valeur |
|---|---|
| Formule à pH 7.4 | **C12H20N3O7⁻** |
| Charge nette | **−1** |
| Atomes | 42 |
| Centres stéréogènes | 3, configuration (R, R, S) |

SMILES isomérique employé, protomère physiologique :

```
C(C[C@@H](C(=O)[O-])[NH3+])C[C@H](C(=O)[O-])NC(=O)CC[C@H](C(=O)[O-])[NH3+]
```

Trois carboxylates déprotonés (pKa 2 à 4), deux amines protonées (pKa 9 à 10),
azote amide neutre. Justification en section 5.1.

---

## 5. Méthodes

Vue d'ensemble :

```
 05  Préparation et docking     modèle AlphaFold, découpe LRR, protonation pH 7.4
                                ligand depuis SMILES, vérification CIP, miroir
                                docking smina, paramètres identiques aux 2 bras
        |
 07  Dynamique moléculaire      report de pose sur référence vérifiée
                                bras miroir = reflet resuperposé du naturel
                                charges AM1-BCC partagées
                                minimisation, NVT 100 ps, NPT 500 ps, 50 ns
        |
 08  MM/GBSA                    retrait du solvant, GBn2, écrantage 150 mM
                                ΔG par image, ΔΔG et incertitude
        |
 09  Analyse                    RMSD, distance de liaison, RMSF, contacts
```

### 5.1 État de protonation

**Ligand.** Le SMILES canonique de PubChem décrit la forme neutre. À pH 7.4,
iE-DAP porte trois carboxylates et deux ammoniums, soit une charge nette de −1.

Cette correction n'est pas cosmétique et ne s'annule pas entre les deux bras :
NOD1 lit le carboxylate du DAP dans une poche basique, les charges AM1-BCC d'un
COOH n'ont rien de commun avec celles d'un COO⁻, et l'observable du projet est
précisément la complémentarité électrostatique différentielle. Une assertion sur
la formule brute et la charge nette bloque le pipeline si la forme neutre
réapparaît.

Les descripteurs de stéréochimie du SMILES sont inchangés par la déprotonation :
l'ordre des voisins de chaque centre reste identique, donc la parité et la
configuration absolue sont préservées.

**Récepteur.** Protonation à pH 7.4 par PDBFixer, qui applique des pKa de
référence par type de résidu et une heuristique de liaison hydrogène pour les
tautomères d'histidine. **Aucun calcul de pKa tenant compte de l'environnement
local n'est effectué.** C'est une limite reconnue, d'autant que le site est
fortement ionique (section 7.5).

### 5.2 Construction du motif miroir

L'énantiomère complet s'obtient par une **opération de symétrie impropre**
appliquée à l'ensemble des coordonnées, en l'occurrence une réflexion à travers
un plan, de déterminant −1. Cette approche est préférée à une inversion centre
par centre : une réflexion globale ne peut pas oublier un centre, par
construction, là où une énumération le peut.

**Vérification systématique.** `is_true_enantiomer()` exige que **chaque** centre
stéréogène ait basculé R vers S ou S vers R. Un seul centre inchangé fait échouer
la chaîne. La configuration est relue depuis la géométrie 3D, pas héritée d'une
étiquette.

**Calibration de la convention L/D.** Le signe du volume chiral signé autour du
carbone alpha n'est pas supposé mais mesuré sur structure réelle : sur le PDB
3V47, il est positif dans 1368 cas sur 1368, soit 100.0000 %.

**Cas iE-DAP.** Le motif porte déjà des centres non-L : un glutamate D et un
méso-DAP associant un centre R et un centre S. Son miroir n'est donc pas « la
version D » mais l'énantiomère de ses trois centres simultanément. Une inversion
manuelle y serait particulièrement exposée à l'erreur.

### 5.3 Docking

smina 2017.11.9, fork d'AutoDock Vina 1.1.2, ligand flexible et récepteur rigide.

| Paramètre | Valeur |
|---|---|
| Exhaustivité | 16 |
| Modes retournés | 20 |
| Graine aléatoire | 42, identique aux deux bras |
| Processeurs | 8 |
| Boîte | cube de 26 Å centré sur le centre de géométrie du LRR |

Le centre de géométrie d'un domaine LRR en fer à cheval ne tombe pas sur les
atomes mais dans la cavité concave, ce qui vise le site présumé. Il s'agit d'un
docking à l'aveugle assumé (section 7.4).

**Le docking ne tranche rien.** Mesure : le score est identique pour les deux
bras, à −5.9 kcal/mol. La fonction de score Vina ne discrimine pas deux
énantiomères, résultat attendu. Le docking sert à produire un point de départ
plausible, pas un résultat.

### 5.4 Report de pose

La sortie de smina **n'est pas la molécule qui lui a été fournie**. Le format
PDBQT d'AutoDock ne conserve que les hydrogènes polaires et ne porte aucune
charge formelle. Mesure effectuée sur la forme neutre alors en usage, 43 atomes :
30 atomes rendus sur 43, atomes lourds réordonnés.

Paramétrer cette sortie telle quelle laisse les règles de valence recompléter les
hydrogènes, ce qui **reconstruit la forme neutre en silence**, sans erreur ni
avertissement.

La pose n'est donc exploitée que pour ce qu'elle apporte réellement, la
**position des atomes lourds**, reportée sur la molécule de référence vérifiée.
L'appariement se fait sur le squelette constitutionnel, insensible à la
protonation. La stéréochimie du résultat est relue depuis la géométrie de la
pose, et non héritée de la référence, afin qu'une inversion introduite par le
docking reste détectable.

### 5.5 Pose de départ du bras miroir

Le bras miroir **n'est pas docké indépendamment**. Il part du reflet de la pose
naturelle, resuperposé par la rotation **propre** (déterminant +1, algorithme de
Kabsch avec correction de déterminant) qui minimise le RMSD des atomes lourds.

**Justification.** Le score de docking étant achiral en pratique et le champ de
force chiral-symétrique, seule la géométrie échantillonnée peut discriminer les
deux bras. Docker le miroir indépendamment reviendrait à choisir sa pose parmi
des poses quasi dégénérées, c'est-à-dire à tirer au sort : le ΔΔG mesurerait ce
tirage autant que la chiralité.

**Garde-fou.** Une rotation impropre défait exactement la réflexion et restitue
le ligand naturel. Le bras miroir simulerait alors le naturel, le ΔΔG vaudrait
zéro, et on conclurait à l'absence de différence sans qu'aucune erreur ne soit
levée. Une assertion `is_true_enantiomer()` après resuperposition interdit ce
scénario.

**Limite assumée.** Ce point de départ favorise le mode de liaison du naturel. On
mesure donc la perte de reconnaissance **dans ce mode de liaison**, et non
l'impossibilité de toute liaison (section 7.3).

### 5.6 Paramétrisation

| Composant | Champ de force |
|---|---|
| Récepteur | **ff14SB** |
| Ligand | GAFF2, version gaff-2.11 |
| Charges du ligand | AM1-BCC, antechamber (AmberTools 23.6) |
| Solvant | TIP3P |

**ff14SB et non ff19SB.** ff19SB introduit des cartes CMAP ajustées sur des
L-acides aminés. Un terme paramétré sur une seule série brise la symétrie chirale
du potentiel : le miroir d'un système n'aurait plus la même énergie que
l'original, ce que la physique interdit.

Mesure sur l'ubiquitine, E(miroir) − E(naturel) :

| Champ de force | Artefact |
|---|---|
| ff19SB | **−356.9 kJ/mol**, entièrement imputable au `CMAPTorsionForce` |
| ff14SB | **0.000000 kJ/mol** |

Ce choix est verrouillé par test. Nuance à ne pas surinterpréter : sur les
couples à ligand petite molécule, seul le ligand est reflété et il relève de
GAFF2, achiral et sans CMAP ; le récepteur reste de série L dans les deux bras,
donc le CMAP se compenserait largement. Le choix est décisif pour TLR5, où le
ligand protéique est reflété, pas pour NOD1.

**Charges AM1-BCC calculées une fois et partagées.** Mécanisme vérifié dans le
code installé : openmmforcefields appelle `assign_partial_charges("am1bcc")` sans
transmettre de conformère, et openff en génère alors un nouveau par ETKDG, sans
graine fixée. Chaque appel tire donc un conformère différent.

| Comparaison | Écart maximal sur les charges |
|---|---|
| Deux conformères de la même molécule | **0.085 e** |
| Naturel contre miroir, réflexion exacte | **< 1e-4 e** |

Facteur d'environ 1000. La réflexion ne change rien, l'hamiltonien AM1 ne
dépendant que de distances interatomiques, invariantes par réflexion ; le
conformère change tout. Les charges sont donc calculées une fois sur un
conformère connu, mises en cache, et appliquées aux deux bras. Le partage est
**exact**, pas approché, et cette exactitude est mesurée.

### 5.7 Dynamique moléculaire

| Paramètre | Valeur |
|---|---|
| Solvatation | TIP3P explicite, boîte dodécaédrique, marge 10 Å |
| Force ionique | 150 mM NaCl, système neutralisé |
| Électrostatique | PME, coupure 1.0 nm |
| Contraintes | HBonds, eau rigide |
| Masses | **HMR à 4 uma** sur les hydrogènes |
| Pas d'intégration | **4 fs** |
| Thermostat | Langevin-Middle, 300 K, 1 ps⁻¹ |
| Barostat | Monte-Carlo, 1 bar, 300 K |
| Protocole | minimisation, NVT 100 ps par paliers, NPT 500 ps, production 50 ns |
| Enregistrement | une image toutes les 10 ps, checkpoint toutes les 10 images |
| Système NOD1 | **54 738 atomes** solvatés, dont 4758 récepteur et 42 ligand |

**Boîte dodécaédrique.** −29.6 % d'atomes par rapport à un cube à marge égale,
soit +33 % de débit, mesuré. Sans contrepartie sur un soluté allongé comme un
domaine LRR.

**Répartition de masse des hydrogènes.** Chaque hydrogène est porté à 4 uma, la
masse étant prélevée sur l'atome lourd porteur. Gain mesuré **×1.99** en NVT
(204.8 vers 407.2 ns/jour). Légitimité : seules les masses changent, jamais le
potentiel. La thermodynamique d'équilibre est inchangée, les masses n'entrant pas
dans la distribution de Boltzmann des positions, et le ΔΔG est une quantité
d'équilibre. La symétrie chirale est intacte, une masse étant un scalaire
invariant par réflexion. La masse totale est conservée, vérifié par test, donc
aucune densité faussée sous barostat. Ce qui change réellement : les modes de
vibration rapides ralentissent, donc les temps de corrélation exprimés en pas
diffèrent, sans effet sur une moyenne d'équilibre.

**Journalisation de l'équilibration.** Énergie, température, volume et densité
sont enregistrés dans un fichier distinct pendant le chauffage et l'équilibration
en pression, afin que la convergence soit démontrable et non supposée.

### 5.8 Énergie libre de liaison

MM/GBSA **mono-trajectoire** : les trois états sont extraits de la même
trajectoire, le solvant explicite est retiré, et l'énergie est réévaluée en
solvant implicite.

```
ΔG_liaison ≈ ⟨E_complexe⟩ − ⟨E_récepteur⟩ − ⟨E_ligand⟩
```

| Paramètre | Valeur |
|---|---|
| Solvant implicite | GBn2, via `implicit/gbn2.xml` |
| Méthode non liée | NoCutoff |
| Écrantage ionique | κ calculé pour 150 mM, cohérent avec la production |
| Images analysées | une sur 10, après rejet des 10 premières ns |
| Précision des points d'énergie | **double**, sur GPU |

La variante mono-trajectoire est adaptée à la comparaison de deux énantiomères du
même ligand : les termes internes et entropiques du ligand se compensent
largement dans la différence. L'entropie de mode normal n'est donc pas calculée,
choix assumé à revisiter si le ΔΔG s'avère petit (section 7.2).

**Précision double, imposée par la mesure.** Le ΔG est une différence de grands
nombres, environ 0.3 % de leur magnitude, régime où toute erreur d'arrondi est
amplifiée. Comparaison sur les mêmes images :

| Précision | ΔG obtenu | Coût par image |
|---|---|---|
| double | **−27.09 kcal/mol** | 102.0 ms |
| mixte | −24.62 kcal/mol | 19.2 ms |
| simple | −24.62 kcal/mol | 17.4 ms |

La précision mixte est 5.3 fois plus rapide pour **2.48 kcal/mol de biais**, soit
l'ordre de grandeur du signal recherché. Elle est donc écartée.

**Incertitude.** Erreur standard par moyennes de blocs, cinq blocs contigus. Les
images d'une trajectoire étant autocorrélées, l'écart-type inter-images
surestimerait l'information. L'incertitude du ΔΔG est la propagation en
quadrature des deux bras. Il s'agit d'une erreur **intra-run**, borne inférieure
de l'incertitude réelle, qui ne remplace pas des répliques indépendantes.

### 5.9 Analyse de trajectoire

| Grandeur | Rôle |
|---|---|
| RMSD du squelette | convergence structurale du récepteur |
| RMSD du ligand | réaménagement de pose |
| **Distance minimale ligand-récepteur** | critère de décrochage |
| RMSF par résidu | flexibilité locale |
| Carte de contacts différentielle | résidus perdus, conservés, **apparus** |

La distance minimale est le critère honnête de décrochage : le RMSD depuis la
pose de départ confond décrochage et simple réaménagement, un petit ligand
flexible pouvant afficher 7 Å de RMSD tout en restant au contact.

Le traitement des conditions périodiques emploie la convention d'image minimale
en boîte **triclinique**, indispensable ici puisque la boîte est dodécaédrique.
Un traitement orthorhombique naïf produisait un RMSD de 90 Å pour un ligand resté
physiquement à 1.5 à 3 Å du site.

Les résidus **apparus** au bras miroir sont d'un intérêt particulier : ils
signalent un point d'accroche résiduel que l'hypothèse de départ n'anticipe pas,
donc une prise possible pour concevoir un détecteur.

---

## 6. Contrôles et validations

### 6.1 Invariant central

Les deux bras doivent différer **uniquement par la géométrie**. Vérification
effectuée en comparant terme à terme les deux systèmes OpenMM du ligand :

| Terme | Naturel contre miroir |
|---|---|
| Charges et Lennard-Jones | identiques |
| Liaisons harmoniques | identiques |
| Angles harmoniques | identiques |
| Torsions périodiques | identiques |
| Paramètres GB | identiques |
| Masses | identiques |
| **Géométrie** | **RMSD 1.388 Å** |

La réflexion est une isométrie exacte : écart maximal des distances internes
**9.08e-05 Å**, centroïdes des atomes lourds distants de **0.000005 Å**, donc
même volume et même site.

### 6.2 Mécanismes silencieux identifiés

Cinq mécanismes violaient l'invariant central ou la validité du calcul sans
lever la moindre erreur. Chacun a été trouvé en exécutant et en mesurant, non en
relisant.

| Mécanisme | Effet | Correction |
|---|---|---|
| Protomère neutre | mauvaise molécule simulée | assertion formule et charge nette |
| Aller-retour PDBQT | 30 atomes rendus sur 43, protomère reconstruit neutre | report de pose sur référence vérifiée |
| Conformère aléatoire pour AM1-BCC | jusqu'à 0.423 e d'écart entre production et scoring | charges calculées une fois et partagées |
| Docking indépendant du miroir | pose tirée au sort parmi des poses dégénérées | reflet resuperposé de la pose naturelle |
| Rotation impropre à la resuperposition | le bras miroir simulerait le naturel, ΔΔG nul | assertion d'énantiomérie après superposition |

### 6.3 Couverture par les tests

**160 tests** au 5 septembre 2026. La couverture est délibérément asymétrique :
le module d'inversion de chiralité est couvert de façon disproportionnée, une
erreur d'énantiomère invalidant silencieusement tout l'aval.

Trois modules restent sans test : `prepare`, `alphafold`, `inventory`. Le premier
produit le récepteur employé par tout le reste, ce qui en fait la lacune la plus
sérieuse.

### 6.4 Vérification de bout en bout

La chaîne complète 07, 08 et 09 est exécutée en mode réduit, dans un répertoire
isolé, avant toute production. Cette vérification a révélé deux défauts
invisibles sur des coordonnées synthétiques : une reconstruction de liaisons
incompatible avec les gabarits du champ de force aux résidus terminaux, et une
sélection qui classait le ligand parmi les atomes protéiques, produisant une
partition complète, disjointe et pourtant absurde.

### 6.5 Robustesse d'exécution

- **Contrôle pré-vol** avant toute production : absence de production
  concurrente, plateforme CUDA, espace disque, fichiers de préparation,
  protomère des deux ligands, cohérence du cache de charges, poses de départ
  constructibles et énantiomérie exacte.
- **Reprise sur checkpoint** : une production interrompue reprend au lieu de
  recommencer. Le système est **relu** et jamais re-solvaté, la solvatation
  plaçant les ions sans graine fixée : une re-solvatation produirait le même
  nombre d'atomes dans un ordre différent, ce que le chargement de checkpoint
  accepterait en appliquant les vitesses aux mauvais atomes.
- **Garde-fou d'instabilité** : la production s'interrompt dès que l'énergie
  cesse d'être finie, au lieu d'écrire des valeurs non numériques pendant des
  heures.

---

## 7. Limites

Cette section est déterminante pour l'interprétation. Aucune correction
éditoriale ne lève ces limites ; seul du calcul supplémentaire le peut.

### 7.1 Absence de calibration (bloquant)

La chaîne n'a jamais reproduit une valeur d'énergie libre de liaison connue
expérimentalement. Sans ce contrôle, un ΔΔG reste **ininterprétable en absolu**.
Un instrument non étalonné peut produire une valeur très précise et fausse.

### 7.2 Absence de répliques et d'expérience nulle (bloquant)

Une seule trajectoire par bras. Le signal attendu est de l'ordre du bruit d'un
MM/GBSA mono-trajectoire. Sans expérience nulle, c'est-à-dire sans ΔΔG calculé
entre deux répliques du **même** bras naturel, le plancher de bruit de la méthode
est inconnu, et aucune valeur ne peut être déclarée significative.

L'entropie de mode normal n'est pas calculée. Le choix est défendable si les
conformères liés des deux énantiomères sont quasi images l'un de l'autre, ce
qu'une poche chirale rend improbable. À revisiter si le ΔΔG s'avère petit.

### 7.3 Portée de la mesure

Par construction, le bras miroir démarre dans le mode de liaison du naturel. Le
ΔΔG mesure donc la **perte de reconnaissance dans ce mode de liaison**, et non
l'impossibilité de toute liaison. Répondre à la question du meilleur mode de
liaison du miroir exigerait plusieurs poses de départ et des répliques.

### 7.4 Site de liaison de NOD1

Aucun site expérimental n'est connu pour le domaine LRR de NOD1. Le site est
déduit de la géométrie du récepteur, et les résidus de contact sont extraits
*a posteriori* de la trajectoire, ce qui est circulaire. Le critère de
recouvrement de la boîte de docking n'écarte qu'un placement grossièrement faux ;
il ne valide pas le site.

C'est le caveat le plus lourd du couple pilote. Le lever demanderait un argument
externe, par homologie ou par la littérature, ou un test de robustesse au choix
de boîte.

### 7.5 Protonation du site

Le protomère du ligand a fait l'objet d'un traitement soigné. Le récepteur, lui,
est protoné par heuristique, sans calcul de pKa tenant compte de l'environnement
local. Or le site est fortement ionique, mesuré autour de la pose de départ :

| Résidu | Distance au ligand |
|---|---|
| LYS 790 | 1.6 Å |
| ARG 734 | 1.9 Å |
| GLU 816 | 3.7 Å |
| ASP 711 | 4.4 Å |
| LYS 793 | 4.5 Å |
| HIS 788 | 5.5 Å |

Dix résidus ionisables dans un rayon de 10 Å, dont deux histidines. L'état de
HIS 788 modifie directement l'électrostatique de la poche, celle-là même dont ce
travail fait son observable. C'est une asymétrie de rigueur à corriger.

### 7.6 Portée de l'argument sur le champ de force

Le choix ff14SB est correct en soi, mais il n'est décisif que pour les couples à
ligand protéique reflété, donc pour TLR5. Le présenter comme le résultat
méthodologique principal des couples petite molécule serait attaquable.

---

## 8. Coûts et durées mesurés

Matériel : NVIDIA RTX 4070, 12 Go, pilote 591.86, sous WSL2. Variance machine
observée d'un run à l'autre : environ 20 %.

### 8.1 Par étape, couple NOD1

| Étape | Rôle | Durée |
|---|---|---|
| 05 préparation et docking | cible, ligands, poses initiales | **≈ 2 min**, dont 19 s de docking par bras |
| 07 dynamique moléculaire | trajectoire des deux bras | **≈ 4.7 h par bras** |
| 08 MM/GBSA | ΔG des deux bras et ΔΔG | **≈ 2 min** |
| 09 analyse | RMSD, contacts, figures | **≈ 10 min** |

**Total d'un couple : environ 9 h 30**, dont 98 % en dynamique moléculaire. Toute
optimisation ne portant pas sur l'étape 07 est cosmétique.

Détail de l'étape 07, par bras :

| Sous-étape | Durée |
|---|---|
| Solvatation et paramétrisation | ≈ 20 s |
| Minimisation | 1 à 2 min |
| Chauffage NVT, 100 ps | ≈ 30 s |
| Équilibration NPT, 500 ps | ≈ 3 min |
| **Production, 50 ns** | **≈ 4.6 h** |

### 8.2 Débit

| Conditions | Débit | Origine |
|---|---|---|
| NVT, sans barostat | 385.8 ns/jour | mesuré, médiane de 3 |
| **NPT, production réelle** | **262.8 ns/jour** | mesuré, médiane de 3, plage 256 à 317 |
| NPT, production en cours | 254 à 266 ns/jour | observé sur le run du 5 septembre |

Le débit de production est celui à retenir : le barostat, actif pendant toute la
production, coûte environ un tiers du débit.

### 8.3 Par couple

| Couple | Système solvaté | Débit | Production 50 ns × 2 bras |
|---|---|---|---|
| NOD1 / iE-DAP | 54 738 atomes, mesuré | 263 ns/j, mesuré | **9.1 h** |
| TLR1/2 / Pam3CSK4 | 187 958 atomes, mesuré | ≈ 127 ns/j, extrapolé | ≈ 19 h |
| TLR4-MD2 / lipide A | non préparé | inconnu | non estimable |
| TLR5 / flagelline | non préparé | inconnu | non estimable |

L'extrapolation pour TLR1/2 applique au débit mesuré sans barostat le rapport
NPT/NVT mesuré sur NOD1, soit 0.68. Pour TLR5, la seule structure disponible
dérive de 3V47, écartée du projet ; aucun chiffre n'est défendable tant que le
modèle humain n'existe pas.

---

## 9. Reproductibilité

### 9.1 Environnement logiciel

| Composant | Version |
|---|---|
| Python | 3.11.15 |
| OpenMM | 8.2.0 |
| openmmforcefields | 0.15.1 |
| OpenFF Toolkit | 0.18.0 |
| RDKit | 2024.03.5 |
| MDTraj | 1.11.1 |
| MDAnalysis | 2.10.0 |
| ParmEd | 4.3.1 |
| NumPy | 1.26.4 |
| PDBFixer | 1.12 |
| AmberTools | 23.6, antechamber 22.0 |
| smina | 2017.11.9, AutoDock Vina 1.1.2 |
| Open Babel | 3.1.1 |
| Système | Linux 6.6.87.2, WSL2, glibc 2.35 |
| GPU | NVIDIA RTX 4070, pilote 591.86 |

Environnement reconstructible depuis `environment.yml`.

### 9.2 Graines aléatoires

| Étape | Graine |
|---|---|
| Construction 3D du ligand, ETKDG | 0xC0FFEE |
| Docking smina | 42 |
| Intégrateur et barostat | 20260714 |

Reproductibilité bit à bit non garantie : `DeterministicForces` est désactivé,
son coût GPU n'étant pas justifié ici. La graine rend le tirage aléatoire
déterministe, ce qui suffit à la comparabilité des bras.

**Point non reproductible identifié** : la solvatation place les ions par tirage
aléatoire sans graine fixée. Deux solvatations du même complexe donnent le même
nombre d'atomes dans un ordre différent. Sans conséquence pour la comparaison des
bras, qui sont de toute façon solvatés séparément, mais cela interdit de
reconstruire un système identique après coup, d'où la relecture du PDB solvaté à
la reprise.

### 9.3 Exécution

Le calcul tourne sous WSL2 : AmberTools n'a pas de distribution Windows et le
calcul de charges en dépend. GPU NVIDIA via passthrough CUDA.

```bash
micromamba env create -f environment.yml
bash scripts/00_check_install.sh
```

Les données lourdes vont dans `~/miroir-data`, sur le système de fichiers natif
WSL, les entrées-sorties sur `/mnt/c` étant deux à trois fois plus lentes.
Surchargeable par la variable `MIROIR_DATA`.

```bash
# Validation technique
bash scripts/run.sh scripts/02_validate_forcefield.py   # symétrie chirale
bash scripts/run.sh scripts/03_benchmark_gpu.py         # débit réel

# Couple NOD1
bash scripts/run.sh scripts/05_dock_nod1.py             # cible, ligands, docking
bash scripts/run.sh scripts/preflight.py --ns 50        # contrôle pré-vol
bash scripts/launch.sh 50                               # production détachée
bash scripts/status.sh                                  # suivi
bash scripts/run.sh scripts/09_analyse_nod1.py          # figures et contacts

# Tests
bash scripts/test.sh                                    # 160 tests
bash scripts/test.sh -m "not slow"                      # sans champ de force
```

La production est lancée détachée par `setsid`, afin de survivre à la fermeture
du terminal. Une production interrompue reprend sur son checkpoint : il suffit de
relancer la même commande.

---

## 10. Organisation du code

```
src/miroir/
  chirality.py   réflexion des coordonnées, classification L/D calibrée
  structure.py   entrées-sorties PDB, audit de chiralité, invariants géométriques
  prepare.py     nettoyage, protonation pH 7.4, sélection de chaînes
  validate.py    test de symétrie chirale d'un champ de force
  inventory.py   audit organisme, résolution, intégrité du ligand
  alphafold.py   récupération AlphaFold DB, découpe de domaine
  ligand.py      construction 3D, miroir, vérification CIP, report de pose
  pocket.py      définition de la boîte de docking
  docking.py     interface smina, paramètres imposés identiques aux deux bras
  complex_md.py  assemblage récepteur et ligand, charges partagées
  md.py          protocole MD complet, reprise, garde-fou d'instabilité
  mmgbsa.py      énergie libre de liaison MM/GBSA
  analysis.py    RMSD, distance de liaison, RMSF, contacts, PBC triclinique
  config.py      registre des couples, décisions encodées

scripts/         points d'entrée numérotés par phase, plus lanceurs et contrôles
tests/           suite pytest, 160 tests
docs/            justification chiffrée de chaque décision, audit critique
```

---

## 11. Traçabilité

Chaque décision méthodologique est accompagnée de la mesure qui la justifie et du
test qui la verrouille.

- **[docs/00-choix-methodologiques.md](docs/00-choix-methodologiques.md)** : seize
  points, chacun avec sa mesure, y compris les corrections d'erreurs commises en
  cours de route et les optimisations rejetées par la mesure.
- **[docs/01-audit-relecture-2026-07-15.md](docs/01-audit-relecture-2026-07-15.md)** :
  relecture critique adversariale, statut vérifié dans le code pour chaque point.
  Un point du rapport de synthèse s'était révélé faux positif à la vérification,
  d'où la règle adoptée depuis : vérifier dans le code, jamais se fier à la
  synthèse.
- **[HANDOFF.md](HANDOFF.md)** : état complet pour reprise ou audit sans contexte.

### Principe de travail

Deux règles ont structuré ce projet, l'une et l'autre issues de l'expérience
plutôt que de principes affichés.

**Toute décision doit être justifiée par une mesure, pas par une intuition.**
Plusieurs raisonnements plausibles se sont révélés faux à la vérification :
l'origine des écarts de charges, l'innocuité de la précision mixte, l'effet de la
fréquence du barostat, et jusqu'au débit de production lui-même, mesuré pendant
des semaines dans une configuration qui n'était pas celle de la production.

**Du code qui n'a jamais été exécuté est présumé cassé.** Le module produisant
l'observable centrale du projet ne s'était jamais exécuté et portait deux erreurs
de construction, alors que 124 tests passaient au vert, tous portant sur
l'arithmétique et aucun sur la physique.
