# Audit de relecture critique — 15 juillet 2026

Relecture adversariale du projet menée par 4 agents parallèles (stéréochimie,
champ de force, Phase 0/MD, ancrage TLR1/2), **puis vérifiée point par point dans
le code réel** avant toute action. Objectif : trouver ce qui casserait en
relecture scientifique par des pairs, pas valider.

Chaque point porte un statut (**CONFIRMÉ** avec la ligne qui le prouve /
**INFIRMÉ** / **INCERTAIN**), une classe, et un coût :

- **A** — nécessaire pour que le pilote NOD1 soit *informatif*
- **B** — nécessaire avant tout ΔΔG *publiable*, pas pour le pilote
- **C** — vrai mais mineur

⚠️ Un point du rapport de synthèse s'est révélé **faux positif** à la
vérification (2.2, SASA) : les agents n'avaient pas ouvert le fichier de champ de
force. D'où la règle : vérifier dans le code, pas se fier à la synthèse.

Aucun ΔΔG n'a encore été produit — `mmgbsa.py` n'a jamais tourné de bout en
bout. Toute la critique porte sur le **protocole tel que codé**.

---

## TIER 1

### 1.1 — Contrôle Phase 0 (reproduire un ΔG connu ±1 kcal/mol) absent
**CONFIRMÉ.** `scripts/` va de `00` à `12` sans script de contrôle antibiotique ;
`results/` ne contient que symétrie-FF, benchmark GPU, inventaire. La chaîne
docking→MD→MM/GBSA n'a jamais été confrontée à une valeur expérimentale. Or
c'est le critère de validité n°1 du CDC (§11).
→ **Classe B** — le pilote peut tourner sans, mais son chiffre reste non
interprétable en absolu. **Coût : vraie refonte** (trouver un complexe à ΔG connu,
passer toute la chaîne dessus — jours).

### 1.2 — Pas de barre d'erreur / 1 réplique / pas de seed / pas d'expérience nulle
**CONFIRMÉ.** `mmgbsa.py:170` renvoie un `float` nu. `md.py:118-122` et
`complex_md.py:99-103` : intégrateur sans `setRandomNumberSeed`. `md.py:80` :
`DeterministicForces=false`. Aucune expérience nulle nat-vs-nat (le « ΔΔG » entre
deux répliques du même énantiomère, qui donne le plancher de bruit). Le signal
attendu (~qq kcal/mol) est de l'ordre du bruit inter-images d'un MM/GBSA
mono-trajectoire.
→ **Classe A** — c'est ce qui distingue un pilote informatif d'un nombre au
hasard. **Coût mixte** : seeds + propagation d'erreur = quelques lignes ;
répliques + expérience nulle = coût *calcul* (×3-4 GPU).

### 1.3 — Charges AM1-BCC non garanties identiques entre bras
**CONFIRMÉ.** `complex_md.py:58` `cache=None` ; `assemble()` appelé par bras
(`07:88`) → antechamber tourne 2× séparément sur 2 géométries. `mmgbsa.py:122-129`
**recrée** les charges 3×/bras au scoring, indépendamment de la production. La
ligne `complex_md.py:9-13` affirme le contraire — faux pour les charges.
Nuance : pour un ligand rigide les charges d'énantiomères sont quasi identiques ;
risque réel surtout pour iE-DAP (flexible, chargé) et Pam3CSK4 (86 rot. bonds).
Ampleur incertaine, mécanisme réel.
→ **Classe A** — assurance peu coûteuse qui élimine tout un doute. **Coût moyen** :
charger une fois sur le naturel, copier bit-à-bit sur le miroir, asserter,
réutiliser en MM/GBSA (~30-40 lignes + réorganisation).

### 1.4 — Chaîne de custody stéréochimique rompue au docking
**CONFIRMÉ.** `07:56-64` prend la 1ʳᵉ pose de `dock_{arm}.sdf` (sortie
smina/OpenBabel, format PDBQT sans stéréo explicite) → `pose_depart_{arm}.sdf` →
`assemble()` (`07:88`). Aucun `cip_codes`/`is_true_enantiomer` entre la sortie
smina et l'entrée en MD. Une inversion/perte de centre par OpenBabel passerait
sans erreur ni traceback.
→ **Classe A** — si la mauvaise molécule entre en MD, le pilote compare autre
chose. **Coût : quelques lignes** — recontrôler les CIP sur `pose_depart` avant
`assemble()`, échouer si divergence.

### 1.5 — ΔΔG déterminé par la pose de départ, et les poses sont biaisées
**CONFIRMÉ.** Puisque le docking est achiral (−5.9 = −5.9) et le FF
chiral-symétrique, seule la géométrie échantillonnée discrimine → tout dépend de
la pose initiale.
- NOD1 : `05:102,107` docke naturel et miroir **indépendamment** ; `07:80` prend
  la meilleure de chaque → confond « placement différent » et « chiralité
  différente ».
- TLR2 (`docs point 9`) : part de la pose cristallo du *naturel* pour les deux
  bras → injecte la réponse ; 50 ns ne ré-équilibrent pas 86 liaisons rotatives.
→ **Classe A pour NOD1** (le pilote). **Coût moyen** — amorcer le miroir depuis la
pose naturelle *réfléchie et resuperposée*, pas un docking indépendant (~30-50
lignes). Volet TLR2 = **B**.

---

## TIER 2

### 2.1 — L'argument ff14SB vs ff19SB ne « mord » pas sur NOD1
**CONFIRMÉ (conceptuel).** Seul le ligand est reflété, il est en GAFF2
(`complex_md.py:31`, achiral, pas de CMAP) ; le récepteur reste L dans les deux
bras. Le CMAP se compense dans chaque bras (E_complexe − E_récepteur) et au
second ordre dans le ΔΔG. Donc ff19SB donnerait quasi le même ΔΔG **et** serait
plus exact sur le récepteur L. Le choix ff14SB est correct mais **décisif
uniquement pour TLR5** (ligand protéique reflété). Le vendre comme « résultat
méthodologique le plus important » est attaquable : il n'affecte pas le chiffre
produit sur les 3 couples petite-molécule.
→ **Classe C** — le ΔΔG reste valide ; c'est la survente qui est attaquable.
**Coût : documentation** (trivial), ou repasser en ff19SB pour les couples
petite-molécule (quelques lignes + rerun).

### 2.2 — « MM/GB, pas MM/GBSA : terme SASA manquant »
**INFIRMÉ.** Le fichier réel
`openmm/app/data/implicit/gbn2.xml` contient `solventArgs = {'SA':'ACE'}` : le
terme non-polaire de surface (méthode ACE) **est inclus par défaut** dans le GBn2
tel qu'utilisé (`mmgbsa.py:73`, `implicitSolvent=GBn2`). Il y a bien un terme
non-polaire. Les agents l'ont supposé absent sans ouvrir le fichier.
→ Reste au plus **C** : l'ACE est une approximation grossière vs un vrai SASA
moléculaire (type LCPO de MMPBSA.py) — à mentionner, pas à corriger.

### 2.3 — MM/GBSA sur pose non convergée + critère de convergence faible
**CONFIRMÉ.** `analysis.py::assess_convergence` = heuristique à 2 seuils codés en
dur, appliqué au **RMSD squelette**, jamais à la série d'énergie. `08_mmgbsa_nod1.py`
calcule le ΔΔG **sans vérifier** la convergence de `per_frame`. Or le RMSD ligand
est « non convergé » (réaménagement de pose).
→ **Classe A** — juger la convergence de l'*énergie* avant de croire le ΔΔG.
**Coût : quelques lignes** (test sur `per_frame`) ; potentiellement plus de
sampling (calcul).

### 2.4 — iE-DAP modélisé neutre (pas pH 7,4)
**CONFIRMÉ.** SMILES iE-DAP dans `config.py` : trois `C(=O)O` protonés + amines
`N` neutres, net = 0. À pH 7,4 c'est un anion/zwitterion (net ≈ −1). Le récepteur,
lui, est protoné à pH 7,4 (`clean_pdb`) → incohérence. Aucune étape de pKa. La
reconnaissance d'iE-DAP par NOD1 repose sur ses carboxylates.
→ **Classe A** — c'est le ligand du pilote. **Coût moyen** — fixer l'état de
protonation, encoder les charges dans le SMILES, passer le net à antechamber,
re-docker/re-run.

### 2.5 — Entropie négligée
**CONFIRMÉ** (`mmgbsa.py:18-21`, choix assumé). Défendable seulement si les
conformères liés des deux énantiomères sont quasi-images — ils ne le sont pas
(poche chirale). Or le ΔΔG sera petit, régime où l'entropie compte.
→ **Classe B** — tolérable pour le pilote si documenté. **Coût réel** (mode normal
/ interaction entropy).

### 2.6 — Équilibration sans restraints et non journalisée
**CONFIRMÉ.** `md.py:186-202` : minimisation → `heat` → `equilibrate_npt` →
production, **aucune contrainte de position**. `attach_reporters` seulement à la
production (`md.py:201`) → aucune trace de l'équilibration. Pose dockée lâchée dès
le chauffage.
→ **Classe A** pour le volet « non journalisée » (prouver la convergence de
l'équilibration). Restraints sur pose dockée = **A/B**. **Coût** : logging trivial ;
restraints ~20-40 lignes.

### 2.7 — GBn2 sans écrantage ionique
**CONFIRMÉ.** `mmgbsa.py:63-76` : aucun `implicitSolventSaltConc` → 0 M en
implicite, alors que la production est à 150 mM (`md.py:39`). Surestime
l'électrostatique de liaison pour un ligand chargé.
→ **Classe B**, mais **coût : quelques lignes** (un kwarg
`implicitSolventSaltConc=0.15*molar`). À faire tôt, c'est gratuit.

### 2.8 — Site de docking NOD1 = pari, `box_covers_selection` faible
**CONFIRMÉ (nuancé).** `pocket.py:43` : fraction ≥ 0.15 n'attrape qu'une boîte
grossièrement mal placée, ne valide **pas** le site. `concave_face_box` centre sur
le barycentre du LRR ; aucun site expérimental pour NOD1. Les 6 résidus de contact
sont extraits *a posteriori* de la MD (circulaire). Nuance : la MD montre une
poche stable — pas absurde, mais rien d'externe ne confirme que c'est LE site.
→ **Classe B** — caveat d'interprétation avant toute affirmation sur NOD1.
**Coût réel** (argument homologie/littérature, ou test de robustesse au choix de
boîte).

### Bonus (latent) — pièges parmed dans `mmgbsa.py`
**CONFIRMÉ latent.** `mmgbsa.py:120-133` : cohérence d'ordre atomique entre
mdtraj / ParmEd / tableau xyz jamais assertie ; code jamais exécuté. Un décalage
d'indices récepteur/ligand donnerait un ΔG faux **sans plantage**.
→ **Classe A** — le pilote calcule via ce code. **Coût moyen** — assertions
(`len(rec)+len(lig)==natoms`, égalité d'ordre) + test sur copie courte avant la
vraie prod.

---

## Récapitulatif

| # | Statut | Classe | Coût |
|---|---|---|---|
| 1.1 contrôle Phase 0 | CONFIRMÉ | B | refonte |
| 1.2 barre d'erreur / seed / nulle | CONFIRMÉ | A | code trivial + calcul |
| 1.3 charges non identiques | CONFIRMÉ | A | moyen |
| 1.4 custody stéréo | CONFIRMÉ | A | quelques lignes |
| 1.5 biais pose départ (NOD1) | CONFIRMÉ | A | moyen |
| 2.1 ff14SB non pertinent NOD1 | CONFIRMÉ | C | doc |
| 2.2 SASA manquant | **INFIRMÉ** | — | — |
| 2.3 convergence énergie | CONFIRMÉ | A | quelques lignes (+calcul) |
| 2.4 protonation iE-DAP | CONFIRMÉ | A | moyen |
| 2.5 entropie | CONFIRMÉ | B | réel |
| 2.6 équilibration (log/restraints) | CONFIRMÉ | A | trivial / moyen |
| 2.7 sel GBn2 | CONFIRMÉ | B | trivial |
| 2.8 site NOD1 | CONFIRMÉ (nuancé) | B | réel |
| parmed (latent) | CONFIRMÉ | A | moyen |

**Classe A (pilote non informatif sans ça)** : 1.2, 1.3, 1.4, 1.5, 2.3, 2.4, 2.6,
parmed.

Fixes bon marché à faire **avant même de relancer une production qu'on prétendra
analyser** : 1.4 (custody, ~5 lignes), 2.6 (logging équilibration), 2.7 (sel, un
kwarg), 1.2-seeds (`setRandomNumberSeed`).

## Ce qui tient (vérifié)
- Réflexion mono-axe = énantiomère exact (opération impropre, det = −1).
- Isométrie du miroir vérifiée par test.
- Termes liés/LJ GAFF identiques entre bras (achiraux).
- Choix ff14SB correct en soi (juste survendu, cf. 2.1).
- Terme non-polaire présent dans GBn2 (cf. 2.2 infirmé).
