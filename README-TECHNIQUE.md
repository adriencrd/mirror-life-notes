# Note technique : carte computationnelle de l'évasion immunitaire par la vie miroir

*Version détaillée. Pour une introduction sans prérequis, voir
[README-SIMPLE.md](README-SIMPLE.md).*

---

## 1. Question scientifique et observable

L'homochiralité du vivant est quasi universelle : traduction ribosomique
exclusivement en L-acides aminés, sucres et acides nucléiques en D. Les
exceptions connues sont non ribosomiques et circonscrites, mais réelles, et
elles concernent directement ce projet : le peptidoglycane bactérien contient du
D-glutamate, de la D-alanine et du méso-DAP (lui-même achiral). La glycine, elle,
n'a pas de centre stéréogène. Une bactérie miroir inverserait aussi ces
briques-là, ce qui rend l'énantiomère complet, et non « la version D », la seule
construction correcte (voir §4.2). Les récepteurs de l'immunité innée (TLR, NLR)
reconnaissent des motifs microbiens conservés (PAMP) par **complémentarité de
forme et de charge dans un site de liaison chiral**. Une bactérie miroir
présenterait l'énantiomère de chaque PAMP.

**Hypothèse testée :** la reconnaissance chute, et le déficit est quantifiable.

**Observable :**

```
ΔΔG = ΔG_liaison(motif miroir) − ΔG_liaison(motif naturel)
```

| Signe | Lecture |
|---|---|
| ΔΔG ≫ 0 | perte de reconnaissance : l'immunité innée serait aveugle |
| ΔΔG ≈ 0 | signal résiduel, résultat **falsifiant**, et piste de contre-mesure |
| ΔΔG < 0 | liaison renforcée : inattendu, à investiguer |

Le design est intentionnellement falsifiable : les trois issues sont
publiables. Recherche strictement *in silico* et défensive.

---

## 2. Couples récepteur / motif

| Couple | Cible | Source | Ligand | État |
|---|---|---|---|---|
| **NOD1 / iE-DAP** | LRR, res 650–953 | AlphaFold Q9Y239 (pLDDT 95.0) | iE-DAP, CID 45480617 | **pilote, prêt à produire** |
| TLR1/2 / Pam3CSK4 | 2Z7X (humain, 2.1 Å) | PDB | Pam3CSK4, CID 130704 | récepteur + ligands prêts |
| TLR4-MD2 / lipide A | 4G8A (humain, 2.4 Å) | PDB | Re-LPS | structure choisie |
| TLR5 / flagelline | ColabFold Q9NR61 | modèle | FliC (295 res) | FASTA préparé |
| Dectine-1 / β-glucane | aucune | aucune | aucun | abandonné (structure murine) |

Ordre de traitement décidé (≠ ordre du cahier des charges) : valider la pipeline
petite-molécule sur le système le plus léger, puis monter en difficulté. Elle est
réutilisée sur 4 couples sur 5 ; TLR5 exige un docking protéine–protéine
(LightDock) et une chaîne distincte.

Le modèle AlphaFold du LRR de NOD1 est solide là où on l'utilise : pLDDT moyen
**95.0** sur les 2358 atomes du domaine découpé, minimum 73.2, et **aucun atome
sous 70**. La confiance du modèle n'est donc pas le maillon faible ; le choix du
site de liaison l'est (§7).

Coût GPU d'une production 50 ns × 2 bras, par couple :

| Couple | Système solvaté | Débit | Production 50 ns × 2 bras |
|---|---|---|---|
| NOD1 / iE-DAP | 54 738 at. (mesuré) | 263 ns/j NPT (mesuré) | **9.1 h** |
| TLR1/2 / Pam3CSK4 | 187 958 at. (mesuré) | ≈ 127 ns/j (extrapolé) | ≈ 19 h |
| TLR4-MD2 / lipide A | non préparé | inconnu | non estimable |
| TLR5 / flagelline | non préparé | inconnu | non estimable |

« Extrapolé » pour TLR1/2 : le débit a été mesuré sans barostat (186 ns/jour) et
ramené aux conditions de production par le rapport NPT/NVT mesuré sur NOD1
(0.68). Pour TLR5, la seule structure disponible est le complexe issu de 3V47,
**écarté du projet** parce que chimérique et non humain ; le modèle ColabFold du
TLR5 humain n'existe pas encore, donc aucun chiffre de durée n'est défendable.

**Aucune structure prescrite à l'origine ne passait tous les critères**
(humain + résolution < 2.5 Å + ligand exploitable). Remplacements motivés par
audit automatisé (`inventory.py`) : 3FXI (3.1 Å, ligand éclaté sur 13 résidus et
4 chaînes) → 4G8A ; 3V47 (poisson-zèbre chimérique) → modèle ColabFold humain.

---

## 3. Pipeline

```
 05  préparation + docking      AlphaFold → LRR → protonation pH 7.4
                                iE-DAP depuis SMILES → miroir → vérif CIP
                                smina, exhaustivité 16, graine fixe, 2 bras
        │
 07  dynamique moléculaire      pose reportée sur la référence vérifiée
                                miroir = reflet resuperposé du naturel
                                charges AM1-BCC partagées (cache)
                                min → NVT 100 ps → NPT 500 ps → 50 ns
        │
 08  MM/GBSA                    solvant retiré, GBn2 + écrantage 150 mM
                                ΔG par image, ΔΔG + incertitude
        │
 09  analyse                    RMSD, distance de liaison, RMSF, contacts
```

### Paramètres de production

| | Valeur |
|---|---|
| Champ de force protéine | **ff14SB**, pas ff19SB (voir §4.1) |
| Ligand | GAFF2 + charges AM1-BCC (antechamber) |
| Solvant | TIP3P explicite, boîte dodécaédrique, marge 10 Å |
| Force ionique | 150 mM NaCl, neutralisé |
| Électrostatique | PME, coupure 1.0 nm |
| Contraintes | HBonds, eau rigide |
| Pas / masses | **4 fs**, HMR à 4 uma (voir §4.7) |
| Thermostat / barostat | Langevin-Middle 300 K, 1 ps⁻¹ / Monte-Carlo 1 bar |
| Production | 50 ns par bras, image toutes les 10 ps |
| Système NOD1 | 54 738 atomes solvatés (4758 + 42 secs) |
| Débit mesuré (NPT) | **263 ns/jour** (médiane de 3, 256-317) |
| Durée, 50 ns × 2 bras | **≈ 9.1 h** (fourchette 7.6-9.4 h) |

### À quoi sert chaque étape, et ce qu'elle coûte

Toutes les durées ci-dessous sont **mesurées** sur le couple NOD1 / iE-DAP
(54 738 atomes solvatés, RTX 4070), sauf mention contraire. La variance machine
d'un run à l'autre est de l'ordre de 20 %.

| Étape | Ce qu'elle produit, et pourquoi | Durée |
|---|---|---|
| **05** préparation et docking | Récupère le modèle AlphaFold, découpe le domaine LRR, protone à pH 7.4. Construit iE-DAP en 3D depuis son SMILES, en vérifiant formule, charge nette et codes CIP, puis son miroir par réflexion globale. Docke les deux bras avec des paramètres strictement identiques. Sert à obtenir un point de départ crédible : le docking ne tranche rien, il place. | **≈ 2 min**, dont 19 s de docking par bras |
| **07** dynamique moléculaire | Reporte la pose sur la molécule de référence vérifiée, réfléchit et resuperpose pour le bras miroir, applique les charges AM1-BCC partagées, solvate, minimise, chauffe puis équilibre en pression, et produit la trajectoire. C'est l'étape qui laisse le complexe se réarranger : sans elle, on ne compare que deux placements rigides. | **≈ 4.7 h par bras**, soit **9.3 h** pour les deux |
| **08** MM/GBSA | Retire le solvant explicite, réévalue complexe, récepteur et ligand en solvant implicite GBn2, et produit le ΔG de chaque bras puis le ΔΔG avec son incertitude. C'est ici que naît l'observable du projet. | **≈ 2 min** pour les deux bras (300 images chacun) |
| **09** analyse | RMSD du squelette et du ligand, distance minimale de liaison, RMSF, et carte de contacts comparée entre bras. Répond à deux questions que le ΔΔG seul ne pose pas : le ligand est-il resté lié, et **quels contacts** ont été perdus, conservés ou gagnés. | **≈ 10 min** pour les deux bras (5000 images chacun) |

**Total d'un couple, de la préparation au résultat : environ 9 h 30**, dont
98 % en dynamique moléculaire. Toute optimisation qui ne porte pas sur l'étape
07 est cosmétique.

Détail de l'étape 07, par bras :

| Sous-étape | Durée | Rôle |
|---|---|---|
| Solvatation et paramétrisation | ≈ 20 s | boîte dodécaédrique, ions à 150 mM ; charges relues du cache |
| Minimisation | 1 à 2 min | résout les contacts trop courts hérités du docking |
| Chauffage NVT, 100 ps | ≈ 30 s | montée en température par paliers, pour ne pas déformer le complexe |
| Équilibration NPT, 500 ps | ≈ 3 min | ajustement de la densité sous barostat, journalisé pour prouver sa convergence |
| **Production, 50 ns** | **≈ 4.6 h** | la trajectoire analysée |


---

## 4. Décisions méthodologiques, chacune adossée à une mesure

### 4.1 ff14SB et non ff19SB

ff19SB introduit des cartes **CMAP** ajustées sur des L-acides aminés. Un terme
paramétré sur une seule main **brise la symétrie chirale du potentiel** : le
miroir d'un système n'a plus la même énergie que l'original, alors que la
physique l'exige.

Mesure sur l'ubiquitine, E(miroir) − E(naturel) :

| Champ de force | Artefact |
|---|---|
| ff19SB | **−356.9 kJ/mol**, entièrement imputable au `CMAPTorsionForce` |
| ff14SB | **0.000000 kJ/mol** |

Verrouillé par test. **Nuance honnête** (audit 2.1) : sur les couples
petite-molécule, seul le *ligand* est reflété et il est en GAFF2 (achiral, sans
CMAP) ; le récepteur reste L dans les deux bras, donc le CMAP se compenserait
largement. Le choix est décisif pour **TLR5** (ligand protéique reflété), pas
pour NOD1. Le survendre serait attaquable.

### 4.2 Miroir par réflexion globale des coordonnées

L'énantiomère complet s'obtient par une **opération de symétrie impropre**
(det = −1) appliquée à toutes les coordonnées, et non par inversion centre par
centre. Une énumération atome par atome peut en oublier un ; une réflexion ne
le peut pas, par construction.

Vérification systématique : `is_true_enantiomer()` exige que **chaque** centre
stéréogène ait basculé R↔S. Un seul centre inchangé fait échouer la chaîne.

La convention L/D n'est pas supposée mais **calibrée sur structure réelle** :
mesurée sur le PDB 3V47, le volume chiral signé autour du Cα est positif dans
**1368 cas sur 1368** (100.0000 %).

Cas particulier iE-DAP : la molécule porte déjà des centres non-L (glutamate D,
méso-DAP R + S). Son miroir n'est donc pas « la version D » mais l'énantiomère
complet des 3 centres. Une inversion manuelle y serait particulièrement
casse-gueule.

### 4.3 Protomère à pH 7.4

Le SMILES canonique de PubChem est la forme **neutre**. À pH 7.4, iE-DAP porte
3 carboxylates (pKa 2–4) et 2 ammoniums (pKa 9–10), l'azote amide restant
neutre : **net −1**, formule `C12H20N3O7⁻`.

L'erreur ne s'annule pas entre bras : NOD1 lit le carboxylate du DAP dans une
poche basique, les charges AM1-BCC d'un COOH n'ont rien à voir avec celles d'un
COO⁻, et l'observable du projet *est* la complémentarité électrostatique
différentielle. Verrouillé par une assertion formule + charge dans le script 05.

### 4.4 La sortie de smina n'est pas la molécule fournie

smina passe par **PDBQT**, format AutoDock qui ne conserve que les hydrogènes
polaires et ne porte aucune charge formelle. Mesure faite sur la forme neutre
alors en usage (43 atomes) : **30 atomes rendus sur 43**, atomes lourds
réordonnés (CCCCOON… → CCOONCO…). Le protomère corrigé en compte 42 (§4.3) ; la
perte est de même nature.

Paramétrer cette sortie telle quelle laisse openff/RDKit recompléter les
hydrogènes par les règles de valence, ce qui **reconstruit la forme neutre en
silence**, sans erreur ni avertissement.

Correction : on ne retient de la pose que ce qu'elle apporte réellement, la
**position des atomes lourds**, reportée sur la molécule de référence vérifiée
(`transfer_pose`), avec appariement sur le squelette constitutionnel
(insensible à la protonation). Stéréochimie relue depuis la géométrie de la
pose, pour qu'une inversion introduite par le docking reste détectable.

### 4.5 Charges AM1-BCC calculées une fois, partagées

Mécanisme vérifié dans le code installé, pas supposé : openmmforcefields appelle
`assign_partial_charges("am1bcc")` **sans** `use_conformers`, et openff
(`ambertools_wrapper.py:175`, `rec_confs=1`) **jette le conformère fourni** pour
en générer un neuf par ETKDG, sans graine. Chaque appel tire donc un conformère
différent.

| Comparaison | Écart max sur les charges |
|---|---|
| Deux conformères de la même molécule | **0.085 e** |
| Naturel vs miroir (réflexion exacte) | **< 1e-4 e** |

Facteur ≈ 1000. La réflexion ne change rien, l'hamiltonien AM1 ne dépendant que
de distances interatomiques, invariantes par réflexion ; le conformère, lui,
change tout.
Partager les charges est donc **exact**, pas approché, et c'est mesuré.

*(NAGL, réseau de neurones insensible au conformère, est essayé en premier par le
registre openff mais **refuse** le nom `am1bcc` ; c'est bien antechamber qui
calcule. L'avertissement NAGL au passage est sans conséquence.)*

### 4.6 Le bras miroir part de la pose naturelle réfléchie

Le score de docking est **achiral en pratique** (mesuré : −5.9 pour les deux
bras) et le champ de force est chiral-symétrique. Seule la **géométrie
échantillonnée** peut donc discriminer les bras : tout dépend de la pose
initiale. Docker les deux bras indépendamment revient à tirer la pose du miroir
à pile ou face parmi des poses quasi dégénérées : le ΔΔG mesurerait alors ce
tirage autant que la chiralité.

Le miroir part donc du reflet du naturel, **resuperposé** par la rotation
**propre** (Kabsch avec correction de déterminant, det = +1) qui minimise le
RMSD des atomes lourds.

⚠ Le piège critique désamorcé ici : une rotation **impropre** défait exactement
la réflexion et rend le ligand naturel. Le bras miroir simulerait le naturel,
ΔΔG = 0, et on conclurait « aucune différence », sans qu'aucune erreur ne soit
levée. D'où l'assertion `is_true_enantiomer()` après resuperposition.

**Limite assumée :** ce point de départ favorise le mode de liaison du naturel.
On mesure donc la perte de reconnaissance *dans ce mode de liaison*, pas
l'impossibilité de toute liaison. Répondre à « quel est le meilleur mode de
liaison du miroir ? » demanderait plusieurs poses de départ et des répliques.

### 4.7 HMR et pas de 4 fs

Chaque hydrogène est porté à 4 uma, la masse étant retirée à l'atome lourd
porteur. Mesure sur le complexe NOD1 solvaté (54 738 atomes) :

| Configuration | Débit (NVT) |
|---|---|
| 2 fs, sans HMR | 204.8 ns/jour |
| 4 fs, avec HMR | **407.2 ns/jour** |

Soit **×1.99**. Ces deux mesures sont prises sans barostat, ce qui isole
l'effet du pas d'intégration. Le débit de PRODUCTION, barostat actif, est plus
bas : **263 ns/jour** (médiane de 3 mesures, 256 à 317). C'est ce dernier
chiffre qui fixe les durées annoncées, soit environ 9.1 h pour 50 ns × 2 bras
au lieu de quelque 18 h sans HMR.

Légitimité : **seules les masses changent, jamais le potentiel.** La
thermodynamique d'équilibre est inchangée (les masses n'entrent pas dans la
distribution de Boltzmann des positions) et le ΔΔG est une quantité d'équilibre.
La symétrie chirale est intacte, une masse étant un scalaire invariant par
réflexion, et le test de symétrie du champ de force passe inchangé. Masse totale
conservée (vérifiée par test), donc pas de densité faussée sous barostat.

Ce qui change réellement : les modes de vibration rapides ralentissent, donc les
temps de corrélation exprimés en *pas* diffèrent. Sans effet sur une moyenne
d'équilibre. Désactivable : `MDConfig(timestep=2.0, hydrogen_mass=0.0)`.

### 4.8 Boîte dodécaédrique

−29.6 % d'atomes vs un cube à marge égale (+33 % de débit), mesuré. Sans
contrepartie sur un soluté globulaire ou allongé comme un fer à cheval LRR.

---

## 5. MM/GBSA

Mono-trajectoire : les trois états (complexe, récepteur, ligand) sont extraits de
la **même** trajectoire, le solvant explicite est retiré, et l'énergie est
évaluée en Generalized Born **GBn2** avec écrantage ionique 150 mM cohérent avec
la production.

```
ΔG_liaison ≈ ⟨E_complexe⟩ − ⟨E_récepteur⟩ − ⟨E_ligand⟩
```

Variante adaptée à la comparaison de deux **énantiomères** du même ligand : les
termes internes et entropiques du ligand se compensent largement dans le ΔΔG.
L'entropie de mode normal n'est donc pas calculée : choix assumé et documenté,
à revisiter si le ΔΔG s'avère petit (audit 2.5).

**Incertitude** : erreur standard par moyennes de blocs (les images sont
autocorrélées ; l'écart-type inter-images surestime l'information), propagée en
quadrature sur le ΔΔG. C'est une erreur **intra-run**, borne inférieure de
l'incertitude vraie, jamais à présenter comme l'incertitude finale.

### Vérifications du système implicite

Le solvant implicite se construit par le **fichier** `implicit/gbn2.xml` chargé
dans le champ de force, et non par un kwarg `implicitSolvent` (celui-ci
appartient à la route `AmberPrmtopFile`). `nonbondedMethod` vit dans
`nonperiodic_forcefield_kwargs`, les trois systèmes étant apériodiques.

La topologie est lue par **OpenMM**, avec les mêmes liaisons que celles ayant
servi à la production, et mdtraj n'est sollicité que pour les **coordonnées** : sa
reconstruction de liaisons ne satisfait pas les gabarits ff14SB aux résidus
terminaux, et sa sélection `"protein"` n'isole pas un ligand organique. Le
ligand est donc identifié par élimination des acides aminés et du solvant, puis
**vérifié** contre la molécule attendue (nombre d'atomes *et* suite des
éléments), et la partition récepteur/ligand est assertée non dégénérée.

Contrôles passés :

| Vérification | Résultat |
|---|---|
| Couverture GB, complexe | 4800 / 4800 particules |
| Couverture GB, récepteur | 4758 / 4758 |
| Couverture GB, ligand (GAFF) | 42 / 42, pas d'évaluation dans le vide |
| Écrantage `implicitSolventKappa` effectif | −0.188 kcal/mol à 150 mM vs sans sel |
| Découpage ParmEd, ordre atomique | préservé (récepteur et ligand) |
| Chaîne 07 → 08 → 09 sur trajectoire réelle | complète, 4758 + 42 pour les deux bras |

### Performance

| | Coût par image | 300 images × 2 bras |
|---|---|---|
| CPU | 1260.6 ms | 12.6 min |
| CUDA, double précision | **107.8 ms** | **1.1 min** |

×11.7, et ×23 par rapport au code d'origine qui évaluait chaque énergie deux
fois. Double précision imposée sur GPU : accord CPU/CUDA mesuré à **0.001
kcal/mol** sur le ΔG, très en dessous du signal recherché.

---

## 6. Invariant central, vérifié terme à terme

Les deux bras doivent différer **uniquement par la géométrie**. Vérifié en
comparant les deux systèmes OpenMM du ligand paramètre par paramètre :

| Terme | Naturel vs miroir |
|---|---|
| Charges + Lennard-Jones | identiques |
| Liaisons harmoniques | identiques |
| Angles harmoniques | identiques |
| Torsions périodiques | identiques |
| Paramètres GB | identiques |
| Masses | identiques |
| **Géométrie** | **RMSD 1.388 Å** |

Et la réflexion est une isométrie exacte : écart maximal des distances internes
**9.08e-05 Å**, centroïdes des atomes lourds distants de **0.000005 Å** (même
site).

---

## 7. Limites connues : rien de publiable en l'état

| # | Limite | Classe |
|---|---|---|
| 1.1 | **Aucun contrôle de calibration** : la chaîne n'a jamais reproduit un ΔG expérimental connu | bloquant |
| 1.2 | **Une seule réplique, pas d'expérience nulle** nat-vs-nat (plancher de bruit) | bloquant |
| 2.3 | Convergence jugée sur le RMSD, jamais sur la **série d'énergie** | à faire |
| 2.5 | Entropie négligée, alors que le ΔΔG attendu est petit | à documenter |
| 2.8 | **Le site de liaison de NOD1 est un pari** : aucun site expérimental, docking à l'aveugle sur la face concave, résidus de contact extraits *a posteriori* de la MD (circulaire) | caveat majeur |
| 4.6 | Le ΔΔG mesure la perte dans le **mode de liaison du naturel**, pas l'impossibilité de toute liaison | par construction |

Le signal attendu est de l'ordre du bruit d'un MM/GBSA mono-trajectoire. Sans
1.1 et 1.2, un chiffre resterait ininterprétable en absolu.

---

## 8. Ce que la pipeline permet

**Un ΔΔG avec son incertitude**, par couple récepteur / motif, et non une
appréciation qualitative. L'erreur est propagée en quadrature depuis les deux
bras, par moyennes de blocs.

**Une lecture mécanistique, pas seulement un scalaire.** `contact_map_diff()`
classe chaque résidu du site en *conservé*, *perdu* ou **apparu**. Les contacts
apparus au bras miroir sont d'un intérêt particulier : ils signalent un point
d'accroche résiduel que l'hypothèse de départ n'anticipe pas, donc une prise
exploitable pour concevoir un détecteur.

**Une comparaison inter-récepteurs.** La chaîne petite-molécule est réutilisée
telle quelle sur NOD1, TLR4-MD2 et TLR1/2. Les ΔΔG deviennent alors comparables
entre eux : une hiérarchie de vulnérabilité, qui indique par quel récepteur une
contre-mesure aurait le plus d'effet.

**Une garantie de contraste propre.** L'invariant du §6 assure que les deux bras
ne diffèrent que par la chiralité : le ΔΔG ne peut pas absorber une différence
de paramétrisation. C'est ce qui distingue une comparaison d'énantiomères d'une
comparaison de deux systèmes vaguement similaires.

**Un budget de calcul compatible avec la rigueur.** Un couple coûte ≈ 9.1 h de
GPU en production et ≈ 1 min de scoring. Ce coût rend abordables les contrôles
qui manquent encore (répliques et expérience nulle nat-vs-nat), là où un
protocole trois fois plus lent les rendrait hors de portée. La performance n'est
pas ici un confort : elle est la condition de la validité statistique.

**Une anticipation sans création de risque.** Aucune bactérie miroir n'existe et
il ne faut pas en fabriquer. La simulation est le seul moyen d'estimer le
déficit de détection sans produire l'objet dont on redoute les effets.

**Un point d'entrée pour l'expérience.** Le classement des couples et la carte
des contacts indiquent quels systèmes justifieraient une validation en
laboratoire sur peptides D synthétiques, travail lent et coûteux qu'il vaut
mieux cibler.

---

## 9. Installation et exécution

Le calcul tourne dans **WSL2 / Ubuntu**, pas sous Windows natif : AmberTools
(`antechamber`) n'a pas de build Windows et le calcul de charges en dépend. GPU
NVIDIA via passthrough CUDA.

```bash
micromamba env create -f environment.yml
bash scripts/00_check_install.sh      # plateforme CUDA + AmberTools/smina/lightdock
```

Données lourdes dans `~/miroir-data` (filesystem natif WSL ; `/mnt/c` passe par
le relais 9p, 2 à 3× plus lent). Surchargeable par `MIROIR_DATA`.

```bash
# validation technique
bash scripts/run.sh scripts/02_validate_forcefield.py   # symétrie chirale du FF
bash scripts/run.sh scripts/03_benchmark_gpu.py         # débit MD réel

# pilote NOD1 / iE-DAP
bash scripts/run.sh scripts/05_dock_nod1.py             # cible, ligands, docking
bash scripts/launch.sh 50                               # production détachée (~9 h)
bash scripts/status.sh                                  # suivi
bash scripts/run.sh scripts/08_mmgbsa_nod1.py           # ΔΔG
bash scripts/run.sh scripts/09_analyse_nod1.py          # figures + contacts

# tests
bash scripts/test.sh                                    # tout
bash scripts/test.sh -m "not slow"                      # rapide
```

Production lancée en `setsid` : elle survit à la fermeture du terminal. Un run
antérieur était mort une minute après la fermeture de l'onglet, `nohup` ne
protégeant que de SIGHUP.

### Robustesse d'une production longue

`launch.sh` passe d'abord par un **contrôle pré-vol** (`scripts/preflight.py`,
lançable seul) : production concurrente, plateforme CUDA, espace disque,
fichiers de préparation, protomère des deux ligands, cohérence du cache de
charges, poses de départ constructibles et miroir énantiomère exact. Trente
secondes qui évitent de découvrir un problème après six heures de GPU.

Une production interrompue **reprend sur son checkpoint** au lieu de repartir de
zéro : le script détecte `prod_<bras>.chk`, relit le système solvaté, tronque la
trajectoire au dernier état sauvegardé et prolonge les sorties. Le système est
relu et jamais re-solvaté, `addSolvent` plaçant les ions sans graine fixée : une
re-solvatation donnerait le même nombre d'atomes dans un ordre différent, que
`loadCheckpoint` accepterait sans broncher (voir `docs/00`, point 15).

Un `InstabilityGuard` interrompt la production dès que l'énergie cesse d'être
finie, plutôt que d'écrire des NaN pendant des heures en occupant le GPU.

---

## 10. Carte du code

```
src/miroir/
  chirality.py   réflexion des coordonnées + classification L/D calibrée
  structure.py   E/S PDB, audit de chiralité, invariants géométriques
  prepare.py     nettoyage, protonation pH 7.4, sélection de chaînes
  validate.py    test de symétrie chirale d'un champ de force
  inventory.py   audit organisme / résolution / intégrité ligand d'un PDB
  alphafold.py   récupération AlphaFold DB + découpe de domaine
  ligand.py      build 3D, miroir, vérif CIP, report de pose, protomère
  pocket.py      définition de la boîte de docking
  docking.py     wrapper smina (paramètres identiques imposés aux 2 bras)
  complex_md.py  assemblage protéine + ligand GAFF2, charges partagées
  md.py          protocole MD complet sur GPU (min → NVT → NPT → prod)
  mmgbsa.py      ΔG de liaison MM/GBSA (GBn2)  ← observable centrale
  analysis.py    RMSD, distance de liaison, RMSF, contacts (PBC triclinique)
  config.py      registre des couples + décisions encodées
```

Le module de chiralité est couvert de façon disproportionnée : l'erreur
d'énantiomère est classée à impact **critique**, car elle invaliderait
silencieusement tout l'aval.

---

## 11. Traçabilité

- [docs/00-choix-methodologiques.md](docs/00-choix-methodologiques.md) : chaque
  décision, sa justification chiffrée, le test qui la verrouille (14 points).
- [docs/01-audit-relecture-2026-07-15.md](docs/01-audit-relecture-2026-07-15.md),
  relecture adversariale dont le statut est vérifié **dans le code** pour chaque point
  (un point du rapport de synthèse s'est révélé faux positif à la vérification :
  d'où la règle « vérifier dans le code, pas se fier à la synthèse »).
- [HANDOFF.md](HANDOFF.md) : état complet pour reprise ou audit sans contexte.
