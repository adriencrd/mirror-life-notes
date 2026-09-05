# Choix méthodologiques et écarts au cahier des charges

Ce document recense les points où l'implémentation **s'écarte** du CDC v1.0, avec
la mesure ou le fait qui justifie chaque écart. Il est destiné à être relu avant
rédaction du preprint : plusieurs de ces points sont eux-mêmes des résultats
méthodologiques publiables.

---

## 1. ff19SB est inutilisable sur les systèmes miroir — utiliser ff14SB

**Le CDC dit** (section 7, Phase 3) : « Champ de force : AMBER ff19SB (protéines) ».

**Problème.** Toute l'étude repose sur une comparaison d'énergies entre un système
naturel et son miroir. Cela n'a de sens que si le potentiel traite les deux
chiralités à égalité. Or ff19SB = ff14SB **+ cartes de correction CMAP**, ajustées
sur des L-acides aminés. Un terme CMAP encode une surface E(φ,ψ) ; sous réflexion,
(φ,ψ) → (−φ,−ψ), et la grille de ff19SB n'a pas cette symétrie. Elle pénalise donc
un D-résidu pour le seul fait d'être D.

Tous les autres termes d'Amber sont, eux, rigoureusement symétriques :

| Terme | Comportement sous réflexion | Symétrique ? |
|---|---|---|
| Liaisons, angles | fonctions de distances/angles | oui |
| Torsions propres et impropres | `k(1 + cos(nφ − δ))`, φ → −φ | oui, **car** toutes les phases δ d'Amber valent 0° ou 180° |
| Lennard-Jones, Coulomb | fonctions de distances | oui |
| **CMAP** (ff19SB seulement) | grille E(φ,ψ) → E(−φ,−ψ) | **non** |

**Mesure** (`scripts/02_validate_forcefield.py`, ubiquitine 1UBQ protonée, énergie
d'un système unique évaluée sur coordonnées naturelles puis réfléchies) :

| Champ de force | ΔE = E(miroir) − E(naturel) | Verdict |
|---|---|---|
| ff19SB | **−356.878 kJ/mol** | brise la symétrie |
| ff14SB | **0.000000 kJ/mol** | symétrique |

L'écart de ff19SB est imputable **au seul `CMAPTorsionForce`** : le terme fautif
vaut −356.877957 kJ/mol sur un ΔE total de −356.877957 kJ/mol. Tous les autres
termes sont identiques jusqu'au dernier chiffre significatif. Le diagnostic est
donc établi, pas seulement constaté.

**Ordre de grandeur.** ≈ 1.1 kcal/mol par résidu de biais systématique, sur une
protéine de 76 résidus. Le ΔΔG recherché se chiffre en quelques kcal/mol.

**Nuance importante — ne pas surinterpréter ce chiffre.** Ce ΔE est un artefact
*intramoléculaire*. En MM/PBSA **mono-trajectoire** (le cas usuel : récepteur et
ligand sont extraits des images du complexe, sans relaxation), l'énergie interne
du ligand apparaît à l'identique dans `G_complexe` et dans `G_ligand` : elle se
soustrait **exactement**. Affirmer « ff19SB fausse le ΔΔG de 357 kJ/mol » serait
donc faux, et un relecteur le relèverait aussitôt.

L'argument décisif est ailleurs : **ff19SB biaise l'échantillonnage lui-même.**
Le CMAP n'est pas qu'un terme d'énergie, il exerce une force. Pendant 100 ns de
production, il pousse continûment le squelette du D-ligand vers des
conformations de type L — il déforme activement l'objet qu'on prétend mesurer.
Aucune soustraction ne rattrape une trajectoire échantillonnée sur le mauvais
paysage énergétique. (En approche 3-trajectoires, où les partenaires sont
relaxés séparément, l'artefact ne se compense même plus.)

**Décision.** `PROTEIN_FF = "amber14/protein.ff14SB.xml"` (`src/miroir/md.py`).
ff14SB est de la même famille, largement validé, et sans CMAP. C'est déjà le choix
usuel de la littérature sur les D-peptides.

Contrepartie assumée : ff19SB est plus exact que ff14SB **sur les protéines
naturelles**, et on y renonce. C'est le bon arbitrage ici, car l'étude est
**comparative** : les deux bras doivent tourner sous le même potentiel, et ce
potentiel doit être équitable envers les deux chiralités. Un ΔΔG issu de deux
champs de force différents ne voudrait rien dire ; un ΔΔG issu d'un champ de
force qui favorise un bras non plus. La symétrie prime sur l'exactitude absolue
d'un seul bras.

**Non-régression.** `tests/test_validate.py` échoue si quelqu'un revient à ff19SB.

**Piste si la précision de ff19SB devient nécessaire** : symétriser les grilles
CMAP en construisant, pour chaque type de résidu, une carte D telle que
E_D(φ,ψ) = E_L(−φ,−ψ), et l'appliquer aux résidus miroir. C'est faisable (ffxml
custom) et ce serait une contribution méthodologique à part entière — mais hors
du chemin critique du premier résultat.

---

## 2. AutoDock Vina / smina ne peuvent pas traiter le couple TLR5 / flagelline

**Le CDC dit** (section 6) : « Démarrage impératif par TLR5/flagelline [...]
ligand protéique donc inversion triviale en coordonnées ». Et (section 7,
Phase 2) : « Outil : AutoDock Vina ou smina ».

**Le raisonnement est juste sur l'inversion, mais ces deux prémisses sont
incompatibles.** L'inversion est effectivement triviale (c'est le cas le plus
simple du projet — une réflexion globale, sans aucune chimie à décider). Mais
Vina et smina sont des dockers **petite molécule** : ils traitent un ligand
flexible de quelques dizaines d'atomes lourds via ses torsions rotatives. Le
ligand ici est la flagelline — **295 résidus**. Vina ne convergera pas ; l'espace
de torsion est hors de son domaine.

Le docking protéine-protéine est un problème algorithmiquement différent (corps
rigides + rééchantillonnage de surface, fonctions de score dédiées).

**Décision.** Deux voies de docking selon le type de ligand
(`ligand_type` dans `src/miroir/config.py`) :

| Couple | Ligand | Taille | Outil |
|---|---|---|---|
| 1. TLR5 / flagelline | protéine | 295 rés. | **LightDock** (protéine-protéine, scriptable en Python) |
| 2. TLR1-TLR2 / Pam3CSK4 | lipopeptide | ~1.5 kDa | smina |
| 3. NOD1 / iE-DAP | petite molécule | ~0.4 kDa | smina |
| 4. TLR4-MD2 / lipide A | glycolipide | ~1.8 kDa | smina |
| 5. Dectine-1 / β-glucane | oligosaccharide | variable | smina |

**Conséquence sur le planning.** Le CDC prévoit « Phase 2 — Docking TLR5 :
1 semaine ». Le docking protéine-protéine est plus coûteux et plus délicat à
valider que le docking petite molécule. Cette estimation est à revoir.

**Conséquence sur l'ordre.** L'argument « TLR5 d'abord car c'est le plus simple »
ne tient que pour l'étape d'inversion. Sur le docking, TLR5 est le **plus
difficile** des cinq, et le seul à demander un outil différent. Deux options :

- **(a)** garder TLR5 en premier — meilleure structure, meilleur alignement avec
  la littérature, résultat le plus attendu ; accepter le surcoût de docking.
- **(b)** commencer par NOD1/iE-DAP ou TLR2/Pam3CSK4 — la voie smina valide la
  pipeline complète (docking → MD → MM/PBSA) plus vite et donne un premier ΔΔG
  plus tôt ; TLR5 ensuite, en terrain connu.

À arbitrer. L'argument pour (b) : la Phase 0 exige de valider la pipeline avant
les motifs miroir, et il est plus sain de la valider sur la voie qu'on
réutilisera 4 fois sur 5.

---

## 3. 3V47 n'est pas du TLR5 humain

**Le CDC dit** (section 3) : « reconnaissance des motifs microbiens miroirs par
les récepteurs de l'immunité innée **humaine** ».

**Fait.** L'en-tête de 3V47 :

```
TITLE   CRYSTAL STRUCTURE OF THE N-TERMINAL FRAGMENT OF ZEBRAFISH TLR5 IN
TITLE 2 COMPLEX WITH SALMONELLA FLAGELLIN
COMPND  MOLECULE: TOLL-LIKE RECEPTOR 5B AND VARIABLE LYMPHOCYTE RECEPTOR B.61
COMPND  CHIMERIC PROTEIN;  ENGINEERED: YES;  MUTATION: YES
```

C'est du **TLR5b de poisson-zèbre**, fusionné à un **VLR de myxine** (échafaudage
de cristallisation), en complexe avec la flagelline de *Salmonella*. Résolution
2.47 Å — le CDC exige < 2.5 Å (section 12) : ça passe, de justesse.

Ce n'est pas une erreur du CDC : **aucune co-cristallographie TLR5
humain/flagelline n'existe**. 3V47 est la référence universelle du domaine faute
de mieux.

**Conséquences.**
1. Le titre et les conclusions ne peuvent pas dire « humain » sans réserve. Un
   relecteur le verra immédiatement.
2. La partie VLR est un artefact de cristallisation, à retirer avant MD (elle
   n'a aucun rôle biologique et alourdit le système).
3. Option à considérer : modèle par homologie / AlphaFold du TLR5 humain aligné
   sur 3V47. Le CDC prévoit déjà ColabFold « si PDB absent » (section 8). Cela
   déplace le risque vers la qualité du modèle — à documenter honnêtement.

**Statut : non tranché.** Décision à prendre avant de lancer la MD de production.

---

## 4. Windows seul ne suffit pas — le calcul tourne dans WSL2

**Le CDC dit** (section 8) : tout en local.

**Fait.** AmberTools (MMPBSA.py, antechamber, tleap) n'a **pas de build Windows**
sur conda-forge : Linux et macOS uniquement. Or la Phase 4 (calcul d'énergie
libre) en dépend, et c'est le résultat principal du projet.

**Décision.** Tout le calcul tourne dans **WSL2 / Ubuntu 22.04**, avec passthrough
CUDA vers la RTX 4070 (vérifié : OpenMM détecte la plateforme CUDA et calcule les
forces dans les tolérances). Le code reste sur le disque Windows, éditable
normalement.

Les données lourdes vont sur le filesystem **natif WSL** (`~/miroir-data`) et non
sur `/mnt/c` : les E/S traversant le relais 9p sont 2 à 3× plus lentes, ce qui est
sensible sur une trajectoire de 100 ns. Surchargeable via `MIROIR_DATA`.

Le CDC mentionne aussi un « Mac M5 » pour RDKit (section 8) : sans objet, tout est
sur la même machine.

---

## 5. Le miroir s'applique au ligand seul, par réflexion globale

**Choix d'implémentation**, cohérent avec le CDC mais qui mérite d'être explicite.

- **Réflexion globale plutôt qu'inversion centre par centre.** Le CDC (Phase 1.4)
  propose d'« inverser les centres stéréogènes concernés via RDKit ». Pour un
  ligand peptidique, on applique à la place une réflexion de **toutes** les
  coordonnées (négation d'un axe). C'est plus sûr : une réflexion globale inverse
  *tous* les centres d'un coup — y compris ceux des chaînes latérales (Ile, Thr)
  et l'hélicité des structures secondaires — sans qu'on puisse en oublier un. Le
  CDC classe justement l'erreur d'énantiomère comme risque d'impact **critique**
  (section 12). RDKit reste nécessaire pour les ligands non peptidiques (couples
  2 à 5), où la molécule doit être reconstruite atome par atome.

- **Réflexion et non inversion centrale.** Les deux sont des opérations impropres
  et donnent le même énantiomère, à une rotation propre près. La pose étant de
  toute façon redéterminée par docking, le choix est sans conséquence.

- **Seul le ligand est réfléchi ; le récepteur reste naturel.** C'est la question
  posée. Réfléchir les *deux* donnerait, par symétrie exacte du hamiltonien, très
  exactement l'énergie du complexe naturel — un résultat vrai mais vide.

- **Garde-fou automatique.** `structure.audit_chirality` classe chaque carbone α
  par le signe du produit mixte autour de CA. La convention L a été **calibrée sur
  3V47** (1368 centres, 100.0000 % homogènes), pas supposée. Vérifié sur le couple
  préparé : complexe naturel 100 % L, ligand miroir 100 % D, récepteur 100 % L.

---

## 6. Coût GPU réel — le planning MD tient, mais sans marge

**Le CDC dit** (section 9) : « Phase 3 — MD TLR5 : 3 semaines », et exige (§11)
« MD ≥ 100 ns par système ».

**Mesure** (`scripts/03_benchmark_gpu.py`, RTX 4070, OpenMM/CUDA, précision mixte,
PME, HBonds, pas de 2 fs) :

| Système | Atomes solvatés | Débit | 100 ns |
|---|---|---|---|
| Ubiquitine (référence) | 17 144 | 216.2 ns/j | 11.1 h |
| TLR5/flagelline, boîte **cubique** | 486 061 | 38.4 ns/j | 62.6 h (2.6 j) |
| TLR5/flagelline, boîte **dodécaédrique** | **342 179** | **51.1 ns/j** | **47.0 h (2.0 j)** |

**Enseignement contre-intuitif : le débit ne s'extrapole pas linéairement.** Un
facteur 28 sur le nombre d'atomes ne coûte qu'un facteur 5.6 sur le débit. Le
petit système est limité par la latence de lancement des noyaux CUDA et
sous-exploite le GPU ; la 4070 ne donne sa pleine mesure qu'au-delà de ~100 k
atoms. **Ne jamais dimensionner ce projet par extrapolation depuis un système
test** — mesurer sur la vraie structure.

**Verdict planning.** 2.0 j par système × 2 (naturel + miroir) ≈ **3.9 jours** de
production. Les 3 semaines du CDC tiennent. Mais la marge est plus courte qu'elle
n'en a l'air :

- la Phase 4 (MM/PBSA sur 10 000 images d'un système de ~500 k atomes) est
  elle-même coûteuse, et le CDC ne lui alloue qu'une semaine ;
- le CDC exige la « convergence RMSD vérifiée » (§11). En pratique cela demande
  des **répliques** (3 × 2 systèmes ≈ 12 jours), ce que le planning ne prévoit
  pas ;
- une production de 2 jours est longue devant un plantage : d'où les checkpoints
  dans `md.py`.

**Optimisation appliquée et mesurée.** Boîte **dodécaédrique** au lieu de cubique
(`MDConfig.box_shape`, défaut). Le coût d'une MD est dominé par le nombre de
molécules d'eau, et un dodécaèdre rhombique pave l'espace avec ~29 % de volume en
moins à marge égale. Le fer à cheval LRR de TLR5 étant très allongé, une boîte
cubique le noie dans du solvant inutile.

Gain mesuré, pas supposé : **−29.6 % d'atomes** (486 061 → 342 179), **+33 % de
débit** (38.4 → 51.1 ns/j), soit **15.6 h économisées par système** et 1.3 jour
sur le seul couple TLR5. Aucune contrepartie : la marge de solvatation autour du
soluté reste de 10 Å dans toutes les directions, conformément au CDC.

---

## 7. NOD1/iE-DAP : cible AlphaFold, et le docking ne discrimine pas la chiralité

**Décisions prises** (couple n°1, premier traité) :

**Cible = modèle AlphaFold du NOD1 humain (Q9Y239), domaine LRR seul.** Vérifié :
*aucune* structure expérimentale du LRR de NOD1 n'existe, dans aucun organisme.
Les seules structures déposées (2DBD, 2NZ7, 4JQW) couvrent le domaine **CARD** —
le mauvais domaine, celui de la signalisation en aval, pas de la reconnaissance
du ligand. Le modèle AlphaFold AF-Q9Y239-F1 (v6) prédit le LRR à **pLDDT 95**
(mesuré sur notre copie), soit la région la mieux prédite de toute la protéine.
C'est un cas où ColabFold/AlphaFold n'est pas un pis-aller mais la seule voie.

**Ligand = iE-DAP, PubChem CID 45480617.** Piège évité : le CID 194426, qui
remonte sous le même nom, a une **structure erronée** (pont cétone au lieu de
l'amide isopeptidique). iE-DAP porte déjà des centres non-L (glutamate D,
méso-DAP 2S/6R) ; son miroir est donc l'**énantiomère complet** — les trois
centres basculent (`{S,R,R}` → `{R,S,S}`, vérifié). Inverser les seuls centres
DAP régénérerait le même méso-DAP : seule la réflexion globale donne le vrai
miroir, ce que `ligand.is_true_enantiomer` contrôle.

**Résultat de docking : ~0, et c'est le résultat attendu.** Docking à l'aveugle
(pas de pose de référence) sur la face concave du LRR, paramètres identiques aux
deux bras (exhaustivité 16, graine 42) :

| Bras | Meilleure pose | Moyenne top-5 |
|---|---|---|
| iE-DAP naturel | −5.9 kcal/mol | −5.6 kcal/mol |
| iE-DAP miroir | −5.9 kcal/mol | −5.3 kcal/mol |

Les meilleures poses coïncident à la décimale affichée ; les distributions
complètes diffèrent, mais de l'ordre du **bruit de docking** (~1 kcal/mol). Ce
n'est pas un bug : la fonction de score de Vina/smina est une somme de termes de
paires **achiraux** (elle n'a aucun terme de chiralité), donc elle sépare mal
deux énantiomères. C'est une démonstration, sur notre propre système, de
*pourquoi* le CDC (section 12) interdit de publier un résultat de docking seul et
impose la MD. Le docking ne sert ici qu'à produire la pose de départ ; le vrai
ΔΔG viendra du MM/PBSA sur trajectoire.

**Système MD validé.** Complexe NOD1-LRR + iE-DAP paramétré (GAFF2 + charges
AM1-BCC via antechamber, 40 s), solvaté à **54 605 atomes** — six fois plus léger
que TLR5. Minimisation, chauffage et intégration GPU confirmés stables. La
production 100 ns par bras sera courte (quelques heures), ce qui fait bien de
NOD1 le bon cas de rodage de la pipeline.

---

## 8. Aucune structure du CDC ne passe tous les critères — audit et remplacements

Le CDC (section 6) donne un tableau de 5 structures PDB comme s'il était acquis.
`scripts/10_inventaire_structures.py` les a auditées contre trois critères :
100 % humain, résolution < 2.5 Å (§12), et **ligand chimiquement exploitable**
(ce troisième critère n'est pas dans le CDC, mais 3FXI montre qu'il est décisif).

**Résultat : aucune des 4 ne passe.**

| PDB | Récepteur humain ? | Résolution | Ligand | Sort |
|---|---|---|---|---|
| 3V47 (TLR5) | ❌ poisson-zèbre | 2.47 Å ✅ | flagelline (protéine) | **écarté** → ColabFold Q9NR61 |
| 3FXI (TLR4) | ✅ | 3.1 Å ❌ | éclaté : 13 résidus, 4 chaînes | **remplacé** → 4G8A |
| 2Z7X (TLR1/2) | ✅ (VLR retirable) | **2.1 Å** ✅ | Pam3CSK4, synthétique | **meilleur du lot** |
| 2CL8 (Dectine-1) | ❌ **souris** | 2.8 Å ❌ | 1 seul glucose | **à revoir entièrement** |

**Distinction critique — deux chimères VLR de myxine, deux verdicts opposés.**
3V47 et 2Z7X portent toutes deux un échafaudage de cristallisation VLR
(*Eptatretus burgeri*). Mais dans 2Z7X, `OTHER_DETAILS` dit `TLR2, UNP RESIDUES
27-508 (HUMAN)` et `TLR1, UNP RESIDUES 25-476 (HUMAN)` : les récepteurs **sont**
humains, le VLR n'est qu'une greffe qu'on découpe. Dans 3V47, c'est le récepteur
lui-même (TLR5b) qui est du poisson-zèbre : rien à sauver. Ne pas confondre
« contient une séquence non humaine » et « n'est pas humain ».

**Dectine-1 (2CL8) est de la souris**, et son « β-glucane » se réduit à un seul
glucose (`BGC`, 34 atomes). Le couple 5 est à reconstruire entièrement, ou à
abandonner (le CDC le donne déjà comme « si temps »).

### 3FXI → 4G8A pour le couple TLR4

Le CDC prescrit 3FXI (Park et al. 2009). Il existe mieux depuis :

| | 3FXI | **4G8A** (Ohto et al. 2012) |
|---|---|---|
| Résolution | 3.1 Å ❌ | **2.4 Å ✅** |
| Ligand | Ra-LPS, cœur sucré complet | **Re-LPS** (lipide A + 1 KDO) |
| Modélisation | 13 résidus, **4 chaînes** ; acyles mêlés aux chaînes protéiques | 10 résidus, **2 chaînes** ; tout le ligand dans la chaîne C |
| Composition | `FTT`×4, `MYR`, `DAO`, `PO4`, `KDO`, `GMH`, `PA1`, `GCS` | `LP4`, `LP5`, `MYR`, `DAO`, `KDO` |

`LP4`/`LP5` sont les deux glucosamines **déjà acylées** (C68N2O23P2, 14 centres
stéréogènes) : 4G8A modélise le lipide A en blocs cohérents au lieu de le
pulvériser. Les liaisons acyles↔sucres, absentes de 3FXI, cessent d'être un
problème pour l'essentiel du ligand.

**Réserve levée par la mesure.** 4G8A porte les polymorphismes D299G/T399I, pas
le TLR4 sauvage. Plutôt que de s'en remettre à « les auteurs disent que c'est
identique », on a mesuré : **D299G est à 15.2 Å du ligand, T399I à 17.2 Å**. Les
résidus TLR4 réellement en contact (< 5 Å) sont 264, 296, 341, 389, 414, 415,
436, 439-441, 463 — ni 299 ni 399. Les SNP sont hors du site : sans effet pour
notre docking, et défendable par un chiffre en relecture.

**Écartés après vérification** : 5IJD (souris + chimère VLR, 2.7 Å) ; 8WTA
(cryo-EM humaine mais 2.9 Å, et son ligand DLAM3 est un **mimétique
synthétique**, pas le MAMP naturel).

**Reste à faire pour TLR4** : les 5 résidus du ligand doivent quand même être
reliés en une molécule valide pour GAFF2 — ce n'est pas un SMILES propre comme
iE-DAP. Deux voies : reconstruire la connectivité LP4-LP5-MYR-DAO-KDO, ou bâtir
le lipide A d'*E. coli* depuis un SMILES de référence et l'aligner sur la pose
cristallographique (bien définie à 2.4 Å).

---

## 9. TLR1/2 (2Z7X) : récepteur prêt, un piège de stéréochimie sur le ligand

Couple retenu comme n°2 effectif (avant TLR4), car c'est le chemin le plus
propre : récepteurs humains, 2.1 Å, ligand synthétique.

**Récepteur préparé et validé.** 2Z7X est un hétérodimère TLR1-TLR2 humain, mais
chaque chaîne porte un VLR de myxine en bout (échafaudage de cristallisation).
Bornes lues dans `OTHER_DETAILS` et vérifiées sur la structure : TLR2 = chaîne A
27-508, TLR1 = chaîne B 25-476. Après découpe du VLR : 1069 → 934 résidus,
protoné, **100 % L**. `scripts/11_prepare_tlr2.py`.

**Site de liaison bien défini.** Pam3CSK4 ponte l'interface TLR1-TLR2 : contacts
étendus sur TLR2 (résidus 266-376) et TLR1 (258-339). La pose cristallographique
(2.1 Å) donne directement la boîte de docking, centrée sur (-11.6, -13.4, 7.1) —
pas de docking à l'aveugle ici, contrairement à NOD1.

**Le piège — le glycérol de Pam3CSK4 est stéréochimiquement indéfini dans
PubChem.** CID 130704, C₈₁H₁₅₆N₁₀O₁₃S. La forme biologique est
N-palmitoyl-S-[**(R)**-2,3-bis(palmitoyloxy)propyl]-(R)-Cys-(S)-Ser-(S)-Lys₄ :
7 centres stéréogènes. Mais le SMILES PubChem laisse le **carbone central du
glycérol sans configuration** (`...CSCC(COC(=O)...)OC(=O)...`, ni `@` ni `@@`),
parce que le produit commercial est un mélange à ce centre.

Conséquences :
- On ne peut pas refléter un centre indéfini — il n'y a rien à inverser.
- `openff.toolkit.Molecule.from_file(..., allow_undefined_stereo=False)`
  (utilisé dans `complex_md.load_ligand`) **refusera** la molécule.

Il faut donc **fixer le glycérol à (R)** avant de construire, puis vérifier par
CIP que les 7 centres sont bien définis et que le miroir les inverse tous les 7.
C'est le même garde-fou que pour iE-DAP (`ligand.is_true_enantiomer`), mais avec
une étape d'assignation en amont. À ne pas bâcler : un centre oublié donnerait un
diastéréoisomère, pas l'énantiomère.

**Stéréochimie résolue et vérifiée.** `scripts/12_build_pam3csk4.py` :
`ligand.build_defined_isomer` détecte le glycérol comme seul centre indéfini
(atome 22), le fixe à (R) par énumération + vérification CIP, et confirme les 7
centres : `{Cys R, glycérol R, Ser S, Lys×4 S}`. Le miroir les inverse tous les
7 (`is_true_enantiomer` OK). Conforme à la forme biologique documentée.

**Le vrai obstacle, mesuré : 86 liaisons rotatives** (MM 1510 g/mol). Vina/smina
sont fiables jusqu'à ~15-20 liaisons rotatives ; à 86 (les 3 chaînes C16), un
docking global produirait du bruit. **On ne dockera pas Pam3CSK4 globalement.**

Stratégie retenue (à implémenter) : la pose cristallographique existe et est
fiable (2.1 Å). On l'utilise comme point de départ plutôt que de la chercher.
- Bras naturel : ancrer le ligand construit sur les positions cristallographiques
  (embedding contraint RDKit sur les atomes lourds de la chaîne C de 2Z7X), puis
  relaxer par MD.
- Bras miroir : le reflet de la pose naturelle irait dans un récepteur miroir,
  pas dans le nôtre. On superpose donc le miroir sur le squelette de la pose
  naturelle (ajustement RMSD des atomes communs) et on laisse la MD résoudre les
  chocs stériques résiduels par minimisation.

C'est plus rigoureux qu'un docking non convergé, et ça exploite la donnée
expérimentale au lieu de la jeter. Le docking smina reste la voie pour les
petits ligands (NOD1 fait, lipide A à venir).

---

## 10. iE-DAP doit être simulé chargé (−1), pas sous la forme neutre de PubChem

**Découvert le 2026-07-15**, en vérifiant le point 2.4 de l'audit après un reboot
qui avait tué la production à 37.6 ns / 50.

**Problème.** Le SMILES canonique de PubChem (CID 45480617) décrit la forme
**neutre** : trois `C(=O)O` et deux amines `N`. C'est ce qui a été construit,
docké et simulé — `iedap_naturel.sdf` faisait 43 atomes sans aucun bloc `M  CHG`.
Or à pH 7.4 :

| Fonction | Nombre | pKa | État à pH 7.4 |
|---|---|---|---|
| Acide carboxylique | 3 | ~2–4 | COO⁻ |
| Amine libre | 2 | ~9–10 | NH₃⁺ |
| Azote amide | 1 | — | neutre |

Soit **C12H20N3O7⁻, net −1** (42 atomes), et non C12H21N3O7 neutre (43 atomes).

**Pourquoi ça ne s'annule pas dans le ΔΔG.** L'objection naturelle est que les
deux bras portent la même erreur, donc qu'elle se compense. Elle ne se compense
pas :

- NOD1 reconnaît le **carboxylate du DAP** dans une poche basique. Neutraliser
  les acides supprime les ponts salins qui *sont* la reconnaissance : le complexe
  échantillonné n'est pas le bon, dans les deux bras.
- Les charges AM1-BCC d'un COOH et d'un COO⁻ n'ont rien à voir. L'observable du
  projet est la **complémentarité électrostatique différentielle** entre
  énantiomères — elle serait calculée sur la mauvaise distribution de charge.
- Le docking lui-même est parti de ce protomère.

**Décision.** `config.py` porte désormais le protomère physiologique, et non le
SMILES PubChem tel quel :

```
C(C[C@@H](C(=O)[O-])[NH3+])C[C@H](C(=O)[O-])NC(=O)CC[C@H](C(=O)[O-])[NH3+]
```

Les tags `@`/`@@` sont inchangés : l'ordre des voisins de chaque centre est
identique, donc la parité l'est aussi. **Vérifié** — les codes CIP restent 2 R +
1 S, comme attendu (les priorités CIP N > COO > CH₂ > H ne dépendent pas de la
protonation du carboxyle).

**Verrouillage.** Trois garde-fous, parce que rien dans la chaîne ne signalait
l'erreur :

1. `05_dock_nod1.py` échoue si la formule ≠ `C12H20N3O7-` ou la charge ≠ −1 ;
2. `ligand.identity_preserved()` compare formule + charge + signature stéréo — le
   PDBQT interne de smina/OpenBabel ne porte **ni stéréo ni charge formelle**, un
   aller-retour peut donc reprotoner un carboxylate sans lever d'erreur ;
3. `07_run_md_nod1.py` l'applique à la pose entrant en MD (custody, audits 1.4 et
   2.4 réunis).

**Coût.** Les 37.6 ns du bras naturel sont jetées. Elles ont servi : elles ont
validé la chaîne (155 ns/day, T et densité stables, checkpoints fonctionnels) et
fourni une trajectoire réelle pour valider `analysis.py` (c'est là qu'on a trouvé
le bug PBC triclinique).

---

## 11. La sortie de smina n'est pas la molécule qu'on lui a donnée

**Mesuré le 2026-07-15**, en cherchant si smina réordonnait les atomes (pour
savoir si les charges AM1-BCC pouvaient être partagées entre bras par indice).

**Le fait.** smina passe par **PDBQT**, le format d'AutoDock. Ce format applique
la convention « hydrogènes polaires seulement » : les H portés par C sont fusionnés
dans leur carbone, et aucune charge formelle n'est stockée. Mesure sur iE-DAP :

| Fichier | Atomes | Composition |
|---|---|---|
| `iedap_naturel.sdf` (référence) | 43 | C12 **H21** N3 O7 |
| `dock_naturel.sdf` (sortie smina) | **30** | C12 **H8** N3 O7 |
| `solv_naturel.pdb` (ce qui a été simulé) | 43 | C12 **H21** N3 O7, **neutre** |

Les 13 hydrogènes manquants ont été **recomplétés par les règles de valence**
d'openff/RDKit. Le protomère n'est donc pas choisi : il est *reconstruit*, et les
règles de valence rendent toujours la forme neutre — un `COO⁻` privé de sa charge
formelle redevient `COOH`. Silencieusement, sans avertissement.

Conséquence : la correction du point 10 **ne suffisait pas**. Un SMILES chargé
serait revenu neutre au passage par smina.

**Deuxième fait.** smina **réordonne les atomes lourds** :

```
référence    CCCCOON CCCOON COCCCCOON
sortie smina CCOONCO CCCNCO OCCCCNCOO
```

Toute correspondance par indice entre les deux bras est donc invalide — ce qui
condamnait l'implémentation naïve du point 1.3 de l'audit (« calculer les charges
une fois, les copier »).

**Décision.** La sortie de smina n'entre plus dans la paramétrisation. On n'en
garde que ce qu'elle apporte réellement — **la position des atomes lourds** — et
on la reporte sur la molécule de référence (`ligand.transfer_pose`), qui est la
seule vérifiée. L'appariement se fait sur le **squelette constitutionnel**
(éléments + connectivité, charges neutralisées, ordres de liaison uniformisés),
donc insensible à la protonation et à la forme de Kekulé. Les hydrogènes sont
replacés par géométrie, puis relaxés par la minimisation qui ouvre la MD.

**Vérifié** sur les fichiers réels : référence chargée + pose smina → `C12H20N3O7⁻`,
charge −1, stéréo (R,R,S) préservée, atomes lourds à **0.000000 Å** des positions
de smina. Le protomère survit, la pose reste celle du docking.

**Limites assumées.**

- Les deux oxygènes d'un carboxylate deviennent équivalents sur le squelette et
  peuvent être échangés à l'appariement. Sans conséquence : positions quasi
  confondues, relaxées à la minimisation. Aucun **centre stéréogène** ne peut être
  échangé — ses quatre substituants sont distincts par définition.
- Les positions d'hydrogène sortent de règles géométriques, pas du docking. Vina
  ne les optimise pas non plus (elles ne sont pas dans le PDBQT) : rien n'est perdu.

**Effet de bord voulu.** Les deux bras héritent de l'ordre atomique de leur
référence, identique entre naturel et miroir (le miroir n'est qu'une copie des
coordonnées). C'est ce qui rend le partage des charges AM1-BCC légitime — voir
point 12.

---

## 12. Les charges AM1-BCC sont calculées une fois et partagées par les deux bras

**Mesuré le 2026-07-15** (point 1.3 de l'audit, confirmé — mais pas par le
mécanisme que l'audit supposait).

**Le fait.** Les types d'atomes GAFF et les incréments BCC dérivent de la
topologie, qui est achirale : sur ce plan les deux bras sont identiques. Mais les
charges **AM1-BCC se calculent depuis un conformère 3D**, et le chemin par défaut
en tire un **au hasard** :

```
template_generators.py:648   assign_partial_charges("am1bcc")        # sans use_conformers
ambertools_wrapper.py:175    if use_conformers is None:
                                 generate_conformers(n_conformers=1)  # rec_confs = 1
```

openff **jette le conformère fourni** et en génère un neuf par ETKDG, sans graine
fixée, à chaque appel. Chaque bras recevait donc ses charges d'un conformère
aléatoire distinct.

**Les deux mesures qui tranchent** (iE-DAP, 42 atomes) :

| Comparaison | Écart max sur les charges |
|---|---|
| Deux **conformères** de la même molécule | **0.085 e** |
| Naturel vs **miroir** (même conformère, réfléchi) | **< 1e-4 e** |

Facteur ~1000. La réflexion ne change rien ; le conformère change tout. L'écart de
0.085 e n'est pas du bruit d'arrondi : c'est une différence massive sur une
électrostatique de liaison, et elle partait intégralement dans le ΔΔG, sans le
moindre rapport avec la chiralité.

**Pourquoi le partage est *exact* et non une approximation.** L'hamiltonien AM1 ne
dépend que de distances interatomiques, invariantes par réflexion. La structure
électronique du miroir est donc l'image exacte de celle du naturel, et les charges
— des scalaires — sont rigoureusement identiques atome par atome. C'est le même
argument de symétrie que pour ff14SB (point 1), et il est **mesuré** ici
(`test_mirror_gets_identical_charges`, marqué `slow`) plutôt qu'affirmé. Si ce
test cassait, le partage redeviendrait une approximation à documenter.

**Décision.** `ligand_partial_charges()` calcule les charges **une seule fois**,
sur le conformère de référence d'`iedap_naturel.sdf`, avec `use_conformers`
explicite (donc déterministe), les met en cache, et les deux bras les
réutilisent. La signature atomique (suite des éléments) voyage avec les charges
et interdit de les appliquer à une molécule d'ordre différent — risque réel, smina
réordonnant les atomes (point 11).

Filet supplémentaire : la somme des charges partielles doit rendre la charge
formelle nette (−1 pour iE-DAP). Un écart signalerait un protomère mal transmis
— ce qui relie ce point au point 10.

**Note sur NAGL.** Le registre openff essaie `NAGLToolkitWrapper` en premier — un
réseau de neurones sur graphe, insensible au conformère, qui aurait rendu le
problème inexistant. Mais il **refuse** le nom `am1bcc` (il n'accepte que ses
propres modèles `.pt`) et passe la main. C'est bien antechamber qui calcule.
L'avertissement NAGL sur `use_conformers` au passage est sans conséquence.

**Limite assumée.** AM1-BCC sur **un seul conformère**. La méthode robuste pour
les molécules chargées est `am1bccelf10` (moyenne sur conformères ELF), qui exige
le toolkit OpenEye, sous licence commerciale, non disponible ici. Le choix du
conformère de référence influence donc les charges en valeur absolue — mais
**identiquement pour les deux bras**, donc sans effet sur le ΔΔG, qui est
l'observable.

---

## 13. Le bras miroir part de la pose naturelle réfléchie, pas d'un docking indépendant

**Décidé le 2026-07-15** (point 1.5 de l'audit, classe A pour le pilote).

**Problème.** Le score de docking est achiral en pratique — mesuré : **−5.9 kcal/mol
pour les deux bras** — et le champ de force est chiral-symétrique (ff14SB point 1,
GAFF2 achiral). Seule la **géométrie échantillonnée** peut donc discriminer les
deux bras : autrement dit, tout dépend de la pose initiale.

Or `05` dockait naturel et miroir **indépendamment**, et `07` prenait la meilleure
pose de chacun. Comme smina ne sait pas classer les poses du miroir (elles sont
quasi dégénérées pour son score), sa « meilleure » revient à un tirage à pile ou
face. Le ΔΔG mesurait donc ce tirage autant que la chiralité — il confondait
« placement différent » et « chiralité différente ».

**Décision.** Le bras miroir part du **reflet resuperposé** de la pose naturelle
(`ligand.mirror_superposed`) :

1. réflexion des coordonnées de la pose naturelle → l'énantiomère exact ;
2. resuperposition sur la pose naturelle par la rotation **propre** qui minimise le
   RMSD des atomes lourds.

L'étape 2 est nécessaire : la réflexion seule envoie la molécule ailleurs dans
l'espace (l'image par un plan arbitraire). Après resuperposition, le miroir occupe
le même volume, dans le même site, et ne diffère du naturel que par sa chiralité —
exactement la variable étudiée. Toute rotation propre supplémentaire donnerait un
autre point de départ valide ; on retient celui qui recouvre le mieux la pose
naturelle, seul choix canonique et sans paramètre libre.

**Le piège, et pourquoi il est testé.** Si la superposition autorisait une rotation
**impropre** (det = −1), Kabsch défairait exactement la réflexion et rendrait le
ligand **naturel**. Le bras miroir simulerait alors le naturel, le ΔΔG vaudrait 0,
et on conclurait « aucune différence » — l'erreur d'impact critique du CDC
section 12, sous une forme que rien ne signalerait. La correction de déterminant
force det = +1, et trois tests verrouillent :

- le résultat est toujours un énantiomère exact (tous les centres basculent) ;
- le RMSD au naturel est > 0.1 Å (un énantiomère ne peut pas se superposer à
  lui-même par rotation propre ; un RMSD nul signalerait la réflexion annulée) ;
- les distances internes sont conservées (isométrie).

**Limite assumée.** Ce point de départ favorise le mode de liaison du naturel : le
miroir n'y fait pas les mêmes contacts et part d'une pose contrainte. C'est voulu
— on isole la chiralité — et la MD le laisse se réarranger sur 50 ns. Mais cela ne
répond **pas** à « quel est le meilleur mode de liaison du miroir ? », qui est la
question biologiquement pertinente pour l'évasion immunitaire. Y répondre demande
plusieurs poses de départ et des répliques (classe B). À dire explicitement dans le
preprint : le ΔΔG produit ici mesure la perte de reconnaissance **dans le mode de
liaison du motif naturel**, pas l'impossibilité de toute liaison.

**Ce qui ne change pas.** `05` continue de docker le miroir : le ΔΔG de docking
documente justement que le docking ne discrimine pas. Il n'amorce simplement plus
la MD.

---

## 14. La Phase 4 ne s'était jamais exécutée — réparation et accélération

**Décidé le 2026-09-04.** Relecture du code avant de relancer une production.

### Le constat

`mmgbsa.py` produit l'observable centrale du projet — le ΔΔG. Il **plantait dès
la première ligne de travail réel**, et personne ne le savait : aucun test ne
couvrait `binding_energy`, et son unique appelant, `scripts/08`, n'avait jamais
tourné.

C'est le trou de couverture qui compte, plus que le bug : la suite comptait 124
tests, dont 94 lignes sur `mmgbsa.py` — **toutes sur l'arithmétique des
dataclasses** (propagation d'incertitude, critère de significativité). Les tests
étaient épais là où le code est facile et absents là où il est décisif. Une
production de 15 h aurait tourné, puis se serait arrêtée sur une `ValueError` à
l'étape finale.

### Deux erreurs de construction (mesurées, pas déduites)

1. **`implicitSolvent=GBn2` n'est pas un paramètre de
   `ForceField.createSystem`.** Il appartient à la route
   `AmberPrmtopFile.createSystem`. Sur la route ForceField, c'est le *fichier*
   `implicit/gbn2.xml` qui ajoute la force GB. Le kwarg lève
   `The argument 'implicitSolvent' was specified to createSystem() but was never used`.
   Le code portait déjà, pour le sel, la remarque que les deux routes diffèrent —
   sans en tirer la conséquence pour le solvant lui-même.
2. **`nonbondedMethod` ne peut pas vivre dans `forcefield_kwargs`** :
   openmmforcefields l'exige dans `periodic_` ou `nonperiodic_forcefield_kwargs`.

Vérifié après correction : `CustomGBForce` couvre **100 % des particules des
trois systèmes** (complexe 4800/4800, récepteur 4758/4758, ligand 42/42). Le
ligand GAFF reçoit donc bien des paramètres GB — il n'est pas évalué dans le
vide, ce qui aurait faussé le bilan sans planter. Et `implicitSolventKappa` est
réellement consommé : E(ligand) change de −0.188 kcal/mol entre 150 mM et sans
sel.

### Un troisième piège, introduit puis mesuré

La première correction retirait les vecteurs de boîte sur une `copy.deepcopy` de
la topologie, pour ne pas muter l'objet de l'appelant. Résultat :
`No template found for residue 1 (UNK)`. **`deepcopy` d'une topologie OpenMM
duplique aussi les objets `Element`**, et l'appariement de templates, qui les
compare, échoue sur la copie. Corrigé par un retrait/restauration sur place
(gestionnaire de contexte `_without_box`), et verrouillé par un test qui vérifie
qu'une topologie périodique passe *et* ressort intacte.

### Le correctif 1.3 n'avait été appliqué qu'à moitié

L'audit 1.3 (charges AM1-BCC partagées) avait été corrigé côté **production**, pas
côté **scoring** — c'est-à-dire pas là où le ΔΔG se calcule. `scripts/08`
appelait `load_ligand()` sans charges, et openmmforcefields en recalculait depuis
un conformère ETKDG tiré au hasard, trois fois par bras.

Mesures sur iE-DAP (2026-09-04) :

| Comparaison | Écart max sur les charges |
|---|---|
| Deux passages du **même fichier** dans le chemin de scoring | **0.092 e** |
| Bras naturel vs bras miroir, au scoring | **0.390 e** |
| Charges de **production** vs charges de **scoring** | **0.423 e** |

Ordre atomique et graphe de liaisons vérifiés identiques entre les deux bras : ce
n'est donc pas un décalage d'indices, c'est bien le tirage du conformère. Le
ΔΔG aurait été calculé sur une paramétrisation différente de celle qui avait été
simulée. Les charges du cache de production sont désormais passées explicitement
à `binding_energy`.

### Travail inutile supprimé

`binding_energy` évaluait **chaque énergie deux fois** : une boucle
`_mean_energy` pour les moyennes, puis une seconde boucle identique pour la série
par image. Un seul passage suffit — exactement moitié moins de points d'énergie.

De même, `_implicit_generator` était appelé **trois fois par bras** (complexe,
récepteur, ligand), chaque appel relançant antechamber. Un seul générateur sert
désormais les trois systèmes, et le ligand portant déjà ses charges partagées,
antechamber ne tourne plus du tout en Phase 4.

### Accélérations, mesurées

**GPU pour les points d'énergie.** Le GB `NoCutoff` est en O(N²) ; sur 4800 atomes
il domine toute la Phase 4. Mesure sur 40 images, *setup exclu* (une première
mesure sur 5 images concluait à tort « aucun gain » : elle était dominée par la
construction des systèmes) :

| Plateforme | Coût par image | 300 images × 2 bras |
|---|---|---|
| CPU | 1260.6 ms | 12.6 min |
| CUDA (double) | **107.8 ms** | **1.1 min** |

Soit **×11.7**, et ×23 par rapport au code d'origine qui doublait le travail. La
double précision est imposée sur GPU pour que les énergies restent comparables :
accord mesuré CPU vs CUDA à **0.001 kcal/mol** sur le ΔG, très en dessous du
signal recherché.

**HMR : pas d'intégration de 4 fs.** Chaque hydrogène est porté à 4 uma, la masse
étant retirée à l'atome lourd qui le porte. Mesure sur le complexe NOD1 solvaté
(54 738 atomes) :

| Configuration | Débit | Production 50 ns × 2 bras |
|---|---|---|
| 2 fs, sans HMR | 204.8 ns/jour | — |
| 4 fs, avec HMR | **407.2 ns/jour** | — |

**×1.99.** Pourquoi c'est légitime ici, et pas un raccourci :

- **seules les masses changent, jamais le potentiel.** La thermodynamique
  d'équilibre est inchangée — les masses ne figurent pas dans la distribution de
  Boltzmann des positions — et le ΔΔG est une quantité d'équilibre ;
- **la symétrie chirale est intacte** : une masse est un scalaire, invariante par
  réflexion. Le test de symétrie du champ de force passe inchangé ;
- **les deux bras reçoivent exactement le même traitement** ;
- la masse totale est conservée (vérifiée par test) : le barostat NPT ne voit
  aucune densité faussée.

Ce qui change réellement : les modes de vibration les plus rapides ralentissent,
donc les temps de corrélation en *pas* diffèrent. Sans conséquence sur une
moyenne d'équilibre, à noter si l'on exploitait un jour la dynamique elle-même.

Écart au CDC (section 7), assumé et documenté ici. Désactivable :
`MDConfig(timestep=2.0, hydrogen_mass=0.0)`.

### Ce que ça ne règle pas

La Phase 4 s'exécute maintenant de bout en bout et vite. Elle ne devient pas pour
autant publiable : il manque toujours le contrôle de calibration (audit 1.1), les
répliques et l'expérience nulle naturel-vs-naturel (audit 1.2), le critère de
convergence sur la série d'énergie (audit 2.3) et l'ancrage externe du site NOD1
(audit 2.8). Un ΔΔG rapide reste un ΔΔG non calibré.

---

## 15. Optimisation : deux pistes rejetées par la mesure, la fiabilité prise au sérieux

**Décidé le 2026-09-04**, après le point 14.

### Deux optimisations mesurées, deux refusées

Les deux paraissaient évidentes. Les deux sont mauvaises, et seule la mesure le
dit.

**Précision mixte pour les points d'énergie MM/GBSA.** La RTX 4070 calcule en
FP64 à 1/64 du débit FP32 : imposer la double précision paraissait un luxe.

| Précision | ΔG (mêmes images) | Coût par image |
|---|---|---|
| double | **−27.0933** kcal/mol | 102.0 ms |
| mixed | −24.6152 kcal/mol | 19.2 ms |
| single | −24.6152 kcal/mol | 17.4 ms |

×5.3 de gain, pour **2.48 kcal/mol de biais** — l'ordre de grandeur du signal
recherché. Rejeté. La raison est structurelle : le ΔG est une différence de
grands nombres (−7271 − (−7170) = −24.6, soit 0.3 % des magnitudes), régime où
l'annulation catastrophique amplifie toute erreur d'arrondi. La double précision
n'était pas de la prudence excessive, elle est nécessaire — et c'est désormais
mesuré, plus supposé.

Décomposition du coût restant : `setPositions` 0.16 ms, `getState(getEnergy)`
43.4 ms. **Le surcoût Python est de 0 %** : les 102 ms sont du calcul. Ce module
est à son plancher.

**Fréquence du barostat Monte-Carlo.** Chaque tentative évalue l'énergie du
système entier ; la valeur par défaut (25 pas) paraissait coûteuse.

| Fréquence | Débit |
|---|---|
| sans barostat (plafond) | 343.5 ns/jour |
| tous les 25 pas | 291.2 ns/jour |
| tous les 100 pas | 245.3 ns/jour |
| tous les 250 pas | 277.3 ns/jour |

Non monotone : espacer les tentatives devrait accélérer, on observe l'inverse
puis le contraire. Le signal est **noyé dans la variance machine**, de l'ordre de
±20 % d'un run à l'autre (horloges GPU, thermique). Aucun gain reproductible à
réclamer, donc aucun changement. Corollaire à retenir pour les estimations de
durée : le « 407 ns/jour » du point 14 porte cette même incertitude, et il a de
surcroît été mesuré **sans barostat** — erreur corrigée au point 16.

### La fiabilité, elle, avait un vrai gisement

Cette machine a redémarré **deux fois** pendant le projet. Les deux fois, la
production est repartie de zéro — dont une à 37.6 ns sur 50. Le checkpoint
existait à chaque fois ; rien ne savait s'en servir.

**Reprise sur checkpoint** (`run_full_protocol(..., resume=True)`). Trois
précautions, chacune nécessaire :

1. *Le barostat est ajouté avant le chargement.* `loadCheckpoint` exige un
   système identique à celui qui l'a écrit ; la production tournant en NPT,
   l'oublier ferait échouer le chargement.
2. *La trajectoire est tronquée au dernier état sauvegardé.* Le DCD et le
   checkpoint sont écrits par deux reporters distincts et le processus peut
   mourir entre les deux. Conserver les images orphelines décalerait toute la
   série temporelle.
3. *Les sorties sont ouvertes en prolongement.* Une reprise qui écraserait le
   DCD détruirait exactement ce qu'elle vient sauver.

**Le système est relu, jamais re-solvaté** — et ce point n'est pas une
optimisation. Vérifié dans le code installé (`openmm/app/modeller.py:353`) :
`addSolvent` remplace des molécules d'eau par des ions via
`random.choice(replaceableList)`, **sans graine fixée**. Deux solvatations du
même complexe donnent le même *nombre* d'atomes dans un *ordre* différent. Or
`loadCheckpoint` ne vérifie que le nombre de particules : il accepterait, et
appliquerait positions et vitesses aux mauvais atomes. Pas d'erreur, pas
d'avertissement, trajectoire corrompue qui continue. On relit donc le PDB
solvaté écrit au lancement, qui *est* le système du checkpoint.

Vérifié de bout en bout sur le système NOD1 réel : production courte, salissure
délibérée du DCD par une image orpheline, reprise à durée doublée. Image
orpheline supprimée, images accumulées sans doublon, série temporelle
strictement croissante à intervalle constant, et une reprise sur production déjà
complète ne fait rien.

**Fréquence du checkpoint : mon propre arbitrage corrigé par la mesure.** J'avais
synchronisé le checkpoint sur chaque image, pour une reprise « exacte ». Mesure
(54 738 atomes) : une écriture coûte **54 ms**, synchronisation GPU → CPU
comprise.

| Fréquence | Écritures sur 50 ns | Coût | Part de la production |
|---|---|---|---|
| chaque image | 5000 | 272 s | **2.56 %** |
| toutes les 10 images | 500 | 27 s | **0.26 %** |

Payer 2.3 % sur *chaque* production pour économiser 90 ps lors des *rares*
reprises est un mauvais échange. Retour à 10, la troncature rendant de toute
façon n'importe quel intervalle exact.

**Garde-fou d'instabilité.** Un système qui explose à 30 ns continuait d'écrire
des NaN pendant des heures : trajectoire perdue, GPU occupé jusqu'au bout.
`InstabilityGuard` lève dès que l'énergie cesse d'être finie. Coût nul :
l'énergie est déjà calculée pour le journal, on l'échantillonne au même rythme.

**Contrôle pré-vol** (`scripts/preflight.py`, appelé par `launch.sh`). Trente
secondes qui vérifient, dans l'ordre où les choses cassent réellement : aucune
production concurrente, plateforme CUDA, espace disque suffisant pour les
trajectoires, fichiers de préparation présents, protomère des deux ligands,
cohérence du cache de charges, poses de départ constructibles et miroir
énantiomère exact. Chaque échec y coûte trente secondes ; le même échec
découvert en production coûte la production.

Observation au passage, notée pendant la mesure des autres couples : le
récepteur TLR1/2 **produit des NaN dès le premier pas** si la minimisation est
écourtée à 100 itérations. La minimisation complète (105 s) est nécessaire.

### Ce que ça change

Rien sur la physique, rien sur le ΔΔG. Une production interrompue reprend au
lieu de recommencer, une production qui diverge s'arrête au lieu de tourner à
vide, et une erreur de préparation se voit avant l'engagement du GPU et non
après. Sur une machine qui redémarre, c'est le poste où le temps se perdait
réellement.

---

## 16. Relecture de crédibilité : cinq affirmations fausses ou fragiles, corrigées

**Décidé le 2026-09-05.** Relecture des documents comme le ferait un rapporteur :
en recalculant les chiffres et en cherchant les contradictions internes, plutôt
qu'en relisant le texte.

### 16.1 — Le débit publié n'était pas celui de la production (grave)

Le chiffre affiché partout, **407 ns/jour**, a été mesuré sans barostat. Or
`run_full_protocol` appelle `equilibrate_npt`, qui ajoute un `MonteCarloBarostat`
**conservé pendant toute la production** : la production tourne en NPT.

Mesure refaite dans la configuration exacte de production, 5000 pas chronométrés,
trois répétitions :

| Conditions | Débit (médiane de 3) | Plage |
|---|---|---|
| NVT, sans barostat | 385.8 ns/jour | 339–389 |
| **NPT, production réelle** | **262.8 ns/jour** | 256–317 |

Le chiffre publié était donc **1.55× trop optimiste**, et toutes les durées qui
en découlaient avec lui : une production 50 ns × 2 bras prend **≈ 9.1 h**
(fourchette 7.6–9.4 h), et non 5.9 h.

Le gain HMR de ×1.99 (point 15) reste valide : ses deux termes ont été mesurés
dans les mêmes conditions NVT, et le rapport est ce qui était en jeu. Seules les
valeurs absolues étaient inutilisables comme durées de production.

Corrigé dans `README.md`, `README-TECHNIQUE.md`, `HANDOFF.md` et
`scripts/preflight.py`, qui annonçait la durée estimée au lancement.

**Leçon** : une mesure de performance doit reproduire la configuration réelle, y
compris ses composants « accessoires ». Le barostat n'était pas un détail, il
coûte un tiers du débit.

### 16.2 — Une estimation reposait sur une structure écartée

Le coût GPU annoncé pour TLR5 (342 179 atomes, 23.2 h) a été mesuré sur
`complexe_naturel.pdb`, qui provient de **3V47** : la structure chimérique de
poisson-zèbre explicitement écartée du projet (point 8). Le modèle ColabFold du
TLR5 humain n'a pas encore été produit.

Le chiffre ne mesurait donc pas le système qui sera simulé. Retiré. TLR5 est
désormais marqué **non estimable** jusqu'à ce que le modèle humain existe.

Le tableau des coûts distingue maintenant explicitement *mesuré*, *extrapolé* et
*non estimable*. TLR1/2 est extrapolé : débit mesuré sans barostat (186 ns/jour)
ramené aux conditions de production par le rapport NPT/NVT mesuré sur NOD1
(0.68), soit ≈ 127 ns/jour et ≈ 19 h.

### 16.3 — « Sans exception » était faux, et contredit par le ligand du projet

`README-SIMPLE` affirmait que tout le vivant est homochiral « sans exception » ;
`README-TECHNIQUE` parlait d'homochiralité « universelle ». Deux problèmes :

- la **glycine** n'a pas de centre stéréogène : elle n'est ni L ni D ;
- les **D-acides aminés existent** dans le vivant, par des voies non
  ribosomiques : D-alanine et D-glutamate du peptidoglycane, méso-DAP (lui-même
  achiral), plusieurs antibiotiques peptidiques.

Surtout, l'affirmation était **contredite deux pages plus loin par le projet
lui-même** : iE-DAP, le ligand du couple pilote, contient un γ-D-glutamyl et un
méso-DAP. C'est précisément pour cela que son miroir n'est pas « la version D »
mais l'énantiomère complet de ses trois centres (point 4.2). Un rapporteur aurait
relevé la contradiction immédiatement.

Reformulé : l'homochiralité vaut pour la **traduction ribosomique**, exclusivement
en L. Les exceptions sont réelles, circonscrites, et une bactérie miroir les
inverserait aussi. L'argument du projet en sort renforcé, pas affaibli.

### 16.4 — pLDDT annoncé à 96, mesuré à 95.0

Recalculé sur les 2358 atomes du domaine LRR découpé : moyenne **95.0**, médiane
97.0, minimum 73.2, et **aucun atome sous 70**. L'écart au chiffre annoncé est
mineur, mais c'est exactement le genre de valeur qu'un rapporteur recalcule.

La formulation gagne au passage : « aucun atome sous 70 » est un argument
vérifiable, là où « pLDDT ≈ 96 » n'était qu'une moyenne arrondie.

### 16.5 — Deux incohérences de moindre portée

**Renvoi interne faux.** Le tableau des paramètres de production renvoyait au
§4.5 pour le HMR, qui se trouve au §4.7 ; §4.5 traite des charges AM1-BCC.

**Comptage d'atomes incohérent.** Le §4.4 annonçait « 30 atomes rendus sur 43 »
pour iE-DAP, à côté d'un protomère `C12H20N3O7⁻` qui en compte 42. La mesure
avait été faite sur la forme neutre alors en usage, avant la correction du
protomère (point 10). Précisé.

**Raison d'abandon absente du code.** `config.py` conservait Dectine-1 avec la
note « si temps », alors que les README annoncent son abandon pour cause de
structure murine. Un lecteur confrontant code et documentation y aurait vu une
divergence. La raison est désormais dans le code.

### Ce que cette relecture ne corrige pas

Les limites de fond restent entières et sont listées au §7 du README technique :
absence de contrôle de calibration, absence de répliques et d'expérience nulle,
convergence jugée sur le RMSD plutôt que sur l'énergie, entropie négligée, et
surtout le site de liaison de NOD1 qui reste déduit de la forme du récepteur
plutôt qu'observé. Aucune correction éditoriale ne les lève.

---

## Récapitulatif des points ouverts

| # | Question | Statut |
|---|---|---|
| 2 | Ordre des couples | **tranché : NOD1 → TLR4 → TLR5 → TLR1/2** |
| 3 | TLR5 poisson-zèbre ou humain | **tranché : ColabFold sur Q9NR61, 100% humain** |
| 7 | Cible NOD1 | **tranché : AlphaFold Q9Y239, LRR** |
| 1 | Symétriser les CMAP de ff19SB ? | écarté du chemin critique |

**Décidés le 2026-07-14** (points 2, 3, 7). L'ordre privilégie la validation de
la pipeline petite-molécule (réutilisée 4 fois sur 5) sur un système léger, puis
TLR4 pour l'impact clinique (seule structure déjà humaine), TLR5 en dernier des
lourds car il demande un outil de docking distinct.
