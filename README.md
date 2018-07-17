# Hassan Kachalo Ali


Stage	de	Master 1	de	Bioinformatique	&	Génomique
---------------------------------------------------------


Période	de	[Avril ~ Juillet]	2017	à	l'Inria	de	Rennes
--------------------------------------------------------------



### Les dépendances
-------------------

[![Dependencies Paloma](https://img.shields.io/badge/Paloma-v0.1.0-green.svg)](http://protomata-learner.genouest.org/)
[![Dependencies Pymol](https://img.shields.io/badge/Pymol-v0.1.0-blue.svg)](https://www.pymol.org/)
[![Dependencies Clingo v5.2](https://img.shields.io/badge/Clingo-v5.2-0077ea.svg)](https://github.com/potassco/clingo/releases/tag/v5.2.0)
[![Dependencies blastp](https://img.shields.io/badge/Dependencies-blastp-cc9900.svg)](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
[![Dependencies Bio](https://img.shields.io/badge/Dependencies-Bio-cc9900.svg)](http://biopython.org/wiki/Download)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

>> ### Le module Bio est à insérer dans le dossier 3D


Mots clés résumant le stage: 
============================

| **Proteines** | **Alignement des séquences(Paloma)** | **Classification** | **Structures 3D** | **Visualization** |

Principe:
=========

On part d'un ensemble de séquences protéiques connues et/ou non (fichier .fasta), de la même super-famille (souhaité).

1. On éffectue un alignement multiple, local et partiel de séquencces (Paloma outil de F. Coste)
2. On détermine la famille des séquences inconnues.
3. Visualisation des blocs des familles en 3D.


Fichiers d'entrées: un fichier **fasta** et un fichier **d'attributs**

Explication du code
====================
On prend pour exemple de fichiers **foo.fasta** et **foo.attr** pour la suite de description

Étape 1: Paloma
===============

L'outil Paloma est utilisé pour effectuer l'alignement. Ce programme prend en paramètre un fichier 
fasta des séquencces protéiques afin de produire un fichier PLMA (**P**artial **L**ocal **M**ultiple **A**lignment).
Il s'exécute en ligne de commande de la manière suivante:
```bash

shell> paloma -i foo.fasta

```

--------------------

On produit un fichier .plma. Soit foo.plma un exemple de fichier produit à ce niveau.
Ce fichier sera utiliser par le programme qui va suivre dans la suite de l'explication.


Étape 2: plma2dot
=================

On utilise ce programme pour convertir le fichier .plma précédemment produit en un fichier dot.
Le programme prend en paramètre un fichier PLMA afin de produire un fichier dot.
Il s'exécute en ligne de commande de la manière suivante:
En prenant comme fichier PLMA le fichier foo.plma, 
```bash

shell> plma2dot -i foo.plma

```

------------------

Le fichier dot de sortie sera du nom: foo_plma.dot


c) colordot2.py: 
On utilise ce programme pour affecter les couleurs aux séquences selon les familles distinguées
au niveau du fichier d'attributs d'entrée (foo.attr).
Il s'exécute en ligne de commande de la manière suivante:
(shell> python colordot2.py -d foo_plma.dot -a foo.attr -o foo_col.dot )
où foo_col.dot est le le nom du fichier de sortie que va produire le programme.


d) plmadot2logic2.py:
Comme son nom l'indique, ce programme permet de convertir le fichier dot du PLMA en un fichier au
format logique.
Il s'exécute en ligne de commande de la manière suivante:
(shell> python plmadot2logic2.py -f foo.fasta -d foo_col.dot -c black -a foo.attr > foo_logic.lp)
où 
foo_logic.lp: est le nom du fichier au format logique que le programme va produire.

La stucture de ce fichier (logique) est la suivante:
sequence(Ident seq, Class).
block(N° block, Ident seq).

e) blockclassconcept3.lp (outil de J. Nicolas):
C'est un programme codé en ASP(Answer Set Programming) pour determiner les concepts formels contenus
dans un PLMA.
Il s'exécute en ligne de commande de la manière suivante:
(shell> ./clingo foo_logic.lp blockclassconcept3.lp --heu=Vsids --stats -n 0)
Soit foo_concepts.lp le fichier de sortie après l'exécution du programme blockclassconcept3.lp.

clingo: est un solveur pour les programmes logiques (http://potassco.sourceforge.net/).
Il faut utiliser la version 5.2 de clingo (supportant python).
Le programme doit etre compilé sur la machine qui va l'utiliser.
S'il est compilé sur une machine et transporté sur autre, son binaire ne marche pas.


f) fullclassify7.lp (outil de J. Nicolas):
C'est un programme codé en ASP qui permet de faire la classification des séquences.
Il s'exécute en ligne de commande de la manière suivante:
(shell> ./clingo foo_concepts.lp blockclassconcept7.lp --heu=Vsids --stats -n 0 > foo_classif.lp)
où
foo_classif.lp: est le fichier de sortie après l'exécution du programme


g) ecrire2.py:
C'est un programme qui prend en paramètre le fichier de classification afin d'extraire et d'écrire
le résultat dans un fichier tabulé & en pdf



*****************************************************************************************************
Partie: Visualisation en 3D des blocs PLMA
*****************************************************************************************************


L'objectif est de pouvoir associer des structures pdb au blocs PLMA et de les rendre cliquable afin
qu'ils puissent etre visualisés en 3D. Les étapes sont les suivantes:

a) sp2pdb.py:
La 1ère étape est celle qui va consiter à établir un fichier d'attributs pdb que nous allons nommer:
foo_pdb.attr
Ce fichier est à générer automatiquement. Pour cela nous disposons de deux fichiers qui coniennent
les identifiants des séquences protéiques (identifiants swissprot/uniprot) associés à des identifiants
pdb des structures.
Ces deux fichiers vont nous servir de base de données pour rechercher les identifiants pdb associés à
un identifiants swissprot.
Ces fichiers sont contenus dans le dossier database du reperpoire "Programs"
-> pdbtosp.tab: est un fichier que l'on trouve au niveau de (http://www.uniprot.org/docs/pdbtosp)
-> uniprotData.tab: est un fichier que l'on peut trouver au nieveau de 
(ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/)
Les scripts qui permettent l'extraction des 2 fichiers fichiers précédents sont dans le meme dossier
que les fichiers cités (il s'agit de: pdbtosp.py et uniprotData.py).

Une fois que tous ces 2 fichiers sont disponibles, le scripts sp2pdb.py permet d'établir le fichier
d'attributs pdb sous le principe qui suit:
- on établit la liste des identifiants des protéines de test
- pour chaque identifiant, on récherche les identifiants pdb qui lui sont associés
- on télécharge un à un de manière automatique les identifiants pdb associés
- on extrait les séquences protéiques de ces structures associées (le script: pdb2fasta.sh)
- on éffectue un blastp (alignement pairwise) entre la séquence de test et toutes les séquences des 
strucures associées.
On récupères les e-values, pour enfin sélectionner la meilleure structure (ayant la plus petite e-value)
- les structures choisies sont stockées dans un dossier qui s'appelle PDB.
Souvent, certaines stuctures viennent avec l'extention .cif, on se sert du script (remane.sh) pour les
renommer en leur donnant une extention pdb
-En fin, on établit le fichier d'attributs pdb pour les séquences de test

b) color3d2.py
Ce programme prend en paramètre le fichier fasta, le fichier dot des couleurs des séquences et le 
fichier d'attributs pbd de ces séquences dont le but est de rendre les blocs de PLMA cliquable et 
d'associer des structures pdb à ces blocs afin de pouvoir etre visualiser.

Il utilise le programme alinpdbtofas.py qui est un programme qui réaligne la structure en fonction de 
la séquence par un blast, crée un fichier PDB, renumérote la structure en fonction de la séquence fasta.

Les blocs qui sont associés à des structures pdb sont cliquables, ceux qui ne sont pas associés à des 
structures pdb, ne sont pas cliquables.

Le code couleur pour les résultats:
Couleur verte: le bloc de séquence sur lequel on clique, 
Couleur rouge: tous les blocs où il y a exactement les mêmes séquences qui passent autrement dit, les
blocs qui partagent le meme concept que le bloc sur lequel on a cliqué.
Couleur bleu: tous les blocs où il y a au moins les mêmes séquences qui passent.
Couleur orange: les sites actifs documentés dans fichier pdb

Comme on vient de le voir, les programmes sont dans deux dossiers differents:
- le dossier classification qui contient tous les programmes nécéssaires pour l'exécution de la partie
classification et 
- le dossier 3D qui contient tous les programmes nécéssaires pour la partie 3D.
Ces deux dossiers sont dans le dossiers "Programs"



******************************************************************************************************
Descriptions des jeux de données: Explication des séquences par Thierry dans le dossier Fwd-Listes-HADs
---------------------------------

*************
Arabidopsis *
*************
152 séquences + 14 séquences = 166 séquences total

Le jeu de données actuel est de 153 séquences au total: certaines séquences n’ont pas été récupérés car 
il y’ a des identifiants des séquences surlignés en rouge que Thierry a recommande d’enlever de l’analyse. 
Il pense que ces protéines en rouge ne sont pas des HADs. Les 14 séquences ont été récupérées sur le site 
de ncbi et les autres sur La base de données TAIR.


*************
Ectocarpus  *
*************
Le jeu de données actuel est de 178 séquences récupérées à partir du site de Orcae, dans les archives 
(https://bioinformatics.psb.ugent.be/gdb/ectocarpus/Archive/) au niveau du fichier 
(Ectsi_prot_2014Sep.tfa.gz). Les 178 identifiants ont été donnés par  Thierry 
(colonne A de la première feuille dans le fichier A_thaliana_HADs_v2.xls)



******************
Escherichia coli *
******************
Le jeu de données actuel est de 23 séquences récupérées à partir de UniProt/SwissProt. Les 23 identifiants 
des séquences ont été donnés par  Thierry 
(colonne G de la première feuille Kuznetsova et al. 2006 dans le fichier E. coli_23.xlsx).



*******
Human *
*******
40 séquences

Le jeu de données actuel est de 39 séquences récupérées à partir de UniProt /SwissProt. Dans  les identifiants 
des séquences donnés par  Thierry, il y a deux identifiants des protéines qui sont les mêmes mais les gènes 
codants pour ces deux protéines sont différents. J’ai pas pu retrouver les deux séquences protéiques, j’ai 
récupéré une dans   UniProt. c’est pourquoi il y a 39 séquences au lieu de 40. les identifiants des protéines 
(colonne C de la première feuille du fichier 40_human_HAD.xlsx).


**************************************************************************************************************** 
Remarques & Perspectives
****************************************************************************************************************

Pour la visualisation des blocs PLMAs:
- Il sera bien de supprimer les blocs non spécifiques et d'afficher uniquement les blocs caractéristiques
- Le Cadre grisé de la superfamille doit etre non visible (?)
- Les labels des séquences non connue, doivent-ils rester black après la classification (?)
- Il faut utiliser xpdf pour la visualisation des pdf
