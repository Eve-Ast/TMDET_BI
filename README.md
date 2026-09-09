# PDB Transmembrane Protein - Membrane Orientation Scanner

Outil en ligne de commande pour prédire l'orientation et la position de la membrane lipidique par rapport à une protéine transmembranaire, à partir d'une structure PDB. Le programme centre la protéine, génère une grille de scan sur une sphère de directions, et identifie l'axe et la position membranaire les plus favorables au niveau hydrophobe selon plusieurs méthodes de calcul. Le résultat peut être visualisé directement dans PyMOL.

## Fonctionnalités

- Chargement d'une structure PDB de protéine transmembranaire et des données de surface accessible au solvant (SASA, fichier `.rsa`)
- Centrage automatique de la protéine à l'origine
- Génération d'une grille de points pour scanner les orientations possibles de la membrane autour de la protéine
- 4 méthodes de scan disponibles (vectorisées et non vectorisées, 2 approches différentes)
- Estimation de l'axe membranaire optimal, du score hydrophobe et de la position du centre de la membrane le long de cet axe
- Visualisation de la membrane prédite (plans hydrophobes) autour de la protéine dans PyMOL
- Script de benchmark pour comparer les performances des 4 méthodes

## Prérequis

- Python 3.x
- [PyMOL](https://pymol.org/) (pour la visualisation)
- [NACCESS](http://www.bioinf.manchester.ac.uk/naccess/) (pour générer les fichiers `.rsa` de SASA, si vous ne les avez pas déjà)

## Installation

```bash
git clone <url-du-repo>
cd <nom-du-repo>
```

Aucune dépendance Python externe n'est requise en dehors de la bibliothèque standard (`argparse`, `math`, `os`, `time`) et des modules internes du projet (`Protein`, `Grid`).

## Structure du projet

```
.
├── src/
│   ├── main.py          # Script principal
│   ├── benchmark.py     # Script de benchmark des méthodes
│   ├── Protein.py        # Classe Protein
│   └── Grid.py            # Classe Grid
├── data/                 # Fichiers PDB en entrée
└── result/               # Fichiers RSA (SASA) en entrée
```

> ⚠️ Adaptez cette arborescence si elle diffère de votre organisation réelle.

## Utilisation

### Scanner une protéine

```bash
python src/main.py <pdb_file> [--method <method>]
```

**Arguments :**
- `pdb_file` : nom ou chemin du fichier PDB de la protéine transmembranaire à analyser (l'extension `.pdb` peut être omise)
- `--method` : méthode de scan à utiliser parmi :
  - `method1_vec`
  - `method2_vec` *(par défaut)*
  - `method1_novec`
  - `method2_novec`

**Exemples :**

```bash
python src/main.py 1PRN.pdb
python src/main.py 1PRN.pdb --method method1_vec
```

**Sortie :**
- Meilleur axe d'orientation membranaire détecté
- Score hydrophobe associé
- Position estimée du centre de la membrane
- Ouverture d'une fenêtre PyMOL avec la visualisation

### Comparer les performances des méthodes (benchmark)

```bash
python src/benchmark.py <pdb_file> <rsa_file> <n_points>
```

**Exemple :**

```bash
python src/benchmark.py 1PRN.pdb 1PRN.rsa 1000
```

Le script affiche, pour chacune des 4 méthodes, le temps d'exécution et le score hydrophobe obtenu.

## Auteur

Eve-Angeline STEPHEN

## Licence

*(à préciser)*