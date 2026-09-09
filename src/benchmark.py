"""
Benchmarking Script for Protein Membrane Orientation Methods

This script compares the execution time of the four membrane orientation
scanning methods implemented in the Grid class.

The protein structure is loaded from a PDB file and the corresponding
solvent accessible surface area (SASA) data are loaded from an RSA file.
The protein is centered before the four scanning methods are executed.

Methods compared:
- Method 1 - non-vectorized
- Method 1 - vectorized
- Method 2 - non-vectorized
- Method 2 - vectorized

For each method, the execution time and the best hydrophobicity score
are displayed in the terminal.

Usage:
python src/benchmark.py <pdb_file> <rsa_file> <n_points>

Example:
python src/benchmark.py 1PRN.pdb 1PRN.rsa 1000

Arguments:
pdb_file:
Name or path of the PDB file containing the protein structure.

```
rsa_file:
    Name or path of the RSA file containing the SASA data.

n_points:
    Number of points used to generate the scanning grid.
```

Author:
Eve-Ang

"""

import argparse
import os
import time

from Protein import Protein
from Grid import Grid

def benchmark(pdb_path, rsa_path, n_points):
    """
    Benchmark the four protein scanning methods.

    ```
    Parameters
    ----------
    pdb_path : str
        Path to the PDB file.

    rsa_path : str
        Path to the RSA file.

    n_points : int
        Number of points in the scanning grid.
    """

    print("=" * 65)
    print("       BENCHMARKING DES MÉTHODES DE SCAN")
    print("=" * 65)

    # ---------------------------------------------------------
    # 1. Chargement de la protéine
    # ---------------------------------------------------------

    print("\n[1/3] Chargement de la protéine...")
    print(f"  → Fichier PDB : {pdb_path}")
    print(f"  → Fichier RSA : {rsa_path}")

    proteine = Protein("prot")

    proteine.extract_calpha_coords_and_sasa(
        pdb_path,
        rsa_path
    )

    print("  ✓ Coordonnées C-alpha et SASA chargées.")

    # ---------------------------------------------------------
    # 2. Préparation de la protéine et de la grille
    # ---------------------------------------------------------

    print("\n[2/3] Préparation de la protéine...")

    proteine.center_protein()

    print("  ✓ Protéine centrée à l'origine.")

    # La même grille est utilisée pour les quatre méthodes
    # afin de garantir une comparaison équitable.
    print(f"  → Création de la grille ({n_points} points)...")

    grid_test = Grid(n_points=n_points)

    # ---------------------------------------------------------
    # 3. Benchmark des quatre méthodes
    # ---------------------------------------------------------

    variants = {
        "Méth. 1 - non vectorisée":
            lambda: grid_test.scan_protein_meth_1_non_vectorized(
                proteine
            ),

        "Méth. 1 - vectorisée":
            lambda: grid_test.scan_protein_meth_1_vectorized(
                proteine
            ),

        "Méth. 2 - non vectorisée":
            lambda: grid_test.scan_protein_meth_2_non_vectorized(
                proteine
            ),

        "Méth. 2 - vectorisée":
            lambda: grid_test.scan_protein_meth_2_vectorized(
                proteine
            ),
    }

    print("\n[3/3] Benchmark des quatre méthodes...")
    print("-" * 65)

    print(
        f"{'Méthode':30s} | "
        f"{'Temps':>10s} | "
        f"{'Score':>10s}"
    )

    print("-" * 65)

    # Exécution successive des quatre méthodes.
    for name, func in variants.items():

        print(f"  → Exécution : {name}")

        # Mesure du temps d'exécution.
        t0 = time.perf_counter()

        result = func()

        t1 = time.perf_counter()

        elapsed_time = t1 - t0

        print(
            f"{name:30s} | "
            f"{elapsed_time:8.4f} s | "
            f"score = {result['best_score']:.3f}"
        )

    print("-" * 65)
    print("✓ Benchmark terminé.")
    print("=" * 65)


if __name__ == "__main__":


    # ---------------------------------------------------------
    # Récupération des arguments en ligne de commande
    # ---------------------------------------------------------

    parser = argparse.ArgumentParser(
        description=(
            "Benchmark des quatre méthodes de scan "
            "de l'orientation membranaire."
        )
    )

    parser.add_argument(
        "pdb_file",
        type=str,
        help="Nom ou chemin du fichier PDB (ex: 1PRN.pdb)"
    )

    parser.add_argument(
        "rsa_file",
        type=str,
        help="Nom ou chemin du fichier RSA (ex: 1PRN.rsa)"
    )

    parser.add_argument(
        "n_points",
        type=int,
        help="Nombre de points utilisés pour la grille (ex: 1000)"
    )

    args = parser.parse_args()

    # ---------------------------------------------------------
    # Construction des chemins des fichiers
    # ---------------------------------------------------------

    pdb_filename = args.pdb_file
    rsa_filename = args.rsa_file

    # Si aucun chemin n'est fourni pour le PDB,
    # le fichier est recherché dans le dossier data/.
    if not os.path.isabs(pdb_filename) and not os.path.exists(
        pdb_filename
    ):
        pdb_path = os.path.join("data", pdb_filename)
    else:
        pdb_path = pdb_filename

    # Si aucun chemin n'est fourni pour le RSA,
    # le fichier est recherché dans le dossier result/.
    if not os.path.isabs(rsa_filename) and not os.path.exists(
        rsa_filename
    ):
        rsa_path = os.path.join("result", rsa_filename)
    else:
        rsa_path = rsa_filename

    # ---------------------------------------------------------
    # Lancement du benchmark
    # ---------------------------------------------------------

    benchmark(
        pdb_path,
        rsa_path,
        args.n_points
    )
