"""
PDB Protein Membrane Orientation Scanner

This script provides a command-line interface for analyzing the membrane
orientation of a protein from a PDB structure. It loads the protein
structure and solvent accessible surface area (SASA) data, centers the
protein, generates a scanning grid, and identifies the most favorable
membrane orientation using one of several available scanning methods.

The program also allows visualization of the predicted membrane
orientation in PyMOL.

Usage:
python src/main.py <pdb_file> [--method <method>]

Arguments:
pdb_file:
Name or path of the PDB file to analyze. The ".pdb" extension
can be omitted.

--method:
    Scanning method used to determine the optimal membrane
    orientation. Available methods are:
        - method1_vec
        - method2_vec
        - method1_novec
        - method2_novec

    The default method is "method2_vec".

Input:
- PDB file containing the protein structure.
- RSA file containing solvent accessible surface area (SASA) data.

Output:
- Best membrane orientation axis.
- Best hydrophobicity score.
- Estimated position of the membrane center.
- PyMOL visualization of the predicted membrane orientation.

Example:
python src/main.py 1PRN.pdb

python src/main.py 1PRN.pdb --method method1_vec

Dependencies:
- argparse
- math
- os
- Protein
- Grid
- PyMOL (for visualization)

Author:
Eve-Angeline STEPHEN
"""

import argparse
import math
from Protein import Protein, Residue
from Grid import Grid
import os


if __name__ == "__main__":

    # 1. Configuration des arguments en ligne de commande
    print("\n" + "=" * 60)
    print("       PDB PROTEIN MEMBRANE ORIENTATION SCANNER")
    print("=" * 60)
    print("\n[1/5] Configuration des arguments...")

    parser = argparse.ArgumentParser(
        description="PDB Protein Membrane Orientation Scanner"
    )

    # Argument obligatoire (positionnel)
    parser.add_argument(
        "pdb_file",
        type=str,
        help="Nom du fichier PDB ou chemin complet (ex: 1PRN ou 1PRN.pdb)",
    )

    # Argument optionnel avec valeur par défaut
    parser.add_argument(
        "--method",
        type=str,
        choices=["method1_vec", "method2_vec", "method1_novec", "method2_novec"],
        default="method2_vec",
        help="Méthode de scan à utiliser (par défaut: method2)",
    )

    args = parser.parse_args()

    # 2. Traitement du nom de fichier et des chemins
    method = args.method

    print(f"  → Fichier PDB demandé : {args.pdb_file}")
    print(f"  → Méthode sélectionnée : {method}")

    # Traitement de l'argument (ajoute .pdb si non présent)
    pdb_input = args.pdb_file

    if not pdb_input.endswith(".pdb"):
        pdb_name = pdb_input
        pdb_filename = f"{pdb_input}.pdb"
    else:
        pdb_name = os.path.splitext(os.path.basename(pdb_input))[0]
        pdb_filename = pdb_input

    # Construction du chemin du fichier PDB
    pdb_path = (
        pdb_filename
        if os.path.isabs(pdb_filename) or os.path.exists(pdb_filename)
        else os.path.join("data", pdb_filename)
    )

    print(f"  → Chemin du fichier PDB : {pdb_path}")

    # Création de l'objet Protein
    print("\n[2/5] Chargement de la protéine...")
    proteine = Protein("prot")
    print(proteine)

    # Calcul de la surface accessible au solvant avec NACCESS
    rsa_path = proteine.run_naccess(pdb_path)
    # rsa_path = r"result\1PRN.rsa"

    print(f"  → Fichier RSA utilisé : {rsa_path}")
    print("  → Extraction des coordonnées C-alpha et des SASA...")

    proteine.extract_calpha_coords_and_sasa(pdb_path, rsa_path)

    print("  ✓ Coordonnées et SASA extraites.")

    # Centrage de la protéine
    print("\n[3/5] Centrage de la protéine...")

    proteine.center_protein()

    cx, cy, cz = proteine.calc_com()

    assert math.isclose(cx, 0.0, abs_tol=1e-5), f"Erreur sur X: {cx}"
    assert math.isclose(cy, 0.0, abs_tol=1e-5), f"Erreur sur Y: {cy}"
    assert math.isclose(cz, 0.0, abs_tol=1e-5), f"Erreur sur Z: {cz}"

    print("  ✓ Protéine centrée avec succès à l'origine (0, 0, 0) !")

    # Création de la grille
    print("\n[4/5] Création de la grille de scan...")

    grid = Grid(n_points=1000)

    print("  → Nombre de points : 1000")
    print("  → Lancement du scan de la protéine...")

    # 3. La grille scanne la protéine
    if method == "method1_novec":
        result = grid.scan_protein_meth_1_non_vectorized(proteine)

    elif method == "method1_vec":
        result = grid.scan_protein_meth_1_vectorized(proteine)

    elif method == "method2_novec":
        result = grid.scan_protein_meth_2_non_vectorized(proteine)

    elif method == "method2_vec":
        result = grid.scan_protein_meth_2_vectorized(proteine)

    print("  ✓ Scan terminé.")

    # Récupération des meilleurs résultats
    best_axis = result["best_axis"]
    best_score = result["best_score"]
    z_center = result["z_center"]

    print("\n" + "=" * 60)
    print("                 RÉSULTATS DU SCAN")
    print("=" * 60)
    print(f"  Meilleur axe détecté              : {best_axis}")
    print(f"  Score hydrophobe moyen (30 Å)     : {best_score:.3f}")
    print(f"  Position Z du centre de la membrane : {z_center:.2f} Å")
    print("=" * 60)

    # Ouverture de PyMOL
    print("\n[5/5] Visualisation dans PyMOL...")

    if best_axis is not None:
        print("  → Ouverture de la fenêtre PyMOL...")
        print(f"  → Axe utilisé : {best_axis}")
        print(f"  → Centre de la membrane : Z = {z_center:.2f} Å")
        print("  → Épaisseur de la membrane : 30 Å")

        proteine.show_in_pymol(
            best_axis,
            z_center,
            memb_thickness=30
        )

        print("  ✓ Visualisation PyMOL lancée.")

    else:
        print(
            "  ✗ Erreur : Aucun axe valide n'a pu être calculé sur la protéine."
        )

    print("\n" + "=" * 60)
    print("                    FIN DU PROGRAMME")
    print("=" * 60 + "\n")

