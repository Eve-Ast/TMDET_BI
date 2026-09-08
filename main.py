import argparse
import math
from src.Protein import Protein, Residue
from src.Grid import Grid
import os


if __name__ == "__main__" : 
# Setup de l'analyseur d'arguments
    parser = argparse.ArgumentParser(
        description="PDB Protein Membrane Orientation Scanner"
    )
    parser.add_argument(
        "pdb_file",
        type=str,
        help="Nom du fichier PDB ou chemin complet (ex: 1QHJ ou 1QHJ.pdb)",
    )
    args = parser.parse_args()

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

    proteine = Protein("prot")
    print(proteine)
    #rsa_path = proteine.run_naccess(pdb_path)
    rsa_path = r"result\1G90.rsa"
    proteine.extract_calpha_coords_and_sasa(pdb_path, rsa_path)
    # print(proteine)

    proteine.center_protein()

    cx, cy, cz = proteine.calc_com()
    assert math.isclose(cx, 0.0, abs_tol=1e-5), f"Erreur sur X: {cx}"
    assert math.isclose(cy, 0.0, abs_tol=1e-5), f"Erreur sur Y: {cy}"
    assert math.isclose(cz, 0.0, abs_tol=1e-5), f"Erreur sur Z: {cz}"

    print("Protéine centrée avec succès à l'origine (0, 0, 0) !")

    grid = Grid(n_points=1000)
    # print(grid)

    # 3. La grille scannne la protéine
    results = grid.scan_protein(proteine)

    best_axis = results["best_axis"]
    best_score = results["best_score"]
    z_center = results["z_center"]

    print("\n=== RÉSULTATS DU SCAN PAR LA GRILLE ===")
    print(f"Meilleur axe détecté : {best_axis}")
    print(f"Score hydrophobe moyen (sur 30 Å) : {best_score:.3f}")
    print(f"Position Z du centre de la membrane : {z_center:.2f} Å")

    # Ouverture de PyMOL
    if best_axis is not None:
        print("\nOuverture de la fenêtre PyMOL...")
        proteine.show_in_pymol(best_axis, z_center, memb_thickness=30)
        # show_axis_in_pymol(proteine, best_axis, z_center)
    else:
        print(
            "Erreur : Aucun axe valide n'a pu être calculé sur la protéine."
        )