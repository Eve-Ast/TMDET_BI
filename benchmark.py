import time
from src.Protein import Protein
from src.Grid import Grid

pdb_path = r"data\1PRN.pdb"
rsa_path = r"result\1PRN.rsa"

proteine = Protein("prot")
proteine.extract_calpha_coords_and_sasa(pdb_path, rsa_path)
proteine.center_protein()

grid_test = Grid(n_points=700)  # réduit pour comparer les 4 variantes équitablement

variants = {
    "Méth. 1 - non vectorisée": lambda: grid_test.scan_protein_meth_1(proteine),
    "Méth. 1 - vectorisée": lambda: grid_test.scan_protein_meth_1_vectorized(proteine),
    "Méth. 2 - non vectorisée": lambda: grid_test.scan_protein_meth_2_non_vectorized(proteine),
    "Méth. 2 - vectorisée": lambda: grid_test.scan_protein_meth_2(proteine),
}

for name, func in variants.items():
    t0 = time.perf_counter()
    result = func()
    t1 = time.perf_counter()
    print(f"{name:30s} | {t1-t0:8.4f} s | score = {result['best_score']:.3f}")