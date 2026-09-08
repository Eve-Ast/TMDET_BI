import subprocess
import numpy as np
import os
import shutil
import tempfile
import textwrap

KYTE_DOLITTLE = {"ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
    "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3, "PRO": -1.6,
    "HIS": -3.2, "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5, 
    }

class Residue :
    """Resume court 
     Explication plus detaillé si nécessaire
     Args : 
        paramètres : description
    """

    def __init__(self, residue_name, res_num, x, y, z , sasa, hydrophobicity):
        self.residue_name = residue_name
        self.res_num = res_num
        self.x = x
        self.y = y
        self.z = z
        self.sasa = sasa 
        self.hydrophobicity = hydrophobicity


    def __str__(self):
        chaine = f"residue {self.residue_name} {self.res_num}, " \
                    f"coord: ({self.x:.3f}, {self.y:.3f}, {self.z:.3f}), "\
                    f"sasa: ({self.sasa}), " \
                    f"hydrophobicity: ({self.hydrophobicity}) \n"
        return chaine
                     
    

class Protein :

    def __init__(self, protein_name):
        self.protein_name = protein_name
        self.list_res = []

    def add_res(self, res):
        """resume"""
        if isinstance (res, Residue): 
            self.list_res.append(res)

    def run_naccess(self, pdb_path, naccess_bin="naccess", output_dir = "result"):
        """"""
        abs_pdb_path = os.path.abspath(pdb_path)
        # Verification existence du fichier pdb
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"Le fichier PDB {pdb_path} n'existe pas.")

        base_name = os.path.splitext(os.path.basename(abs_pdb_path))[0]

        # Execution de NACCESS
        try:
            subprocess.run(
                [naccess_bin, pdb_path],
                check=True,
                capture_output= True
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Erreur lors de l'exécution de NACCESS : {e.stderr}"
            )

        # Création du dossier de destination
        os.makedirs(output_dir, exist_ok=True)

        # Déplacement des extensions générées (.rsa, .asa, .log)
        for ext in [".rsa", ".asa", ".log"]:
            generated_file = f"{base_name}{ext}"
            if os.path.exists(generated_file):
                shutil.move(generated_file, os.path.join(output_dir, generated_file))

        rsa_path = os.path.join(output_dir, f"{base_name}.rsa")
        if not os.path.exists(rsa_path):
            raise FileNotFoundError(f"Le fichier de sortie {rsa_path} est introuvable.")

        return rsa_path
        

    def extract_calpha_coords_and_sasa (self, pdb_path, rsa_path): 
        """"""
        sasa_dict = {}

        with open(rsa_path, "r") as rsa_file:
            for line in rsa_file:
                # Les lignes 'RES' contiennent les valeurs par résidu
                if line.startswith("RES"):
                    chain_id = line[7:9].strip()

                    if chain_id == "A" : 
                        res_name = line[3:7].strip()
                        res_num = int(line[9:13].strip())
                        # SASA relative de la chaîne latérale (colonnes 55 à 62 selon la spécification NACCESS)
                        total_rel_sasa = float(line[24:28].strip())

                        # Clé unique pour identifier le résidu : (Chaîne, Numéro, Nom)
                        sasa_dict[(res_name, res_num)] = total_rel_sasa


        with open(pdb_path, "r") as pdb_file : 
            for line in pdb_file :
                if (line.startswith("ATOM")) :

                    chain_id = line[21:22].strip()
                    atom_name = line[12:16].strip()

                    if chain_id == "A" and atom_name == "CA" :
                        res_name = line[17:20].strip()
                        res_num = int(line[22:26].strip())
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])

                        sasa_val = sasa_dict.get((res_name, res_num), 0.0)

                        hydro_val = KYTE_DOLITTLE.get(res_name, 0.0)

                        res = Residue(res_name, res_num, x, y, z, sasa=sasa_val, hydrophobicity=hydro_val)

                        self.list_res.append(res)

    def calc_com(self): 
        """"""
        nb_res = len(self.list_res)
        sum_x = 0 
        sum_y = 0 
        sum_z = 0

        for res in self.list_res : 
            sum_x += res.x
            sum_y += res.y
            sum_z += res.z 

        com_x = sum_x / nb_res
        com_y = sum_y / nb_res
        com_z = sum_z / nb_res

        return com_x, com_y, com_z

    def center_protein(self):
        """"""
        com_x, com_y, com_z = self.calc_com()

        for res in self.list_res: 
            res.x -= com_x
            res.y -= com_y
            res.z -= com_z

    def get_arrays(self, min_sasa=15.0):
        """Vectorise les données pour accélérer le calcul avec NumPy."""
        filtered = [r for r in self.list_res if r.sasa >= min_sasa]
        if not filtered:
            filtered = self.list_res  # Secours si SASA n'est pas remplie

        coords = np.array([[r.x, r.y, r.z] for r in filtered])
        hydros = np.array([r.hydrophobicity for r in filtered])
        return coords, hydros

    @staticmethod
    def compute_rotation_matrix(v1, v2=np.array([0.0, 0.0, 1.0])):
        """Calcule la matrice de rotation 3x3 pour aligner le vecteur v1 sur v2 (axe Z canonique)."""
        v1 = v1 / np.linalg.norm(v1)
        v2 = v2 / np.linalg.norm(v2)
        v = np.cross(v1, v2)
        c = np.dot(v1, v2)
        s = np.linalg.norm(v)

        if s == 0:
            return np.identity(3) if c > 0 else -np.identity(3)

        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.identity(3) + vx + np.dot(vx, vx) * ((1 - c) / (s**2))

    

    def show_in_pymol(self, best_axis, z_center, memb_thickness=30):
        """Génère un PDB réorienté et lance PyMOL avec la membrane affichée."""
        R = self.compute_rotation_matrix(best_axis)
        half_thick = memb_thickness / 2.0

        # 1. Écriture du PDB temporaire réorienté
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False
        ) as pdb_temp:
            for res in self.list_res:
                vec = np.array([res.x, res.y, res.z])
                # 1. Alignement sur l'axe Z
                new_vec = np.dot(R, vec)
                # 2. Soustraction du décalage le long de cet axe
                new_vec[2] -= z_center

                pdb_temp.write(
                    f"ATOM  {res.res_num:5d}  CA  {res.residue_name:3s} A{res.res_num:4d}    "
                    f"{new_vec[0]:8.3f}{new_vec[1]:8.3f}{new_vec[2]:8.3f}  1.00  0.00           C\n"
                )

            pdb_path_out = pdb_temp.name

        # 2. Écriture du script PML (nettoyé de l'indentation grâce à textwrap.dedent)
        pml_content = textwrap.dedent(f"""\
            load {pdb_path_out}, prot_obj
            show cartoon, prot_obj
            color gray70, prot_obj

            python
            from pymol.cgo import *
            from pymol import cmd

            z_top = {half_thick}
            z_bot = -{half_thick}

            top_plane = [
                BEGIN, TRIANGLE_FAN,
                COLOR, 1.0, 0.4, 0.2,
                VERTEX, 0.0, 0.0, z_top,
                VERTEX, 30.0, 0.0, z_top,
                VERTEX, 0.0, 30.0, z_top,
                VERTEX, -30.0, 0.0, z_top,
                VERTEX, 0.0, -30.0, z_top,
                VERTEX, 30.0, 0.0, z_top,
                END
            ]
            cmd.load_cgo(top_plane, "top_plane")

            bottom_plane = [
                BEGIN, TRIANGLE_FAN,
                COLOR, 0.2, 0.6, 1.0,
                VERTEX, 0.0, 0.0, z_bot,
                VERTEX, 30.0, 0.0, z_bot,
                VERTEX, 0.0, 30.0, z_bot,
                VERTEX, -30.0, 0.0, z_bot,
                VERTEX, 0.0, -30.0, z_bot,
                VERTEX, 30.0, 0.0, z_bot,
                END
            ]
            cmd.load_cgo(bottom_plane, "bottom_plane")
            python end

            set cgo_transparency, 0.4, top_plane
            set cgo_transparency, 0.4, bottom_plane
            center prot_obj
            zoom prot_obj, 10
        """)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pml", delete=False
        ) as pml_temp:
            pml_temp.write(pml_content)
            pml_path_out = pml_temp.name

        # 3. Lancement de PyMOL
        try:
            subprocess.run(["PyMOLWin", pml_path_out])
        except FileNotFoundError:
            print(
                f"Erreur : Impossible de trouver PyMOL à l'emplacement."
            )
