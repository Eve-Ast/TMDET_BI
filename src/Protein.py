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

class Residue:
    """
    Represent a protein residue and its associated properties.

    A residue is described by its name, residue number, three-dimensional
    coordinates, solvent-accessible surface area (SASA), and hydrophobicity.

    Attributes
    ----------
    residue_name : str
        Name of the residue using its three-letter code.
    res_num : int
        Residue number in the protein sequence.
    x : float
        X-coordinate of the residue.
    y : float
        Y-coordinate of the residue.
    z : float
        Z-coordinate of the residue.
    sasa : float
        Solvent-accessible surface area of the residue in Å².
    hydrophobicity : float
        Hydrophobicity value associated with the residue.
    """

    def __init__(self, residue_name, res_num, x, y, z, sasa, hydrophobicity):
        """
        Initialize a Residue object.
        """
        self.residue_name = residue_name
        self.res_num = res_num
        self.x = x
        self.y = y
        self.z = z
        self.sasa = sasa
        self.hydrophobicity = hydrophobicity

    def __str__(self):
        """
        Return a formatted string describing the residue.

        Returns
        -------
        str
            String containing the residue name, residue number,
            three-dimensional coordinates, SASA, and hydrophobicity.
        """
        chaine = (
            f"residue {self.residue_name} {self.res_num}, "
            f"coord: ({self.x:.3f}, {self.y:.3f}, {self.z:.3f}), "
            f"sasa: ({self.sasa}), "
            f"hydrophobicity: ({self.hydrophobicity})\n"
        )
        return chaine
                     
    

class Protein :
    """Represent a protein and its associated residues.
    
    Attributes
    ----------
    protein_name : str
        Name or identifier of the protein.
    list_res : list of Residue
        List containing the residues belonging to the protein.
    """

    def __init__(self, protein_name):
        """
        Initialize a Protein object.

        Parameters
        ----------
        protein_name : str
            Name or identifier of the protein.
        """

        self.protein_name = protein_name
        self.list_res = []

    def add_res(self, res):
        """
        Add a residue to the protein.

        The residue is added only if it is an instance of the
        ``Residue`` class.

        Parameters
        ----------
        res : Residue
            Residue to add to the protein.
        """

        if isinstance (res, Residue): 
            self.list_res.append(res)

    def run_naccess(self, pdb_path, naccess_bin="naccess", output_dir = "result"):
        """
        Run NACCESS to calculate the solvent-accessible surface area.

        NACCESS is executed on the provided PDB file. The generated
        ``.rsa``, ``.asa`` and ``.log`` files are moved to the specified
        output directory. The path to the generated ``.rsa`` file is
        returned.

        Parameters
        ----------
        pdb_path : str
            Path to the input PDB file.
        naccess_bin : str, optional
            Name or path of the NACCESS executable. Default is
            ``"naccess"``.
        output_dir : str, optional
            Directory where the NACCESS output files are stored.
            Default is ``"result"``.

        Returns
        -------
        str
            Path to the generated ``.rsa`` file.

        Raises
        ------
        FileNotFoundError
            If the PDB file does not exist or the expected ``.rsa`` file
            cannot be found.
        RuntimeError
            If NACCESS fails during execution.
        """
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
        """
        Extract C-alpha coordinates, SASA, and hydrophobicity for each residue.

        The PDB file is parsed to retrieve the three-dimensional coordinates
        of C-alpha atoms from chain A. The corresponding residue SASA values
        are obtained from the NACCESS ``.rsa`` file. A hydrophobicity value
        is assigned to each residue using the Kyte-Doolittle scale.

        The resulting residues are created as ``Residue`` objects and added
        to the protein.

        Parameters
        ----------
        pdb_path : str
            Path to the input PDB file containing the protein structure.
        rsa_path : str
            Path to the NACCESS ``.rsa`` file containing residue SASA values.
        """
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
        """
        Calculate the geometric center of the protein.

        The center is calculated as the mean of the x, y, and z coordinates
        of all residues in the protein.

        Returns
        -------
        tuple of float
            Coordinates of the center of mass as ``(com_x, com_y, com_z)``.
        """
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
        """
        Center the protein coordinates around its geometric center.

        The geometric center of the protein is first calculated using
        ``calc_com``. The center coordinates are then subtracted from
        the coordinates of every residue so that the protein is centered
        at the origin.
        """
        com_x, com_y, com_z = self.calc_com()

        for res in self.list_res: 
            res.x -= com_x
            res.y -= com_y
            res.z -= com_z

    def get_arrays(self, min_sasa=15.0):
        """
        Convert residue data into NumPy arrays for numerical computations.

        Only residues with a SASA greater than or equal to ``min_sasa``
        are selected. If no residue satisfies this threshold, all residues
        are used as a fallback.

        Parameters
        ----------
        min_sasa : float, optional
            Minimum SASA value required for a residue to be included.
            Default is 15.0 Å².

        Returns
        -------
        coords : numpy.ndarray
            Array of residue coordinates with shape ``(n, 3)``.
        hydros : numpy.ndarray
            One-dimensional array containing the hydrophobicity value
            of each selected residue.
        """

        filtered = [r for r in self.list_res if r.sasa >= min_sasa]
        if not filtered:
            filtered = self.list_res  # Secours si SASA n'est pas remplie

        coords = np.array([[r.x, r.y, r.z] for r in filtered])
        hydros = np.array([r.hydrophobicity for r in filtered])
        return coords, hydros

    @staticmethod
    def compute_rotation_matrix(v1, v2=np.array([0.0, 0.0, 1.0])):
        """
        Compute a rotation matrix that aligns one vector with another.

        The input vectors are normalized before calculating the rotation.
        The rotation matrix is constructed using the cross product and
        Rodrigues' rotation formula. By default, the target vector is
        the canonical Z-axis.

        Parameters
        ----------
        v1 : numpy.ndarray
            Three-dimensional vector to be aligned.
        v2 : numpy.ndarray, optional
            Three-dimensional target vector. Default is
            ``[0.0, 0.0, 1.0]``.

        Returns
        -------
        numpy.ndarray
            A 3 × 3 rotation matrix that aligns ``v1`` with ``v2``.
            If the vectors are already aligned, the identity matrix is
            returned. If they are opposite, the negative identity matrix
            is returned.
        """

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
        """
        Generate a reoriented PDB file and visualize the protein in PyMOL.

        The protein is rotated so that the selected membrane axis is aligned
        with the Z-axis. The protein is then shifted along the Z-axis so that
        the membrane center is positioned at zero. A PyMOL script is generated
        to display the protein together with the two planes representing the
        membrane boundaries.

        Parameters
        ----------
        best_axis : numpy.ndarray
            Three-dimensional vector corresponding to the optimal membrane
            orientation.
        z_center : float
            Position of the center of the membrane along the selected axis.
        memb_thickness : float, optional
            Thickness of the membrane in Å. Default is 30 Å.
        """
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
            subprocess.run(["pymol", pml_path_out])
        except FileNotFoundError:
            print(
                f"Erreur : Impossible d'ouvrir PyMOL avec la commande indiquée "
            )
