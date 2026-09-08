import math
import numpy as np


class Grid : 
    """"""
    def __init__(self, n_points = 1000) :
        """"""
        self.n_points = n_points
        self.axis_list = []
        self._generate_fibonnacci_sphere()

    def _generate_fibonnacci_sphere(self):
        """"""
        phi = math.pi * (math.sqrt(5) - 1) # golden angle in radian

        for i in range (self.n_points) : 
            y = 1 - (i / float(self.n_points - 1)) * 2
            radius = math.sqrt(1 - y * y)

            theta = phi * i 

            x = math.cos(theta) * radius
            z = math.sin(theta) * radius


            self.axis_list.append(np.array([x, y, z]))

    def scan_protein_meth_2(self, protein, min_sasa=15.0, memb_thickness=30.0):
        """Scanne la protéine vectoriellement avec un score d'octanol/eau ou de contraste."""
        coords, hydros = protein.get_arrays(min_sasa=min_sasa)

        if len(coords) == 0:
            return {
                "best_score": -float("inf"),
                "best_axis": np.array([0, 0, 1]),
                "z_center": 0.0,
            }

        half_thick = memb_thickness / 2.0

        best_score = -float("inf")
        best_axis = None
        best_z_shift = 0.0

        # Balayage vectorisé de tous les axes
        for axis in self.axis_list:
            axis_norm = axis / np.linalg.norm(axis)

            # Projection de tous les résidus sur l'axe en 1 seule opération matricielle
            projections = np.dot(coords, axis_norm)

            min_p, max_p = np.min(projections), np.max(projections)

            # Balayage par pas de 1 Å du centre de la membrane
            for z_center in np.arange(min_p, max_p, 1.0):
                # Masque booléen : True si dans le bloc membranaire de 30 Å
                in_membrane = (projections >= z_center - half_thick) & (
                    projections <= z_center + half_thick
                )

                # Score de contraste :
                # + Hydrophobicité moyenne DANS la membrane
                # - Hydrophobicité moyenne HORS de la membrane
                hydro_in = hydros[in_membrane]
                hydro_out = hydros[~in_membrane]

                if len(hydro_in) == 0 or len(hydro_out) == 0:
                    continue

                score = np.mean(hydro_in) - np.mean(hydro_out)

                if score > best_score:
                    best_score = score
                    best_axis = axis_norm
                    best_z_shift = z_center

        return {
            "best_score": best_score,
            "best_axis": best_axis,
            "z_center": best_z_shift,
        }
    
    def compute_axis_profile(self, protein, axis_vector, min_sasa=30.0):
        projections = []
        hydrophobicities = []

        # Normalisation du vecteur axe
        axis_norm = axis_vector / np.linalg.norm(axis_vector)

        for res in protein.list_res:
            if res.sasa >= min_sasa:
                # Produit scalaire = projection sur l'axe
                proj = (
                    res.x * axis_norm[0]
                    + res.y * axis_norm[1]
                    + res.z * axis_norm[2]
                )
                projections.append(proj)
                hydrophobicities.append(res.hydrophobicity)

        if not projections:
            return [], 0.0

        projections = np.array(projections)
        hydrophobicities = np.array(hydrophobicities)

        min_p = float(np.min(projections))
        max_p = float(np.max(projections))

        # Tranches de 1 Angstrom
        bins = np.arange(math.floor(min_p), math.ceil(max_p) + 1, 1.0)
        hydro_profile = []

        for i in range(len(bins) - 1):
            mask = (projections >= bins[i]) & (projections < bins[i + 1])
            if np.any(mask):
                hydro_profile.append(np.mean(hydrophobicities[mask]))
            else:
                hydro_profile.append(0.0)

        return hydro_profile, min_p

    def scan_protein_meth_1(self, protein, min_sasa=15.0, memb_thickness=30):
        best_score = -float("inf")
        best_axis = None
        best_z_shift = 0.0

        for axis in self.axis_list:
            hydro_profile, min_p = self.compute_axis_profile(
                protein, axis, min_sasa=min_sasa
            )

            if len(hydro_profile) < memb_thickness:
                continue

            for i in range(len(hydro_profile) - memb_thickness + 1):
                slice_30 = hydro_profile[i : i + memb_thickness]
                score = np.mean(slice_30)

                if score > best_score:
                    best_score = score
                    best_axis = axis
                    # Décalage réel par rapport à l'origine (0,0,0)
                    best_z_shift = min_p + i + (memb_thickness / 2.0)

        return {
            "best_score": best_score,
            "best_axis": best_axis,
            "z_center": best_z_shift,
        }

    def __str__(self): 
        if self.vector : 
            message = f"Grid with {len(self.vector)} direction vectors generated."
        else : 
            message = "No vectors"
        return message
