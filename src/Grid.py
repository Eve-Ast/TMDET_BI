import math
import numpy as np

"""
Grid module for scanning protein membrane orientations.

This module provides the ``Grid`` class, which generates a set of
candidate orientation vectors distributed over a unit sphere and uses
these vectors to analyze the hydrophobicity of a protein along different
axes.

Two approaches are implemented for identifying the optimal membrane
orientation. The first method searches for the region with the highest
mean hydrophobicity, while the second method maximizes the contrast
between the mean hydrophobicity inside and outside the membrane region.

Both non-vectorized and partially vectorized implementations are
provided to compare their computational performance.
"""

class Grid : 
    """
    A class to generate and representes a grid of vectors

    ...

    Attributes
    ----------
    n_points : int
        Number of vector generated.
    axis_list : list
        List of all the vector generated

    Methods
    -------
    _generate_fibonnaci_sphere : 
        Generation of the vector 
    compute_axis_profile : 
        Compute axis profil for the first methode of 
    scan_protein_method1_non_vectorized and scan protein_method1_vectorized : 
        first method of scanning the axes
    scan_protein_method2_non_vectorized and scan protein_method2_vectorized : 
        second method of scanning the axes
    """
    def __init__(self, n_points = 1000) :
        self.n_points = n_points
        self.axis_list = []
        self._generate_fibonnacci_sphere()

    def _generate_fibonnacci_sphere(self):
        """
        Generate a set of evenly distributed points on a unit sphere
        using the Fibonacci sphere algorithm.
        """
        phi = math.pi * (math.sqrt(5) - 1) # golden angle in radian

        for i in range (self.n_points) : 
            y = 1 - (i / float(self.n_points - 1)) * 2
            radius = math.sqrt(1 - y * y)

            theta = phi * i 

            x = math.cos(theta) * radius
            z = math.sin(theta) * radius


            self.axis_list.append(np.array([x, y, z]))

    
    def compute_axis_profile(self, protein, axis_vector, min_sasa=25.0):
        """
        Compute the hydrophobicity profile of a protein along a given axis.
        Parameters 
        ---------- 
        protein : Protein 
            Protein object containing the residues to analyze. 
        axis_vector : numpy.ndarray 
            Three-dimensional vector defining the axis along which the protein is analyzed. 
        min_sasa : float, optional 
            Minimum SASA value required for a residue to be included in the analysis. Default is 25.0 Å. 
        
        Returns 
        ------- 
        hydro_profile : list of float 
            Mean hydrophobicity calculated for each 1 Å interval along the selected axis. 
            Empty intervals are assigned a value of 0.0. 
        min_p : float 
            Minimum projection value along the normalized axis. 
            Returns 0.0 if no residue satisfies the SASA threshold.
        """
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

    def scan_protein_meth_1_non_vectorized(self, protein, min_sasa=25.0, memb_thickness=30):
        """
        Scan a protein over all candidate axes to identify the best membrane orientation using a non-vectorized approach.
        
        Parameters 
        ---------- 
        protein : Protein 
            Protein object containing the residues to analyze. 
        min_sasa : float, optional 
            Minimum SASA value required for a residue to be included in the hydrophobicity profile. Default is 25.0. 
        memb_thickness : int, optional 
            Thickness of the membrane in Å, corresponding to the size of the sliding window used to calculate the hydrophobicity score.
            Default is 30 Å. 
        
        Returns 
        ------- 
        dict 
            Dictionary containing the results of the scan: 
            ``best_score`` Highest mean hydrophobicity score found. 
            ``best_axis`` Axis vector corresponding to the highest score. 
            ``z_center`` Position of the center of the optimal membrane region along the selected axis.
        """
        
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

    def compute_axis_profile_vectorized(self, protein, axis_vector, min_sasa=25.0):
        """
        Compute the hydrophobicity profile of a protein along a given axis using a vectorized NumPy implementation.

        Parameters 
        ----------
        protein : Protein 
            Protein object containing the residue coordinates, SASA values, and hydrophobicity values. 
        axis_vector : numpy.ndarray 
            Three-dimensional vector defining the axis along which the protein is analyzed. 
        min_sasa : float, optional 
            Minimum SASA value required for a residue to be included in the analysis. Default is 25.0. 
            
        Returns 
        ------- 
        profile : numpy.ndarray 
            Hydrophobicity profile calculated for each 1 Å interval along the selected axis. 
            Empty intervals are assigned a value of 0.0. 
        min_p : float 
            Minimum projected coordinate along the normalized axis. Returns 0.0 if no residue satisfies the SASA threshold.
        """
        
        axis_norm = axis_vector / np.linalg.norm(axis_vector)
        coords, hydros = protein.get_arrays(min_sasa=min_sasa)

        if len(coords) == 0:
            return np.array([]), 0.0

        projections = np.dot(coords, axis_norm)
        min_p = math.floor(np.min(projections))
        max_p = math.ceil(np.max(projections))
        bins = np.arange(min_p, max_p + 1, 1.0)
        n_bins = max(len(bins) - 1, 1)

        bin_idx = np.digitize(projections, bins) - 1
        bin_idx = np.clip(bin_idx, 0, n_bins - 1)

        sums = np.bincount(bin_idx, weights=hydros, minlength=n_bins)
        counts = np.bincount(bin_idx, minlength=n_bins)

        with np.errstate(invalid="ignore", divide="ignore"):
            profile = np.where(counts > 0, sums / counts, 0.0)

        return profile, float(min_p)

    def scan_protein_meth_1_vectorized(self, protein, min_sasa=25.0, memb_thickness=30):
        """ 
        Scan a protein over all candidate axes to identify the best membrane orientation using a vectorized approach.

        Parameters 
        ---------- 
        protein : Protein 
            Protein object containing the residue coordinates, SASA values, and hydrophobicity values. 
        min_sasa : float, optional 
            Minimum SASA value required for a residue to be included in the hydrophobicity profile. Default is 25.0. 
        memb_thickness : int, optional 
            Thickness of the membrane in Å, corresponding to the size of the sliding window used to calculate the hydrophobicity score. 
            Default is 30 Å. 
        
        Returns 
        ------- 
        dict
          Dictionary containing the results of the scan:
            ``best_score`` Highest mean hydrophobicity score found among all candidate membrane regions. 
            ``best_axis`` Normalized axis vector corresponding to the highest score. 
            ``z_center`` Position of the center of the optimal membrane region along the selected axis.
        """

        best_score = -float("inf")
        best_axis = None
        best_z_shift = 0.0
        w = memb_thickness

        for axis in self.axis_list:
            profile, min_p = self.compute_axis_profile_vectorized(protein, axis, min_sasa=min_sasa)
            if len(profile) < w:
                continue

            cumsum = np.cumsum(np.insert(profile, 0, 0.0))
            window_sums = cumsum[w:] - cumsum[:-w]
            window_means = window_sums / w

            i = int(np.argmax(window_means))
            score = window_means[i]
            if score > best_score:
                best_score = score
                best_axis = axis / np.linalg.norm(axis)
                best_z_shift = min_p + i + (w / 2.0)

        return {"best_score": best_score, "best_axis": best_axis, "z_center": best_z_shift}

    def scan_protein_meth_2_non_vectorized(self, protein, min_sasa=25.0, memb_thickness=30.0):
        """ 
        Scan a protein over all candidate axes to identify the best membrane orientation using a non-vectorized approach. 
        
        For each candidate axis, the residues satisfying the SASA threshold are projected onto the normalized axis. 
        A sliding membrane region centered at different positions along the axis is then considered. 
        For each position, the mean hydrophobicity inside the membrane region is compared with the mean hydrophobicity 
        outside the membrane region. The difference between these two values is used as the score.

        Parameters 
        ---------- 
        protein : Protein 
            Protein object containing the residue coordinates, SASA values, and hydrophobicity values. 
        min_sasa : float, optional 
            Minimum SASA value required for a residue to be included in the analysis. Default is 25.0. 
        memb_thickness : float, optional 
            Thickness of the membrane in Å. The membrane region is defined as a symmetric interval of half-thickness around its center. 
            Default is 30.0 Å. 
            
        Returns 
        ------- 
        dict 
            Dictionary containing the results of the scan: 
                ``best_score`` Highest difference between the mean hydrophobicity inside and outside the membrane region. 
                ``best_axis`` Normalized axis vector corresponding to the highest score. 
                ``z_center`` Position of the center of the optimal membrane region along the selected axis.
        """
        half_thick = memb_thickness / 2.0
        filtered = [r for r in protein.list_res if r.sasa >= min_sasa]
        if not filtered:
            return {"best_score": -float("inf"), "best_axis": (0.0, 0.0, 1.0), "z_center": 0.0}

        best_score = -float("inf")
        best_axis = None
        best_z_shift = 0.0

        for axis in self.axis_list:
            norm = math.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
            ax, ay, az = axis[0]/norm, axis[1]/norm, axis[2]/norm

            projections = [r.x*ax + r.y*ay + r.z*az for r in filtered]
            min_p, max_p = min(projections), max(projections)

            z_center = min_p
            while z_center <= max_p:
                hydro_in, hydro_out = [], []
                for proj, r in zip(projections, filtered):
                    if z_center - half_thick <= proj <= z_center + half_thick:
                        hydro_in.append(r.hydrophobicity)
                    else:
                        hydro_out.append(r.hydrophobicity)

                if hydro_in and hydro_out:
                    score = (sum(hydro_in)/len(hydro_in)) - (sum(hydro_out)/len(hydro_out))
                    if score > best_score:
                        best_score = score
                        best_axis = (ax, ay, az)
                        best_z_shift = z_center

                z_center += 1.0

        return {"best_score": best_score, "best_axis": best_axis, "z_center": best_z_shift}


    def scan_protein_meth_2_vectorized(self, protein, min_sasa=25.0, memb_thickness=30.0):
        """
        Scan a protein over all candidate axes to identify the best membrane orientation using a partially vectorized approach.

        Parameters 
        ---------- 
        protein : Protein 
            Protein object containing the residue coordinates, SASA values, and hydrophobicity values. 
        min_sasa : float, optional 
            Minimum SASA value required for a residue to be included in the analysis. Default is 25.0. 
        memb_thickness : float, optional 
            Thickness of the membrane in Å. The membrane region is defined as a symmetric interval around its center. Default is 30.0 Å. 
        
        Returns 
        ------- 
        dict 
            Dictionary containing the results of the scan: 
                ``best_score`` Highest difference between the mean hydrophobicity inside and outside the membrane region. 
                ``best_axis`` Normalized axis vector corresponding to the highest score. 
                ``z_center`` Position of the center of the optimal membrane region along the selected axis.
        """
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

    def __str__(self): 
        if self.vector : 
            message = f"Grid with {len(self.vector)} direction vectors generated."
        else : 
            message = "No vectors"
        return message
