"""
Example script for composits laminaties and their homoginization.

This file calculates the properties of simple composit plates.
The stacking sequence should be defined starting from the top of the laminate,
which is the negative :math:`z` direction.

A.J.J. Lagerweij
COHMAS Mechanical Engineering KAUST
2020
"""

# Import external packages.
import numpy as np

# Import local packages.
import homogenization
import abdcal
import deformation
import failure


###############################################################################
# Ply Properties                                                              #
###############################################################################
# List the elastic properties of the ply.
El = 142 * 1e3  # MPa
Et = 13 * 1e3  # MPa
G = 5 * 1e3  # MPa
nult = 0.3  # -

# List the failure properties of the ply.
Xt = 2200  # MPa
Xc = 1850  # MPa
Yt = 55  # MPa
Yc = 200  # MPa
Smax = 120  # MPa

# List the other properties of the ply.
t = 0.16  # mm

# Calculate the ply stiffness matricess matrix.
Q = abdcal.QPlaneStress(El, Et, nult, G)


###############################################################################
# Laminate Properites                                                         #
###############################################################################
# Define the stacking sequence.
angles_deg = [0, 0, 45, 90, -45, -45, 90, 45, 0, 0]
thickness = [t] * len(angles_deg)
Q = [Q] * len(angles_deg)

# Calculate the ABD matrix and its inverse.
abd = abdcal.abd(Q, angles_deg, thickness)
abd_inv = abdcal.matrix_inverse(abd)


###############################################################################
# Applied Running Loads                                                       #
###############################################################################
# Calculate the deformation caused by a given running load.
NM = np.matrix([0, 1, 0, 1, 0, 0]).T  # MPa/mm and MPa*mm/mm
deformed = deformation.load_applied(abd_inv, NM)

# Calculate the stress in each layer caused by the running loads.
strain = deformation.ply_strain(deformed, Q, angles_deg, thickness)
stress = deformation.ply_stress(deformed, Q, angles_deg, thickness, plotting=True)


###############################################################################
# Test Ply Failure with various criteria                                      #
###############################################################################
# Testing whether the failure criterias are violated.
failure.max_stress(stress, Xt, Xc, Yt, Yc, Smax)
failure.tsai_wu(stress, Xt, Xc, Yt, Yc, Smax)
failure.tsai_hill(stress, Xt, Xc, Yt, Yc, Smax)
