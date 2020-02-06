"""
Example script for composits laminaties and their homoginization.

This file calculates the properties of simple composit plates.
The stacking sequence should be defined starting from the top of the laminate,
which is the negative :math:`z` direction.

A.J.J. Lagerweij
COHMAS Mechanical Engineering KAUST
2020
"""

# import packages
import numpy as np

# import local packages
import abdCal
import deformation
import failure


###############################################################################
# Ply Properties                                                              #
###############################################################################
# Elastic properties
El = 142 * 1e3  # MPa
Et = 13 * 1e3  # MPa
G = 5 * 1e3  # MPa
nult = 0.3  # -

# Failure properties
Xt = 2200  # MPa
Xc = 1850  # MPa
Yt = 55  # MPa
Yc = 200  # MPa
Smax = 120  # MPa

# Others
t = 0.16  # mm

# Calculate the ply stiffness matricess matrix
Q = abdCal.QPlaneStress(El, Et, nult, G)


###############################################################################
# Laminate Properites                                                         #
###############################################################################
# Stacking sequence
angles_deg = [0, 90, 90, 0]
thickness = [t] * len(angles_deg)
Q = [Q] * len(angles_deg)

# Calculate the ABD matrix and its inverse
abd = abdCal.abd(Q, angles_deg, thickness)
abd_inv = abdCal.abd_inverse(abd)


###############################################################################
# Applied Running Loads                                                       #
###############################################################################
# Calculate the deformation caused by a given running load
NM = np.matrix([0, 1, 0, 1, 0, 0]).T  # MPa/mm and MPa*mm/mm
deformed = deformation.load_applied(abd_inv, NM)

# Calculate the stress in each layer caused by the running loads
strain = deformation.ply_strain(deformed, angles_deg, thickness, Q)
stress = deformation.ply_stress(deformed, angles_deg, thickness, Q, plotting=True)


###############################################################################
# Test Ply Failure with various criteria                                      #
###############################################################################
# testing if the the criteria is violated
failure.max_stress(stress, Xt, Xc, Yt, Yc, Smax)
failure.tsai_wu(stress, Xt, Xc, Yt, Yc, Smax)
failure.tsai_hill(stress, Xt, Xc, Yt, Yc, Smax)

