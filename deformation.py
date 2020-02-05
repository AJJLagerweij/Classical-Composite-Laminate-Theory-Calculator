"""
Calculate the resulting deformation or stresses for loading conditions.

A.J.J. Lagerweij
COHMAS Mechanical Engineering KAUST
2020
"""

# Importing packages
import numpy as np
import matplotlib.pyplot as plt


# Calculating the strain in the plate and curvature (Kirchhoff plate theory)
def load_applied(abd_inv, load):
    r"""
    Calculate the strain and curvature of the full plate under a given load using Kirchhoff plate theory.

    Parameters
    ----------
    abd : matrix
        The inverse of the ABD matrix.
    load : vector
        The load vector consits of are :math:`(N_x, N_y, N_{xy}, M_x, M_y, M_{xy})^T`

    Returns
    -------
    deformation : vector
        This deformation consists of :math:`(\varepsilon_x, \varepsilon_y
        \varepsilon_{xy},\kappa_x, \kappa_y, \kappa_{xy})^T`
    """
    deformation = abd_inv.dot(load)
    return deformation


# Calculating the running force & moment on the plate (Kirchhoff plate theory)
def deformation_applied(abd, deformation):
    r"""
    Calculate the running load and moment of the plate under a given using Kichhoff plate theory.

    Parameters
    ----------
    abd : matrix
        The ABD matrix.
    deformation : vector
        This deformation consists of :math:`(\varepsilon_x, \varepsilon_y,
        \varepsilon_{xy},\kappa_x, \kappa_y, \kappa_{xy})^T`

    Returns
    -------
    load : vector
        The load vector consits of are :math:`(N_x, N_y, N_{xy}, M_x, M_y, M_{xy})^T`
    """
    load = abd.dot(deformation)
    return load


# Calculate the strain in each ply
def ply_strain(deformed, angles, thickness, Q):
    r"""
    Calculate the strain at the top and bottom of each ply.

    Small and linear deformations are assumed. For each ply two strain states
    are retured, one for the top and bottom of each ply. As bending moments can
    leed to different stresses depending on the :math:`z` location in the ply.
    Top plies should be listed first in the lists of Q, angles and thickness..


    Parameters
    ----------
    deformed : vector
        This deformation consists of :math:`(\varepsilon_x, \varepsilon_y
        \varepsilon_{xy},\kappa_x, \kappa_y, \kappa_{xy})^T`
    angles : list
        The rotation of each ply in degrees.
    thickness : list
        The thickness of each ply.
    Q : list
        The local stiffness matrix of each ply.
    plotting : bool, optional
        Plotting the stress distribution or not.

    Returns
    -------
    strain : list
        The stress vector :math:`(\sigma_{xx}, \sigma_{yy}, \tau_{xy})^T` of
        the top, middle and bottom of each ply in the laminate.
    """
    # Calculating total thickness
    h = np.sum(thickness) / 2

    # Setting up a list for the output
    strain = []

    # iterate over all plies
    for i in range(len(thickness)):
        # Calculating z coordinates around the ply
        z_top = np.sum(thickness[:i]) - h
        z_bot = np.sum(thickness[:i+1]) - h

        # Deformation of the midplane of the plate
        strain_membrane = deformed[:3, 0]
        curvature = deformed[3:, 0]

        # Caluclate strain in plies
        strain_top = strain_membrane + z_top * curvature
        strain_bot = strain_membrane + z_bot * curvature

        # Rotating strain from global to ply axis sytstem
        strain_lt_top = strain_rotation(strain_top, -angles[i])
        strain_lt_bot = strain_rotation(strain_bot, -angles[i])

        # Store the stress and z coordinate results
        strain_ply = [strain_lt_top, strain_lt_bot]
        strain.append(strain_ply)

    return strain


# Calculate the stress in each ply
def ply_stress(deformed, angles, thickness, Q, plotting=False):
    r"""
    Calculate the stresses at the top and bottom of each ply.

    Small and linear deformations are assumed. For each ply two stress states
    are retured, one for the top and bottom of each ply. As bending moments can
    leed to different stresses depending on the :math:`z` location in the ply.
    Top plies should be listed first in the lists of Q, angles and thickness.

    Parameters
    ----------
    deformed : vector
        This deformation consists of :math:`(\varepsilon_x, \varepsilon_y
        \varepsilon_{xy},\kappa_x, \kappa_y, \kappa_{xy})^T`
    angles : list
        The rotation of each ply in degrees.
    thickness : list
        The thickness of each ply.
    Q : list
        The local stiffness matrix of each ply.
    plotting : bool, optional
        Plotting the through thickness stress distribution or not.

    Returns
    -------
    stress : list
        The stress vector :math:`(\sigma_{xx}, \sigma_{yy}, \tau_{xy})^T` of
        the top, middle and bottom of each ply in the laminate.
    """
    # Calculating the strain in local ply axes system
    strain = ply_strain(deformed, angles, thickness, Q)

    # Setting up a list for the output
    stress = []

    # empty stress list for when plotting is required
    if plotting is True:
        h = np.sum(thickness)/2
        z = []
        sigma_xx = []
        sigma_yy = []
        tau_xy = []

    # iterate over all plies
    for i in range(len(thickness)):
        # Rotating strain from global to ply axis sytstem
        strain_lt_top = strain[i][0]
        strain_lt_bot = strain[i][1]

        # Now this is converted into stress
        stress_lt_top = Q[i].dot(strain_lt_top)
        stress_lt_bot = Q[i].dot(strain_lt_bot)

        # store the stress and z coordinate results
        stress_ply = [stress_lt_top, stress_lt_bot]
        stress.append(stress_ply)

        # for the plotting
        if plotting is True:
            z_top = np.sum(thickness[:i]) - h
            z_bot = np.sum(thickness[:i+1]) - h
            z += [z_top, z_bot]
            sigma_xx += [stress_lt_top.item(0), stress_lt_bot.item(0)]
            sigma_yy += [stress_lt_top.item(1), stress_lt_bot.item(1)]
            tau_xy += [stress_lt_top.item(2), stress_lt_bot.item(2)]

    if plotting is True:
        plt.subplot(1, 3, 1)
        plt.plot(sigma_xx, z)
        plt.xlabel(r"Stress $\sigma_{xx}$")
        plt.ylabel(r"Location along $z$")
        plt.ylim((min(z), max(z)))
        ax = plt.gca()
        ax.invert_yaxis()

        plt.subplot(1, 3, 2)
        plt.plot(sigma_yy, z)
        plt.xlabel(r"Stress $\sigma_{yy}$")
        plt.ylim((min(z), max(z)))
        ax = plt.gca()
        ax.invert_yaxis()

        plt.subplot(1, 3, 3)
        plt.plot(tau_xy, z)
        plt.xlabel(r"Stress $\tau_{xy}$")
        plt.ylim((min(z), max(z)))
        ax = plt.gca()
        ax.invert_yaxis()

        plt.tight_layout()
        plt.show()

    return stress


# Stress vertor rotation
def stress_rotation(stress, angle):
    """
    Rotates a stress vector over a given angle.

    The stress vector must be in Voigt notation and engineering stress is used.

    Parameters
    ----------
    stress : vector
        The matrix that must be rotated.
    angle : float
        The rotation angle in degrees.

    Returns
    -------
    stress_rot : vector
        A rotated version of the matrix.
    """
    angle = angle * np.pi/180  # convert to radians
    m = np.cos(angle)
    n = np.sin(angle)
    T1 = np.matrix([[m**2, n**2, 2*m*n],
                    [n**2, m**2, -2*m*n],
                    [-m*n, m*n, m**2-n**2]])
    stress_rot = T1 * stress
    return stress_rot


# Strain vertor rotation
def strain_rotation(strain, angle):
    """
    Rotates a strain vector over a given angle.

    The strain vector must be in Voigt notation and engineering strain is used.

    Parameters
    ----------
    strain : vector
        The matrix that must be rotated.
    angle : float
        The rotation angle in degrees.

    Returns
    -------
    strain_rot : vector
        A rotated version of the matrix.
    """
    angle = angle * np.pi/180  # convert to radians
    m = np.cos(angle)
    n = np.sin(angle)
    T2 = np.matrix([[m**2, n**2, m*n],
                    [n**2, m**2, -m*n],
                    [-2*m*n, 2*m*n, m**2-n**2]])
    strain_rot = T2 * strain
    return strain_rot
