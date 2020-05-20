"""
Calculate the ply properties from the properties of the constituents.

A.J.J. Lagerweij
COHMAS Mechanical Engineering KAUST
2020
"""

# Import external packages.
import numpy as np


# Convert mass fractions into volume fractions.
def massfrac_to_volfrac(fm1, rho1, rho2):
    r"""
    Caluculate the volume fraction from the mass fraction.

    Parameters
    ----------
    fm1 : float
        The mass fraction of material 1 defined as, fm1 = massa 1 / massa
        total.
    rho1 : float
        The density of material 1.
    rho2 : float
        The density of material 2.

    Return
    ------
    fv1 : float
        The fraction material 1 is used in volume.
    fv2 : float
        The fraction material 2 is used in volume.
    """
    v1 = fm1 / rho1         # Volume occupied by material 1.
    v2 = (1 - fm1) / rho2   # Volume orrupied by material 2.
    fv1 = v1 / (v1 + v2)
    fv2 = v2 / (v1 + v2)
    return fv1, fv2


# Convert volume fraction into mass fractions.
def volfrac_to_massfrac(fv1, rho1, rho2):
    r"""
    Caluculate the mass fraction from the volume fraction.

    Parameters
    ----------
    fv1 : float
        The volume fraction of material 1 defined as, fv1 = volume 1 / volume
        total.
    rho1 : float
        The density of material 1.
    rho2 : float
        The density of material 2.

    Return
    ------
    fm1 : float
        The fraction material 1 is used in mass
    fm2 : float
        The fraction material 2 is used in mass.
    """
    fm1, fm2 = massfrac_to_volfrac(fv1, 1/rho1, 1/rho2)
    return fm1, fm2


# Orthotropic constitutive equation.
def orthotropic3D(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23):
    r"""
    Determine stiffness & compliance matrix of a 3D orthotropic material.

    This notation is in Voigt notation with engineering strain
    :math:`\gamma_{12}=2\varepsilon_{12}`.

    Parameters
    ----------
    E1 : float
        Young's modulus in 1 direction.
    E2 : float
        Young's modulus in 2 direction.
    E3 : float
        Young's modulus in 3 direction.
    nu12 : float
        Poisson's ratio over 12.
    nu13 : float
        Poisson's ratio over 13.
    nu23 : float
        Poisson's ratio over 23.
    G12 : float
        Shear modulus over 12.
    G13 : float
        Shear modulus over 13.
    G23 : float
        Shear modulus over 23.

    Returns
    -------
    C : matrix
        2D stiffness matrix in Voigt notaion (6x6).
    S : matrix
        2D compliance matrix in Voigt notation (6x6).
    """
    S = np.matrix([[1/E1, -nu12/E1, -nu13/E1, 0, 0, 0],
                   [-nu12/E1, 1/E2, -nu23/E2, 0, 0, 0],
                   [-nu13/E1, -nu23/E2, 1/E3, 0, 0, 0],
                   [0, 0, 0, 1/G12, 0, 0],
                   [0, 0, 0, 0, 1/G13, 0],
                   [0, 0, 0, 0, 0, 1/G23]])
    C = np.linalg.inv(S)
    return C, S


# Transversly isotropic constitutive equation.
def trans_isotropic3D(E1, E2, nu12, nu23, G12):
    r"""
    Determine stiffness & compliance matrix of 3D transverse isotropic material.

    This notation is in Voigt notation with engineering strain
    :math:`\gamma_{12}=2\varepsilon_{12}`.

    Parameters
    ----------
    E1 : float
        Young's modulus in 1 direction.
    E3 : float
        Young's modulus in 3 direction.
    nu12 : float
        Poisson's ratio over 12.
    nu23 : float
        Poisson's ratio over 23.
    G13 : float
        Shear modulus over 13.

    Returns
    -------
    C : matrix
        2D stiffness matrix in Voigt notaion (6x6).
    S : matrix
        2D compliance matrix in Voigt notation (6x6).
    """
    E3 = E2
    nu13 = nu12
    G23 = E2/(2*(1 + nu23))
    G13 = G12
    C, S = orthotropic3D(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23)
    return C, S


# Transversly isotropic constitutive equation in plane stress.
def trans_isotropic2D(E1, E2, nu12, G12):
    r"""
    Determine stiffness & compliance matrix of 2D plane stress transverse isotropic material.

    This notation is in Voigt notation with engineering strain
    :math:`\gamma_{12}=2\varepsilon_{12}`.

    Parameters
    ----------
    E1 : float
        Young's modulus in 1 direction.
    E2 : float
        Young's modulus in 3 direction.
    nu12 : float
        Poisson's ratio over 12.
    G12 : float
        Shear modulus over 13.

    Returns
    -------
    C : matrix
        2D stiffness matrix in Voigt notaion (3x3).
    S : matrix
        2D compliance matrix in Voigt notation (3x3).
    """
    S = np.matrix([[1/E1, -nu12/E1, 0],
                   [-nu12/E1, 1/E2, 0],
                   [0, 0, 1/(2*G12)]])
    C = np.linalg.inv(S)
    return C, S


# Isotropic constitutive equation.
def isotropic3D(E, nu):
    """
    Determine stiffness & compliance matrix of 3D isotropic material.

    Parameters
    ----------
    E : float
        Young's modulus.
    nu : float
        Poisson's ratio.

    Returns
    -------
    C : matrix
        3D stiffness matrix in Voigt notation (6x6).
    S : matrix
        3D compliance matrix in Voigt notation (6x6).
    """
    G = E/(2*(1 + nu))
    C, S = orthotropic3D(E, E, E, nu, nu, nu, G, G, G)
    return C, S


# Isotropic constitutive equation in plane stress.
def isotropic2D(E, nu):
    """
    Determine stiffness & compliance matrix of 2D isotropic material.

    Parameters
    ----------
    E : float
        Young's modulus.
    nu : float
        Poisson's ratio.

    Returns
    -------
    C : matrix
        3D stiffness matrix in Voigt notation (6x6).
    S : matrix
        3D compliance matrix in Voigt notation (6x6).
    """
    G = E/(2*(1 + nu))
    C, S = trans_isotropic2D(E, E, nu, G)
    return C, S


# Voigt limit of the rule of mixtures.
def voigt(C1, C2, volf):
    """
    Calculate the mixed stiffness matrix with the Voigt.

    This upper limit of the rule of mixtures is generaly used for the
    longitudional stiffness :math:`E_l`.

    Parameters
    ----------
    C1 : matrix
        The stiffness matrix of material 1.
    C2 : float
        The stiffness matrix of material 2.
    volf : matrix
        The volume fraction of material 1.

    Returns
    -------
    C_hat : matrix
        The stiffness matrix of the mixed material.
    S_hat : matrix
        The compliance matrix of the mixed material.
    """
    C_hat = volf*C1 + (1-volf)*C2
    S_hat = np.linalg.inv(C_hat)
    return C_hat, S_hat


# Reuss limit of the rule of mixtures.
def reuss(S1, S2, volf):
    """
    Calculate the compliance matrix with the Reuss limit.

    This lower limit of the rule of mixtures is generaly used for the
    transverse stiffness :math:`E_l`.

    Parameters
    ----------
    S1 : matrix
        The compliance matrix of material 1.
    S2 : float
        The compliance matrix of material 2.
    volf : matrix
        The volume fraction of material 1.

    Returns
    -------
    C_hat : matrix
        The stiffess matrix of the mixed material.
    S_hat : matrix
        The compliance matrix of the mixed material.
    """
    S_hat = volf*S1 + (1-volf)*S2
    C_hat = np.linalg.inv(S_hat)
    return C_hat, S_hat
