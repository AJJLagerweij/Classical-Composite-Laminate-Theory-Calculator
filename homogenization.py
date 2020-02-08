"""
Calculate the ply properties from the properties of the constituents.

A.J.J. Lagerweij
COHMAS Mechanical Engineering KAUST
2020
"""


# Convert mass fractions into volume fractions.
def massfrac_to_volfrac(fm1, rho1, rho2):
    """
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
    """
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


# Voigt limit of the rule of mixtures.
def rule_of_mixtures(prop_f, prop_m, vol_f):
    """
    Calculate the mixed property with the Voight limit of the rule of mixtures.

    This upper limite of the rule of mixtures is generaly used for the
    longitudional stiffness :math:`E_l`.

    Parameters
    ----------
    prop_f : float
        The property of the fiberous material.
    prop_m : float
        The property of the matrix material.
    vol_f : float
        The volume fraction of the fiberous material.

    Returns
    -------
    prop : float
        The property of the mixed material.
    """
    prop = False

    return prop
