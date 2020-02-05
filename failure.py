"""
Calculate the failure (critera) of the composite material.

A.J.J. Lagerweij
COHMAS Mechanical Engineering KAUST
2020
"""

import numpy as np


# Max stress failure (No interaction)
def max_stress(stress, sl_max, sl_min, st_max, st_min, tlt_max):
    """
    Compare stresses to the max stress criteria, returns which layers failed.

    Failure stresses and ply stresses must be in the same orientation.

    .. warning::
        This method does not take interaction of stresses in different
        directions in account.

    Parameters
    ----------
    stress : list
        A list containing the stress vector of the bottom, middle
        and top of each layer.
    sl_max : float
        Maximum tensile stress in longitudional direction.
    sl_min : float
        Maximum compressive stress in longitudional direction.
    st_max : float
        Maximum tensile stress in transverse direction.
    st_min : float
        Maximum compressive stress in transverse direction.
    tlt_max : float
        Maximum shear stress in material axis system.

    Returns
    -------
    pass : bool
        True if the load was below the maximum allowables.
    """
    # Location inside each ply
    loc = ['top', 'bottom']

    # Make arrays with only stress in a certain direction
    sl = np.array(stress)[:, :, 0, 0]
    st = np.array(stress)[:, :, 1, 0]
    tlt = np.array(stress)[:, :, 2, 0]

    # Check if anything fails, value is True if it is outside admissible
    sl_tens = np.any(sl >= sl_max)
    sl_comp = np.any(-sl_min >= sl)
    st_tens = np.any(st >= st_max)
    st_comp = np.any(-st_min >= st)
    shear = np.any(np.abs(tlt) >= tlt_max)

    if any([sl_tens, sl_comp, st_tens, st_comp, shear]):
        print("The max stress criteria was violated at:")

        if np.max(sl) >= sl_max:
            error = np.unravel_index(np.argmax(sl, axis=None), sl.shape)
            print('    - Max longitudional tension at the', loc[error[1]], 'of Layer', error[0]+1)
        if np.min(sl) <= -sl_min:
            error = np.unravel_index(np.argmin(sl, axis=None), sl.shape)
            print('    - Max longitudional compression at the', loc[error[1]], 'of Layer', error[0]+1)
        if np.max(st) >= st_max:
            error = np.unravel_index(np.argmax(st, axis=None), st.shape)
            print('    - Max transverse tension at the', loc[error[1]], 'of Layer', error[0]+1)
        if np.min(st) <= -st_min:
            error = np.unravel_index(np.argmin(st, axis=None), st.shape)
            print('    - Max transverse compression at the', loc[error[1]], 'of Layer', error[0]+1)
        if np.max(np.abs(tlt)) >= tlt_max:
            error = np.unravel_index(np.argmax(tlt, axis=None), tlt.shape)
            print('    - Max Shear failure at the', loc[error[1]], 'Layer', error[0]+1)

        return False

    return True


# Tsai-Wu Criteria (Bad compression approximation)
def tsai_wu(stress, sl_max, sl_min, st_max, st_min, tlt_max):
    """
    Test wether the stresses are outside the plane stress Tsai-Wu criteria.

    Failure stresses and ply stresses must be in the same orientation.

    .. warning::
        This criteria is a bad approximation for compression failure.

    Parameters
    ----------
    stress : list
        A list containing the stress vector of the bottom, middle
        and top of each layer.
    sl_max : float
        Maximum tensile stress in longitudional direction.
    sl_min : float
        Maximum compressive stress in longitudional direction.
    st_max : float
        Maximum tensile stress in transverse direction.
    st_min : float
        Maximum compressive stress in transverse direction.
    tlt_max : float
        Maximum shear stress in material axis system.

    Returns
    -------
    pass : bool
        True if the load was below the maximum allowables.
    """
    # Location inside each ply
    loc = ['top', 'bottom']

    # Make arrays with only stress in a certain direction
    sl = np.array(stress)[:, :, 0, 0]
    st = np.array(stress)[:, :, 1, 0]
    tlt = np.array(stress)[:, :, 2, 0]

    # Calculating constants for Tsai-Wu Criteria
    f1 = 1 / sl_max - 1 / sl_min
    f11 = 1 / (sl_min * sl_min)
    f2 = 1 / st_max - 1 / st_min
    f22 = 1 / (st_max * st_min)
    f66 = 1 / tlt_max**2
    f12 = -1/2 * np.sqrt(f11 * f22)  # This is an approximation

    # Appliing the criteria itself
    Criteria = -f1*sl + f2*st + f11*sl**2 + f22*st**2 + f66*tlt**2 + 2*f12*sl*st

    # Returning list with plies failing
    if np.max(Criteria) > 1:
        error = np.where(Criteria > 1)
        print("The Tsai-Wu criteria was violated at:")
        for i in range(len(error[0])):
            print("    - The", loc[error[1][i]], "of layer", error[0][i]+1)
        return False

    return True


# Tsai-Hill Criteria for plane stress(Probramming should be speed up...)
def tsai_hill(stress, sl_max, sl_min, st_max, st_min, tlt_max):
    """
    Test wether the stresses are outside the plane stress Tsai-Hill criteria.

    Failure stresses and ply stresses must be in the same orientation.
    The method is an extension of the von Mieses stress criteria.

    Parameters
    ----------
    stress : list
        A list containing the stress vector of the bottom, middle
        and top of each layer.
    sl_max : float
        Maximum tensile stress in longitudional direction.
    sl_min : float
        Maximum compressive stress in longitudional direction.
    st_max : float
        Maximum tensile stress in transverse direction.
    st_min : float
        Maximum compressive stress in transverse direction.
    tlt_max : float
        Maximum shear stress in material axis system.

    Returns
    -------
    pass : bool
        True if the load was below the maximum allowables.
    """
    # Location inside each ply
    loc = ['top', 'bottom']

    # Make arrays with only stress in a certain direction
    sl = np.array(stress)[:, :, 0, 0]
    st = np.array(stress)[:, :, 1, 0]
    tlt = np.array(stress)[:, :, 2, 0]

    # Calculating compression and tension maxima:
    slm = ((np.sign(sl)+1)/(2))*sl_max + ((np.sign(sl)-1)/(2))*sl_min
    stm = ((np.sign(st)+1)/(2))*st_max + ((np.sign(st)-1)/(2))*st_min

    # Applying the criteria itelf
    Criteria = (sl/slm)**2 + (st/stm)**2 - (sl*st)/(slm**2)+(tlt/tlt_max)**2

    # Returning list with plies failing
    if np.max(Criteria) > 1:
        error = np.where(Criteria > 1)
        print("The Tsai-Hill criteria was violated at:")
        for i in range(len(error[0])):
            print("    - The", loc[error[1][i]], "of layer", error[0][i]+1)
        return False

    return True
