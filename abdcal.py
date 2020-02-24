"""
ABD matrix calculator for a given stacking sequence.

The methods in this file will call create a ABD matrix of a composit for given
ply properties and stacking sequence.

A.J.J. Lagerweij
COHMAS Mechanical Engineering KAUST
2020
"""

# Import external packages.
import numpy as np


# Determine the stiffness matrix for an in plane stress ply.
def QPlaneStress(El, Et, nult, G):
    r"""
    Generate the plane stress local stiffness matrix.

    Parameters
    ----------
    El : float
        Elastic modulus in the longitudional direction.
    Et : float
        Elastic modulus in the transverse direction.
    nult : float
        Poisson ratio in longitudional-transverse direction.
    G : float
        Shear modulus in longitudional-transverse directions.

    Returns
    -------
    Q : matrix
        Stiffness matrix in longitudional-transverse directions.
    """
    nutl = nult*Et/El
    Q = np.matrix([[El/(1-nult*nutl), (nult*Et)/(1-nult*nutl), 0],
                   [(nutl*El)/(1-nult*nutl), Et/(1-nult*nutl), 0],
                   [0, 0, G]])
    return Q


# Determine the stiffness matrix for an in plane strain ply.
def QPlaneStrain(El, Et, nult, G):
    r"""
    Generate the plane strain local stiffness matrix.

    .. warning::
        Not yet implemented.
        It raises a `NotImplementedError`.

    Parameters
    ----------
    El : float
        Elastic modulus in the longitudional direction.
    Et : float
        Elastic modulus in the transverse direction.
    nult : float
        Poisson ratio in longitudional-transverse direction.
    G : float
        Shear modulus in longitudional-transverse directions.

    Returns
    -------
    Q : matrix
        Stiffness matrix in longitudional-transverse directions.
    """
    raise NotImplementedError(
        'QPlaneStrain is undefined. Use plane stress version QPlaneStress.')
    return 0


# Rotate the stiffness matrix over a given angle.
def stiffness_rotation(stiffness, angle):
    r"""
    Rotate the stiffness matrix over a given angle.

    Parameters
    ----------
    stiffness : matrix
        The matrix that must be rotated.
    angle : float
        The rotation angle in degrees.

    Returns
    -------
    stiffness_rot : matrix
        A rotated version of the matrix.
    """
    angle = angle * np.pi/180  # convert to radians
    m = np.cos(angle)
    n = np.sin(angle)
    T1 = np.matrix([[m**2, n**2, 2*m*n],
                    [n**2, m**2, -2*m*n],
                    [-m*n, m*n, m**2-n**2]])
    T2 = np.matrix([[m**2, n**2, m*n],
                    [n**2, m**2, -m*n],
                    [-2*m*n, 2*m*n, m**2-n**2]])
    stiffness_rot = np.linalg.inv(T1) * stiffness * T2
    return stiffness_rot


# Rotate the compliance matrix over a given angle.
def compliance_rotation(compliance, angle):
    r"""
    Rotate the compliance matrix over a given angle.

    Parameters
    ----------
    compliance : matrix
        The matrix that must be rotated.
    angle : float
        The rotation angle in degrees.

    Returns
    -------
    stiffness_rot : matrix
        A rotated version of the matrix.
    """
    angle = angle * np.pi/180  # convert to radians
    m = np.cos(angle)
    n = np.sin(angle)
    T1 = np.matrix([[m**2, n**2, 2*m*n],
                    [n**2, m**2, -2*m*n],
                    [-m*n, m*n, m**2-n**2]])
    T2 = np.matrix([[m**2, n**2, m*n],
                    [n**2, m**2, -m*n],
                    [-2*m*n, 2*m*n, m**2-n**2]])
    compliance_rot = np.linalg.inv(T2) * compliance * T1
    return compliance_rot


# Calculate the stiffness matrix of a thin (no benbing) laminate.
def abdthin(Q, angles, thickness):
    r"""
    ABD matrix calculator for a thin laminate.

    In the thin laminate theroy it is assumed that the out of plane stiffness
    is negligible. Hence only the membrane (A part) of the ABD matrix remains.
    Top plies should be listed first in the lists of Q, angles and thickness.

    Parameters
    ----------
    Q : list
        The stiffness matrix of each ply in its l-t axis system.
    angles : list
        The rotation of each ply in degrees.
    thickness : list
        The thickness of each ply.

    Returns
    -------
    C : matrix
        The stiffness matrix of the thin laminate.
    """
    # Create an empty stiffness matrix.
    C = np.zeros((3, 3))

    # Loop over all plies.
    for i in range(len(angles)):
        Q_bar = striffness_rotation(Q[i], angles[i]) * thickness[i]
        C += Q_bar

    # Go from the stresses to running loads.
    C = np.matrix(C/(np.sum(thickness)))
    return C


# Calculate the ABD matrix of a given laminate.
def abd(Q, angles, thickness):
    r"""
    Calculate the full ABD matrix of a laminate.

    Top plies should be listed first in the lists of Q, angles and thickness.

    Parameters
    ----------
    Q : list
        The stiffness matrix of each ply in its l-t axis system.
    angles : list
        The rotation of each ply in degrees.
    thickness : list
        The thickness of each ply.

    Returns
    -------
    ABD : matrix
        The stiffness matrix of the thin laminate.
    """
    # Check if the size of the input lists are equal.
    if len(angles) != len(thickness) != len(Q):
        error_message = 'Error occeured, imput arrays not same length' + \
                        'length angle list = ' + len(angles) + \
                        'length thickness list = ' + len(thickness) + \
                        'length ply stiffness matrix list = ', + len(Q)
        raise ValueError(error_message)

    # Calculate the total thickness.
    h = np.sum(thickness) / 2

    # Create empty matricces for A B en D.
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))

    # Loop over all plies
    for i in range(len(angles)):
        # Calculate the z coordinates of the top and bottom of the ply.
        z_top = np.sum(thickness[:i]) - h
        z_bot = np.sum(thickness[:i+1]) - h

        # Rotate the local stiffenss matrix.
        Q_bar = stiffness_rotation(Q[i], angles[i])

        # Calculate the contribution to the A, B and D matrix of this layer.
        Ai = Q_bar * (z_bot - z_top)
        Bi = 1/2 * Q_bar * (z_bot**2 - z_top**2)
        Di = 1/3 * Q_bar * (z_bot**3 - z_top**3)

        # Summ this layer to the previous ones.
        A = A + Ai
        B = B + Bi
        D = D + Di

    # Compile the entirety of the ABD matrix.
    ABD = np.matrix([[A[0, 0], A[0, 1], A[0, 2], B[0, 0], B[0, 1], B[0, 2]],
                     [A[1, 0], A[1, 1], A[1, 2], B[1, 0], B[1, 1], B[1, 2]],
                     [A[2, 0], A[2, 1], A[2, 2], B[2, 0], B[2, 1], B[2, 2]],
                     [B[0, 0], B[0, 1], B[0, 2], D[0, 0], D[0, 1], D[0, 2]],
                     [B[1, 0], B[1, 1], B[1, 2], D[1, 0], D[1, 1], D[1, 2]],
                     [B[2, 0], B[2, 1], B[2, 2], D[2, 0], D[2, 1], D[2, 2]]])

    return ABD


# Invert the ABD matrix.
def abd_inverse(ABD):
    r"""
    Invert a ABD matrix.

    Parameters
    ----------
    ABD : matrix
        ABD matrix from system (N) = [ABD](e).

    Returns
    -------
    abd : matrix
        abd matrix is inverse of ABD, aka (e) = [abd](N)
    """
    # Check if this ABD matrix is invertible.
    if np.linalg.det(ABD) == 0:
        error_message = 'Determinant is equal to zero inverting not posible'
        raise ValueError(error_message)

    # Invert the matrix.
    abd = np.linalg.inv(ABD)

    return abd
