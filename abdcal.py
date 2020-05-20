"""
ABD matrix calculator for a given stacking sequence.

The methods in this file will call create a ABD matrix of a composit for given
ply properties and stacking sequence.

A.J.J. Lagerweij
COHMAS Mechanical Engineering KAUST
2020
"""

# Import external package.
import numpy as np

# Import local package.
import deformation


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
    Rotate the stiffness matrix against given angle.

    This rotates the stiffness matrix from local to the global axis sytem.
    Use a negative angle to rotate from global to local system.

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

    This rotates the complianc matrix from local to the global axis sytem.
    Use a negative angle to rotate from global to local system.

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


# Calculate the stiffness matrix of a thin (no bending) laminate.
def abdthin(Q, angles, thickness, truncate=False):
    r"""
    ABD matrix calculator for a thin laminate.

    In the thin laminate theroy it is assumed that the out of plane stiffness
    is negligible and that the layup is symmetric. Hence only the membrane
    (A part) of the ABD matrix remains. Top plies should be listed first in the
    lists of Q, angles and thickness.

    Parameters
    ----------
    Q : list
        The stiffness matrix of each ply in its l-t axis system.
    angles : list
        The rotation of each ply in degrees.
    thickness : list
        The thickness of each ply.
    truncate : bool
        Truncates very small numbers when true.

    Returns
    -------
    C : matrix
        The stiffness matrix of the thin laminate.
    """
    # Create an empty stiffness matrix.
    C = np.zeros((3, 3))

    # Loop over all plies.
    for i in range(len(angles)):
        Q_bar = stiffness_rotation(Q[i], angles[i])
        C += Q_bar * thickness[i]

    # Go from the stresses * length^2 to running loads.
    C = np.matrix(C) / np.sum(thickness)

    # Truncate very small values.
    if truncate is True:
        C = np.matrix(np.where(np.abs(C) < np.max(C)*1e-6, 0, C))
    return C


# Calculate the ABD matrix of a given laminate.
def abd(Q, angles, thickness, truncate=False):
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
    truncate : bool
        Truncates very small numbers when true.

    Returns
    -------
    ABD : matrix
        The stiffness matrix of the thin laminate.
    """
    # Calculate the total thickness.
    h = sum(thickness) / 2

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

    # Truncate very small values.
    if truncate is True:
        ABD = np.matrix(np.where(np.abs(ABD) < np.max(ABD)*1e-6, 0, ABD))
    return ABD


# Calculate the global coefficent of thermal expansion of a thin laminate.
def ctethin(Q, angles, thickness, alpha):
    r"""
    Thin laminate Coefficient of Thermal Expansion calculator.

    This funcion calculates the CTE in the x-y axis sytem. It summs up the
    rotated CTE of each layer and weights them by layer stiffness and
    thickness. Here it is assumed that there is no bending behaviour this is
    only true for thin and symmetric layups.
    Top plies should be listed first in the lists of Q, angles, thickness and
    alpha.

    Parameters
    ----------
    Q : list
        The stiffness matrix of each ply in its l-t axis system.
    angles : list
        The rotation of each ply in degrees.
    thickness : list
        The thickness of each ply.
    alpha : list
        The coeficcient of termal expansion of each ply in l-t axis system.

    Returns
    -------
    cte : vector
        The coefficient of termal expansion of the laminate, in x-y axis sytem.
    """
    # Calculate the laminate compliance matrix.
    C = abdthin(Q, angles, thickness)
    S = matrix_inverse(C)

    # Create an empty CTE vector.
    cte = np.zeros((3, 1))

    # Loop ver all plies.
    for i in range(len(angles)):
        alpha_bar = deformation.strain_rotation(alpha[i], angles[i])
        Q_bar = stiffness_rotation(Q[i], angles[i])
        cte += Q_bar * alpha_bar * thickness[i]

    # Normalize the cte to be in strain per temperature unit.
    cte = cte / sum(thickness)
    cte = S * cte
    return cte


# Calculate the global coefficent of thermal expansion.
def cte(Q, angles, thickness, alpha):
    r"""
    Coefficient of Thermal Expansion calculator.

    This funcion calculates the CTE in the x-y axis sytem. It summs up the
    rotated CTE of each layer and multiplies them by layer stiffness and
    thickness. The resulting CTE relates thermal change (:math:`\Delta T`) to
    the deformation vector (strains and curvatures).
    Top plies should be listed first in the lists of Q, angles,
    thickness and alpha.

    Parameters
    ----------
    Q : list
        The stiffness matrix of each ply in its l-t axis system.
    angles : list
        The rotation of each ply in degrees.
    thickness : list
        The thickness of each ply.
    alpha : list
        The coeficcient of termal expansion of each ply in l-t axis system.

    Returns
    -------
    cte : vector
        The coefficient of termal expansion of the laminate, in x-y axis sytem.
    """
    # Calculate the total thickness.
    h = sum(thickness) / 2

    # Create an empty cte * stiffness vector vector.
    alphaEt = np.zeros((6, 1))

    # Loop ver all plies.
    for i in range(len(angles)):
        # Calculate the z coordinates of the top and bottom of the ply.
        z_top = np.sum(thickness[:i]) - h
        z_bot = np.sum(thickness[:i+1]) - h

        # rotate the ply cte and Q f
        alpha_bar = deformation.strain_rotation(alpha[i], angles[i])
        Q_bar = stiffness_rotation(Q[i], angles[i])

        # Calculate the contribution if this ply to the total thermal load.
        alphaEt[:3] += Q_bar * alpha_bar * (z_bot - z_top)
        alphaEt[3:] += 1/2 * Q_bar * alpha_bar * (z_bot**2 - z_top**2)

    # Convert alphaEt which links Delta T and running loads to strain and
    # curvature.
    ABD = abd(Q, angles, thickness)
    ABD_inv = matrix_inverse(ABD)
    cte = deformation.load_applied(ABD_inv, alphaEt)  # stren per unit change T
    return cte


# Invert the matrix.
def matrix_inverse(matrix):
    r"""
    Inverts a matrix.

    Parameters
    ----------
    matrix : matrix
        The matrix to be inverted.

    Returns
    -------
    inverse : matrix
        The inverse of the matrix.
    """
    # Check if this ABD matrix is invertible.
    if np.linalg.det(matrix) == 0:
        error_message = 'Determinant is equal to zero inverting not posible'
        raise ValueError(error_message)

    # Invert the matrix.
    inverse = np.linalg.inv(matrix)
    return inverse
