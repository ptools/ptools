
"""
Superposition methods.
"""

from pyptools import spatial
from pyptools.rigidbody import RigidBody

import numpy as np
from scipy.spatial.transform import Rotation


def kabsch_matrix(P, Q):
    """Calculates a rotation matrice using Kabsch algorithm.

    Calculates the rotations matrix between two sets of points.
    The points are expected to be centered on the origin.

    Args:
        P (numpy.ndarray (N x D))
        Q (numpy.ndarray (N x D))

    Returns
        numpy.ndarray: D x D rotation matrix
    """
    C = np.dot(np.transpose(P), Q)

    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    U = np.dot(V, W)

    return np.transpose(U)


def fit_matrix(mobile, target):
    """Return the fit matrix between two RigidBody."""
    t0 = target.get_center()
    t1 = mobile.get_center()

    # Center to origin.
    coords_target = target.coords - t0
    coords_mobile = mobile.coords - t0

    # Calculate rotation matrix.
    kabsch = kabsch_matrix(coords_mobile, coords_target)
    rotation = np.identity(4)
    rotation[:3, :3] = kabsch

    # Calculate translation component.
    result = spatial.translation_matrix(-t1)
    result = np.matmul(rotation, result)
    rotation = spatial.translation_matrix(t0)
    result = np.matmul(rotation, result)

    return result

def fit(mobile, target):
    """Fit two RigidBody."""
    matrix = fit_matrix(mobile, target)
    mobile.move(matrix)
    mobile.writepdb("coucou.pdb")


def main():
    target = RigidBody("/Users/benoist/Desktop/1BTA-new.pdb")
    mobile = RigidBody("/Users/benoist/Desktop/1BTA-2.pdb")

    matrix = fit_matrix(mobile, target)
    print(matrix)


if __name__ == '__main__':
    main()
