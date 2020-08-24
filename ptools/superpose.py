
"""
Superposition methods.
"""

import math

from ptools import spatial
from ptools.rigidbody import RigidBody

import numpy as np


class Screw:
    def __init__(self):
        self.unit = np.zeros(3)
        self.normtranslation = 0.0
        self.point = np.zeros(3)
        self.angle = 0.0

    def __str__(self):  # pragma: no cover
        return (f"Screw(unit={self.unit}, "
                f"normtranslation={self.normtranslation}), "
                f"point={self.point}), "
                f"angle={self.angle})")

    def copy(self):
        s = Screw()
        s.unit = self.unit
        s.normtranslation = self.normtranslation
        s.point = self.point
        s.angle = self.angle
        return s


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
    t0 = target.center()
    t1 = mobile.center()

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


def mat_trans_2_screw(matrix):
    assert matrix.shape == (4, 4)

    EPSILON = 1e-5

    trans = matrix[:3, 3]
    rotmatrix = matrix[:3, :3]

    x, y, z = rotmatrix[0], rotmatrix[1], rotmatrix[2]

    a = rotmatrix[0][0]
    b = rotmatrix[1][1]
    c = rotmatrix[2][2]

    eigenvect = np.zeros(3)
    screw = Screw()

    if abs(1 + a - b - c) > EPSILON:
        eigenvect[0] = x[0] + 1 - b - c
        eigenvect[1] = x[1] + y[0]
        eigenvect[2] = x[2] + z[0]
        screw.unit = eigenvect / np.linalg.norm(eigenvect)
        screw.normtranslation = np.dot(screw.unit, trans)

        s = trans - screw.normtranslation * screw.unit
        screw.point[0] = 0
        screw.point[1] = s[2] * z[1] + s[1] * (1 - z[2])
        screw.point[2] = s[1] * y[2] + s[2] * (1 - y[1])
        screw.point = screw.point / (1 + x[0] - y[1] - z[2])

    elif abs(1 - a + b - c) > EPSILON:
        eigenvect[0] = y[0] + x[1]
        eigenvect[1] = y[1] + 1 - x[0] - z[2]
        eigenvect[2] = y[2] + z[1]

        screw.unit = eigenvect / np.linalg.norm(eigenvect)
        screw.normtranslation = np.dot(screw.unit, trans)

        s = trans - screw.normtranslation * screw.unit
        screw.point[0] =  s.z * z[0] + s[0] * (1 - z[2])
        screw.point[1] =  0
        screw.point[2] =  s[0] * x[2] + s[2] * (1 - x[0])
        screw.point = screw.point / (1 - x[0] + y[1] - z[2])

    elif abs(1 - a - b + c) > EPSILON:
        eigenvect[0] = z[0] + x[2]
        eigenvect[1] = z[1] + y[2]
        eigenvect[2] = z[2] + 1 - x[0] - y[1]

        screw.unitVector = eigenvect / np.linalg.norm(eigenvect)
        screw.normtranslation = np.dot(screw.unitVector, trans)

        s = trans - screw.normtranslation * screw.unitVector
        screw.point[0] = s[1] * y[0] + s[0] * (1 - y[1])
        screw.point[1] = s[0] * x[1] + s[1] * (1 - x[0])
        screw.point[2] =  0
        screw.point = screw.point/(1 - x[0] - y[1] + z[2])

    else:     # angle=0
        screw.point = spatial.coord3d(0, 0, 0)
        norm_trans = np.linalg.norm(trans)
        if norm_trans != 0:
            screw.unit = trans / norm_trans
        else:
            screw.unit =  spatial.coord3d(0, 0, 1)
        screw.normtranslation = np.linalg.norm(trans)
        screw.angle = 0
        return screw

    v = spatial.coord3d((1, 0, 0))
    if abs(spatial.angle(screw.unit, v)) < 0.1:
        v = spatial.coord3d((0, 0, 1))

    u = v - np.dot(v, screw.unit) * screw.unit
    u /= np.linalg.norm(u)

    uprime = np.dot(rotmatrix, u)
    cost = np.dot(u, uprime)
    usec = np.cross(screw.unit, u)
    sint = np.dot(usec, uprime)

    if cost < -1:
        cost = -1
    if cost > 1:
        cost = 1

    screw.angle = math.acos(cost)
    if sint < 0:
        screw.angle = -screw.angle

    return screw


def main():  # pragma: no cover
    target = RigidBody("/Users/benoist/Desktop/1BTA-new.pdb")
    mobile = RigidBody("/Users/benoist/Desktop/1BTA-2.pdb")

    matrix = fit_matrix(mobile, target)
    mat_trans_2_screw(matrix)


if __name__ == '__main__':  # pragma: no cover
    main()
