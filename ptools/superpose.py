"""Superposition methods."""


from dataclasses import dataclass
import math

import numpy as np
from scipy.spatial.transform import Rotation

from ptools import spatial


@dataclass
class Screw:
    """A screw."""

    unit: np.ndarray = np.zeros(3)
    point: np.ndarray = np.zeros(3)
    normtranslation: float = 0.0
    angle: float = 0.0

    def copy(self):
        """Returns a copy of itself."""
        s = Screw(**self.__dict__)
        # Makes sure numpy arrays are actually copied.
        for key, value in s.__dict__.items():
            if isinstance(value, np.ndarray):
                setattr(s, key, value.copy())
        return s


def kabsch_matrix(mobile, target):
    """Calculates a rotation to optimally align two sets of coordinates.

    Uses Kabsch algorithm.
    The coordinates are expected to be centered on the origin.

    This is an alias for `scipy.spatial.transform.Rotation.align_vectors`.

    Args:
        mobile (numpy.ndarray (N, 3))
        target (numpy.ndarray (N, 3))

    Returns
        numpy.ndarray: (3, 3) matrix
    """
    rotation, _ = Rotation.align_vectors(target, mobile)
    return rotation.as_matrix()


def fit_matrix(mobile, target):
    """Return the fit matrix between two RigidBody."""
    t0 = target.center()
    t1 = mobile.center()

    # Center to origin.
    coords_target = target.coords - t0
    coords_mobile = mobile.coords - t1

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


def rmsd(mobile, target, do_fit=False):
    """Returns the Root Mean Square Deviation between two groups of atoms."""
    assert len(mobile) == len(target)
    if do_fit:
        mobile = mobile.copy()
        fit(mobile, target)
    e = np.power((mobile.coords - target.coords), 2).sum(axis=1)
    return np.sqrt(np.mean(e))


def mat_trans_2_screw(matrix):
    """Converts a transformation matrix to a Screw

    Args:
        matrix (np.ndarray): 4 x 4 matrix

    Returns:
        Screw
    """
    assert matrix.shape == (4, 4)

    EPSILON = 1e-5

    trans = matrix[:3, 3]
    rotmatrix = matrix[:3, :3]

    x, y, z = rotmatrix[:, 0], rotmatrix[:, 1], rotmatrix[:, 2]

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
        screw.point[0] =  s[2] * z[0] + s[0] * (1 - z[2])
        screw.point[1] =  0
        screw.point[2] =  s[0] * x[2] + s[2] * (1 - x[0])
        screw.point = screw.point / (1 - x[0] + y[1] - z[2])

    elif abs(1 - a - b + c) > EPSILON:
        eigenvect[0] = z[0] + x[2]
        eigenvect[1] = z[1] + y[2]
        eigenvect[2] = z[2] + 1 - x[0] - y[1]

        screw.unit = eigenvect / np.linalg.norm(eigenvect)
        screw.normtranslation = np.dot(screw.unit, trans)

        s = trans - screw.normtranslation * screw.unit
        screw.point[0] = s[1] * y[0] + s[0] * (1 - y[1])
        screw.point[1] = s[0] * x[1] + s[1] * (1 - x[0])
        screw.point[2] =  0
        screw.point = screw.point/(1 - x[0] - y[1] + z[2])

    else:     # angle=0
        screw.point = np.zeros(3)
        norm_trans = np.linalg.norm(trans)
        if norm_trans != 0:
            screw.unit = trans / norm_trans
        else:
            screw.unit =  np.array((0, 0, 1))
        screw.normtranslation = np.linalg.norm(trans)
        screw.angle = 0
        return screw

    v = np.array((1, 0, 0))
    if abs(spatial.angle(screw.unit, v)) < 0.1:
        v = np.array((0, 0, 1))

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
