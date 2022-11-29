"""Superposition methods."""


from dataclasses import dataclass, field
import math

import numpy as np
from scipy.spatial.transform import Rotation

from .atomcollection import AtomCollection
from . import linalg


def zeros3f():
    return np.zeros(3, dtype="float64")


@dataclass
class Screw:
    """A screw."""

    unit: np.ndarray = field(default_factory=zeros3f)
    point: np.ndarray = field(default_factory=zeros3f)
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


def kabsch_matrix(mobile: AtomCollection, target: AtomCollection) -> np.ndarray:
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


def fit_matrix(mobile: AtomCollection, target: AtomCollection) -> np.ndarray:
    """Returns the fit matrix between two ``AtomCollection`` instances."""
    t0 = target.centroid()
    t1 = mobile.centroid()

    # Center to origin.
    coords_target = target.coords - t0
    coords_mobile = mobile.coords - t1

    # Calculate rotation matrix.
    kabsch = kabsch_matrix(coords_mobile, coords_target)
    rotation = np.identity(4)
    rotation[:3, :3] = kabsch

    # Calculate translation component.
    result = linalg.translation_matrix(-t1)
    result = np.matmul(rotation, result)
    rotation = linalg.translation_matrix(t0)
    result = np.matmul(rotation, result)

    return result


def fit(mobile: AtomCollection, target: AtomCollection):
    """Fits two ``mobile`` onto ``target``."""
    matrix = fit_matrix(mobile, target)
    mobile.move(matrix)


def rmsd(mobile: AtomCollection, target: AtomCollection, do_fit: bool = False) -> float:
    """Returns the Root Mean Square Deviation between two groups of atoms."""
    assert len(mobile) == len(target)
    if do_fit:
        mobile = mobile.copy()
        fit(mobile, target)
    e = np.power((mobile.coords - target.coords), 2).sum(axis=1)
    return np.sqrt(np.mean(e))


# pylint: disable=R0914,R0915,R1730
# This function is actually too long, has too many variables, and cannot be
# understood.
def mat_trans_2_screw(matrix: np.ndarray) -> Screw:
    """Converts a transformation matrix to a Screw

    Args:
        matrix (np.ndarray): 4 x 4 matrix

    Returns:
        Screw

    """

    def initialize_screw(eigenvect):
        screw = Screw()
        screw.unit = eigenvect / np.linalg.norm(eigenvect)
        screw.normtranslation = np.dot(screw.unit, translation)
        return screw

    assert matrix.shape == (4, 4)

    EPSILON = 1e-5

    translation = matrix[:3, 3]
    rotation = matrix[:3, :3]

    x, y, z = rotation[:, 0], rotation[:, 1], rotation[:, 2]
    a, b, c = np.diag(rotation)

    if abs(1 + a - b - c) > EPSILON:
        eigenvect = np.array([x[0] + 1 - b - c, x[1] + y[0], x[2] + z[0]])
        screw = initialize_screw(eigenvect)

        s = translation - screw.normtranslation * screw.unit
        screw.point[0] = 0
        screw.point[1] = s[2] * z[1] + s[1] * (1 - z[2])
        screw.point[2] = s[1] * y[2] + s[2] * (1 - y[1])
        screw.point = screw.point / (1 + x[0] - y[1] - z[2])

    elif abs(1 - a + b - c) > EPSILON:
        eigenvect = np.array([y[0] + x[1], y[1] + 1 - x[0] - z[2], y[2] + z[1]])
        screw = initialize_screw(eigenvect)

        s = translation - screw.normtranslation * screw.unit
        screw.point[0] = s[2] * z[0] + s[0] * (1 - z[2])
        screw.point[1] = 0
        screw.point[2] = s[0] * x[2] + s[2] * (1 - x[0])
        screw.point = screw.point / (1 - x[0] + y[1] - z[2])

    elif abs(1 - a - b + c) > EPSILON:
        eigenvect = np.array([z[0] + x[2], z[1] + y[2], z[2] + 1 - x[0] - y[1]])
        screw = initialize_screw(eigenvect)

        s = translation - screw.normtranslation * screw.unit
        screw.point[0] = s[1] * y[0] + s[0] * (1 - y[1])
        screw.point[1] = s[0] * x[1] + s[1] * (1 - x[0])
        screw.point[2] = 0
        screw.point = screw.point / (1 - x[0] - y[1] + z[2])

    else:  # angle=0
        screw = Screw()
        norm_trans = np.linalg.norm(translation)
        if norm_trans != 0:
            screw.unit = translation / norm_trans
        else:
            screw.unit = np.array((0, 0, 1))
        screw.normtranslation = np.linalg.norm(translation)
        screw.angle = 0
        return screw

    v = np.array((1, 0, 0))
    if abs(linalg.angle(screw.unit, v)) < 0.1:
        v = np.array((0, 0, 1))

    u = v - np.dot(v, screw.unit) * screw.unit
    u /= np.linalg.norm(u)

    uprime = np.dot(rotation, u)
    cost = np.dot(u, uprime)
    usec = np.cross(screw.unit, u)
    sint = np.dot(usec, uprime)

    if cost < -1:
        cost = -1
    elif cost > 1:
        cost = 1

    screw.angle = math.acos(cost)
    if sint < 0:
        screw.angle = -screw.angle

    return screw
