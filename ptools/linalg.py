"""Ptools linear algebra library."""

# Python core libraries.
import math

# Scientific libraries.
import numpy as np
import scipy.linalg as L


def angle(u: np.ndarray, v: np.ndarray) -> float:
    """Returns the angle between two vectors in radians."""
    return math.acos(np.dot(u, v) / (L.norm(u) * L.norm(v)))


def center_of_mass(x: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Return the center of mass (barycenter)."""
    N = x.shape[0]
    if N == 0:
        raise ZeroDivisionError("cannot compute center of mass of empty array")
    if weights.shape[0] != N:
        raise ValueError(
            f"input array and weights should be the same size ({N} != {weights.shape[0]})"
        )
    weights = weights.reshape(N, 1)
    return (x * weights).sum(axis=0) / weights.sum()


def centroid(x: np.ndarray) -> np.ndarray:
    """Return x centroid.

    Centroid is the average coordinate along axis 0.
    """
    return np.mean(x, axis=0)


def inertia_tensor(coords: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Returns the inertia tensor of a set of coordinates.

    Only works on N x 3 arrays.
    """
    if not coords.shape[1] == 3:
        raise ValueError(
            f"inertia tensor can only be calculated on a N x 3 array, (found {coords.shape})"
        )

    com = center_of_mass(coords, weights)
    X = coords - com

    x, y, z = X.T

    Ixx = np.sum(weights * (y**2 + z**2))
    Iyy = np.sum(weights * (x**2 + z**2))
    Izz = np.sum(weights * (x**2 + y**2))
    Ixy = -np.sum(weights * x * y)
    Iyz = -np.sum(weights * y * z)
    Ixz = -np.sum(weights * x * z)

    I = np.array(
        [
            [Ixx, Ixy, Ixz],
            [Ixy, Iyy, Iyz],
            [Ixz, Iyz, Izz],
        ]
    )

    return I


def principal_axes(tensor: np.ndarray, sort: bool = True) -> np.ndarray:
    """Return the principal axes given the tensor of inertia.

    Computes the eigenvalues and right eigenvectors of the tensor of inertia.

    Args:
        tensor (np.ndarray): 3 x 3 matrix
        sort (bool): sort axes by amplitude

    Returns:
        np.ndarray: 3 x 3 matrix
    """
    evalues, evectors = np.linalg.eig(tensor)
    if sort:
        order = np.argsort(evalues)[::-1]  # descending order
        evectors = evectors[:, order].transpose()
    return evectors
