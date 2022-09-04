"""Ptools linear algebra library."""

# Python core libraries.
import math

# Scientific libraries.
import numpy as np
import scipy.linalg as L

def angle(u: np.ndarray, v: np.ndarray) -> float:
    """Returns the angle between two vectors in radians."""
    return math.acos(np.dot(u, v) / (L.norm(u) * L.norm(v)))


def center_of_mass(coords: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Return the center of mass (barycenter)."""
    weights = weights.reshape(-1, 1)
    return (coords * weights).sum(axis=0) / weights.sum()


def centroid(x: np.ndarray) -> np.ndarray:
    """Return x centroid.

    Centroid is the average coordinate along axis 0.
    """
    return np.mean(x, axis=0)


def _tensor_of_inertia_accurate(coords: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Returns the inertia tensors of a set of coordinates."""
    com = center_of_mass(coords, weights)
    X = coords - com

    Ixx = np.sum(weights * (X[:, 1] ** 2 + X[:, 2] ** 2))
    Iyy = np.sum(weights * (X[:, 0] ** 2 + X[:, 2] ** 2))
    Izz = np.sum(weights * (X[:, 0] ** 2 + X[:, 1] ** 2))
    Ixy = -np.sum(weights * X[:, 0] * X[:, 1])
    Iyz = -np.sum(weights * X[:, 1] * X[:, 2])
    Ixz = -np.sum(weights * X[:, 0] * X[:, 2])

    I = np.array(
        [
            [Ixx, Ixy, Ixz],
            [Ixy, Iyy, Iyz],
            [Ixz, Iyz, Izz],
        ]
    )

    return I


def _tensor_of_inertia_fast(coords: np.ndarray) -> np.ndarray:
    """Returns the inertia tensors of a set of atoms.

    Fast: does not take into account the weights
    """
    com = centroid(coords)
    X = coords - com
    tensors = np.dot(X.transpose(), X)
    return tensors


def tensor_of_inertia(
    coords: np.ndarray, weights: np.ndarray = None, method: str = "accurate"
) -> np.ndarray:
    """Returns the inertia tensors of a set of coordinates.

    Allows to choose between "accurate" (weighted) and "fast" methods.
    """
    if method not in ("accurate", "fast"):
        raise ValueError(
            "parameter 'method' should be 'accurate' or 'fast' " f"(found {method=!r})"
        )
    if method == "accurate":
        if weights is None:
            raise ValueError("need weights to compute accurate tensor of " "inertia")
        return _tensor_of_inertia_accurate(coords, weights)
    return _tensor_of_inertia_fast(coords)  # method is "fast"


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

