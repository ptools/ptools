"""ptools.measure - Measure geometric properties."""

# Python core modules.
from typing import Optional

# Scientific libraries.
import numpy as np

# Ptools modules.
from . import linalg as L
from ._typing import ArrayLike, HasCoordinates, Topology


def distance(lhs: HasCoordinates, rhs: HasCoordinates) -> float:
    """Returns the euclidean distance between two objects."""
    return L.distance(lhs.coordinates, rhs.coordinates)


def distance_to_axis(obj: HasCoordinates, axis: ArrayLike) -> float:
    """Returns the distance between an object and an arbitrary axis."""
    return L.distance_to_axis(obj.coordinates, axis)


def centroid(obj: HasCoordinates) -> np.ndarray:
    """Returns an object geometric center, i.e. isobarycenter."""
    return L.centroid(obj.coordinates)


def center_of_mass(obj: Topology) -> np.ndarray:
    """Returns an object center of mass, i.e. barycenter."""
    return L.center_of_mass(obj.coordinates, obj.masses)


def inertia_tensor(obj: HasCoordinates | Topology, weights: Optional[ArrayLike] = None):
    """Returns the inertia tensors of a set of atoms."""
    if weights is None:
        weights = getattr(obj, "masses", np.ones(len(obj.coordinates)))
    return L.inertia_tensor(obj.coordinates, weights)


def principal_axes(obj: Topology | HasCoordinates, sort: bool = True) -> np.ndarray:
    """Returns an object principal axes.

    Args:
        sort (bool): sort axes by importance
    """
    return L.principal_axes(inertia_tensor(obj), sort)


def radius_of_gyration(obj: HasCoordinates) -> float:
    """Returns the isometric radius of gyration (atom mass is not taken
    into account)."""
    centered = obj.coordinates - centroid(obj)
    rgyr2 = np.sum(centered**2) / len(obj.coordinates)
    return rgyr2**0.5
