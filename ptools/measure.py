"""ptools.measure - Measure geometric properties."""

# Python core modules.
from typing import Optional

# Scientific libraries.
import numpy as np

# Ptools modules.
from . import linalg as L
from .pairlist import PairList
from ._typing import ArrayLike, HasCoordinatesType, TopologyType
from .particlecollection import ParticleCollection


def distance(lhs: HasCoordinatesType, rhs: HasCoordinatesType) -> float:
    """Returns the euclidean distance between two objects."""
    return L.distance(lhs.coordinates, rhs.coordinates)


def distance_to_axis(obj: HasCoordinatesType, axis: ArrayLike) -> float:
    """Returns the distance between an object and an arbitrary axis."""
    return L.distance_to_axis(obj.coordinates, axis)


def centroid(obj: HasCoordinatesType) -> np.ndarray:
    """Returns an object geometric center, i.e. isobarycenter."""
    return L.centroid(obj.coordinates)


def center_of_mass(obj: TopologyType) -> np.ndarray:
    """Returns an object center of mass, i.e. barycenter."""
    return L.center_of_mass(obj.coordinates, obj.masses)


def inertia_tensor(
    obj: HasCoordinatesType | TopologyType, weights: Optional[ArrayLike] = None
):
    """Returns the inertia tensors of a set of atoms."""
    if weights is None:
        weights = np.asarray(getattr(obj, "masses", np.ones(len(obj.coordinates))))
    return L.inertia_tensor(obj.coordinates, weights)


def principal_axes(
    obj: TopologyType | HasCoordinatesType, sort: bool = True
) -> np.ndarray:
    """Returns an object principal axes.

    Args:
        sort (bool): sort axes by importance
    """
    return L.principal_axes(inertia_tensor(obj), sort)


def radius_of_gyration(obj: HasCoordinatesType) -> float:
    """Returns the isometric radius of gyration (atom mass is not taken
    into account)."""
    centered = obj.coordinates - centroid(obj)
    rgyr2 = np.sum(centered**2) / len(obj.coordinates)
    return rgyr2**0.5


def contacts(
    obj1: HasCoordinatesType, obj2: HasCoordinatesType, cutoff: float
) -> list[tuple[int, int]]:
    """Returns the indexes of the atoms of `obj1` that are in contact with
    the atoms of `obj2`.
    """
    return PairList(obj1, obj2, cutoff).contacts()


def contacts_by_atom(
    obj1: HasCoordinatesType, obj2: HasCoordinatesType, cutoff: float
) -> list[tuple[int, int]]:
    """Returns the indexes of the atoms of `obj1` that are in contact with
    the atoms of `obj2`.

    This is an alias for `contacts`.
    """
    return contacts(obj1, obj2, cutoff)


def contacts_by_residue(
    lhs: ParticleCollection, rhs: ParticleCollection, cutoff: float
) -> set[tuple[int, int]]:
    """Returns the indexes of the residues of `obj1` that are in contact with
    the residues of `obj2`.
    """
    by_atom = contacts(lhs, rhs, cutoff)
    res_indexes = set(
        (lhs[index1].residue_index, rhs[index2].residue_index)
        for index1, index2 in by_atom
    )
    return res_indexes
