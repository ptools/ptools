"""ptools.measure - Measure geometric properties."""

# Python core modules.

# Scientific libraries.
import numpy as np

# Ptools modules.
from . import linalg as L
from ._typing import ArrayLike, HasCoordinatesType, TopologyType
from .pairlist import PairList
from .particlecollection import ParticleCollection


def bounding_box(obj: HasCoordinatesType) -> np.ndarray:
    """Returns the bounding box of object with coordinates."""
    return np.array((np.min(obj.coordinates, axis=0), np.max(obj.coordinates, axis=0)))


def distance(lhs: HasCoordinatesType, rhs: HasCoordinatesType) -> float:
    """Returns the euclidean distance between two objects."""
    return L.distance(lhs.coordinates, rhs.coordinates)


def distance_to_axis(
    obj: HasCoordinatesType, axis: ArrayLike, center: bool | ArrayLike = False
) -> float:
    """Returns the distance between an object and an arbitrary axis."""
    return L.distance_to_axis(obj.coordinates, axis, center)


def minmax_distance_to_axis(
    obj: HasCoordinatesType, axis: ArrayLike, center: bool | ArrayLike = False
) -> tuple[float, float]:
    """Returns the minimum and maximum distance between an object and an
    arbitrary axis."""
    return L.minmax_distance_to_axis(obj.coordinates, axis, center)


def centroid(obj: HasCoordinatesType) -> np.ndarray:
    """Returns an object geometric center, i.e. isobarycenter."""
    return L.centroid(obj.coordinates)


def center(obj: HasCoordinatesType) -> np.ndarray:
    """Returns an object geometric center, i.e. isobarycenter.

    This is an alias for `centroid`.
    """
    return centroid(obj)


def center_of_mass(obj: TopologyType) -> np.ndarray:
    """Returns an object center of mass, i.e. barycenter."""
    return L.center_of_mass(obj.coordinates, obj.masses)


def inertia_tensor(
    obj: HasCoordinatesType | TopologyType, weights: ArrayLike | None = None
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


def raw_contacts(
    obj1: HasCoordinatesType, obj2: HasCoordinatesType, cutoff: float
) -> tuple[list[int], list[int]]:
    """Returns the indexes of the atoms of `obj1` that are in contact with
    the atoms of `obj2`.

    Returns two lists of atom indexes, one for each object.
    """
    return PairList(obj1, obj2, cutoff).raw_contacts()  # type: ignore[return-value]


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
    """Returns the indexes of the residues of `obj1` that are in contact with the
    residues of `obj2`."""
    lhs_atom_ids, rhs_atom_ids = PairList(lhs, rhs, cutoff).raw_contacts()
    lhs_residue_ids = lhs.residue_indices[lhs_atom_ids]
    rhs_residue_ids = rhs.residue_indices[rhs_atom_ids]
    by_residue = set(zip(lhs_residue_ids, rhs_residue_ids, strict=False))
    return by_residue


def fnat(
    receptor1: ParticleCollection,
    lig1: ParticleCollection,
    receptor2: ParticleCollection,
    lig2: ParticleCollection,
    cutoff: float,
) -> float:
    """Computes the fraction of native contacts between 2 pairs of `RigidBody` instances."""

    res_pair1 = contacts_by_residue(receptor1, lig1, cutoff)
    if len(res_pair1) == 0:
        return 0

    res_pair2 = contacts_by_residue(receptor2, lig2, cutoff)
    intersect = res_pair1 & res_pair2
    return len(intersect) / len(res_pair1)
