"""Provides functions to write reduces topologies.

Reduced PDB files and PDB files are identical for the first part of the line.
After the bead coodinates, the reduced PDB file contains

- the bead type,
- the bead charge,
- 2 zeros (because why not).

"""

from typing import Any, Iterable, Protocol

from ..._typing import FilePath
from ...particlecollection import ParticleCollection
from ...namedarray import NamedArrayContainer
from .pdb import (
    format_chain_token,
    format_atom_name_token,
    format_atom_index_token,
    format_residue_name_token,
    format_residue_index_token,
)


class ReducedPDBConvertible(Protocol):
    """Protocol for objects that can be converted to Reduced PDB format.

    A ``ReducedPDBConvertible`` object has the same attributes as a ``pdb.PDBConvertible``
    object with the addition of the following attributes:

    - typeid (int), bead type id,
    - element (float), bead charge.
    """

    @property
    def coordinates(self) -> tuple[float, float, float]:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def index(self) -> int:
        ...

    @property
    def residue_name(self) -> str:
        ...

    @property
    def residue_index(self) -> int:
        ...

    @property
    def chain(self) -> str:
        ...

    @property
    def typeid(self) -> int:
        ...

    @property
    def charge(self) -> float:
        ...


class ParticleBuffer:
    """Returns a Particle properties."""

    def __init__(self, properties: NamedArrayContainer, i: int):
        self.coordinates = properties["coordinates"][i]
        self.name = properties["names"][i]
        self.index = properties["indices"][i]
        self.residue_name = properties["residue_names"][i]
        self.residue_index = properties["residue_indices"][i]
        self.chain = properties["chains"][i]
        self.element = properties["elements"][i]
        self.typeid = properties["typeids"][i]
        self.charge = properties["charges"][i]


REDUCED_PDB_FORMAT = (
    "{record:<6s}{atom_index:5s} "
    "{atom_name:4s}{altloc}{residue_name:<4s}{chain:s}{residue_index:>4s}{insertion}   "
    "{x:8.3f}{y:8.3f}{z:8.3f}{typeid:5d}{charge:8.3f} 0 0"
)


def write_reduced_pdb(
    atom_or_collection: ReducedPDBConvertible | Iterable[ReducedPDBConvertible],
    path: FilePath,
):
    """Write a PDB file from an atom or a collection of atoms."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(to_reduced_pdb(atom_or_collection))


def to_reduced_pdb(
    atom_or_collection: ReducedPDBConvertible | Iterable[ReducedPDBConvertible],
) -> str:
    """Converts an atom or a collection of atoms to Reduced PDB format."""

    if isinstance(atom_or_collection, ParticleCollection):
        properties = atom_or_collection.atom_properties
        return "\n".join(format_atom(ParticleBuffer(properties, i)) for i in range(len(atom_or_collection)))  # type: ignore

    if isinstance(atom_or_collection, Iterable):
        return "\n".join(format_atom(atom) for atom in atom_or_collection)
    return format_atom(atom_or_collection)


def format_atom(atom: ReducedPDBConvertible) -> str:
    fmt_args: dict[str, Any] = {
        "record": "ATOM",
        "altloc": " ",
        "insertion": " ",
        "occupancy": 1.0,
        "bfactor": 0.0,
        "chain": format_chain_token(atom.chain),
        "atom_name": format_atom_name_token(atom.name),
        "atom_index": format_atom_index_token(atom.index),
        "residue_name": format_residue_name_token(atom.residue_name),
        "residue_index": format_residue_index_token(atom.residue_index),
        "x": atom.coordinates[0],
        "y": atom.coordinates[1],
        "z": atom.coordinates[2],
        "typeid": atom.typeid,
        "charge": atom.charge,
    }
    return REDUCED_PDB_FORMAT.format(**fmt_args)
