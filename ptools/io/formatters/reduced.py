"""Provides functions to write reduces topologies.

Reduced PDB files and PDB files are identical for the first part of the line.
After the bead coodinates, the reduced PDB file contains

- the bead type,
- the bead charge,
- 2 zeros (because why not).

"""

from typing import Any, Iterable, Protocol
from ..._typing import FilePath

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



REDUCED_PDB_FORMAT = (
    "{record:<6s}{atom_index:5s} "
    "{atom_name:4s}{altloc}{residue_name:<4s}{chain:s}{residue_index:>4s}{insertion}   "
    "{x:8.3f}{y:8.3f}{z:8.3f}{typeid:5d}{charge:8.3f} 0 0"
)


def write_reduced_pdb(atom_or_collection: ReducedPDBConvertible | Iterable[ReducedPDBConvertible], path: FilePath):
    """Write a PDB file from an atom or a collection of atoms."""
    with open(path, "w") as f:
        f.write(to_reduced_pdb(atom_or_collection))


def to_reduced_pdb(atom_or_collection: ReducedPDBConvertible | Iterable[ReducedPDBConvertible]) -> str:
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
        "chain": _format_chain(atom.chain),
        "atom_name": _format_atom_name(atom.name),
        "atom_index": _format_atom_index(atom.index),
        "residue_name": _format_residue_name(atom.residue_name),
        "residue_index": _format_residue_index(atom.residue_index),
        "x": atom.coordinates[0],
        "y": atom.coordinates[1],
        "z": atom.coordinates[2],
        "typeid": atom.typeid,
        "charge": atom.charge,
    }
    return REDUCED_PDB_FORMAT.format(**fmt_args)


def _format_chain(chain: str) -> str:
    assert len(chain) < 2
    return " " if not chain else chain


def _format_atom_name(name: str) -> str:
    return f"{name:>4s}" if len(name) > 2 else name.center(4)


def _format_atom_index(index: int) -> str:
    return f"{index:5d}" if index < 100000 else f"{index:05x}"


def _format_residue_name(name: str) -> str:
    return f"{name:<4s}" if len(name) > 3 else name.center(4)


def _format_residue_index(index: int) -> str:
    return f"{index:4d}" if index < 10000 else f"{index:04x}"
