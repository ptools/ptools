from collections.abc import Iterable
from typing import Any, Protocol

from ...particlecollection import ParticleCollection

PDB_FORMAT = (
    "{record:<6s}{atom_index:5s} "
    "{atom_name:4s}{altloc}{residue_name:<4s}{chain:s}{residue_index:>4s}{insertion}   "
    "{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          "
    "{element:>2s}"
)


class PDBConvertible(Protocol):
    """Protocol for objects that can be converted to PDB format.

    Those objects must have the following attributes:

    - coordinates: a 3-tuple of floats
    - name: a string, atom name
    - index: an integer, atom index
    - residue_name: a string, residue name
    - residue_index: an integer, residue index
    - chain: a string, chain name
    - element: a string, element name
    """

    @property
    def coordinates(self) -> tuple[float, float, float]: ...

    @property
    def name(self) -> str: ...

    @property
    def index(self) -> int: ...

    @property
    def residue_name(self) -> str: ...

    @property
    def residue_index(self) -> int: ...

    @property
    def chain(self) -> str: ...

    @property
    def occupancy(self) -> float | None: ...

    @property
    def bfactor(self) -> float | None: ...

    @property
    def element(self) -> str | None: ...


def to_pdb(atom_or_collection: PDBConvertible | Iterable[PDBConvertible]) -> str:
    """Converts an atom or a collection of atoms to PDB format."""

    if isinstance(atom_or_collection, ParticleCollection):
        atoms = atom_or_collection.dump()
        return "\n".join(format_atom(atom) for atom in atoms)  # type: ignore

    if isinstance(atom_or_collection, Iterable):
        return "\n".join(format_atom(atom) for atom in atom_or_collection)
    return format_atom(atom_or_collection)


def format_atom(atom: PDBConvertible) -> str:
    bfactor = getattr(atom, "bfactor", 0.0)
    occupancy = getattr(atom, "occupancy", 1.0)
    element = getattr(atom, "element", "")
    fmt_args: dict[str, Any] = {
        "record": "ATOM",
        "altloc": " ",
        "insertion": " ",
        "occupancy": occupancy,
        "bfactor": bfactor,
        "element": element,
        "chain": format_chain_token(atom.chain),
        "atom_name": format_atom_name_token(atom.name),
        "atom_index": format_atom_index_token(atom.index),
        "residue_name": format_residue_name_token(atom.residue_name),
        "residue_index": format_residue_index_token(atom.residue_index),
        "x": atom.coordinates[0],
        "y": atom.coordinates[1],
        "z": atom.coordinates[2],
    }
    return PDB_FORMAT.format(**fmt_args)


def format_chain_token(chain: str) -> str:
    return " " if not chain else chain[0]


def format_atom_name_token(name: str) -> str:
    return f"{name:>4s}" if len(name) > 2 else name.center(4)


def format_atom_index_token(index: int) -> str:
    return f"{index:5d}" if index < 100000 else f"{index:05x}"


def format_residue_name_token(name: str) -> str:
    return f"{name:<4s}" if len(name) > 3 else name.center(4)


def format_residue_index_token(index: int) -> str:
    return f"{index:4d}" if index < 10000 else f"{index:04x}"
