from typing import Any, Iterable, Protocol
from ..._typing import FilePath

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
    def element(self) -> str:
        ...


def write_pdb(
    atom_or_collection: PDBConvertible | Iterable[PDBConvertible], path: FilePath
):
    """Write a PDB file from an atom or a collection of atoms."""
    with open(path, "w") as f:
        f.write(to_pdb(atom_or_collection))


def to_pdb(atom_or_collection: PDBConvertible | Iterable[PDBConvertible]) -> str:
    if isinstance(atom_or_collection, Iterable):
        return "\n".join(format_atom(atom) for atom in atom_or_collection)
    return format_atom(atom_or_collection)


def format_atom(atom: PDBConvertible) -> str:
    fmt_args: dict[str, Any] = {
        "record": "ATOM",
        "altloc": " ",
        "insertion": " ",
        "occupancy": 1.0,
        "bfactor": 0.0,
        "element": atom.element,
        "chain": _format_chain(atom.chain),
        "atom_name": _format_atom_name(atom.name),
        "atom_index": _format_atom_index(atom.index),
        "residue_name": _format_residue_name(atom.residue_name),
        "residue_index": _format_residue_index(atom.residue_index),
        "x": atom.coordinates[0],
        "y": atom.coordinates[1],
        "z": atom.coordinates[2],
    }
    return PDB_FORMAT.format(**fmt_args)


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
