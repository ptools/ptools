from typing import Any

from ...atomattrs import AtomAttrs
from ...atomcollection import AtomCollection


PDB_FORMAT = (
    "{record:<6s}{atom_index:5s} "
    "{atom_name:4s}{altloc}{residue_name:<4s}{chain:s}{residue_index:>4s}{insertion}   "
    "{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          "
    "{element:>2s}"
)


def to_pdb(atom_or_collection: AtomAttrs | AtomCollection) -> str:
    if isinstance(atom_or_collection, AtomCollection):
        return "\n".join(format_atom(atom) for atom in atom_or_collection)
    return format_atom(atom_or_collection)


def format_atom(atom: "AtomAttrs", **kwargs) -> str:
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
        "x": atom.coords[0],
        "y": atom.coords[1],
        "z": atom.coords[2],
    }
    fmt_args.update(kwargs)
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
