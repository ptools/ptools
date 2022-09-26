from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from ...atom import BaseAtom


class PDBFormatter:
    fmt = (
        "{record:<6s}{atom_index:5s} "
        "{atom_name:4s}{altloc}{residue_name:<4s}{chain:s}{residue_index:>4s}{insertion}   "
        "{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          "
        "{element:>2s}"
    )

    @classmethod
    def format_atom(cls, atom: "BaseAtom", **kwargs) -> str:
        fmt_args: dict[str, Any] = {
            "record": "ATOM",
            "altloc": " ",
            "insertion": " ",
            "occupancy": 1.0,
            "bfactor": 0.0,
            "element": atom.element,
            "chain": cls.format_chain(atom.chain),
            "atom_name": cls.format_atom_name(atom.name),
            "atom_index": cls.format_atom_index(atom.index),
            "residue_name": cls.format_residue_name(atom.resname),
            "residue_index": cls.format_residue_index(atom.resid),
            "x": atom.coords[0],
            "y": atom.coords[1],
            "z": atom.coords[2],
        }
        fmt_args.update(kwargs)
        return cls.fmt.format(**fmt_args)

    @classmethod
    def format_chain(self, chain: str) -> str:
        assert len(chain) < 2
        return " " if not chain else chain

    @classmethod
    def format_atom_name(self, name: str) -> str:
        return f"{name:>4s}" if len(name) > 2 else name.center(4)

    @classmethod
    def format_atom_index(self, index: int) -> str:
        return f"{index:5d}" if index < 100000 else f"{index:05x}"

    @classmethod
    def format_residue_name(self, name: str) -> str:
        return f"{name:<4s}" if len(name) > 3 else name.center(4)

    @classmethod
    def format_residue_index(self, index: int) -> str:
        return f"{index:4d}" if index < 10000 else f"{index:04x}"
