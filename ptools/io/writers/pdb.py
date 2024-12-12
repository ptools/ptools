from collections.abc import Iterable

from ..._typing import FilePath
from ..formatters.pdb import PDBConvertible, to_pdb


def write_pdb(
    atom_or_collection: PDBConvertible | Iterable[PDBConvertible], path: FilePath
):
    """Write a PDB file from an atom or a collection of atoms."""
    with open(path, "w") as f:
        f.write(to_pdb(atom_or_collection))
