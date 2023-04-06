from typing import Iterable

from ..formatters.pdb import to_pdb, PDBConvertible
from ..._typing import FilePath


def write_pdb(
    atom_or_collection: PDBConvertible | Iterable[PDBConvertible], path: FilePath
):
    """Write a PDB file from an atom or a collection of atoms."""
    with open(path, "w") as f:
        f.write(to_pdb(atom_or_collection))
