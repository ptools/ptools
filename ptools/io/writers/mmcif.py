from collections.abc import Iterable

from ..._typing import FilePath
from ..formatters.mmcif import mmCIFConvertible, to_mmCIF


def write_mmCIF(atom_or_collection: mmCIFConvertible | Iterable[mmCIFConvertible], path: FilePath):
    """Write a PDB file from an atom or a collection of atoms."""
    with open(path, "w") as f:
        f.write(to_mmCIF(atom_or_collection))
