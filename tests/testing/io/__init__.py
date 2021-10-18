"""testing.io - I/O testing utilities."""

import tempfile

from .temporaryfile import mk_empty_file, mk_tmp_file
from .pdb import mk_pdb_10_atoms, mk_pdb_no_model, mk_pdb_models, TestPDBBuilder
from .red import mk_red

__all__ = [
    "mk_empty_file",
    "mk_tmp_file",
    "mk_pdb_10_atoms",
    "mk_pdb_no_model",
    "mk_pdb_models",
    "TestPDBBuilder",

    "mk_red"
]


def random_filename() -> str:
    """Returns a random file name."""
    with tempfile.NamedTemporaryFile() as tmpfile:
        return tmpfile.name
