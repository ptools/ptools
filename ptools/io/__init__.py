"""linalg - Defines functions to read/write files."""

import pathlib

from .._typing import FilePath
from .formatters.pdb import to_pdb as to_pdb
from .formatters.mmcif import to_mmCIF as to_mmCIF
from .formatters.reduced import (
    to_reduced_pdb as to_reduced_pdb,
    write_reduced_pdb as write_reduced_pdb,
)
from .readers.pdb import read_pdb as read_pdb
from .readers.attract import read_docking_parameters as read_attract_docking_parameters
from .readers.attract import read_topology as read_attract_topology
from .writers.pdb import write_pdb as write_pdb
from .writers.mmcif import write_mmCIF as write_mmCIF


def check_file_exists(path: FilePath, message: bool | str = False) -> bool:
    """Checks that a file exists.

    Returns a boolean and optionaly prints an error message if file does
    not exist.

    Args:
        path (FilePath): path to file
        message (bool | str): if True, prints default a message if file does not exist.
            If a string, print the string if file does not exists.

    Returns:
        bool: True if file exists, False either.
    """
    path = pathlib.Path(path)
    if not path.exists() and message:
        if message is True:
            print(f"ERROR: file not found: '{path}'")
        else:
            print(message.format(path))
    return path.exists()


def assert_file_exists(path: FilePath, message: str = ""):
    """Makes sure a file does exists.

    Args:
        path (str): path to file that supposedly exists
        message (str): exception message (default: file path)

    Raises:
        FileNotFoundError: if file does not exists
        IsADirectoryError: if ``path`` is a directory
    """
    path = pathlib.Path(path)
    if not message:
        message = f"{path}"
    if not path.exists():
        raise FileNotFoundError(message)
    if path.is_dir():
        raise IsADirectoryError(message)


def backup_if_exists(source: FilePath):
    """Creates a backup of a file if source already exists.

    Files will be renamed file.1, file.2, etc.
    """
    source = pathlib.Path(source)
    if source.exists():
        idx = 1
        target = source.with_suffix(source.suffix + f".{idx}")
        while target.exists():
            idx += 1
            target = source.with_suffix(source.suffix + f".{idx}")
        source.rename(target)
