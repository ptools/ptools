
"""linalg - Defines functions to read/write files."""

import os
import sys

from .pdb import read_pdb

def check_file_exists(path, message=False):
    """Checks that a file exists.

    Returns a boolean and optionaly prints an error message if file does
    not exist.

    Args:
        path (str): path to file
        message (bool/str): if True, prints default a message if file does not exist.
            If a string, print the string if file does not exists.
        exitstatus (int): exit status
        prefix (str): message prefix

    Returns:
        bool: True if file exists, False either.
    """
    exists = os.path.exists(path)
    if not exists and message:
        if message is True:
            print(f"ERROR: file not found: '{path}'")
        else:
            print(message.format(path))
    return exists


def assert_file_exists(path, message=""):
    """Makes sure a file does exists.

    Args:
        path (str): path to file that supposedly exists
        message (str): exception message (default: file path)

    Raises:
        FileNotFoundError: if file does not exists
        IsADirectoryError: if `path` is a directory
    """
    if not message:
        message = path
    if not os.path.exists(path):
        raise FileNotFoundError(message)
    if os.path.isdir(path):
        raise IsADirectoryError(message)


def backup_if_exists(source):
    """Creates a backup of a file if source already exists.
    
    Files will be renamed file.1, file.2, etc.
    """
    def new_backup_name():
        return os.path.join(dirname, f"{basename}.{idx}")

    if os.path.exists(source):
        dirname = os.path.dirname(source)
        basename = os.path.basename(source)
        idx = 1
        target = new_backup_name()

        while os.path.exists(target):
            idx += 1
            target = new_backup_name()
        os.rename(source, target)
