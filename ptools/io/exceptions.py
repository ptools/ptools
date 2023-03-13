"""Defines the exceptions used in ``ptoos.io``."""


class InvalidPDBFormatError(IOError):
    """Raised when the PDB format is incorrect."""


class InvalidPDBAtomLineError(ValueError):
    """Raised when the a PDB line does not describe an atom."""

    def __init__(self, header, *args):
        message = f'Invalid atom line: expected "ATOM" or "HETATM" as header, found "{header}"'
        super().__init__(message, *args)


class InvalidREDFormatError(IOError):
    """Raised when the RED format is incorrect."""
