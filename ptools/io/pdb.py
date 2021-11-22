"""Protein Data Bank format I/O."""

import abc
from dataclasses import dataclass, field
import os
import re
from typing import Protocol

from typing import Sequence, Tuple, Union
from ..atom import BaseAtom, AtomCollection


class _FromPDBBase(Protocol):
    """Abstract Base class for classes which can be initialized from PDB files.

    Ensures that child classes implement the `init_from_pdb` method.
    """

    @abc.abstractmethod
    def init_from_pdb(self, path: Union[str, bytes, os.PathLike]):
        """Initializes internal components from data read in PDB file."""


class FromPDB(_FromPDBBase):
    """Base class for classes which can be initialized from PDB files.

    Implements FromPDB.from_pdb which returns a new instance initialized
    from a PDB file, read with `cls.init_from_pdb`.
    """

    @classmethod
    def from_pdb(cls, path: Union[str, bytes, os.PathLike]):
        """Returns a new instance of the class initialized using `cls.init_from_pdb`."""
        rigid = cls()
        rigid.init_from_pdb(path)
        return rigid


class InvalidPDBFormatError(IOError):
    """Raised when the PDB format is incorrect."""


class PDBLine(str):
    """Generic PDB formatted line.

    Parses the line header and provides helper methods to determine the type of line
    (i.e. atom, new model, end model, etc.).
    """

    def __init__(self, *args) -> None:
        super().__init__()
        self.header = self[0:6].strip()

    def is_atom(self) -> bool:
        """Returns True is the line describes an atom."""
        return self.header in ("ATOM", "HETATM")

    def is_new_model(self) -> bool:
        """Returns True if the line describes a new model."""
        return self.header == "MODEL"

    def is_end_model(self) -> bool:
        """Returns True if the line describes the end of the current model."""
        return self.header == "ENDMDL"


class ModelLine(PDBLine):
    """Helper methods for a model line."""

    @property
    def model_id(self) -> str:
        """Model identifier."""
        return self[10:].strip()


class AtomLine(PDBLine):
    """Helper methods for an atom line."""

    @property
    def atom_index(self) -> int:
        """Atom index.

        Alias for AtomLine.atom_id.
        """
        return self.atom_id

    @property
    def atom_id(self):
        """Atom identifier.

        See also:
            AtomLine.atom_index
        """
        return int(self[6:11])

    @property
    def name(self) -> str:
        """Atom name."""
        return self[12:16].strip().upper()

    @property
    def resname(self) -> str:
        """Atom residue name."""
        return self[17:20].strip().upper()

    @property
    def chain(self) -> str:
        """Chain identifier."""
        return self[21].strip()

    @property
    def resid(self) -> int:
        """Residue identifer.

        Alias for AtomLine.residue_id.
        """
        return self.residue_id

    @property
    def residue_id(self) -> int:
        """Residue identifer.

        See also:
            AtomLine.resid
        """
        return int(self[22:26])

    @property
    def x(self) -> float:
        """X-coordinate"""
        return float(self[30:38])

    @property
    def y(self) -> float:
        """Y-coordinate"""
        return float(self[38:46])

    @property
    def z(self) -> float:
        """Z-coordinate"""
        return float(self[46:54])

    @property
    def coordinates(self) -> Tuple[float, float, float]:
        """Coordinates."""
        return (self.x, self.y, self.z)

    @property
    def extra(self) -> str:
        """Extra properties"""
        return self[54:].strip()

    def to_atom(self) -> BaseAtom:
        """Creates a BaseAtom from an atom line."""
        return BaseAtom(
            index=self.atom_id,
            name=self.name,
            resname=self.resname,
            resid=self.residue_id,
            chain=self.chain,
            coords=self.coordinates,
            meta={"extra": self.extra},
        )


def parse_atom_line(buffer: str) -> BaseAtom:
    """Returns an `atom.BaseAtom` initialized with data read from line."""
    line = AtomLine(buffer)
    if not line.is_atom():
        raise ValueError(f"Not a valid atom line: header is {line.header}")
    return line.to_atom()


def read_pdb(
    path: str, as_dict=False
) -> Union[dict[str, AtomCollection], Sequence[AtomCollection], AtomCollection]:
    """Read a Protein Data Bank file.

    Args:
        path (str): path to file.
        as_dict (bool): if True, returns models in a dictionnary.

    Returns:
        AtomCollection: collection of Atoms
    """

    def register_model(atom_list: Sequence[BaseAtom]):
        """Stores `atom_list` into `models` as an AtomCollection."""
        models.append(AtomCollection(atom_list))

    def register_model_id(line: str):
        """Extracts model id from model header line and stores `model_id` into `model_id_list`."""
        model_id_list.append(ModelLine(line).model_id)

    def register_new_atom(line: str):
        """Parses an ATOM line and stores the atom into the current model."""
        current_model.append(AtomLine(line).to_atom())

    models = []
    model_id_list = []
    current_model = []
    with open(path, "rt", encoding="utf-8") as f:
        for buffer in f:
            line = PDBLine(buffer)
            if line.is_new_model():
                register_model_id(line)
            elif line.is_end_model():
                register_model(current_model)
                current_model.clear()
            elif line.is_atom():
                register_new_atom(line)

    # No "ENDMDL" flag for last model.
    if current_model:
        register_model(current_model)

    # Returns models into a dictionary.
    if as_dict:
        if len(model_id_list) == 0:
            raise InvalidPDBFormatError(
                "can't initialize dictionary without model identifier (no MODEL found)"
            )
        return dict(zip(model_id_list, models))

    # Only 1 model: returns a simple AtomCollection
    if len(models) == 1:
        return AtomCollection(models[0])

    # Multiple models: returns a list of AtomCollection instances.
    return models


class RobustPDBBuilder:

    class _RecordBase(abc.ABC):
        """Generic Record Handler Class."""

        @property
        @abc.abstractmethod
        def pdb(self) -> str:
            """Returns data as a PDB formatted string."""

    @dataclass(frozen=True)
    class RecordBase(_RecordBase):
        """Base class for a Record formatter class.

        Implements the `check_formats` methods which calls `self.check_format_<member>()`
        for each member variable.
        """

        def __post_init__(self):
            """Post initialization method.

            Mandatory call to `self.check_formats()`.
            """
            # Probably a better way to do this using a metaclass.
            for member_name in self.__dict__:
                method = f"check_format_{member_name}"
                if not hasattr(self, method):
                    raise NotImplementedError(f"class {self.__class__.__name__} has no member '{method}'")
            self.check_formats()

        def check_formats(self):
            """Makes sure data is well formatted.

            Calls `check_format_<member>` for each member variable.
            """
            for member_name in self.__dict__:
                getattr(self, f"check_format_{member_name}")()



    @dataclass(frozen=True)
    class HeaderRecord(RecordBase):
        """Header Record.

        Record Format
        =============

        COLUMNS       DATA  TYPE     FIELD             DEFINITION
        ------------------------------------------------------------------------------------
        1 -  6       Record name    "HEADER"
        11 - 50       String(40)     classification    Classifies the molecule(s).
        51 - 59       Date           depDate           Deposition date. This is the date the coordinates were received at the PDB.
        63 - 66       IDcode         idCode            This identifier is unique within the PDB.
        """

        name: str = field(init=False, default="HEADER")
        classification: str
        date: str
        idcode: str

        classification_pattern = re.compile(r".{0,40}$")
        date_pattern = re.compile(r"[0-9]{2}-[A-Z]{3}-[0-9]{2}$")
        idcode_pattern = re.compile(r"\w{4}$")

        def check_format_classification(self):
            """Checks the format for member `classification`."""
            assert self.classification_pattern.match(self.classification) is not None

        def check_format_date(self):
            """Checks the format for member `date`."""
            assert self.date_pattern.match(self.date) is not None

        def check_format_idcode(self):
            """Checks the format for member `idcode`."""
            assert self.idcode_pattern.match(self.idcode) is not None


        @property
        def pdb(self) -> str:
            """Returns the record formatted as a valid PDB string."""
            return f"{self.name:<6}    {self.classification:<40}{self.date:<9}   {self.idcode:<3}"


