"""Protein Data Bank format I/O."""

from typing import Sequence, Tuple, Union

from ...atomattrs import AtomAttrs
from ...particlecollection import ParticleCollection
from ..._typing import FilePath
from ..exceptions import InvalidPDBFormatError, InvalidPDBAtomLineError


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

    def __init__(self, *args) -> None:
        super().__init__(*args)
        if not self.is_atom():
            raise InvalidPDBAtomLineError(self.header)

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
        return self[12:16].strip()

    @property
    def chain(self) -> str:
        """Chain identifier."""
        return self[21].strip()

    @property
    def residue_name(self) -> str:
        """Atom residue name."""
        return self[17:20].strip()

    @property
    def residue_index(self) -> int:
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

    def to_atom(self) -> AtomAttrs:
        """Creates a AtomAttrs from an atom line."""
        return AtomAttrs(
            index=self.atom_id,
            name=self.name,
            residue_name=self.residue_name,
            residue_index=self.residue_index,
            chain=self.chain,
            coordinates=self.coordinates,
            meta={"extra": self.extra},
        )


def parse_atom_line(buffer: str) -> AtomAttrs:
    """Returns an `atom.AtomAttrs` initialized with data read from line."""
    line = AtomLine(buffer)
    if not line.is_atom():
        raise ValueError(f"Not a valid atom line: header is {line.header}")
    return line.to_atom()


def read_single_model_pdb(path: FilePath) -> ParticleCollection:
    """Read a Protein Data Bank file containing a single model.

    Args:
        path (FilePath): path to file.

    Returns:
        ParticleCollection: collection of Atoms
    """
    topology = read_pdb(path)
    if not isinstance(topology, ParticleCollection):
        raise ValueError("Topology file must contain only one model.")
    return topology


def read_pdb(
    path: FilePath, as_dict=False
) -> Union[
    dict[str, ParticleCollection], Sequence[ParticleCollection], ParticleCollection
]:
    """Read a Protein Data Bank file.

    Args:
        path (FilePath): path to file.
        as_dict (bool): if True, returns models in a dictionnary.

    Returns:
        ParticleCollection: collection of Atoms
    """

    def register_model(atom_list: Sequence[AtomAttrs]):
        """Stores `atom_list` into `models` as an ParticleCollection."""
        models.append(ParticleCollection(atom_list))

    def register_model_id(line: str):
        """Extracts model id from model header line and stores `model_id` into `model_id_list`."""
        model_id_list.append(ModelLine(line).model_id)

    def register_new_atom(line: str):
        """Parses an ATOM line and stores the atom into the current model."""
        current_model.append(AtomLine(line).to_atom())

    models: list[ParticleCollection] = []
    model_id_list: list[str] = []
    current_model: list[AtomAttrs] = []

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

    # Only 1 model: returns a simple ParticleCollection
    if len(models) == 1:
        return models[0]

    # Multiple models: returns a list of ParticleCollection instances.
    return models
