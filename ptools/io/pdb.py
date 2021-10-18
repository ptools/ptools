"""Protein Data Bank format I/O."""

from typing import List, Union
from ..atom import BaseAtom, AtomCollection


class InvalidPDBFormatError(IOError):
    """Raised when the PDB format is incorrect."""


def get_header(line: str) -> str:
    """Returns PDB line header i.e. first 6 characters stripped from white spaces."""
    return line[:6].strip()


def is_atom_line(line: str) -> bool:
    """Returns True is the line header is "ATOM  " or "HETATM"."""
    return get_header(line) in ("ATOM", "HETATM")


def parse_atom_line(line: str) -> BaseAtom:
    """Returns an `atom.BaseAtom` initialized with data read from line."""
    index = int(line[6:11])
    name = line[12:16].strip().upper()
    resname = line[17:20].strip().upper()
    chain = line[21].strip()
    resid = int(line[22:26])
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    coords = (x, y, z)
    extra = line[54:].strip()

    atom = BaseAtom(
        index=index,
        name=name,
        resname=resname,
        resid=resid,
        chain=chain,
        coords=coords,
        meta={"extra": extra},
    )
    return atom


def read_pdb(path: str, as_dict=False) -> Union[List[AtomCollection], AtomCollection]:
    """Read a Protein Data Bank file.

    Args:
        path (str): path to file.
        as_dict (bool): if True, returns models in a dictionnary.

    Returns:
        AtomCollection: collection of Atoms
    """
    def register_model(atom_list: list[BaseAtom]):
        """Stores `atom_list` into `models` as an AtomCollection."""
        models.append(AtomCollection(atom_list))

    def register_model_id(line: str):
        """Extracts model id from model header line and stores `model_id` into `model_id_list`."""
        model_id_list.append(line[10:].strip())

    def register_new_atom(line: str):
        """Parses an ATOM line and stores the atom into the current model."""
        current_model.append(parse_atom_line(line))

    models = []
    model_id_list = []
    current_model = []
    with open(path, "rt", encoding="utf-8") as f:
        for line in f:
            header = get_header(line)
            if header == "MODEL":
                register_model_id(line)
            if header == "ENDMDL":
                register_model(current_model)
                current_model.clear()
            elif is_atom_line(line):
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
