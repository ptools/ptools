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
    models = []
    model_id_list = []
    current_model = []
    with open(path, "rt", encoding="utf-8") as f:
        for line in f:
            header = get_header(line)
            if header == "MODEL":
                model_id_list.append(line[10:].strip())
            if header == "ENDMDL":
                models.append(AtomCollection(current_model))
                current_model = []
            elif is_atom_line(line):
                current_model.append(parse_atom_line(line))
    if models:
        if current_model:  # No "ENDMDL" flag for last model.
            models.append(AtomCollection(current_model))
        if len(models) == 1:
            return AtomCollection(models[0])
        if as_dict:
            return dict(zip(model_id_list, models))
        return models
    return AtomCollection(current_model)
