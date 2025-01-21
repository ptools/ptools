from collections.abc import Iterable
from typing import Any, Protocol


class mmCIFConvertible(Protocol):
    """Protocol for objects that can be converted to mmCIF format.

    Will convert only the coordinates section.

    Those objects must have the following attributes:

    - coordinates: a 3-tuple of floats
    - name: a string, atom name
    - index: an integer, atom index
    - residue_name: a string, residue name
    - residue_index: an integer, residue index
    - chain: a string, chain name
    - element: a string, element name
    """

    @property
    def coordinates(self) -> tuple[float, float, float]: ...

    @property
    def name(self) -> str: ...

    @property
    def index(self) -> int: ...

    @property
    def residue_name(self) -> str: ...

    @property
    def residue_index(self) -> int: ...

    @property
    def chain(self) -> str: ...

    @property
    def typeid(self) -> int: ...

    @property
    def charge(self) -> float: ...


# The format is just for aesthetic purposes
MMCIF_FORMAT = (
    "{record:<6s} {atom_index:6d} {element:1s} {atom_name:6s} {altloc:1s} "
    "{residue_name:4s} {chain:3s} {entity_id:1s} {residue_index:6d} {icode:1s} "
    "{x:8.3f} {y:8.3f} {z:8.3f} {occupancy:6.2f} {bfactor:6.2f} {charge:1s} "
    "{residue_index:6d} {residue_name:4s} {chain:3s} {atom_name:6s} {num_model:1d}"
)


def to_mmCIF(
    atom_or_collection: mmCIFConvertible | Iterable[mmCIFConvertible],
    mmCIF_name: str = "Ptools_output",
) -> str:
    """Converts an atom or a collection of atoms to mmCIF format.

    Will convert only the coordinates section.

    Parameters
    ----------
    atom_or_collection : mmCIFConvertible | Iterable[mmCIFConvertible]
        an atom or a collection of atoms.
    mmCIF_name : str, optional
        name of mmCIF data block , by default "Ptools_output"

    Returns
    -------
    str
        the content of the mmCIF data
    """

    mmCIF_header = f"data_{mmCIF_name}" + "\n" + format_header() + "\n"

    if isinstance(atom_or_collection, Iterable):
        return mmCIF_header + "\n".join(format_atom(atom) for atom in atom_or_collection)
    return mmCIF_header + format_atom(atom_or_collection)  # atom_or_collection is a single atom


def format_header() -> str:
    """Returns the mmCIF coordinates headers

    Returns
    -------
    str
        the headers
    """
    headers = (
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_formal_charge",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_model_num",
    )

    return "\n".join(headers)


def format_atom(atom: mmCIFConvertible) -> str:
    """Format one atom to the mmCIF format

    Parameters
    ----------
    atom : mmCIFConvertible
        an atom

    Returns
    -------
    str
        data's atom in the mmCIF format
    """
    bfactor = getattr(atom, "bfactor", 0.0)
    occupancy = getattr(atom, "occupancy", 1.0)
    element = getattr(atom, "element", "")
    fmt_args: dict[str, Any] = {
        "record": "ATOM",
        "atom_index": atom.index,
        "element": element,
        "atom_name": atom.name,
        "altloc": "?",
        "residue_name": atom.residue_name,
        "chain": atom.chain,
        "entity_id": "?",
        "residue_index": atom.residue_index,
        "icode": "?",
        "x": atom.coordinates[0],
        "y": atom.coordinates[1],
        "z": atom.coordinates[2],
        "occupancy": occupancy,
        "bfactor": bfactor,
        "charge": "?",
        "num_model": 1,
    }
    return MMCIF_FORMAT.format(**fmt_args)
