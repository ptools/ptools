from pathlib import Path
from typing import Any


import yaml


from ..io import assert_file_exists
from ..io.pdb import read_single_model_pdb
from ..atomcollection import AtomCollection

from .bead import Bead
from .exceptions import NoReductionRulesError
from .residue import Residue

PathLike = str | Path


# TODO: define this path somewhere else.
PTOOLS_DATA_PATH = Path(__file__).parent.parent / "data"

PTOOLS_ATOM_RENAME_RULES_PATH = PTOOLS_DATA_PATH / "atom_rename_rules.yml"
PTOOLS_REDUCTION_PARAMETERS_PATH = PTOOLS_DATA_PATH / "reduction_parameters"
DEFAULT_ATTRACT1_PROTEIN_REDUCTION_PARAMETERS_PATH = (
    PTOOLS_REDUCTION_PARAMETERS_PATH / "attract1_protein_reduction_parameters.yml"
)


class Reducer:

    all_atoms: AtomCollection
    reduction_parameters: dict[str, Any]
    atom_rename_map: dict[str, dict[str, str]]
    beads: list[Bead]

    def __init__(self, topology_file: PathLike, reduction_parameters_file: PathLike):
        self.all_atoms = read_single_model_pdb(topology_file)
        self.reduction_parameters = read_reduction_parameters(
            reduction_parameters_file
        )["beads"]
        self.atom_rename_map = read_name_conversion_rules(PTOOLS_ATOM_RENAME_RULES_PATH)
        self.beads: list[Bead] = []

    def reduce(self) -> None:
        self._rename_atoms_and_residues()
        residue_list = list(self.all_atoms.iter_residues())

        for residue in residue_list:
            try:
                coarse_residue = self._reduce_residue(residue)
                coarse_residue.check_composition()
            except Exception as e:
                print(e)
            else:
                self.beads.extend(coarse_residue.beads)

    def number_of_atoms(self) -> int:
        return len(self.all_atoms)

    def number_of_beads(self) -> int:
        return len(self.beads)

    def _reduce_residue(self, residue: AtomCollection) -> Residue:
        """Reduces a single residue."""
        residue_name = residue[0].residue_name
        residue_index = residue[0].residue_index

        # Reduction parameters for this residue.
        reduction_parameters = self.reduction_parameters.get(residue_name, {})

        # No reduction parameters: this is an error.
        if not reduction_parameters:
            raise NoReductionRulesError(residue_name, residue_index)

        return Residue(residue, reduction_parameters)

    def _rename_atoms_and_residues(self):
        """Renames atoms and residues according to the rules defined ``self.atom_rename_map``."""
        residue_rename_map = self.atom_rename_map["residues"]
        atom_rename_map = self.atom_rename_map["atoms"]

        for atom in self.all_atoms:
            # Residue renaming.
            if atom.residue_name in residue_rename_map:
                atom.residue_name = residue_rename_map[atom.residue_name]

            # Atom renaming.
            if "*" in atom_rename_map and atom.name in atom_rename_map["*"]:
                atom.name = atom_rename_map["*"][atom.name]

            atom_rename_ = atom_rename_map.get(atom.residue_name, {})
            if atom.name in atom_rename_:
                atom.name = atom_rename_[atom.name]


def read_reduction_parameters(path: PathLike) -> dict[str, Any]:
    """Reads reduction parameters from a YAML file."""
    assert_file_exists(path)
    with open(path) as f:
        reduction_parameters = yaml.safe_load(f)
    return reduction_parameters


def read_topology(path: PathLike) -> AtomCollection:
    """Reads a PDB topology file."""
    assert_file_exists(path)
    return read_single_model_pdb(path)


def read_name_conversion_rules(path: PathLike) -> dict[str, dict[str, str]]:
    """Reads a YAML file containing name conversion rules."""
    assert_file_exists(path)
    with open(path) as f:
        name_conversion_rules = yaml.safe_load(f)
    return name_conversion_rules
