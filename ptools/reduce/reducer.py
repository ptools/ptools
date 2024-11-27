import logging
from pathlib import Path
from typing import Any, Container, Optional, Type

import yaml

from ..io import assert_file_exists
from ..io.readers.pdb import read_single_model_pdb
from ..particlecollection import ParticleCollection
from .bead import Bead
from .exceptions import EmptyModelError, NoReductionRulesError
from .residue import Residue

PathLike = str | Path
ExceptionTypeContainer = Container[Type[Exception]]


# TODO: define this path somewhere else.
PTOOLS_DATA_PATH = Path(__file__).parent.parent / "data"

DEFAULT_ATOM_RENAME_RULES_PATH = PTOOLS_DATA_PATH / "atom_rename_rules.yml"
PTOOLS_REDUCTION_PARAMETERS_DIR = PTOOLS_DATA_PATH / "reduction_parameters"


# Maps force field name to the default reduction parameters file.
FORCEFIELDS = {
    "attract1": PTOOLS_REDUCTION_PARAMETERS_DIR / "attract1_reduction_parameters.yml",
    "attract2": PTOOLS_REDUCTION_PARAMETERS_DIR / "attract2_reduction_parameters.yml",
    "scorpion": PTOOLS_REDUCTION_PARAMETERS_DIR / "scorpion_reduction_parameters.yml",
}


logger = logging.getLogger(__name__)


class Reducer:
    all_atoms: ParticleCollection
    reduction_parameters: dict[str, Any]
    atom_rename_map: dict[str, dict[str, str]]
    beads: list[Bead]

    def __init__(
        self,
        topology_file_or_atoms: PathLike | ParticleCollection,
        reduction_parameters_file: PathLike,
        name_conversion_file: PathLike,
    ):
        if isinstance(topology_file_or_atoms, Path):
            self.all_atoms = read_single_model_pdb(topology_file_or_atoms)
        else:
            self.all_atoms = topology_file_or_atoms

        self.reduction_parameters = read_reduction_parameters(
            reduction_parameters_file
        )["beads"]
        self.atom_rename_map = read_name_conversion_rules(name_conversion_file)
        self.beads: list[Bead] = []

    def number_of_atoms(self) -> int:
        return len(self.all_atoms)

    def number_of_beads(self) -> int:
        return len(self.beads)

    def reduce(
        self,
        ignore_exceptions: Optional[ExceptionTypeContainer] = None,
        warn_exceptions: Optional[ExceptionTypeContainer] = None,
    ) -> None:
        """Actual reduction routine."""
        warn_exceptions = warn_exceptions or []
        ignore_exceptions = ignore_exceptions or []

        self._rename_atoms_and_residues()

        # Reduces each residue.
        for residue in self.all_atoms.iter_residues():  # type: ignore[var-annotated]
            coarse_residue = self._reduce_residue(residue)
            try:
                coarse_residue.check_composition()
            except Exception as error:
                if type(error) in warn_exceptions:
                    logger.warning("%s", error)
                elif type(error) not in ignore_exceptions:
                    raise error
            self.beads.extend(coarse_residue.beads)

        # Properly sets the bead indices.
        for i, bead in enumerate(self.beads):
            bead.index = i + 1

        if not self.beads:
            raise EmptyModelError()

    def _reduce_residue(self, residue: ParticleCollection) -> Residue:
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

    def get_reduced_model(self) -> ParticleCollection:
        """Returns bead atoms concatenated into a single ParticleCollection."""
        atoms = self.beads[0].atoms
        print(atoms)
        pc = ParticleCollection(atoms)
        print(pc)

        exit()
        atoms = []
        for bead in self.beads:
            atoms.extend(bead.atoms)
        return ParticleCollection(atoms)


def read_reduction_parameters(path: PathLike) -> dict[str, Any]:
    """Reads reduction parameters from a YAML file."""
    assert_file_exists(path)
    with open(path) as f:
        reduction_parameters = yaml.safe_load(f)
    return reduction_parameters


def read_topology(path: PathLike) -> ParticleCollection:
    """Reads a PDB topology file."""
    assert_file_exists(path)
    return read_single_model_pdb(path)


def read_name_conversion_rules(path: PathLike) -> dict[str, dict[str, str]]:
    """Reads a YAML file containing name conversion rules."""
    assert_file_exists(path)
    with open(path) as f:
        name_conversion_rules = yaml.safe_load(f)
    return name_conversion_rules
