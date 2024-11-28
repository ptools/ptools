import itertools
import logging

from pathlib import Path
from typing import Any, Container, Optional, Type, Protocol, Iterator, Iterable

import yaml

from ..io import assert_file_exists
from ..io.readers.pdb import read_single_model_pdb
from ..particlecollection import ParticleCollection
from .bead import Bead, BeadDump
from .exceptions import EmptyModelError, NoReductionRulesError, all_exceptions
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


class AtomType(Protocol):
    """Protocol for an atom type."""

    name: str
    residue_name: str
    residue_index: int
    chain: str


def iter_residues(atoms: Iterable[AtomType]) -> Iterator[list[AtomType]]:
    """Iterates over residues in a list of atoms."""
    sorted_atoms = sorted(atoms, key=lambda atom: (atom.residue_index, atom.chain))
    grouped = itertools.groupby(sorted_atoms, key=lambda atom: (atom.residue_index, atom.chain))
    for _, group in grouped:
        yield list(group)


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
            self.all_atoms = read_single_model_pdb(topology_file_or_atoms).dump()
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
        for residue in iter_residues(self.all_atoms):
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

    @property
    def reduced_model(self) -> list[BeadDump]:
        """Returns bead atoms concatenated into a single ParticleCollection."""
        beads = [bead.dump() for bead in self.beads]
        return beads


def reduce(
        atoms: ParticleCollection,
        forcefield: str = "attract1",
        renaming_rules: Path = DEFAULT_ATOM_RENAME_RULES_PATH,
        warn_exceptions: list[Exception] = [],
        ignore_exceptions: list[Exception] = all_exceptions(),
) -> list[BeadDump]:
    """Reduces an all-atom model to a coarse-grained model."""
    if forcefield not in FORCEFIELDS:
        raise ValueError(f"Force field {forcefield!r} is not supported "
                         f"(choices are {', '.join(f'{name!r}' for name in FORCEFIELDS.keys())}).")
    reduction_parameters_file = FORCEFIELDS[forcefield]
    reducer = Reducer(atoms, reduction_parameters_file, renaming_rules)
    reducer.reduce(ignore_exceptions, warn_exceptions)
    return reducer.reduced_model


def read_reduction_parameters(path: PathLike) -> dict[str, Any]:
    """Reads reduction parameters from a YAML file."""
    assert_file_exists(path)
    with open(path, "rt", encoding="utf-8") as f:
        reduction_parameters = yaml.safe_load(f)
    return reduction_parameters


def read_topology(path: PathLike) -> ParticleCollection:
    """Reads a PDB topology file."""
    assert_file_exists(path)
    return read_single_model_pdb(path)


def read_name_conversion_rules(path: PathLike) -> dict[str, dict[str, str]]:
    """Reads a YAML file containing name conversion rules."""
    assert_file_exists(path)
    with open(path, "rt", encoding="utf-8") as f:
        name_conversion_rules = yaml.safe_load(f)
    return name_conversion_rules
