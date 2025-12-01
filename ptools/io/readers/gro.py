"""Protein Data Bank format I/O."""

from typing import TextIO

import tqdm

from ..._typing import FilePath, PathLike
from ...atomattrs import AtomAttrs
from ...particlecollection import ParticleCollection


class GroFormatError(Exception):
    pass


class AtomLine:
    def __init__(self, line: str, line_number: int):
        self.line = line[:-1]  # strip carriage return
        self.line_number = line_number

        if len(line) < 44:
            raise GroFormatError(
                f"{self.line_number}: line too short: expected 44, found {len(line)}"
            )

    @property
    def residue_id(self) -> int:
        token = self.line[0:5]
        try:
            return int(token)
        except ValueError as e:
            raise ValueError(f"{self.line_number}: expected residue id, found '{token}'") from e

    @property
    def residue_name(self) -> str:
        return self.line[5:10].strip()

    @property
    def atom_name(self) -> str:
        return self.line[10:15].strip()

    @property
    def atom_id(self) -> int:
        token = self.line[15:20]
        try:
            return int(token)
        except ValueError as e:
            raise ValueError(f"{self.line_number}: expected atom id, found '{token}'") from e

    @property
    def x(self) -> float:
        """X-coordinate"""
        return float(self.line[20:28])

    @property
    def y(self) -> float:
        """Y-coordinate"""
        return float(self.line[28:36])

    @property
    def z(self) -> float:
        """Z-coordinate"""
        return float(self.line[36:44])

    @property
    def coordinates(self) -> tuple[float, float, float]:
        """Coordinates."""
        return (self.x, self.y, self.z)

    def to_atom(self) -> AtomAttrs:
        return AtomAttrs(
            name=self.atom_name,
            index=self.atom_id,
            residue_name=self.residue_name,
            residue_index=self.residue_id,
            coordinates=self.coordinates,
        )


def _parse_number_of_atoms(fileobj: TextIO, path: PathLike) -> int:
    try:
        line = next(fileobj)
    except StopIteration as e:
        raise GroFormatError(f"{path}: line 2: excepted the number of atoms, found nothing") from e
    try:
        return int(line)
    except ValueError as e:
        raise GroFormatError(
            f"{path}: line 2: excepted the number of atoms, found '{line[:-1]}'"
        ) from e


def read_gro(path: FilePath, show_progress: bool = False) -> ParticleCollection:
    """Read a Protein Data Bank file.

    Args:
        path (FilePath): path to file.
        show_progress (bool): displays a progress bar

    Returns:
        ParticleCollection: collection of Atoms
    """

    model: list[AtomAttrs] = []

    with open(path, encoding="utf-8") as f:
        next(f)  # header
        n_atoms = _parse_number_of_atoms(f, path)
        iterator = range(1, n_atoms + 1)
        if show_progress:
            iterator = tqdm.tqdm(iterator)
        for i in iterator:
            line_number = i + 2

            # Should be able to catch as many lines as the number of atoms we read
            try:
                buffer = next(f)
            except StopIteration as e:
                raise GroFormatError(
                    f"{path}: line {line_number}: expected an atom line, reached end of file"
                ) from e

            line = AtomLine(buffer, line_number)
            model.append(line.to_atom())

    return ParticleCollection(model)
