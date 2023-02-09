from __future__ import annotations
import itertools
from typing import Any, Callable, Iterable, Optional, TypeVar, Type

import numpy as np
from numpy.typing import ArrayLike

from .atomattrs import guess_atom_element, guess_atom_mass
from .namedarray import NamedArrayContainer

ParticleType = TypeVar("ParticleType", bound="Particle")


class Particle:
    """Represents a single particle in a collection."""

    def __init__(self, collection, index):
        self._collection = collection
        self._index = index
        self._singular_to_plural = {
            prop.singular: prop.plural for prop in self._collection.atom_properties
        }

    def properties(self) -> dict[str, Any]:
        """ "Returns a dictionary of all particle properties."""
        return {k: getattr(self, k) for k in self._singular_to_plural.keys()}

    def copy(self: ParticleType) -> ParticleType:
        return self.__class__(self._collection, self._index)

    def __eq__(self, __o: object) -> bool:
        """Compares two particles using their properties read from the ParticleCollection."""
        if isinstance(__o, Particle):
            self_properties = self._collection.atom_properties[self._index]
            other_properties = __o._collection.atom_properties[__o._index]
            return self_properties == other_properties

        rhs = ParticleCollection(atoms=[__o])
        self_properties = self._collection.atom_properties[self._index]
        other_properties = rhs.atom_properties[0]
        return self_properties == other_properties

    def __getattr__(self, name):
        if name in self._singular_to_plural:
            return self._collection.atom_properties.get(self._singular_to_plural[name])[
                self._index
            ]
        raise AttributeError(f"No such attribute: {name!r}")

    def __setattr__(self, name, value):
        # Setting the attributes of the class itself.
        if name in ("_collection", "_index", "_singular_to_plural"):
            super().__setattr__(name, value)
            return

        # Setting the attributes of the particle via the collection.
        if name in self._singular_to_plural:
            self._collection.atom_properties.get(self._singular_to_plural[name])[
                self._index
            ] = value
        else:
            raise KeyError(f"No such property: {name!r}")

    def __str__(self) -> str:
        attrs = ", ".join(
            f"{k}={getattr(self, k)!r}" for k in self._singular_to_plural.keys()
        )
        return f"Particle({attrs})"


ParticleCollectionType = TypeVar("ParticleCollectionType", bound="ParticleCollection")


class ParticleCollection:
    """Represents a collection of particles."""

    class Selection:
        """A class that can be used to set the parent of a collection."""

        parent: ParticleCollection
        indices: np.ndarray

        def __init__(self, parent: ParticleCollection, indices: ArrayLike | slice):
            self.parent = parent
            if isinstance(indices, slice):
                indices = np.arange(*indices.indices(len(parent)), dtype=np.int64)
            self.indices = np.array(indices, dtype=np.int64)

    _atom_properties: NamedArrayContainer
    _selection: Selection | None

    # == Initialization =================================================================

    def __init__(
        self,
        atoms: Optional[Iterable[Any]] = None,
        selection: Optional[Selection] = None,
    ):
        """Initializes a new collection of particles from a list of particles."""
        self._selection = None
        self._atom_properties = NamedArrayContainer()

        if selection is not None and atoms is not None:
            raise ValueError("Cannot specify both selection and atoms.")

        if selection is not None:
            self._init_from_selection(selection)
        elif atoms is not None:
            self._init_from_atoms(atoms)

    def _init_from_selection(self, selection: Selection):
        assert selection is not None

        if len(selection.indices) != len(np.unique(selection.indices)):
            raise ValueError("Indices must be unique.")

        self._selection = selection

    def _init_from_atoms(self, atoms: Iterable[Any]):
        assert atoms is not None
        self._atom_properties = NamedArrayContainer.from_objects(atoms)

    # ===================================================================================
    def has_parent(self):
        """Returns whether the collection has a parent (i.e. is a sub-collection)."""
        return self._selection is not None

    @property
    def atom_properties(self) -> NamedArrayContainer:
        """Returns the properties of the atoms."""
        if self.has_parent():
            return self._selection.parent.atom_properties[self._selection.indices]  # type: ignore
        return self._atom_properties

    @atom_properties.setter
    def atom_properties(self, value: NamedArrayContainer):
        if self.has_parent():
            raise NotImplementedError("Cannot set atom properties on a sub-collection.")
        self._atom_properties = value

    def __eq__(self, __o: object) -> bool:
        """Compares two collections using their properties."""
        if not isinstance(__o, ParticleCollection):
            return NotImplemented
        return self.atom_properties == __o.atom_properties

    def __len__(self) -> int:
        """Returns the number of atoms in the collection.

        This is an alias for ``ParticleCollection.size()``.
        """
        return self.size()

    def __getitem__(self, key: int | slice) -> Particle | ParticleCollection:
        """Returns a new collection with the selected atoms."""
        if isinstance(key, (int, np.integer)):
            return Particle(self, key)
        return ParticleCollection(selection=self.__class__.Selection(self, key))

    def __iter__(self):
        """Iterates over the atoms."""
        return iter(self.particles)

    def __contains__(self, particle: object) -> bool:
        """Returns whether the given particle is part of the collection."""
        return particle in self.particles

    def __add__(
        self: ParticleCollectionType, other: ParticleCollection
    ) -> ParticleCollectionType:
        """Adds two collections together."""
        return self.__class__.from_properties(
            self.atom_properties + other.atom_properties
        )

    def __getattr__(self, name):
        if name in ("_atom_properties", "atom_properties"):
            return super().__getattribute__(name)

        atom_properties = super().__getattribute__("atom_properties")
        if name in atom_properties:
            return self.atom_properties.get(name).values
        raise AttributeError(f"{self.__class__.__name__} has no attribute: {name!r}")

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} with {self.size()} particles>"

    # == Factory methods ================================================================

    @classmethod
    def from_properties(
        cls: type[ParticleCollectionType], properties: NamedArrayContainer
    ) -> ParticleCollectionType:
        """Creates a new collection from a list of properties."""
        obj = cls()
        obj.atom_properties = properties
        return obj

    # ===================================================================================

    def size(self) -> int:
        """Returns the number of atoms in the collection."""
        if self.has_parent():
            return len(self._selection.indices)  # type: ignore
        if self.atom_properties:
            return self.atom_properties.number_of_elements()
        return 0

    # == Accessors ======================================================================

    @property
    def coordinates(self):
        """Returns the coordinates of the atoms."""
        return self.atom_properties.get("coordinates").values

    @coordinates.setter
    def coordinates(self, value: ArrayLike):
        """Sets the coordinates of the atoms."""
        self.atom_properties.set("coordinates", value)

    @property
    def particles(self):
        """Returns the list of particles in the collection."""
        return list(Particle(self, i) for i in range(self.size()))

    # ===================================================================================

    def copy(self: ParticleCollectionType) -> ParticleCollectionType:
        """Returns a copy of the collection."""
        return self.__class__.from_properties(self.atom_properties.copy())

    def guess_elements(self):
        """Guesses the elements of the atoms."""
        values = [guess_atom_element(a.name) for a in self.particles]
        if "elements" not in self.atom_properties:
            self.atom_properties.add_array("element", "elements", values)
        else:
            self.atom_properties.get("elements").values = values

    def guess_masses(self, guess_element: bool = True):
        """Guesses the masses of the atoms."""
        if guess_element:
            self.guess_elements()

        values = np.array([guess_atom_mass(a.element) for a in self.particles])
        if "masses" not in self.atom_properties:
            self.atom_properties.add_array("mass", "masses", values)
        else:
            self.atom_properties.get("masses").values = values

    # == Selection methods ==============================================================

    def groupby(
        self: ParticleCollectionType, key: Callable
    ) -> dict[Any, ParticleCollectionType]:
        data = sorted(self, key=key)
        grouped = itertools.groupby(data, key=key)
        return {key: self.__class__(list(group)) for key, group in grouped}

    def select_atom_type(
        self: ParticleCollectionType, atom_type: str
    ) -> ParticleCollectionType:
        """Returns a new collection with the selected atom type."""
        raise NotImplementedError

    def select_atom_types(
        self: ParticleCollectionType, atom_types: Iterable[str]
    ) -> ParticleCollectionType:
        """Returns a new collection with the selected atom types."""
        raise NotImplementedError

    def select_residue_range(
        self: ParticleCollectionType, start: int, end: int
    ) -> ParticleCollectionType:
        """Returns a new collection with the selected residue range."""
        raise NotImplementedError

    def select_chain(
        self: ParticleCollectionType, chain: str
    ) -> ParticleCollectionType:
        """Returns a new collection with the selected chain."""
        raise NotImplementedError

    # == Iteration methods ==============================================================

    # def iter_particles(self) -> Iterator[Particle]:
    #     return iter(self)

    # def iter_residues(self) -> Iterator[ParticleCollectionType]:
    #     by_residue = self.groupby(lambda atom: (atom.residue_index, atom.chain))
    #     return iter(by_residue.values())

    # def iter_chains(self) -> Iterator[ParticleCollectionType]:
    #     by_chain = self.groupby(lambda atom: atom.chain)
    #     return iter(by_chain.values())