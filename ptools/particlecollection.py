from __future__ import annotations
import itertools
from typing import Any, Callable, Iterable, Iterator, Optional, Sequence, TypeVar, overload

import numpy as np
from numpy.typing import ArrayLike

from .atomattrs import guess_atom_element, guess_atom_mass
from .namedarray import NamedArrayContainer

ParticleType = TypeVar("ParticleType", bound="Particle")
ParticleCollectionType = TypeVar("ParticleCollectionType", bound="ParticleCollection")
ParticleCollectionKeyType = int | slice | Sequence[int] | np.ndarray


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
        if name not in self._singular_to_plural:
            raise KeyError(f"No such property: {name!r}")

        property_name = self._singular_to_plural[name]
        self._collection._set_particle_property(property_name, self._index, value)
        # self._collection.atom_properties.set_at(property_name, self._index, value)

    def __str__(self) -> str:
        attrs = ", ".join(
            f"{k}={getattr(self, k)!r}" for k in self._singular_to_plural.keys()
        )
        return f"Particle({attrs})"


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

    def _set_particle_property(self, name: str, index: int, value: Any):
        """Sets the value of a property of a particle."""
        if self.has_parent():
            self._selection.parent.atom_properties.set_at(  # type: ignore[union-attr]
                name, self._selection.indices[index], value  # type: ignore[union-attr]
            )
            return
        self._atom_properties.set_at(name, index, value)


    @property
    def serial(self) -> np.ndarray:
        """Returns the serial numbers of the atoms, i.e. location in the parent collection."""
        if self.has_parent():
            return self._selection.indices
        return np.arange(self.size())


    # ===================================================================================
    def has_parent(self):
        """Returns whether the collection has a parent (i.e. is a sub-collection)."""
        return self._selection is not None

    @property
    def parent(self) -> ParticleCollection | None:
        """Returns the parent collection of the sub-collection."""
        if self.has_parent():
            return self._selection.parent
        return None

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

    @overload
    def __getitem__(
        self: ParticleCollectionType, key: int
    ) -> Particle:
        ...

    @overload
    def __getitem__(
        self: ParticleCollectionType, key: slice | Sequence[int] | np.ndarray
    ) -> ParticleCollectionType:
        ...

    def __getitem__(
        self: ParticleCollectionType, key: ParticleCollectionKeyType
    ) -> Particle | ParticleCollectionType:
        """Returns a new collection with the selected atoms."""
        if isinstance(key, (int, np.integer)):
            return Particle(self, key)
        return self.__class__(selection=self.__class__.Selection(self, key))

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

    # == Properties manipulation ========================================================
    def add_atom_property(
        self, singular: str, plural: str, value: Sequence | np.ndarray
    ):
        """Adds a property to the collection."""
        if self.has_parent():
            raise NotImplementedError("Cannot add properties to a sub-collection.")
        self.atom_properties.add_array(singular, plural, value)

    def remove_atom_property(self, plural: str):
        """Removes a property from the collection."""
        if self.has_parent():
            raise NotImplementedError("Cannot remove properties from a sub-collection.")
        self.atom_properties.remove_array(plural)

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
            self.atom_properties.set("masses", values)

    # == Selection methods ==============================================================

    def select_atom_type(self: ParticleCollectionType, atom_type: str
                         ) -> ParticleCollectionType:
        """Returns a new collection with the selected atom type."""
        indices = np.where(self.atom_properties.get("names").values == atom_type)[0]
        return self[indices]

    def select_atom_types(
        self: ParticleCollectionType, atom_types: Iterable[str]
    ) -> ParticleCollectionType:
        """Returns a new collection with the selected atom types."""
        indices = np.where(
            np.isin(self.atom_properties.get("names").values, atom_types)  # type: ignore[arg-type]
        )[0]
        return self[indices]

    def select_residue_indices(self: ParticleCollectionType, residues: ArrayLike) -> ParticleCollectionType:
        """Returns a new collection with the selected residues."""
        indices = np.where(np.isin(self.atom_properties.get("residue_indices").values, residues))[0]
        return self[indices]

    def select_residue_range(
        self: ParticleCollectionType, start: int, end: int
    ) -> ParticleCollectionType:
        """Returns a new collection with the selected residue range."""
        indices = np.where(
            np.logical_and(
                self.atom_properties.get("residue_indices").values >= start,
                self.atom_properties.get("residue_indices").values <= end,
            )
        )[0]
        return self[indices]

    def select_chain(
        self: ParticleCollectionType, chain: str
    ) -> ParticleCollectionType:
        """Returns a new collection with the selected chain."""
        indices = np.where(self.atom_properties.get("chains").values == chain)[0]
        return self[indices]

    # == Grouping methods ===============================================================
    def groupby(
        self: ParticleCollectionType, key: Callable
    ) -> dict[Any, ParticleCollectionType]:
        """Groups the atoms by the given key."""
        sorted_atoms = sorted(enumerate(self), key=lambda i_atom: key(i_atom[1]))
        grouped = itertools.groupby(sorted_atoms, key=lambda i_atom: key(i_atom[1]))
        grouped_collections = {
            dict_key: self[[i for i, _ in group]] for dict_key, group in grouped
        }
        return grouped_collections

    # == Iteration methods ==============================================================

    def iter_particles(self) -> Iterator[Particle]:
        """Iterates over the particles in the collection."""
        return iter(self)

    def iter_residues(self) -> Iterator[ParticleCollectionType]:
        """Iterates over the residues in the collection."""
        by_residue = self.groupby(lambda atom: (atom.residue_index, atom.chain))
        return iter(by_residue.values())  # type: ignore[arg-type]

    def iter_chains(self) -> Iterator[ParticleCollectionType]:
        """Iterates over the chains in the collection."""
        by_chain = self.groupby(lambda atom: atom.chain)
        return iter(by_chain.values())  # type: ignore[arg-type]
