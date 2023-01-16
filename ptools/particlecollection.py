from __future__ import annotations
from typing import Any, Iterable, Optional, TypeVar

from .namedarray import NamedArrayContainer
from . import spelling
from .atomattrs import AtomAttrs



class Particle:
    """Represents a single particle in a collection."""

    def __init__(self, collection, index):
        self._collection = collection
        self._index = index
        self._singular_to_plural = {
            prop.singular: prop.plural for prop in self._collection.atom_properties
        }

    def __eq__(self, __o: object) -> bool:
        """Compares two particles using their properties read from the ParticleCollection."""
        if isinstance(__o, Particle):
            self_properties = self._collection.atom_properties[self._index]
            other_properties = __o._collection.atom_properties[__o._index]
            return self_properties == other_properties

        rhs = ParticleCollection.from_objects([__o])
        self_properties = self._collection.atom_properties[self._index]
        other_properties = rhs.atom_properties[0]
        return self_properties == other_properties


    def __getattr__(self, name):
        if name in self._singular_to_plural:
            return self._collection.atom_properties.get(self._singular_to_plural[name])[
                self._index
            ]
        raise KeyError(f"No such property: {name!r}")

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

    def properties(self) -> list[str]:
        return [
            prop.singular for prop in self.self._collection.atom_properties.values()
        ]

    def __str__(self) -> str:
        attrs = ", ".join(f"{k}={getattr(self, k)!r}" for k in self._singular_to_plural.keys())
        return f"Particle({attrs})"



ParticleCollectionType = TypeVar("ParticleCollectionType", bound="ParticleCollection")


class ParticleCollection:
    """Represents a collection of particles."""

    atom_properties: NamedArrayContainer

    def __init__(self, atom_properties: Optional[NamedArrayContainer] = None):
        if atom_properties is None:
            atom_properties = NamedArrayContainer()

        if not isinstance(atom_properties, NamedArrayContainer):
            raise TypeError(
                f"Expected NamedArrayContainer, got {type(atom_properties)}"
            )

        self.atom_properties = atom_properties

    def size(self) -> int:
        if self.atom_properties:
            return self.atom_properties.number_of_elements()
        return 0

    def __len__(self) -> int:
        return self.size()

    def __eq__(self, __o: object) -> bool:
        if not isinstance(__o, ParticleCollection):
            return NotImplemented
        return self.atom_properties == __o.atom_properties

    @classmethod
    def from_objects(cls: type[ParticleCollectionType], objects: Iterable[Any]) -> ParticleCollectionType:
        """Creates a new collection from a list of objects.

        Reads the first object's attributes and stores all atom
        properties in a new ParticleCollection.
        """
        obj = cls()
        objects = list(objects) if objects is not None else []
        if objects:
            attrs = vars(objects[0]).keys()
            for name in attrs:
                plural = spelling.plural(name)
                values = [getattr(o, name) for o in objects]
                obj.atom_properties.add_array(name, plural, values)
        return obj

    def __getitem__(self, key: int | slice) -> Particle | ParticleCollection:
        """Returns a new collection with the selected atoms."""
        if isinstance(key, int):
            return Particle(self, key)
        return ParticleCollection(self.atom_properties[key])

    def __iter__(self):
        """Iterates over the atoms."""
        return iter(self.particles)

    def __contains__(self, particle: object) -> bool:
        """Returns whether the given particle is part of the collection."""
        return particle in self.particles

    def __add__(self, other: ParticleCollection) -> ParticleCollection:
        """Adds two collections together."""
        return ParticleCollection(self.atom_properties + other.atom_properties)

    @property
    def coordinates(self):
        """Returns the coordinates of the atoms."""
        return self.atom_properties.get("coordinates").values

    @property
    def particles(self):
        """Returns an iterator over the particles."""
        return list(Particle(self, i) for i in range(self.size()))

    def copy(self) -> ParticleCollection:
        """Returns a copy of the collection."""
        return ParticleCollection(self.atom_properties.copy())
