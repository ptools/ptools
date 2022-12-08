from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Iterable, Optional, Self

from .namedarray import NamedArrayContainer
from . import spelling


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
        if not isinstance(__o, Particle):
            return NotImplemented
        self_properties = self._collection.atom_properties[self._index]
        other_properties = __o._collection.atom_properties[__o._index]
        return self_properties == other_properties

    def __getattr__(self, name):
        if name in self._singular_to_plural:
            return self._collection.atom_properties.get(self._singular_to_plural[name])[self._index]
        raise KeyError(f"No such property: {name!r}")

    def properties(self) -> list[str]:
        return [prop.singular for prop in self.self._collection.atom_properties.values()]


@dataclass
class ParticleCollection:
    """Represents a collection of particles."""

    atom_properties: NamedArrayContainer = field(default_factory=NamedArrayContainer)
    _particles: list[Particle] = field(init=False, repr=False, default_factory=list)

    def __post_init__(self):
        if not isinstance(self.atom_properties, NamedArrayContainer):
            raise TypeError(f"Expected NamedArrayContainer, got {type(self.atom_properties)}")
        if self.atom_properties is None:
            self.atom_properties = NamedArrayContainer()
        self._particles = [Particle(self, i) for i in range(self.size())]

    def size(self) -> int:
        if self.atom_properties:
            return self.atom_properties.number_of_elements()
        return 0

    def __len__(self) -> int:
        return self.size()

    @classmethod
    def from_attributes(cls, attrs: NamedArrayContainer) -> Self:
        """Creates a new collection from a list of atoms."""
        obj = cls()
        obj.atom_properties = attrs
        return obj

    @classmethod
    def from_objects(cls, objects: Iterable[Any]) -> Self:
        """Creates a new collection from a list of objects.

        Reads the first object's attributes and stores all atom
        properties in a new ParticleCollection.
        """
        obj = cls()
        objects = list(objects) if objects is not None else []
        if objects:
            attrs = list(vars(objects[0]))
            for name in attrs:
                plural = spelling.plural(name)
                values = [getattr(o, name) for o in objects]
                obj.atom_properties.add_array(name, plural, values)
        return obj

    def __getitem__(self, key: int | slice) -> ParticleCollection:
        """Returns a new collection with the selected atoms."""
        return ParticleCollection(self.atom_properties[key])

    def __iter__(self):
        """Iterates over the atoms."""

        return iter(self._particles)

    def __contains__(self, particle: object) -> bool:
        """Returns whether the given particle is part of the collection."""
        if isinstance(particle, Particle):
            return particle in self._particles
        return False

    @property
    def coordinates(self):
        """Returns the coordinates of the atoms."""
        return self.atom_properties.get("coordinates").values


