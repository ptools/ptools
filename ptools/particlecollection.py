from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Iterable, Optional, Self

from .namedarray import NamedArrayContainer

@dataclass
class ParticleCollection:

    atom_properties: Optional[NamedArrayContainer] = None

    def __post_init__(self):
        if self.atom_properties is None:
            self.atom_properties = NamedArrayContainer()

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
            attrs = list(objects[0].__dict__.keys())
            for name in attrs:
                plural = singular_to_plural(name)
                values = [getattr(o, name) for o in objects]
                obj.atom_properties.add_array(name, plural, values)
        return obj

    def __getitem__(self, key: int | slice) -> ParticleCollection:
        """Returns a new collection with the selected atoms."""
        return ParticleCollection(self.atom_properties[key])

    def __iter__(self):
        """Iterates over the atoms."""
        class Atom:
            def __init__(self, collection, index):
                self._collection = collection
                self._index = index

            def __getattr__(self, name):
                return self._collection.atom_properties.get(name)[self._index]

        return iter(Atom(self, i) for i in range(self.size()))


def singular_to_plural(word: str) -> str:
    """Given a word, returns its plural form."""
    def format_output(word: str) -> str:
        output = word
        if tokens:
            output = "_".join(tokens) + "_" + word
        return output

    irregular = {
        "child": "children",
        "goose": "geese",
        "man": "men",
        "woman": "women",
        "tooth": "teeth",
        "foot": "feet",
        "mouse": "mice",
        "person": "people",
        "index": "indices",  # "indexes" actually is correct as well
    }

    tokens = []

    if "_" in word:
        tokens = word.split("_")
        word = tokens.pop(-1)  # make plural form only on last word

    if len(word) == 1:
        return word + "s"

    if word in irregular:
        return format_output(irregular[word])

    if word in ("sheep", "series", "species", "data", "deer", "coordinates"):
        return format_output(word)

    if word.endswith("us"):
        return format_output(word[:-2] + "i")

    if word.endswith("is"):
        return format_output(word[:-2] + "es")

    if word.endswith("on"):
        return format_output(word[:-2] + "a")

    for suffix in ("s", "ss", "sh", "ch", "x", "z"):
        if word.endswith(suffix):
            return format_output(word + "es")

    for suffix in ("f", "fe"):
        exceptions: tuple[str, ...] = ("roof", "belief", "chef", "chief")
        if word.endswith(suffix) and word not in exceptions:
            return format_output(word[: -len(suffix)] + "ves")

    if word.endswith("y") and word[-2] not in "aeiouy":
        return format_output(word[:-1] + "ies")

    if word.endswith("o"):
        exceptions = ("photo", "piano", "halo")
        if word not in exceptions:
            return format_output(word + "es")

    if word.endswith("ex"):
        return format_output(word[:-2] + "ices")

    word = word + "s"

    return format_output(word)
