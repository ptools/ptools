from __future__ import annotations

import collections
import copy
from typing import Iterable, Optional

import numpy as np
from numpy.typing import ArrayLike


class AtomProperty:
    """Atom property.

    Stores a property for all atoms that belong to an AtomCollection.

    Attributes:
        singular (str): name in its singular form
        plural (str): name in its plural form (used in container types)
        values (ArrayLike): array of values for all atoms in a collection
    """

    def __init__(
        self,
        singular: str,
        plural: str,
        values: ArrayLike,
        *,
        array_comparison_func=None,
    ) -> None:
        self.singular: str = singular
        self.plural: str = plural
        self._values: np.ndarray = np.asarray(values)
        self._array_comparison_func = array_comparison_func

        if self._array_comparison_func is None:
            if np.issubdtype(self._values.dtype, np.floating):
                self._array_comparison_func = np.allclose
            else:
                self._array_comparison_func = np.array_equal

    @property
    def values(self) -> np.ndarray:
        return self._values

    @values.setter
    def values(self, values: ArrayLike):
        self._values = np.array(values)

    def size(self) -> int:
        """Returns the number of elements stored in the property."""
        return self.values.size

    def copy(self):
        return copy.deepcopy(self)

    def __repr__(self) -> str:
        attrs = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"{self.__class__.__qualname__}({attrs})"

    def __eq__(self, other: object):
        if not isinstance(other, self.__class__):
            raise ValueError(f"cannot compare {self.__class__} and {type(other)}")
        return (
            self.singular == other.singular
            and self.plural == other.plural
            and self._array_comparison_func(self.values, other.values)
        )

    def __getitem__(self, key: int | slice) -> float | AtomProperty:
        if isinstance(key, slice):
            return AtomProperty(
                self.singular,
                self.plural,
                self.values[key],
                array_comparison_func=self._array_comparison_func,
            )
        return self.values[key]


class AtomPropertyContainer(collections.abc.Container):
    """Container for AtomProperty instances.

    Attributes:
        properties (dict[str, AtomProperty]): maps property name (plural) with actual
            AtomProperty instance.
    """

    def __init__(self, properties: Optional[Iterable[AtomProperty]] = None) -> None:
        self._properties: dict[str, AtomProperty] = {}
        if properties:
            for prop in properties:
                self.register(prop)

    def __contains__(self, name_or_item: object) -> bool:
        """Returns whether a property is present in the collections."""
        if not isinstance(name_or_item, (str, AtomProperty)):
            raise ValueError(f"expects {str} or {AtomProperty}")
        plural = (
            name_or_item.plural
            if isinstance(name_or_item, AtomProperty)
            else name_or_item
        )
        return plural in self._properties

    def __getitem__(self, key: int | slice) -> AtomPropertyContainer:
        """Returns a new collection with a slice of all properties."""
        if isinstance(key, int):
            key = slice(key, key + 1)
        return AtomPropertyContainer(prop[key] for prop in self._properties.values())

    def get(self, plural: object) -> AtomProperty:
        """Returns the property with the given plural name."""
        if not isinstance(plural, str):
            raise ValueError("expects string to fetch property using their plural name")
        return self._properties[plural]

    def number_of_properties(self) -> int:
        return len(self._properties)

    def number_of_elements(self) -> int:
        if len(self._properties) == 0:
            return 0
        return len(next(iter(self._properties.values())).values)

    def add_property(self, singular: str, plural: str, values: ArrayLike):
        self.register(AtomProperty(singular, plural, np.asarray(values)))

    def register(self, item: AtomProperty):
        """Stores a new property in the collection."""
        if item.plural in self._properties:
            raise KeyError(f"property named {item.plural!r} already exists")
        if self.number_of_elements() != 0 and self.number_of_elements() != len(
            item.values
        ):
            err = (
                f"cannot add property {item.plural!r}: expected {self.number_of_elements()} elements, "
                f"got {len(item.values)}"
            )
            raise ValueError(err)
        self._properties[item.plural] = item.copy()
