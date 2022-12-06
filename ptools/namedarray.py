from __future__ import annotations

import collections
import copy
from typing import Iterable, Optional, Self

import numpy as np
from numpy.typing import ArrayLike


class NamedArray:
    """Holds a numpy array besides a singular and a plural name.

    Attributes:
        singular (str): name in its singular form
        plural (str): name in its plural form (used in container types)
        values (ArrayLike): array of values
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

    def __getitem__(self, key: int | slice) -> float | Self:
        if isinstance(key, slice):
            return self.__class__(
                self.singular,
                self.plural,
                self.values[key],
                array_comparison_func=self._array_comparison_func,
            )
        return self.values[key]


class NamedArrayContainer(collections.abc.Container):
    """Container for NamedArray instances.

    Provides a way to store multiple ``NamedArray`` instances with the same
    number of elements.

    The container is indexed by the number of elements.

    The ``get`` method allows to retrieve a property by its name.

    Attributes:
        array_list (dict[str, NamedArray]): maps array names (plural) with actual
            NamedArray instance.
    """

    def __init__(self, array_list: Optional[Iterable[NamedArray]] = None) -> None:
        self._properties: dict[str, NamedArray] = {}
        if array_list:
            for prop in array_list:
                self.register(prop)

    def __contains__(self, name_or_item: object) -> bool:
        """Returns whether a property is present in the collections."""
        if not isinstance(name_or_item, (str, NamedArray)):
            raise ValueError(f"expects {str} or {NamedArray}")
        plural = (
            name_or_item.plural
            if isinstance(name_or_item, NamedArray)
            else name_or_item
        )
        return plural in self._properties

    def __getitem__(self, key: int | slice) -> NamedArrayContainer:
        """Returns a new collection with a slice of all properties."""
        if isinstance(key, int):
            key = slice(key, key + 1)
        return NamedArrayContainer(prop[key] for prop in self._properties.values())

    def get(self, plural: object) -> NamedArray:
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
        self.register(NamedArray(singular, plural, np.asarray(values)))

    def register(self, item: NamedArray):
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
