from __future__ import annotations

import collections
import copy
from typing import Any, Iterable, Iterator, Optional, Sequence

import numpy as np
from numpy.typing import ArrayLike

from . import spelling
from .array3d import array3d


class NamedArray:
    """Holds a numpy array besides a singular and a plural name.

    Attributes:
        singular (str): name in its singular form
        plural (str): name in its plural form (used in container types)
        values (ArrayLike): array of values
    """

    _values: np.ndarray

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

        if not isinstance(values, (collections.abc.Sequence, np.ndarray)):
            raise ValueError(
                f"values must be a sequence or an array, not {type(values).__qualname__}"
            )

        if len(values) == 0:
            self._values = np.empty((0,), dtype=float)
        else:
            if isinstance(values[0], str):
                self._values = np.asarray(values, dtype="O")
            elif isinstance(values[0], array3d):
                self._values = array3d(values)
            else:
                self._values = np.asarray(values)

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

    def empty(self) -> bool:
        """Returns whether the property is empty."""
        return self.size() == 0

    def tolist(self) -> list:
        """Returns the values as a list."""
        return self.values.tolist()

    def copy(self):
        return copy.deepcopy(self)

    def __repr__(self) -> str:
        attrs = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"{self.__class__.__qualname__}({attrs})"

    def __eq__(self, other: object):
        if isinstance(other, self.__class__):
            return (
                self.singular == other.singular
                and self.plural == other.plural
                and self._array_comparison_func(self.values, other.values)
            )
        if isinstance(other, collections.abc.Iterable):
            other_values = np.asarray(other)
            if other_values.shape != self.values.shape:
                return False
            return self._array_comparison_func(self.values, other_values)
        if not isinstance(other, self.__class__):
            raise ValueError(f"cannot compare {self.__class__} and {type(other)}")
        raise NotImplementedError

    def __iter__(self) -> Any:
        return iter(self.values)

    def __getitem__(self, key: int | slice | Sequence[int]) -> Any | NamedArray:
        if isinstance(key, (int, np.integer)):
            return self.values[key]
        return self.__class__(
            self.singular,
            self.plural,
            self.values[key],
            array_comparison_func=self._array_comparison_func,
        )

    def __setitem__(self, key: int | slice, value: ArrayLike):
        self.values[key] = value


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

    @classmethod
    def from_objects(cls, objects: Iterable[object]) -> NamedArrayContainer:
        """Creates a new container from a list of objects."""
        new_container = cls()
        obj = next(iter(objects))
        attrs = vars(obj).keys()
        for name in attrs:
            plural = spelling.pluralize(name)
            values = [getattr(o, name) for o in objects]
            new_container.add_array(name, plural, values)
        return new_container

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

    def __iter__(self) -> Iterator[NamedArray]:
        return iter(self._properties.values())

    def __getitem__(self, key: str | int | slice) -> NamedArrayContainer:
        if isinstance(key, str):
            return self._get_by_name(key)
        if isinstance(key, int):
            key = slice(key, key + 1)
        # mypy-ignore: ``key`` is a slice, ``prop[key]`` is NamedArray
        return NamedArrayContainer(prop[key] for prop in self._properties.values())  # type: ignore

    def __setitem__(self, key: str, value: NamedArrayContainer):
        if isinstance(key, str):
            self._set_by_name(key, value)
            return
        raise KeyError(f"cannot set item with key {key!r}")

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, self.__class__):
            raise TypeError(
                f"cannot compare {self.__class__.__qualname__} and {type(other)}"
            )
        if not self._properties.keys() == other._properties.keys():
            return False
        for lhs, rhs in zip(self.iter_arrays(), other.iter_arrays()):
            if not lhs == rhs:
                return False
        return True

    def __add__(self, other: object) -> NamedArrayContainer:
        if not isinstance(other, self.__class__):
            raise TypeError(
                f"cannot add {self.__class__.__qualname__} and {type(other)}"
            )

        if self.empty():
            return other.copy()

        if other.empty():
            return self.copy()

        # If not empty, check that both collections have the same properties
        if not self._properties.keys() == other._properties.keys():
            raise ValueError("cannot add two collections with different properties")

        return NamedArrayContainer(
            NamedArray(
                lhs.singular,
                lhs.plural,
                np.concatenate((lhs.values, rhs.values)),
                array_comparison_func=lhs._array_comparison_func,
            )
            for lhs, rhs in zip(self.iter_arrays(), other.iter_arrays())
        )

    def empty(self) -> bool:
        """Returns whether the collection is empty."""
        return len(self._properties) == 0

    def names(self) -> list[str]:
        """Returns the plural names of the properties."""
        return list(self._properties.keys())

    def singular_names(self) -> list[str]:
        """Returns the singular names of the properties."""
        return [prop.singular for prop in self._properties.values()]

    def plural_names(self) -> list[str]:
        """Returns the plural names of the properties.

        This is equivalent to the names() method.
        """
        return [prop.plural for prop in self._properties.values()]

    def iter_arrays(self) -> Iterator[NamedArray]:
        return iter(self._properties.values())

    def _get_by_name(self, plural: object) -> NamedArray:
        """Returns the property with the given plural name."""
        if not isinstance(plural, str):
            raise ValueError("expects string to fetch property using its plural name")
        return self._properties[plural]

    def _set_by_name(self, plural: object, value: ArrayLike):
        """Sets the property with the given plural name."""
        if not isinstance(plural, str):
            raise KeyError("expects string to set property using its plural name")
        if plural not in self._properties:
            raise KeyError(f"property {plural!r} not registered")
        value = np.asarray(value)
        if value.shape != () and value.shape != self._properties[plural].values.shape:
            raise ValueError(
                f"cannot set property {plural!r} with array of shape {value.shape} "
                f"(expected {self._properties[plural].values.shape})"
            )
        if value.shape == ():
            value = np.full(self._properties[plural].values.shape, value)
        self._properties[plural].values = value

    def get_at(self, plural: str, index: int) -> Any:
        """Returns the value of the property with the given plural name at the given index."""
        return self._properties[plural].values[index]

    def set_at(self, plural: str, index: int, value: Any):
        """Sets the value of the property with the given plural name at the given index."""
        self._properties[plural].values[index] = value

    def number_of_properties(self) -> int:
        return len(self._properties)

    def number_of_elements(self) -> int:
        if len(self._properties) == 0:
            return 0
        return len(next(iter(self._properties.values())).values)

    def add_array(self, singular: str, plural: str, values: Sequence | np.ndarray):
        if len(values) > 0 and isinstance(values[0], array3d):
            values = array3d(values)
        else:
            values = np.asarray(values)
        self.register(NamedArray(singular, plural, values))

    def register(self, item: NamedArray):
        """Stores a new property in the collection."""
        if item.plural in self._properties:
            raise KeyError(f"property named {item.plural!r} already exists")
        if self.number_of_properties() != 0 and self.number_of_elements() != len(
            item.values
        ):
            err = (
                f"cannot add property {item.plural!r}: expected {self.number_of_elements()} "
                f"elements, got {len(item.values)}"
            )
            raise ValueError(err)
        self._properties[item.plural] = item.copy()

    def remove_array(self, plural: str):
        """Removes a property from the collection."""
        del self._properties[plural]

    def copy(self) -> NamedArrayContainer:
        return self.__class__(self._properties.values())
