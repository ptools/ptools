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
        self._values = np.asarray(values)

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


class AtomPropertyCollection(collections.abc.Mapping):
    """Container for AtomProperty instances.

    Allows to access properties either by singular or plural name.

    Attributes:
        properties (dict[str, AtomProperty]): maps property name (plural) with actual AtomProperty instance.
    """

    def __init__(self, properties: Optional[Iterable[AtomProperty]] = None) -> None:
        self._properties = {}
        if properties:
            self._properties = {prop.plural: prop for prop in properties}

    def __len__(self) -> int:
        """"Returns the numbers of properties in the collection."""
        return len(self._properties)

    def __iter__(self):
        return iter(self._properties)

    def __getitem__(self, plural: object) -> AtomProperty:
        if not isinstance(plural, str):
            raise ValueError("expects string to fetch property using their plural name")
        return self._properties[plural]

    def __contains__(self, name_or_item: object) -> bool:
        """Returns whether a property is present in the collections."""
        if not isinstance(name_or_item, (str, AtomProperty)):
            raise ValueError(f"expects {str} or {AtomProperty}")
        plural = name_or_item.plural if isinstance(name_or_item, AtomProperty) else name_or_item
        return plural in self._properties

    def add_property(self, singular: str, plural: str, values: ArrayLike):
        self.register(AtomProperty(singular, plural, np.asarray(values)))

    def register(self, item: AtomProperty):
        """Stores a new property in the collection."""
        if item.plural in self._properties:
            raise KeyError(f"property named {item.plural!r} already exists")
        self._properties[item.plural] = item.copy()
