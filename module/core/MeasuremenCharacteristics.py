from typing import Any
import pandas as pd
from module.core.Constants import Characteristic
from dataclasses import dataclass


@dataclass(frozen=True, eq=False)
class MeasurementCharacteristics:
    """
    A flexible and immutable measurement class that:
    - Stores arbitrary characteristics as key-value pairs (like a dictionary).
    - Can be used in pandas DataFrame cells and still allow sorting and grouping.
    - Implements tuple-like ordering, hashing, and equality by leveraging characteristic values.
    - Allows convenient selection of attributes via indexing (e.g., measurement['region']) and .get() method.

    Representation:
    - If single attribute: str(measurement) -> the single attribute value.
    - If multiple attributes: str(measurement) -> a '|' joined string of characteristic string representations.
    """

    def __init__(self, characteristics):
        # Convert pd.Series to a dict if needed
        if isinstance(characteristics, pd.Series):
            characteristics = characteristics.to_dict()
        elif not isinstance(characteristics, dict):
            raise ValueError(
                f"Characteristics must be a dict or pd.Series, got {type(characteristics)}"
            )
        object.__setattr__(
            self,
            "_characteristics",
            frozenset(Characteristic(k, v) for k, v in characteristics.items()),
        )
        object.__setattr__(
            self,
            "_characteristic_order",
            tuple(characteristics.keys()),
        )  # TODO dirty hack?

    def __getitem__(self, key):
        """Get a characteristic by key, raises KeyError if not found."""
        for characteristic in self._characteristics:
            if characteristic.type == key:
                return characteristic.value
        raise KeyError(f"Characteristic '{key}' does not exist")

    def keys(self):
        return self._characteristic_order

    def get(self, key, default=None):
        """Get a characteristic with a default."""
        try:
            return self[key]
        except KeyError:
            return default

    def __str__(self):
        """
        String representation based on characteristics:
        """
        return ", ".join(str(c) for c in self._characteristics)

    def __repr__(self):
        """
        String representation based on characteristics:
        """
        return f"MeasurementCharacteristincs({', '.join(str(c) for c in self._characteristics)})"

    def __eq__(self, other):
        """Equality check."""
        if isinstance(other, MeasurementCharacteristics):
            return self._characteristics == other._characteristics
        elif isinstance(other, dict):
            return all(self[k] == other[k] for k in other)
        return False

    def __lt__(self, other):
        """Less-than comparison for sorting."""
        if isinstance(other, MeasurementCharacteristics):
            if self._characteristic_order != other._characteristic_order:
                return self._characteristic_order < other._characteristic_order
            return tuple(self._characteristic_order) < tuple(
                other[k] for k in self._characteristic_order
            )
        return NotImplementedError

    def __hash__(self):
        """Make the object hashable, enabling groupby and set usage."""
        return hash(self._characteristics)

    # def to_dict(self):
    #     """Convert the characteristics back to a dictionary."""
    #     return {c.type: c.value for c in self._characteristics}
