from dataclasses import dataclass, field
from typing import Any
from module.core.Constants import COMPOUNDS_AND_REGIONS, COMPOUNDS_AND_REGIONS_CLASSES

TYPES = ["compound", "region"]

@dataclass
class StandardizedParameter:
    """
    Represents a standardized parameter that can be either a compound or a region.

    This class is designed to behave like a string, delegating most string operations
    to its `name` attribute. It overrides the equality operator to allow for custom
    comparison logic, including aliases or related classes.

    Attributes:
        name (str): The name of the parameter.
        type (str): The type of the parameter; should be either 'compound' or 'region'.
        classes (list): A list of related classes or aliases associated with the parameter.
        full_name (str): The full name corresponding to the abbreviation.

    Methods:
        includes(type: str) -> bool:
            Static method to check if a given type is valid (i.e., 'compound' or 'region').

        initialize(name: str, type: str) -> Union['StandardizedParameter', str]:
            Class method to create an instance if the name is valid for the given type.
            Returns the name unchanged if it's not valid.

    Behavior:
        - The class is immutable and hashable as it is necessary for serialization + df operations.
        - It overrides equality (`__eq__`) to compare based on the `name`, `classes`, and `full_name` attributes.
        - Supports all standard string operations by delegating to the `name` attribute.

    Example:
        >>> param = StandardizedParameter(name='OF', type='region')
        >>> print(param.upper())
        'OF'
        >>> param == 'cortex'
        >>> param == 'Orbitofrontal cortex'

    """

    name: str
    type: str
    classes: list = field(init=False)
    full_name: str = field(init=False)
    def __post_init__(self):
        classes = COMPOUNDS_AND_REGIONS_CLASSES[self.type]
        self.full_name = COMPOUNDS_AND_REGIONS[self.type][self.name]
        self.classes = classes.get_item_classes(self.name)

    @staticmethod
    def includes(type):
        return type in TYPES
    
    @classmethod
    def initialize(cls, name, type):
        valid_items = COMPOUNDS_AND_REGIONS[type]
        return cls(name=name, type=type) if name in valid_items else name
            
    def __eq__(self, other):
        other_str = str(other)
        return self.name == other_str or other_str in self.classes or other_str == self.full_name

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def __lt__(self, other):
        return self.name < str(other)

    def __le__(self, other):
        return self.name <= str(other)

    def __gt__(self, other):
        return self.name > str(other)

    def __ge__(self, other):
        return self.name >= str(other)

    def __add__(self, other):
        return self.name + str(other)

    def __radd__(self, other):
        return str(other) + self.name

    def __contains__(self, item):
        return item in self.name
 