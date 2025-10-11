from typing import Any, Type
from module.core.JSON import JSONMapping
from module.core.FileSystem import FileSystem
import os
import sys
from dataclasses import dataclass, field
import difflib
from module.core.questions import yes_or_no, select_one
from typing import ClassVar


def singular(name):
    return name[:-1] if name[-1] == "s" else name


class Registry(JSONMapping):
    """A JSON mapping for constants.
    Used for validation of predefined constants (regions, compounds, compound classes..)

    Args:
        filepath (str): The path to the JSON file.

    Returns:
        THe content of the JSON file as a dict.
    """

    _instances: ClassVar[dict] = {}

    from_scratch: bool = field(init=False, default=False)
    filepath: str

    def __new__(cls, filepath: str, *args, **kwargs):
        # Check if an instance for the given filepath already exists
        if filepath not in cls._instances:
            # Create and store the new instance
            instance = super().__new__(cls)
            cls._instances[filepath] = instance
        return cls._instances[filepath]

    def generate(self):
        raise FileNotFoundError(f"No constant registry at {self.filepath}")

    def detect(self, key):
        """Uses difflib to detect the closest match in the list.

        Args:
            key (hasable): The key to detect.
            Used by get_valid_choice to detect the closest match.
            Used for incorrect constants (region name with mistake...)

        Returns:
            str: the closest match
        """
        case_converter = {k.upper(): k for k in self.list}
        lazy_guess_upper = difflib.get_close_matches(
            key.upper(), case_converter.keys(), n=1, cutoff=0.6
        )
        return case_converter[lazy_guess_upper[0]] if lazy_guess_upper else None

    def choose_valid_value(self, value):
        """
        Tries to get a valid choice from the user.
        Uses difflib to detect the closest match in the list.
        If no match is found, it asks the user to select one from the list.

        Args:
            invalid_choice (_type_): Key that was invalid

        Returns:
            str: A valid choice from keys
        """
        if value in self:
            return value
        lazy_guess = self.detect(value)
        if lazy_guess:
            is_correct = yes_or_no(
                f"INVALID: {value}. DETECTED {lazy_guess}: {self[lazy_guess]}. CONFIRM?"
            )
            if is_correct:
                return lazy_guess
        try:
            new_choice = select_one(f"INVALID: {value}, SELECT FROM:", self.dict)
            while new_choice not in self:
                new_choice = select_one(
                    f"UNKNOWN CHOICE: {new_choice}, SELECT FROM:", self.dict
                )
            return new_choice
        except SystemExit:
            print("EDIT REGISTRY AND RETRY")
            sys.exit(1)

    def order(self, iterable):
        ordered = [item for item in self if item in iterable]
        rest = [item for item in iterable if item not in ordered]
        return ordered + rest

    def get_order(self, item):
        return self.list.index(item) if item in self else None

    @staticmethod
    def get_filename_from_element_type(element_type):
        return f"{element_type + 's'}"

    @classmethod
    def get_filename(cls, element_type=None, name=None) -> str:
        """Get the filename for the given element_type or name."""
        if element_type:
            base_name = cls.get_filename_from_element_type(element_type)
        elif name:
            base_name = name
        else:
            raise ValueError("Either element_type or name must be provided.")
        return f"{base_name}.{cls.extension}"

    @classmethod
    def get_filepath(cls, element_type=None, name=None) -> str:
        filename = cls.get_filename(element_type=element_type, name=name)
        return os.path.join(FileSystem.CONSTANTS, filename)

    @classmethod
    def get_registry(cls, element_type=None, name=None) -> Type["Registry"]:
        """Public method to get the loaded registry."""
        filepath = cls.get_filepath(element_type=element_type, name=name)
        return cls(filepath=filepath)

    @classmethod
    def exists(cls, element_type=None, name=None) -> bool:
        """Check if the registry file exists."""
        filepath = cls.get_filepath(element_type=element_type, name=name)
        return os.path.exists(filepath)

    @classmethod
    def list_registries(cls):
        return {
            singular(os.path.splitext(filename)[0]): cls(
                filepath=os.path.join(FileSystem.CONSTANTS, filename)
            )  # type: ignorefilename)
            for filename in os.listdir(FileSystem.CONSTANTS)
            if os.path.splitext(filename)[1] == f".{cls.extension}"
        }


class ClassRegistry(Registry):
    def get_item_classes(self, item):
        return [klass for klass, elements in self.items() if item in elements]

    @staticmethod
    def get_filename_from_element_type(element_type):
        return f"{element_type}_classes"


REGISTRIES = Registry.list_registries()


def string_to_numerical(value):
    return int("".join([str(ord(c)) for c in value]))


@dataclass(frozen=True, eq=False)
class Characteristic:
    type: str
    value: str
    position: int = field(init=False)
    classes: Any = field(init=False, default=None)

    def __post_init__(self):
        # Immutable assignment for calculated attributes
        object.__setattr__(
            self,
            "position",
            (
                REGISTRIES[self.type].get_order(self.value)
                if self.type in REGISTRIES
                else string_to_numerical(self.value)
            ),
        )
        if ClassRegistry.exists(element_type=self.type):
            object.__setattr__(
                self,
                "classes",
                ClassRegistry.get_registry(element_type=self.type).get_item_classes(
                    self.value
                ),
            )

    def __str__(self):
        return f"{self.type}={self.value}"

    def __lt__(self, other):
        if isinstance(other, Characteristic):
            return self.position < other.position
        raise NotImplementedError("Cannot compare with non-Characteristic")

    def __eq__(self, other):
        if isinstance(other, Characteristic):
            return self.type == other.type and self.value == other.value
        elif isinstance(other, dict):
            return all(self.dict.get(k) == v for k, v in other.items())
        return False

    def __hash__(self):
        return hash((self.type, self.value))

    @property
    def dict(self):
        return {self.type: self.value}
