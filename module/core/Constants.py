from module.core.JSON import JSONMapping
from module.core.FileSystem import FileSystem
import os
from dataclasses import dataclass
import difflib
from module.core.questions import yes_or_no, select_one
import json


@dataclass
class ConstantRegistry(JSONMapping):
    """A JSON mapping for constants.
    Used for validation of predefined constants (regions, compounds, compound classes..)

    Args:
        filepath (str): The path to the JSON file.

    Returns:
        THe content of the JSON file as a dict.
    """

    filepath: str

    def generate(self):
        raise NotImplementedError(
            "Contant registries should be handled directly in the JSON"
        )

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

    def get_valid_choice(self, invalid_choice):
        """
        Tries to get a valid choice from the user.
        Uses difflib to detect the closest match in the list.
        If no match is found, it asks the user to select one from the list.

        Args:
            invalid_choice (_type_): Key that was invalid

        Returns:
            str: A valid choice from keys
        """
        lazy_guess = self.detect(invalid_choice)
        if lazy_guess:
            is_correct = yes_or_no(
                f"INVALID: {invalid_choice}. DETECTED {lazy_guess}: {self[lazy_guess]}. CONFIRM?"
            )
            if is_correct:
                return lazy_guess
        try:
            new_choice = select_one(
                f"INVALID: {invalid_choice}, SELECT FROM:", self.dict
            )
            while new_choice not in self:
                new_choice = select_one(
                    f"UNKNOWN CHOICE: {new_choice}, SELECT FROM:", self.dict
                )
            return new_choice
        except SystemExit:
            print("EDIT REGISTRY AND RETRY")
            exit(1)

    def order(self, iterable):
        ordered = [item for item in self if item in iterable]
        rest = [item for item in iterable if item not in ordered]
        return ordered + rest
    
class ClassRegistry(ConstantRegistry):
    
    def get_item_classes(self, item):
        return [klass for klass, constituents in self.items() if item in constituents]


REGIONS = ConstantRegistry(filepath=os.path.join(FileSystem.CONSTANTS, "regions"))
COMPOUNDS = ConstantRegistry(filepath=os.path.join(FileSystem.CONSTANTS, "compounds"))
COMPOUND_CLASSES = ClassRegistry(
    filepath=os.path.join(FileSystem.CONSTANTS, "compound_classes")
)
REGION_CLASSES = ClassRegistry(
    filepath=os.path.join(FileSystem.CONSTANTS, "region_classes")
)
REGION_CLASSES_POSITIONS = ConstantRegistry(
    filepath=os.path.join(FileSystem.CONSTANTS, "region_classes_positions")
)
CIRCUITS = ConstantRegistry(filepath=os.path.join(FileSystem.CONSTANTS, "circuits"))
COMPOUNDS_AND_REGIONS = {"region": REGIONS, "compound": COMPOUNDS}
COMPOUNDS_AND_REGIONS_CLASSES = {"region": REGION_CLASSES, "compound": COMPOUND_CLASSES}

