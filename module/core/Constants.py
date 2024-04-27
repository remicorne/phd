from module.core.JSON import JSONMapping
from module.constants import ROOT
from dataclasses import dataclass
import difflib
from module.core.Question import Question


@dataclass
class ConstantRegistry(JSONMapping):

    def detect(self, invalid_name):
        lazy_dict = {k.upper(): k for k in self.list}
        lazy_guess = difflib.get_close_matches(
            invalid_name.upper(), lazy_dict.keys(), n=1, cutoff=0.6
        )
        return lazy_guess[0] if lazy_guess else None

    def get_valid_choice(self, invalid_choice):
        lazy_guess = self.detect(invalid_choice)
        if lazy_guess:
            is_correct = Question.yes_or_no(
                f"INVALID: {invalid_choice}. DETECTED {lazy_guess}: {self[lazy_guess]}. CONFIRM?"
            )
            if is_correct:
                return lazy_guess
        try:
            new_choice = Question.select_one(
                f"UNKNOWN: {invalid_choice}, SELECT FROM:", self.dict.items()
            )
            while new_choice[0] not in self:
                new_choice = Question.select_one(
                    f"UNKNOWN {invalid_choice}, SELECT FROM:", self.dict.items()
                )
            return new_choice[0]
        except SystemExit:
            print("EDIT REGISTRY AND RETRY")
            exit(1)


constant_registry_location = f"{ROOT}/module/json"

# Rewrite the following code but inverting the two parameters of the function select_one
REGIONS = ConstantRegistry(constant_registry_location, "regions")
COMPOUNDS = ConstantRegistry(constant_registry_location, "compounds")
COMPOUND_CLASSES = ConstantRegistry(constant_registry_location, "compound_classes")
MACRO_REGIONS = ConstantRegistry(constant_registry_location, "macro_regions")
CIRCUITS = ConstantRegistry(constant_registry_location, "circuits")
