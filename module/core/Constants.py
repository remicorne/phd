from module.core.JSON import JSONMapping
from module.core.FileSystem import FileSystem
import os
from dataclasses import dataclass
import difflib
from module.core.questions import yes_or_no, select_one

@dataclass
class ConstantRegistry(JSONMapping):

    def generate(self):
        raise NotImplementedError('Contant registries should be handled directly in the JSON')

    def detect(self, invalid_name):
        lazy_dict = {k.upper(): k for k in self.list}
        lazy_guess = difflib.get_close_matches(
            invalid_name.upper(), lazy_dict.keys(), n=1, cutoff=0.6
        )
        return lazy_guess[0] if lazy_guess else None

    def get_valid_choice(self, invalid_choice):
        lazy_guess = self.detect(invalid_choice)
        if lazy_guess:
            is_correct = yes_or_no(
                f"INVALID: {invalid_choice}. DETECTED {lazy_guess}: {self[lazy_guess]}. CONFIRM?"
            )
            if is_correct:  
                return lazy_guess
        try:
            new_choice = select_one(
                f"UNKNOWN: {invalid_choice}, SELECT FROM:", self.dict
            )
            while new_choice[0] not in self:
                new_choice = select_one(
                    f"UNKNOWN {invalid_choice}, SELECT FROM:", self.dict
                )
            return new_choice[0]
        except SystemExit:
            print("EDIT REGISTRY AND RETRY")
            exit(1)


# Rewrite the following code but inverting the two parameters of the function select_one
REGIONS = ConstantRegistry(filepath=os.path.join(FileSystem.CONSTANTS, "regions"))
COMPOUNDS = ConstantRegistry(filepath=os.path.join(FileSystem.CONSTANTS, "compounds"))
COMPOUND_CLASSES = ConstantRegistry(filepath=os.path.join(FileSystem.CONSTANTS, "compound_classes"))
MACRO_REGIONS = ConstantRegistry(filepath=os.path.join(FileSystem.CONSTANTS, "macro_regions"))
CIRCUITS = ConstantRegistry(filepath=os.path.join(FileSystem.CONSTANTS, "circuits"))
COMPOUNDS_AND_REGIONS = {
    "region": REGIONS.list,
    "compund": COMPOUNDS.list
}