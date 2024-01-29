
import json
from typing import Any
from module.constants import ROOT
from dataclasses import dataclass
import difflib
from module.core.Questions import Questions

@dataclass
class JSONMapping:
    name: str
    _location = f"{ROOT}/module/json"
    
    def __post_init__(self):
        self.path = f"{self._location}/{self.name}.json"
    
    def load(self):
        with open(self.path) as outfile:
            mapping = json.load(outfile)
        return mapping

    def save(self, mapping):
        with open(self.path, "w") as json_file:
            json.dump(mapping, json_file)

    def add(self, key, value):
        mapping = self.load()
        mapping[key] = value
        self.save(mapping)
    
    @property
    def dict(self):
        return self.load()        
    
    @property
    def list(self):
        return self.dict.keys()
    
    def detect(self, invalid_name):
        lazy_dict = {k.upper(): k for k in self.list}
        lazy_guess = difflib.get_close_matches(invalid_name.upper(), lazy_dict.keys(), n=1, cutoff=0.6)
        return lazy_guess[0] if lazy_guess else None
    
    def __contains__(self, key):
        return key in self.list
    
    def __getitem__(self, key):
        return self.dict.get(key)
    
    def __setitem__(self, key, value):
        self.add(key, value)
        
    def get_valid_choice(self, invalid_choice):
        lazy_guess = self.detect(invalid_choice)
        if lazy_guess:
            is_correct = Questions.yes_or_no(
                f"INVALID: {invalid_choice}. DETECTED {lazy_guess}: {self[lazy_guess]}. CONFIRM?"
            )
            if is_correct:
                return lazy_guess
        try:
            new_choice = Questions.select_one(
                f"UNKNOWN: {invalid_choice}, SELECT FROM:", self.dict
            )
            while new_choice not in self:
                new_choice = Questions.select_one(
                    f"UNKNOWN {invalid_choice}, SELECT FROM:", self.dict
                )
            return new_choice
        except SystemExit:
            print("EDIT REGISTRY AND RETRY")
            exit(1)

            

REGIONS = JSONMapping('regions')
COMPOUNDS = JSONMapping('compounds')
COMPOUND_CLASSES = JSONMapping('compound_classes')
MACRO_REGIONS = JSONMapping('macro_regions')
CIRCUITS = JSONMapping('circuits')