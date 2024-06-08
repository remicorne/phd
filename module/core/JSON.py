
import json
from dataclasses import dataclass
from module.core.Cacheable import Cacheable
from module.core.questions import input_escape
from typing import ClassVar

@dataclass
class JSONMapping(Cacheable):
    
    extension: ClassVar[str] = "json"
    
    def load(self):
        with open(self.filepath) as outfile:
            mapping = json.load(outfile)
        return mapping

    def save(self, mapping):
        with open(self.filepath, "w") as json_file:
            json.dump(mapping, json_file)

    def add(self, key, value):
        mapping = self.load()
        mapping[key] = value
        self.save(mapping)
    
    def __contains__(self, key):
        return key in self.list
    
    def __getitem__(self, key):
        return self.dict.get(key)
    
    def __setitem__(self, key, value):
        self.add(key, value)    
        
    def __repr__(self) -> str:
        return "\n".join(f"{k}: {v}" for k, v in self.dict.items())
    
    @property
    def dict(self):
        return self.load()        
    
    @property
    def list(self):
        return self.dict.keys()
