
import json
from dataclasses import dataclass
from module.core.Cacheable import Cacheable

@dataclass
class JSONMapping(Cacheable):
    location = str
    name: str
    
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
    
    @property
    def dict(self):
        return self.load()        
    
    @property
    def list(self):
        return self.dict.keys()
    
    def __contains__(self, key):
        return key in self.list
    
    def __getitem__(self, key):
        return self.dict.get(key)
    
    def __setitem__(self, key, value):
        self.add(key, value)    
        
    @property
    def filepath(self):
        return f"{self.location}/{self.name}.json"
