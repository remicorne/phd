from dataclasses import dataclass
from module.constants import ROOT
import os

PROJECTS = f"{ROOT}/PROJECTS"

@dataclass
class Cacheable:

    location: str

    def __post_init__(self):
        if not self.is_saved:
            self.initialize()
          
    def generate(self):
        raise NotImplementedError('This method should be implemented for all custom Cacheables')
    
    def initialize(self):
        data = self.generate()
        self.save(data)
            
    def load(self):
        raise NotImplementedError('This method should be implemented for all custom Cacheables')
            
    def save(self, data):
        pass

    @property
    def filepath(self):
        raise NotImplementedError
    
    @property
    def is_saved(self):
        return os.path.isfile(self.filepath)
