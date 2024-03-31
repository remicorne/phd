import os
from dataclasses import dataclass
from module.core.Project import Project
from typing import ClassVar

ROOT = os.getcwd()  # This gives terminal location (terminal working dir)

@dataclass
class ProjectMember:

    project: str
    _subfolder: ClassVar[str] = None

    def __post_init__(self):
        self.project = Project(self.project)
        self.location = f"{self.project.location}/{self._subfolder}"
            
    def generate(self):
        pass

    def initialize(self):
        data = self.generate()
        self.save(data)
            
    def load(self):
        pass
            
    def save(self, data):
        pass
    
    def get(self):
        if not self.is_saved:
            self.initialize()
        return self.load()

    @property
    def is_saved(self):
        pass

