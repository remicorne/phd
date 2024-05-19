from dataclasses import dataclass, field
from module.constants import ROOT
import os

PROJECTS = f"{ROOT}/PROJECTS"


@dataclass
class Cacheable:

    location: str
    filepath: str = field(init=False, default=None)

    def __post_init__(self):
        if self.filepath is None:
            raise ValueError(
                "Subclasses must initialize 'filepath' before calling super().__post_init__()"
            )
        if not self.is_saved:
            self.initialize()

    def generate(self):
        raise NotImplementedError(
            "This method should be implemented for all custom Cacheables"
        )

    def initialize(self):
        data = self.generate()
        self.save(data)

    def load(self):
        raise NotImplementedError(
            "This method should be implemented for all custom Cacheables"
        )

    def save(self, data):
        raise NotImplementedError(
            "This method should be implemented for all custom Cacheables"
        )

    @property
    def is_saved(self):
        return os.path.isfile(self.filepath)

