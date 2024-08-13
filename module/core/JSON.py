import json
from dataclasses import dataclass
from module.core.Cacheable import Cacheable
from typing import ClassVar
from module.core.utils import is_array_like


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

    def get(self, key, default=None):
        return self.dict.get(key, default)

    def values(self):
        return self.dict.values()

    def keys(self):
        return self.dict.keys()

    def items(self):
        return self.dict.items()

    def __contains__(self, key):
        if is_array_like(key):
            key = tuple(key)
        return key in self.list

    def __getitem__(self, key):
        return self.dict.get(key)

    def __setitem__(self, key, value):
        self.add(key, value)

    def __repr__(self) -> str:
        return "\n".join(f"{k}: {v}" for k, v in self.items())

    @property
    def reversed(self):
        return {
            tuple(value) if is_array_like(value) else value: key
            for key, value in self.items()
        }

    @property
    def dict(self):
        return self.load()

    @property
    def list(self) -> list:
        return list(self.dict.keys())
