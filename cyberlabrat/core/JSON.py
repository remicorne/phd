import json
from dataclasses import dataclass
from cyberlabrat.core.Cacheable import Cacheable
from typing import ClassVar
from cyberlabrat.core.utils import is_array_like


@dataclass
class JSONMapping(Cacheable):
    """Base class for JSON mappings.
    Abstract and currentty only used for constants.
    Used to manipulate JSON file as runtime dictionnaries while maintaing possibility to manually edit and have changes be reflected at runtime.
    """

    extension: ClassVar[str] = "json"

    def load(self) -> dict:
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
        return self.dict.get(tuple(key) if is_array_like(key) else key, default)

    def values(self) -> list:
        return self.dict.values()

    def keys(self) -> list:
        return self.dict.keys()

    def items(self) -> list:
        return self.dict.items()

    def __contains__(self, key):
        if is_array_like(key):
            key = tuple(key)
        return key in self.list

    def __getitem__(self, key):
        return self.dict.get(key)

    def __setitem__(self, key, value):
        self.add(key, value)


    def __iter__(self):
        for item in self.list:
            yield item
        
    def __repr__(self) -> str:
        return "\n".join(f"{k}: {v}" for k, v in self.items())

    @property
    def reversed(self):
        return {
            tuple(value) if is_array_like(value) else value: key
            for key, value in self.items()
        }

    @property
    def dict(self) -> dict:
        return self.load()

    @property
    def list(self) -> list:
        return list(self.dict.keys())
