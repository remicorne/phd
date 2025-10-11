import json
from dataclasses import dataclass
from module.core.Cacheable import Cacheable
from typing import ClassVar
from module.core.utils import is_array_like


@dataclass
class JSONMapping(Cacheable):
    """Base class for JSON mappings.
    Abstract and currentty only used for constants.
    Used to manipulate JSON file as runtime dictionnaries while maintaing possibility to manually edit and have changes be reflected at runtime.
    """

    extension: ClassVar[str] = "json"

    def __post_init__(self):
        super().__post_init__()
        self.dict = self.load()

    def load(self):
        with open(self.filepath, "r", encoding="utf-8") as outfile:
            return json.load(outfile)

    def save(self, data):
        with open(self.filepath, "w", encoding="utf-8") as json_file:
            json.dump(data, json_file)

    def add(self, key, value):
        self.dict[key] = value
        self.save(self.dict)

    def get(self, key, default=None):
        return self.dict.get(tuple(key) if is_array_like(key) else key, default)

    def get_many(self, key, default=None):
        key = [key] if not is_array_like(key) else key
        values = []
        for subkey in key:
            if subkey in self:
                values.extend(self[subkey])
        return values or default

    def values(self) -> list:
        return list(self.dict.values())

    def keys(self) -> list:
        return list(self.dict.keys())

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

    def get_key(self, value_to_find, default=None):
        return {
            tuple(value) if is_array_like(value) else value: key
            for key, value in self.items()
        }.get(
            tuple(value_to_find) if is_array_like(value_to_find) else value_to_find,
            default,
        )

    @property
    def list(self) -> list:
        return list(self.keys())
