from dataclasses import dataclass, field
from typing import ClassVar
from module.core.Constants import COMPOUND_CLASSES, REGION_CLASSES

TYPES = ["compound", "region"]

@dataclass
class StandardizedParameter:
    
    TYPES: ClassVar = field(init=False, default=["compound", "region"])
    
    name: str = field(kw_only=True)
    type: str = field(kw_only=True)
    
    def __new__(cls, name, type):
        instance = super().__new__(cls)
        return instance
    
    def __post_init__(self):
        classes = {"region": REGION_CLASSES, "compound": COMPOUND_CLASSES}[self.type]
        self.classes = classes.get_item_classes(self.name)

    @staticmethod
    def includes(type):
        return type in TYPES
    
    def __eq__(self, value: str) -> bool:
        return self.name == value or value in self.classes
    
    def __str__(self) -> str:
        return self.name
    
    def __str__(self) -> str:
        return str(self)

