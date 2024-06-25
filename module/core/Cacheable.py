from dataclasses import dataclass, field
from typing import ClassVar
import os, platform, subprocess
from module.core.FileSystem import FileSystem
import re

def sanitize_filename(filename):
    illegal_chars = r'[<>:"/\\|?*]'
    sanitized = re.sub(illegal_chars, '_', filename)
    return sanitized


@dataclass
class Cacheable:

    extension: ClassVar[str] = None
    filename: ClassVar[str] = None
    filepath: str = field(default=None, kw_only=True)
    from_scratch: bool = field(default=False, kw_only=True)

    def __post_init__(self):
        if not self.filepath:
            if not self.filename:
                raise ValueError("Child classes must define filename")
            self.filename = sanitize_filename(self.filename)
            # Automatically extrat relevant params for laction building
            self.filepath = os.path.join(
                FileSystem.get_location(**self.__dict__), self.filename
            )
        if not self.extension:
            raise ValueError("Child classes must define extension") 
        # Remove extension if it has already been added
        filepath, _ = os.path.splitext(self.filepath)
        self.filepath = f"{filepath}.{self.extension}"
        if not self.is_saved or self.from_scratch:
            self.initialize()

    def generate(self):
        raise NotImplementedError(
            "This method should be implemented for all custom Cacheables"
        )

    def initialize(self):
        data = self.generate()
        self.save(data) if data is not None else self.save()

    def load(self):
        raise NotImplementedError(
            "This method should be implemented for all custom Cacheables"
        )

    def save(self, data):
        raise NotImplementedError(
            "This method should be implemented for all custom Cacheables"
        )
        
    def validate(self, data):
        pass
        
    def delete(self):
        os.remove(self.filepath)
                
    def open(self):
        if self.is_saved:
            if platform.system() == "Windows":
                os.startfile(self.filepath)
            elif platform.system() == "Darwin":
                subprocess.call(("open", self.filepath))
            elif platform.system() == "Linux":
                print("Can't handle Linux")
            else:
                raise OSError("Unknown operating system")
        else:
            raise FileNotFoundError(self.filepath)

    @property
    def is_saved(self):
        return os.path.isfile(self.filepath)
