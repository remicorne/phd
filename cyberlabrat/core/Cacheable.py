from dataclasses import dataclass, field
from typing import ClassVar
import os, platform, subprocess
from cyberlabrat.core.FileSystem import FileSystem
from cyberlabrat.core.questions import input_escape
import re


def sanitize_filename(filename):
    illegal_chars = r'[<>:"/\\|?*]'
    sanitized = re.sub(illegal_chars, "_", filename)
    return sanitized


@dataclass
class Cacheable:
    """
    Base Abstract class for all Cacheable objects.
    The point of cacheables is to have runtime-like variables that reflect the state of what is stored to disk.
    Especially useful for configuration excel files such as Palette.
    Handles loading, saving and filesystem architecture (see FileSystem.py)
    Unless filepath is provided, it will be generated from filename and extension.
    If filepath is provided, it will be used as-is.
    All child classes must implement the generate(), load() and save() methods.
    Generate must either return something that can be saved or the save method should know how to handle it.

    Args:
        filepath (str, optional): The path to the file. Defaults to None.

    Returns:
        object: What the child class returns
    """

    extension: ClassVar[str] = None
    filename: ClassVar[str] = None
    filepath: str = field(default=None, kw_only=True)
    from_scratch: bool = field(default=False, kw_only=True)

    def __post_init__(self):
        """
        Generates the filepath from the filename and extension.
        If the Cacheable is not saved, it will be initialized.

        """
        if not self.filepath:
            if not self.filename:
                raise ValueError("Child classes must define filename")
            self.filename = sanitize_filename(self.filename)
            # Automatically extrat relevant params for laction building
            self.filepath = os.path.join(FileSystem.get_location(**self.__dict__))
        if not self.extension:
            raise ValueError("Child classes must define extension")
        # Remove extension if it has already been added
        filepath, _ = os.path.splitext(self.filepath)
        self.filepath = f"{filepath}.{self.extension}"
        if not self.is_saved or self.from_scratch:
            self.initialize()

    def get_location(self):
        return FileSystem.get_location(**self.__dict__)

    def get_filename(self):
        return f"{sanitize_filename(self.filename or input_escape('Enter filename'))}.{self.extension}"

    def get_filepath(self):
        return self.filepath or os.path.join(self.get_location(), self.get_filename())

    def generate(self):
        raise NotImplementedError(
            "This method should be implemented for all custom Cacheables"
        )

    def initialize(self):
        """
        Calls the generate method and saves what it return to disk.
        """
        data = self.generate()
        self.save(data) if data is not None else self.save()
        print(f"CREATED AND CACHED {self.filepath}")

    def load(self):
        raise NotImplementedError(
            "This method should be implemented for all custom Cacheables"
        )

    def save(self, data):
        raise NotImplementedError(
            "This method should be implemented for all custom Cacheables"
        )

    def delete(self):
        os.remove(self.get_filepath())

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
