import os, platform, subprocess
import pandas as pd
from dataclasses import dataclass
from typing import ClassVar
from module.core.Cacheable import Cacheable

ROOT = os.getcwd()  # This gives terminal location (terminal working dir)


def mask(df, mask_conditions):
    selected = True
    for column, value in mask_conditions.items():
        if isinstance(value, list):
            sub_selection = df[column].isin(value)
        else:
            sub_selection = df[column] == value
        selected &= sub_selection
    return selected


def sub_select(df, selector):
    df = df[mask(df, selector)]
    return df


@dataclass
class Dataset(Cacheable):

    _name: ClassVar[str] = None
    _template: ClassVar[pd.DataFrame] = None

    def validate(self, data):
        if data.empty:
            raise ValueError("NO DATA")

    def initialize(self):
        super().initialize()
        print(f"CREATED AND CACHED {self.filepath}")
            
    def save(self, data):
        print(f"Caching {self._name} dataframe")
        data.to_pickle(self.pkl_path)
        print(f"Cached {self._name} datarame to {self.pkl_path}")
        print(f"Saving {self._name} datarame (may take a minute or two)")
        data.to_excel(self.excel_path)
        print(f"Saved {self._name} datarame to {self.excel_path}")

    def load(self):
        return pd.read_pickle(self.pkl_path)
            
    def select(self, selector):
        sub_selection = sub_select(self.get(), selector)
        if sub_selection.empty:
            raise ValueError(f"EMPTY SELECTION: {selector}")
        return sub_selection

    def _open_excel(self):
        if self.is_exceled:
            if platform.system() == "Windows":
                os.startfile(self.excel_path)
            elif platform.system() == "Darwin":
                subprocess.call(("open", self.excel_path))
            elif platform.system() == "Linux":
                subprocess.call(("xdg-open", self.excel_path))
            else:
                raise OSError("Unknown operating system")
        else:
            raise FileNotFoundError(self.excel_path)
        
    def clear_cache(self):
        if self.is_pickeled:
            os.remove(self.pkl_path)
            
    def delete(self):
        self.clear_cache()
        if self. is_exceled:
            os.remove(self.excel_path)
    
    @property
    def filepath(self):
        return f"{self.location}/{self._name}"

    # Probably useless to have all these properties if were only saving proper stuff
    @property
    def pkl_path(self):
        return f"{self.filepath}.pkl"

    @property
    def excel_path(self):
        return f"{self.filepath}.xlsx"
    
    @property
    def is_saved(self):
        return self.is_pickeled

    @property
    def is_pickeled(self):
        return os.path.isfile(self.pkl_path)

    @property
    def is_exceled(self):
        return os.path.isfile(self.excel_path)

    @property
    def df(self):
        return self.get()
