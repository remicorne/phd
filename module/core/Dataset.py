import os, platform, subprocess, re, sys
import pandas as pd
from dataclasses import dataclass
from module.core.Questions import Questions
from typing import ClassVar

ROOT = os.getcwd()  # This gives terminal location (terminal working dir)
CACHE = "cache"


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
class Dataset:

    project: str
    _type: ClassVar[str] = None

    def __post_init__(self):
        self.location = f"{ROOT}/{self.project}/{CACHE}"
        file_path = f"{self.location}/{self._type}"
        self.pkl_path = f"{file_path}.pkl"
        self.excel_path = f"{file_path}.xlsx"
        if not os.path.exists(self.location):
            os.mkdir(self.location)
        if not self.is_pickeled:
            self.initialize()
            
    def generate_data(self):
        return pd.DataFrame()

    def initialize(self):
        self.save(self.generate_data())
        print(f"CREATED AND CACHED {self.project}/{self._type}")
            
    def save(self, dataset):
        dataset.to_pickle(self.pkl_path)
        try:
            dataset.to_excel(self.excel_path)
        except ValueError as e:
            if "This sheet is too large" in str(e):
                print(f"COULD NOT SAVE TO EXCEL: {self._type} TOO LARGE")

    def get(self):
        if not self.is_pickeled:
            self.initialize()
        return pd.read_pickle(self.pkl_path)
            
    def edit_excel(self):
        self.open_excel()
        user_finished = Questions.yes_or_no("Have you finished editing and saved the file?")
        updated_dataset = pd.read_excel(self.excel_path)
        self.save(updated_dataset)

    def select(self, selector):
        return sub_select(self.get(), selector)

    def open_excel(self):
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

    @property
    def is_pickeled(self):
        return os.path.isfile(self.pkl_path)

    @property
    def is_exceled(self):
        return os.path.isfile(self.pkl_path)

