import os, platform, subprocess
import pandas as pd
from dataclasses import dataclass
from module.core.Questions import Questions
from typing import ClassVar
from module.core.ProjectMember import ProjectMember

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
class Dataset(ProjectMember):

    _name: ClassVar[str] = None
    _subfolder: ClassVar[str] = 'cache'
    _with_validation: ClassVar[str] = False
    _template: pd.DataFrame = None

    def __post_init__(self):
        file_path = f"{self.location}/{self._name}"
        self.pkl_path = f"{file_path}.pkl"
        self.excel_path = f"{file_path}.xlsx"
        if not self.is_saved:
            self.initialize()
            
    def generate(self):
        data = pd.DataFrame(self._template)
        return data
    
    def validate(self, data):
        if data.empty:
            raise ValueError("NO DATA")

    def initialize(self):
        super().initialize()
        print(f"CREATED AND CACHED {self.project}/{self._name}")
            
    def save(self, dataset):
        if self._with_validation:
            self.validate(dataset)
        dataset.to_pickle(self.pkl_path)
        try:
            dataset.to_excel(self.excel_path)
        except ValueError as e:
            if "This sheet is too large" in str(e):
                print(f"COULD NOT SAVE TO EXCEL: {self._name} TOO LARGE")

    def load(self):
        return pd.read_pickle(self.pkl_path)
            
    def edit_excel_and_save(self):
        self._open_excel()
        user_finished = Questions.yes_or_no("Have you finished editing and saved the file?")
        while not user_finished:
            user_finished = Questions.yes_or_no("Have you finished editing and saved the file?")
        updated_dataset = pd.read_excel(self.excel_path)
        self.save(updated_dataset)

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
        
    def delete(self):
        if self.is_pickeled:
            os.remove(self.pkl_path)
        if self. is_exceled:
            os.remove(self.pkl_path)
        
    @property
    def is_saved(self):
        return self.is_pickeled

    @property
    def is_pickeled(self):
        return os.path.isfile(self.pkl_path)

    @property
    def is_exceled(self):
        return os.path.isfile(self.excel_path)

