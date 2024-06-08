import os
import pandas as pd
from dataclasses import dataclass, field
from typing import ClassVar
from module.core.Cacheable import Cacheable
import pandas as pd


ROOT = os.getcwd()  # This gives terminal location (terminal working dir)


def mask(df: pd.DataFrame, mask_conditions: dict):
    selected = df.index != None # Select all
    nan = mask_conditions.pop('nan', None)
    absent_columns = [col for col in mask_conditions if col not in df]
    if absent_columns:
        raise ValueError(f"Unknown columns: {absent_columns}, possible columns are {df.columns}")
    if nan is not None:
        selected &= df.value.isna() if nan else df.value.notna()
    for column, value in mask_conditions.items(): # Refine selection
        if isinstance(value, list):
            sub_selection = df[column].isin(value)
        else:
            sub_selection = df[column] == value
        selected &= sub_selection
    return selected


def sub_select(df, selector):
    df = df[mask(df, selector)].copy()
    return df

class SelectableDataFrame(pd.DataFrame):
    
    @property
    def _constructor(self):
        return SelectableDataFrame

    def select(self, **selector):
        """
        Filter the DataFrame based on a selector.

        Args:
            selector (dict): A dictionary of column conditions to filter by.

        Returns:
            CustomDataFrame: Filtered DataFrame that also includes the select method.
        """
        sub_selection = sub_select(self, selector)
        # if sub_selection.empty:
        #     raise ValueError(f"NO DATA FOR SELECTION: {selector}")
        return sub_selection
    
    def extend(self, other):
        if isinstance(other, Dataset):
            df = other.df
        common_columns = self.columns.intersection(df.columns).tolist()
        return self.merge(df, on=common_columns, how="left")

@dataclass
class Dataset(Cacheable):
    
    _loader: ClassVar = None
    _saver: ClassVar = None
    
    def initialize(self):
        super().initialize()
        print(f"CREATED AND CACHED {self.filepath}")

    def validate(self, data):
        if data.empty:
            raise ValueError("NO DATA")

    def select(self, **selector) -> SelectableDataFrame:
        return self.df.select(**selector)
        
    @property
    def list(self):
        return self.df.to_dict(orient="index").values()

    @property
    def df(self) -> SelectableDataFrame:
        return SelectableDataFrame(self.load())
    
    def extend(self, other) -> SelectableDataFrame:
        return self.df.extend(other)
    
    def __contains__(self, column):
        return column in self.df
        
    def __repr__(self) -> str:
        """Called by terminal to display the dataframe (pretty)

        Returns:
            str: Pretty representation of the df
        """
        return repr(self.df)
    
    def _repr_html_(self) -> str:
        """Called by jupyter notebook to display the dataframe as html (pretty)

        Returns:
            str: Pretty representation of the df
        """
        return self.df._repr_html_()
    
@dataclass
class PickleDataset(Dataset):

    extension: ClassVar[str] = "pkl"
    
    def save(self, data):
        data.to_pickle(self.filepath)

    def load(self):
        return pd.read_pickle(self.filepath)


@dataclass
class ExcelDataset(Dataset):

    extension: ClassVar[str] = "xlsx"
    
    def save(self, data):
        data.to_excel(self.filepath, ignore_index=True)

    def load(self):
        return pd.read_excel(self.filepath)
