import os
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import ClassVar
from module.core.Cacheable import Cacheable
import pandas as pd
from module.core.utils import is_array_like

ROOT = os.getcwd()  # This gives terminal location (terminal working dir)


def mask(df: pd.DataFrame, mask_conditions: dict):
    selected = df.index != None  # Select all
    absent_columns = set(mask_conditions) - set(df.columns)
    if absent_columns:
        raise ValueError(
            f"Unknown columns: {absent_columns}, possible columns are {df.columns}"
        )
    for column, value in mask_conditions.items():  # Refine selection
        if value is None:
            print(f"Skipping {column}, .select() ignores None for practical purposes, use 'nan' (str) instead.")
            continue
        elif is_array_like(value):
            sub_selection = df[column].isin(value)
        else:
            if value in ["nan", "notna"]:
                sub_selection = df[column].isna() if value is None else df[column].notna()
            else:
                sub_selection = df[column] == value
        selected &= sub_selection
    return selected


def sub_select(df, selector, copy=True):
    df = df.loc[mask(df, selector)]
    return df.iloc[0] if len(df) == 1 else df


class SelectableDataFrame(pd.DataFrame):

    @property
    def _constructor(self):
        return SelectableDataFrame

    def select(self, **selector) -> "SelectableDataFrame":
        """
        Filter the DataFrame based on a selector.

        Args:
            selector (dict): A dictionary of column conditions to filter by.

        Returns:
            CustomDataFrame: Filtered DataFrame that also includes the select method.
        """
        sub_selection = sub_select(self, selector)
        return sub_selection
    
    def extend(self, df) -> "SelectableDataFrame":
        if isinstance(df, Dataset):
            df = df.df
        common_columns = self.columns.intersection(df.columns).to_list()
        return self.merge(df, on=common_columns, how="left")

@dataclass
class Dataset(Cacheable):

    _loader: ClassVar = None
    _saver: ClassVar = None

    def initialize(self):
        super().initialize()
        print(f"CREATED AND CACHED {self.filepath}")

    def select(self, **selector) -> SelectableDataFrame:
        return self.df.select(**selector)

    @property
    def list(self):
        return list(self.df.to_dict(orient="index").values())

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
    
    def __iter__(self):
        return iter([row for _, row in self.df.iterrows()])


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
        data.to_excel(self.filepath, index=False)

    def load(self):
        return pd.read_excel(self.filepath)