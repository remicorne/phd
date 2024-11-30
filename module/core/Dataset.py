import os
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import ClassVar
from module.core.Constants import COMPOUNDS_AND_REGIONS_CLASSES
from module.core.Cacheable import Cacheable
import pandas as pd
from module.core.utils import is_array_like

ROOT = os.getcwd()  # This gives terminal location (terminal working dir)


def handle_class_selectors(classes, select_value):

    if is_array_like(select_value):
        values = []
        for item in select_value:
            values.extend(classes.get(item, [item]))
        return values

    return classes.get(select_value, select_value)


def mask(df: pd.DataFrame, mask_conditions: dict):
    selected = df.index != None  # Select all
    absent_columns = set(mask_conditions) - set([*df.columns, "index"])
    if absent_columns:
        raise ValueError(
            f"Unknown columns: {absent_columns}, possible columns are {df.columns}"
        )
    for key, value in mask_conditions.items():  # Refine selection
        column = pd.Series(df.index) if key == "index" else df[key]
        if value is None:
            print(
                f"Skipping {column.name}, .select() ignores None for practical purpose s, use 'nan' (str) instead."
            )
        else:

            if callable(value):
                sub_selection = column.apply(value)
            else:
                if key in COMPOUNDS_AND_REGIONS_CLASSES:
                    value = handle_class_selectors(
                        COMPOUNDS_AND_REGIONS_CLASSES[key], value
                    )
                if is_array_like(value):
                    sub_selection = column.isin(value)
                else:
                    if value in ["na", "notna"]:
                        sub_selection = (
                            column.isna() if value == "na" else column.notna()
                        )
                    else:
                        sub_selection = column == value
            selected &= sub_selection
    return selected


def sub_select(df, selector):
    df = df.loc[mask(df, selector)]
    return df


class SelectableDataFrame(pd.DataFrame):

    @property
    def _constructor(self):
        return SelectableDataFrame

    def select(self, **selector) -> "SelectableDataFrame":
        """
        Filter the DataFrame based on a selector.

        Args:
            selector (dict): A dictionary of column conditions to filter by.
            'nan' and 'notna' are supported using strings.
            None is ignored for dict unpacking purposes and because it is not a valid value.

        Returns:
            SelectableDataFrame: Filtered DataFrame that also includes the select method.
            Series: if selection conditions result in a single row
        """
        sub_selection = sub_select(self, selector)
        return sub_selection

    def extend(
        self, df: "Dataset|SelectableDataFrame|pd.DataFrame"
    ) -> "SelectableDataFrame":
        """
        Extend the DataFrame with another DataFrame. Automatically selects common columns.

        Args:
            df (_type_): the df to left join to self

        Returns:
            SelectableDataFrame:  Resulting DataFrame of left join
        """
        if isinstance(df, Dataset):
            df = df.df
        common_columns = self.columns.intersection(df.columns).to_list()
        return self.merge(df, on=common_columns)


@dataclass
class Dataset(Cacheable):
    """
    Base class for datasets ie dataframes stored in Excel or Pickle files.
    Similar to JSONmapping interface for json/dict.
    Actual dataframe is accessed through the df property and read directly from the file.

    Returns:
        Dataset: Wrapper for dataframes
    """

    def select(self, **selector) -> SelectableDataFrame:
        return self.df.select(**selector)

    @property
    def list(self):
        return list(self.df.to_dict(orient="index").values())

    @property
    def df(self) -> SelectableDataFrame:
        return SelectableDataFrame(self.load())

    def extend(self, other) -> SelectableDataFrame:
        """
        Extend the DataFrame with another DataFrame. Automatically selects common columns.

        Args:
            df (_type_): the df to left join to self

        Returns:
            SelectableDataFrame:  Resulting DataFrame of left join
        """
        return self.df.extend(other)

    def replace(self, column, mapping):
        data = self.df
        self.save(data.replace(mapping))

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
        if self.is_saved:
            return self.df._repr_html_()
        else:
            return repr(self)


@dataclass
class PickleDataset(Dataset):
    """
    Dataset wrapper for pickle files

    """

    extension: ClassVar[str] = "pkl"

    def save(self, data: pd.DataFrame, filepath=None):
        data.to_pickle(filepath or self.filepath)

    def load(self) -> SelectableDataFrame:
        return SelectableDataFrame(pd.read_pickle(self.filepath))


@dataclass
class ExcelDataset(Dataset):
    """
    Dataset wrapper for excel files

    """

    extension: ClassVar[str] = "xlsx"

    def save(self, data: pd.DataFrame):
        data.to_excel(self.filepath, index=False)

    def load(self) -> SelectableDataFrame:
        return SelectableDataFrame(pd.read_excel(self.filepath))
