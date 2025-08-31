import os
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import ClassVar
from abc import ABC, abstractmethod
from module.core.Registry import Registry
from module.core.Cacheable import Cacheable
import pandas as pd
from module.core.utils import is_array_like

ROOT = os.getcwd()  # This gives terminal location (terminal working dir)


class SelectionError(Exception):
    pass


def mask(df: pd.DataFrame, mask_conditions: dict):
    selected = df.index.notna()  # Select all
    absent_columns = set(mask_conditions) - set([*df.columns, "index"])
    if absent_columns:
        print("mask", df, df.columns)
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


def sub_select(df: "SelectableDataFrame", selector: dict) -> "SelectableDataFrame":
    df = df.loc[mask(df, selector)]
    return df.copy()


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
        if unknown_cols := set(selector.keys()) - set(self.columns):
            raise ValueError(f"Unknown columns: {unknown_cols}")
        catgorical_cols = [col for col in selector if self[col].dtype == "category"]
        sub_selection = sub_select(self, selector)
        for col in catgorical_cols:
            sub_selection.loc[:, col] = sub_selection[
                col
            ].cat.remove_unused_categories()
        return sub_selection

    def select_one(self, **selector):
        sub_selection = self.select(**selector)
        if sub_selection.empty:
            raise SelectionError(f"No rows found for selector: {selector}")
        if len(sub_selection) > 1:
            raise SelectionError(f"Multiple rows found for selector: {selector}")
        return sub_selection.iloc[0]

    def extend(
        self, other: "CachedDataFrame|SelectableDataFrame|pd.DataFrame"
    ) -> "SelectableDataFrame":
        """
        Extend the DataFrame with another DataFrame. Automatically selects common columns.

        Args:
            other (_type_): the other to left join to self

        Returns:
            SelectableDataFrame:  Resulting DataFrame of left join
        """
        if isinstance(other, CachedDataFrame):
            other = other.df
        common_columns = self.columns.intersection(other.columns).to_list()
        if not common_columns:
            return self.merge(other, how="cross")
        return self.merge(other, on=common_columns)


class DataframeWrapperMixin(ABC):
    """
    Mixin for datasets ie dataframes stored in Excel or Pickle files.
    Similar to JSONmapping interface for json/dict.
    Actual dataframe is accessed through the df property and read directly from the file.

    Returns:
        CachedDataFrame: Wrapper for dataframes
    """

    def select(self, **selector) -> SelectableDataFrame:
        return self.df.select(**selector)

    def select_one(self, **selector):
        return self.df.select_one(**selector)

    @property
    def list(self):
        return list(self.df.to_dict(orient="index").values())

    @property
    @abstractmethod
    def df(self) -> SelectableDataFrame:
        """The dataframe property that must be implemented by subclasses."""

    def extend(self, other) -> SelectableDataFrame:
        """
        Extend the DataFrame with another DataFrame. Automatically selects common columns.

        Args:
            df (_type_): the df to left join to self

        Returns:
            SelectableDataFrame:  Resulting DataFrame of left join
        """
        return self.df.extend(other)

    def __contains__(self, column):
        return column in self.df

    def __setitem__(self, key, value):
        self.df[key] = value


class CachedDataFrame(Cacheable, DataframeWrapperMixin):
    @abstractmethod
    def save(self, data: pd.DataFrame):
        pass

    @abstractmethod
    def load(self, **kwargs) -> SelectableDataFrame:
        pass

    def select_one(self, **selector) -> SelectableDataFrame:
        try:
            return super().select_one(**selector)
        except SelectionError as e:
            raise SelectionError(f"{e} for {self.filename}")

    @property
    def df(self) -> SelectableDataFrame:
        return self.load()


@dataclass
class PickleCachedDataFrame(CachedDataFrame):
    """
    CachedDataFrame wrapper for pickle files

    """

    extension: ClassVar[str] = "pkl"

    def save(self, data: pd.DataFrame, filepath=None):
        data.to_pickle(filepath or self.filepath)

    def load(self, **kwargs) -> SelectableDataFrame:
        return SelectableDataFrame(pd.read_pickle(self.filepath, **kwargs))


@dataclass
class ExcelCachedDataFrame(CachedDataFrame):
    """
    CachedDataFrame wrapper for excel files

    """

    sheet_name: ClassVar[str] = None
    extension: ClassVar[str] = "xlsx"

    def save(self, data: pd.DataFrame):
        data.to_excel(self.filepath, index=False)

    def load(self, **kwargs) -> SelectableDataFrame:
        data = pd.read_excel(self.filepath, sheet_name=self.sheet_name, **kwargs)
        return SelectableDataFrame(data) if isinstance(data, pd.DataFrame) else data

    def select_one(self, **selector) -> SelectableDataFrame:
        try:
            return super().select_one(**selector)
        except SelectionError as e:
            raise SelectionError(
                str(e) + f" - {self.sheet_name}" if self.sheet_name else e
            )
