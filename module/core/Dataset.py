import os
import pandas as pd
from dataclasses import dataclass
from typing import ClassVar, Iterable, Callable
from abc import ABC, abstractmethod
from module.core.Cacheable import Cacheable
from module.core.utils import is_array_like

ROOT = os.getcwd()  # This gives terminal location (terminal working dir)


class SelectionError(Exception):
    pass


def mask(df: pd.DataFrame, mask_conditions: dict[str, str | Iterable | Callable]):
    """
    Select rows in df based on mask_conditions.
    None is considered a wildcard and selects everything
    'na' and 'notna' are supported using strings.
    """
    selected = df.index.notna()  # Select all
    absent_columns = set(mask_conditions) - set(df.columns)
    if absent_columns:
        raise ValueError(
            f"Unknown columns: {absent_columns}, possible columns are {df.columns}"
        )
    for key, value in mask_conditions.items():  # Refine selection
        column = df[key]
        if value is None:
            continue
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


def sub_select(df: "CustomDataFrame", selector: dict) -> "CustomDataFrame":
    df = df.loc[mask(df, selector)]
    return df.copy()


class CustomDataFrame(pd.DataFrame):
    @property
    def _constructor(self):
        return CustomDataFrame

    def select(
        self, *, select_one=False, **selector: dict[str, str | Iterable | Callable]
    ) -> "CustomDataFrame":
        """
        Filter rows using a column → condition mapping.

        Applies one or more **selectors** to return only the rows that satisfy *all*
        conditions (logical AND). Each selector targets a column (or the special key
        `"index"`) and supports several matching modes:

        - **Literal**: `col=value` keeps rows where `col == value`.
        - **Iterable**: `col=[v1, v2, ...]` uses membership (`.isin(...)`).
        - **Null checks**: `col="na"` keeps `NA/NaN`; `col="notna"` keeps non-nulls.
        - **Wildcard**: `col=None` is ignored (selects everything for that key).
        - **Predicate**: `col=lambda x: <bool>` keeps rows where the function returns `True`.

        Args:
            select_one (bool, optional): If True, returns only one row. Defaults to False.
                Raises SelectionError if no rows are found or if multiple rows are found.
            **selector: dict[str, Any]
                Mapping of column names (or `"index"`) to selector values as described above.

        Returns:
            SelectableDataFrame
                A new DataFrame of the same type as the input containing rows that match
                all selectors. For categorical columns included in `selector`, unused
                categories are removed.

        Raises:
            ValueError: If any key in `selector` does not refer to an existing column.

        Examples:
            >>> df.select(status="open", priority=[1, 2, 3])
            >>> df.select(score=lambda s: s > 0)
            >>> df.select(category="na")  # keep rows where 'category' is missing
            >>> df.select(index=lambda i: i.str.startswith("A"))
        """
        if unknown_cols := set(selector.keys()) - set(self.columns):
            raise ValueError(f"Unknown columns: {unknown_cols}")
        catgorical_cols = [col for col in selector if self[col].dtype == "category"]
        sub_selection = sub_select(self, selector)
        if sub_selection.empty:
            raise SelectionError(f"No rows found for selector: {selector}")
        if select_one and len(sub_selection) > 1:
            raise SelectionError(f"Multiple rows found for selector: {selector}")
        for col in catgorical_cols:
            sub_selection.loc[:, col] = sub_selection[
                col
            ].cat.remove_unused_categories()
        return sub_selection.iloc[0] if select_one else sub_selection

    def extend(
        self, other: "CachedDataFrame|CustomDataFrame|pd.DataFrame"
    ) -> "CustomDataFrame":
        """
        Extend the DataFrame with another DataFrame. Automatically selects common columns.

        Args:
            other (_type_): the other to left join to self

        Returns:
            CustomDataFrame:  Resulting DataFrame of left join
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

    def select(self, *, select_one=False, **selector) -> CustomDataFrame:
        """
        Filter rows using a column → condition mapping.

        Applies one or more **selectors** to return only the rows that satisfy *all*
        conditions (logical AND). Each selector targets a column (or the special key
        `"index"`) and supports several matching modes:

        - **Literal**: `col=value` keeps rows where `col == value`.
        - **Iterable**: `col=[v1, v2, ...]` uses membership (`.isin(...)`).
        - **Null checks**: `col="na"` keeps `NA/NaN`; `col="notna"` keeps non-nulls.
        - **Wildcard**: `col=None` is ignored (selects everything for that key).
        - **Predicate**: `col=lambda x: <bool>` keeps rows where the function returns `True`.

        Args:
            select_one (bool, optional): If True, returns only one row. Defaults to False.
                Raises SelectionError if no rows are found or if multiple rows are found.
            **selector: dict[str, Any]
                Mapping of column names (or `"index"`) to selector values as described above.

        Returns:
            SelectableDataFrame
                A new DataFrame of the same type as the input containing rows that match
                all selectors. For categorical columns included in `selector`, unused
                categories are removed.

        Raises:
            ValueError: If any key in `selector` does not refer to an existing column.

        Examples:
            >>> df.select(status="open", priority=[1, 2, 3])
            >>> df.select(score=lambda s: s > 0)
            >>> df.select(category="na")  # keep rows where 'category' is missing
            >>> df.select(index=lambda i: i.str.startswith("A"))
        """
        return self.df.select(select_one=select_one, **selector)

    @property
    def list(self):
        return list(self.df.to_dict(orient="index").values())

    @property
    @abstractmethod
    def df(self) -> CustomDataFrame:
        """The dataframe property that must be implemented by subclasses."""

    def extend(self, other) -> CustomDataFrame:
        """
        Extend the DataFrame with another DataFrame. Automatically selects common columns.

        Args:
            df (_type_): the df to left join to self

        Returns:
            CustomDataFrame:  Resulting DataFrame of left join
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
    def load(self, **kwargs) -> CustomDataFrame:
        pass

    def select(self, *, select_one=False, **selector) -> CustomDataFrame:
        try:
            return super().select(select_one=select_one, **selector)
        except SelectionError as e:
            raise SelectionError(f"{e} for {self.filename}")

    @property
    def df(self) -> CustomDataFrame:
        return self.load()


@dataclass
class PickleCachedDataFrame(CachedDataFrame):
    """
    CachedDataFrame wrapper for pickle files

    """

    extension: ClassVar[str] = "pkl"

    def save(self, data: pd.DataFrame, filepath=None):
        data.to_pickle(filepath or self.filepath)

    def load(self, **kwargs) -> CustomDataFrame:
        return CustomDataFrame(pd.read_pickle(self.filepath, **kwargs))


@dataclass
class ExcelCachedDataFrame(CachedDataFrame):
    """
    CachedDataFrame wrapper for excel files

    """

    sheet_name: ClassVar[str] = None
    extension: ClassVar[str] = "xlsx"

    def save(self, data: pd.DataFrame):
        data.to_excel(self.filepath, index=False)

    def load(self, **kwargs) -> CustomDataFrame:
        data = pd.read_excel(self.filepath, sheet_name=self.sheet_name, **kwargs)
        return CustomDataFrame(data) if isinstance(data, pd.DataFrame) else data

    def select(self, *, select_one=False, **selector) -> CustomDataFrame:
        try:
            return super().select(select_one=select_one, **selector)
        except SelectionError as e:
            raise SelectionError(
                str(e) + f" - {self.sheet_name}" if self.sheet_name else e
            )
