import pandas as pd
from typing import Dict, List
from module.core.constants import DatasetColumn
from module.core.Dataset import SelectableDataFrame


class Ratio:
    def __init__(self, column: str, value: str):
        """
        Initialize with columnumn=ratio_string, e.g. compound="DA/5HT" for ratio in all regions
        or compound="DA/5HT", region="VTA" for ratio in VTA or compound="DA/5HT", region="VTA/OF" for ratio in VTA and OF
        """
        if not Ratio.is_ratio(value):
            raise ValueError(f"{value} is not a valid ratio")
        self.column = column
        self.value = value
        self.num, self.den = value.split("/")

    def is_ratio(value: str) -> bool:
        components = value.split("/")
        return len(components) == len(set(components)) == 2

    def compute(self, data: SelectableDataFrame) -> SelectableDataFrame:
        """
        Compute the ratio DataFrame using the ratio, assuming:
        - 'value' column holds the numerical values
        - 'unit' column holds the units
        - All non-specified columns are matched exactly
        """
        if self.column not in data.columns:
            raise ValueError(f"Data does not have: {self.column}")

        num_df, den_df = (
            data[data[self.column] == self.num],
            data[data[self.column] == self.den],
        )

        ratio_columns = [self.column, DatasetColumn.VALUE, DatasetColumn.UNIT]
        merge_cols = list(set(data.columns) - set(ratio_columns))
        merged = pd.merge(num_df, den_df, on=merge_cols, suffixes=("_num", "_den"))
        merged[DatasetColumn.VALUE] = (
            merged[f"{DatasetColumn.VALUE}_num"] / merged[f"{DatasetColumn.VALUE}_den"]
        )

        def resolve_unit(row):
            unit_num, unit_den = (
                f"{DatasetColumn.UNIT}_num",
                f"{DatasetColumn.UNIT}_den",
            )
            return (
                f"{row[unit_num]}/{row[unit_den]}"
                if row[unit_num] != row[unit_den]
                else ""
            )

        merged[DatasetColumn.UNIT] = merged.apply(resolve_unit, axis=1)
        merged[self.column] = self.value
        return merged[data.columns]
