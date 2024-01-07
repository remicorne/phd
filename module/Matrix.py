from dataclasses import dataclass, field
import scipy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from module.utils import subselectDf
from module.utils import (
    parallel_process,
)


def correlate(method, pvalue_threshold):
    def executor(x, y):
        match method:
            case "pearson":
                result = scipy.stats.pearsonr(x, y)
                correlation, pvalue = result.statistic, result.pvalue
            case "spearman":
                result = scipy.stats.spearmanr(x, y)
                correlation, pvalue = result.correlation, result.pvalue
            case "kendall":
                result = scipy.stats.kendalltau(x, y)
                correlation, pvalue = result.correlation, result.pvalue
        if pvalue < pvalue_threshold:  # Return correlation score if significant
            return correlation
        else:  # Return 0 otherwiss, this is done to avoid having to generate a mask which takes time
            return 0

    return executor


@dataclass
class Matrix:

    """Creates a reusable matrix class
    selects the relevant data from the larer set
    creates a matrix with it
    eliminates date where n < n_minimum
    Args:
        data (pd.DataFrame): the original df
        treatment (str): the treatment group
        between (str): variables to correlate (compound or region)
        variables (str): list of variables to correlate from 'between'. If only one: self correlation
        accross (str): the column that will constitute the rows/cols of the matrix
        n_minimum (int): minumum occurnces of overlapping var1 and var2. Default to 5
        columns (list[str]): List to select from the 'across' column. Defaults to None
        method (str): one of 'pearson', 'kendall', 'spearman'. Defaults to "pearson"
        pvalue_threshold: float = 0.05

    Returns:
        filtered_data (pd.DataFrame): Data after subselection and replacing zeros with NaN.
        pivot (pd.DataFrame): Result of the data pivot operation.
        corr (pd.DataFrame): correlations value or 0.0 if not significant or nan if insuficient_overlapp
        missing_values (list): List of missing values (< n_minimum).
        missing_overlap (list): List of pairs with insufficient overlap.

    """

    data: pd.DataFrame
    grouping: str
    between: str
    var1: str
    var2: str
    accross: str
    columns: list[str] = None
    n_minimum: int = 5
    method: str = "pearson"
    pvalue_threshold: float = 0.05

    filtered_data: pd.DataFrame = field(init=False)
    pivot: pd.DataFrame = field(init=False)
    corr: pd.DataFrame = field(init=False)
    missing_values: list = field(init=False)
    missing_overlap: list = field(init=False)

    def __post_init__(self):
        self.filter_missing_values()
        self.pivot_data()
        self.order_columns()
        self.calculate_correlations()
        self.find_missing_overlap()
        self.process_triangle_correlogram()
        self.replace_insignificant_correlations()

    def filter_missing_values(self):
        self.missing_values = []
        missing_indices = []
        for col, df in self.data.groupby(by=[self.between, self.accross]):
            if df.value.notna().sum() < self.n_minimum:
                self.missing_values.append(col)
                missing_indices.extend(df.index)
        self.filtered_data = self.data.drop(missing_indices)
        if self.missing_values:
            print(
                f"{self.grouping} missing data for {self.missing_values}, deleted from analysis"
            )

    def pivot_data(self):
        self.pivot = self.filtered_data.pivot_table(
            values="value",
            index=self.filtered_data["mouse_id"],
            columns=[self.between, self.accross],
        )

    def order_columns(self):
        columns = (
            sorted(
                self.pivot.columns,
                key=lambda x: self.columns.index(x[1])
                if x[1] in self.columns
                else float("inf"),
            )
            if self.columns
            else self.pivot.columns
        )
        self.pivot = self.pivot[columns]

    def calculate_correlations(self):
        method = correlate(self.method, self.pvalue_threshold)
        self.corr = self.pivot.corr(method=method, min_periods=self.n_minimum + 1).loc[
            tuple([self.var1, self.var2])
        ]

    def find_missing_overlap(self):
        self.missing_overlap = [
            ((self.var1, index), (self.var2, column))
            for (index, column), is_na in self.corr.T.isna().stack().items()
            if is_na
        ]
        if self.missing_overlap:
            print(
                f"{self.treatment} insuficient overlapp for {self.missing_ovelapp} pairs"
            )
            print("Inspect with self.corr to adjust {columns} and redo analysis")

    def process_triangle_correlogram(self):
        if not self.is_square:
            mask = np.triu(np.ones(self.corr.shape, dtype=bool), k=1)
            self.corr[mask] = np.nan
            np.fill_diagonal(self.corr.values, 1)

    def replace_insignificant_correlations(self):
        self.corr = self.corr.where(self.corr != 0, np.nan)

    @property
    def is_square(self):
        return len(self.var1) != len(self.var2)

    def get_title(self):
        if self.is_square:
            return f"{'-'.join([self.var1, self.var2])} in {self.grouping}"
        return f"{self.var1} in {self.grouping}"


@dataclass
class Matrices:
    """
    Creates a collection of Matrix objects from a larger dataset. This class
    helps in processing and analyzing data by grouping, selecting relevant variables,
    and building individual matrices for further analysis.

    Attributes:
        data (pd.DataFrame): The original dataframe containing the data.
        group_by (str): The column name in 'data' to group by.
        between (str): The first variable to correlate (e.g., compound or region).
        variables (str): The specific variables to correlate from 'between'.
        accross (str): The second variable to correlate against 'between'.
        sub_selector (str): Additional filtering criteria for sub-selecting the data.
        columns (list[str]): Columns to select from the 'accross' column. Defaults to None.
        n_minimum (int): Minimum number of occurrences for a valid correlation. Defaults to 5.
        method (str): Correlation method, one of 'pearson', 'kendall', 'spearman'. Defaults to "pearson".
        pvalue_threshold (float): P-value threshold for significance. Defaults to 0.05.

    Returns:
        matrices (list[Matrix]): A list of Matrix objects created from the grouped data.
        var1 (str): The first variable derived from 'variables'.
        var2 (str): The second variable derived from 'variables'.
    """
    data: pd.DataFrame
    group_by: str
    between: str
    variables: str
    accross: str
    sub_selector: str = None
    columns: list[str] = None
    n_minimum: int = 5
    method: str = "pearson"
    pvalue_threshold: float = 0.05

    matrices: list[Matrix] = field(init=False)
    var1: str = field(init=False)
    var2: str = field(init=False)

    def __post_init__(self):
        self.define_variables()
        self.sub_select()
        self.build_matrices()
        self.homogenize_datasets()

    def define_variables(self):
        self.variables = (
            self.variables if len(self.variables) == 2 else self.variables * 2
        )
        self.var1, self.var2 = self.variables

    def sub_select(self):
        self.sub_selector[self.between] = self.variables
        if self.columns:
            self.sub_selector[self.accross] = self.columns
        self.data = subselectDf(self.data, self.sub_selector)

    def build_matrices(self):
        cases = [
            (
                group_df,
                grouping[0],
                self.between,
                self.var1,
                self.var2,
                "compound" if self.between == "region" else "region",
                self.columns,
                self.n_minimum,
                self.method,
                self.pvalue_threshold,
            
            ) for grouping, group_df in self.data.groupby(by=[self.group_by], sort=False)
        ]  # Setup multiprocessing pool
        self.matrices = parallel_process(Matrix, cases)

    def homogenize_datasets(self):
        conserved_rows = set.intersection(
            *(set(matrix.corr.index) for matrix in self.matrices)
        )
        conserved_cols = set.intersection(
            *(set(matrix.corr.columns) for matrix in self.matrices)
        )

        for matrix in self.matrices:
            rows_to_drop = [
                row for row in matrix.corr.index if row not in conserved_rows
            ]
            cols_to_drop = [
                col for col in matrix.corr.columns if col not in conserved_cols
            ]
            matrix.corr = matrix.corr.drop(index=rows_to_drop, columns=cols_to_drop)
