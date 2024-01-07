from dataclasses import dataclass, field
import scipy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from module.utils import subselectDf
from scipy.stats._warnings_errors import ConstantInputWarning
import warnings


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
        experiment (str): the experiment group
        treatment (str): the treatment group
        between (str): variables to correlate (compound or region)
        variables (str): list of variables to correlate from 'between'. If only one: self correlation
        accross (str): the column that will constitute the rows/cols of the matrix
        n_minimum (int): minumum occurnces of overlapping var1 and var2. Default to 5
        columns (list[str]): List to select from the 'across' column. Defaults to None
        method (str): one of 'pearson', 'kendall', 'spearman'. Defaults to "pearson"
        pvalue_threshold: float = 0.05

    Returns:
        self.pivot: result of self.data.pivot()
        self.corr: correlations value or 0.0 if not significant or nan if insuficient_overlapp
        self.variables: string with the variables being correlated
        self.missing: missing rows/columns due to lack of data or pvalue

    """

    data: pd.DataFrame
    experiment: str
    treatment: str
    between: str
    variables: list[str]
    accross: str
    n_minimum: int = 5
    columns: list[str] = None
    method: str = "pearson"
    pvalue_threshold: float = 0.05

    def __post_init__(self):
        self.filtered_data = self.sub_select().replace(0, np.nan)
        # Handle missing data (n goup < n minimum)
        self.missing_values = []
        missing_indices = []
        for col, df in self.filtered_data.groupby(by=[self.between, self.accross]):
            if df.value.notna().sum() < self.n_minimum:
                self.missing_values.append(col)
                missing_indices.extend(df.index)
        self.filtered_data = self.filtered_data.drop(missing_indices)
        if self.missing_values:
            print(f"{self.treatment} missing data for {self.missing_values}, deleted from analysis")
        # Define variables to correlate
        self.variables = (
            self.variables if len(self.variables) == 2 else self.variables * 2
        )
        self.var1 = self.variables[0]
        self.var2 = self.variables[1]
        # Pivot df and order columns
        self.pivot = self.filtered_data.pivot_table(
            values="value",
            index=self.filtered_data["mouse_id"],
            columns=[self.between, self.accross],
        )
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
        # Calculate correlations
        method = correlate(self.method, self.pvalue_threshold)
        self.corr = (
            self.pivot.corr(method=method, min_periods=self.n_minimum + 1)
            .loc[tuple(self.variables)]
        )
        # Delete rows where cells are nan
        self.missing_ovelapp = [((self.var1, index), (self.var2, column)) for (index, column), is_na in self.corr.T.isna().stack().items() if is_na]
        if self.missing_ovelapp:
            print(f"{self.treatment} insuficient overlapp for {self.missing_ovelapp} pairs")
            print("Inspect with self.corr to adjust {columns} and redo analysis")
        # TRIANGLE CORRELOGRAMS remove duplicate data and make diagonal correlations of 1 visible
        if not self.is_square:  
            mask = np.triu(np.ones(self.corr.shape, dtype=bool), k=1)
            self.corr[mask] = np.nan
            np.fill_diagonal(
                self.corr.values, 1
            )
        # Replace non-significant correlations (= 0.0) by nan
        self.corr = self.corr.where(self.corr != 0, np.nan)

    @property
    def is_square(self):
        return len(self.var1) != len(self.var2)

    def get_title(self):
        if self.is_square:
            return f"{'-'.join(self.variables)} in {self.treatment}"
        return f"{self.var1} in {self.treatment}"

    def sub_select(self):
        selector = {
            "treatment": self.treatment,
            "experiment": self.experiment,
            self.between: self.variables,
        }
        if self.columns:
            selector[self.accross] = self.columns
        return subselectDf(self.data, selector)

