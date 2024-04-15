from dataclasses import dataclass, field
import scipy
import pandas as pd
import numpy as np
from module.utils import parallel_process, subselectDf


def calculate_correlation(method, x, y):
    """Calculate the correlation and p-value based on the specified method."""
    if method == "pearson":
        result = scipy.stats.pearsonr(x, y)
        return result.statistic, result.pvalue
    elif method == "spearman":
        result = scipy.stats.spearmanr(x, y)
        return result.correlation, result.pvalue
    elif method == "kendall":
        result = scipy.stats.kendalltau(x, y)
        return result.correlation, result.pvalue
    else:
        raise ValueError(f"Unknown method: {method}")


def correlate(method, return_type):
    """Return a correlation function based on the specified method, p-value threshold, and return type."""

    def executor(x, y):
        correlation, pvalue = calculate_correlation(method, x, y)
        if return_type == "pvalues":
            return pvalue
        elif return_type == "correlations":
            return correlation
        else:
            raise ValueError(f"Unknown return type: {return_type}")

    return executor


@dataclass
class Matrix:

    """
    Creates a reusable matrix class for a given eperimnet. 
    Args:
        data (pd.DataFrame):    The original dataset.
        treatment (str):        Identifier for the treatment group.
        between (str):          Variable type for correlation (e.g., 'compound' or 'region').
        variables (str):         'var1-var2' to correlate from type 'between'. If only one: self correlation.
        accross (str): The column that will constitute the rows/cols of the matrix.
        columns (list[str]): Columns to include in the analysis. If None, all columns are included.
        n_minimum (int): Minumum occurnces of overlapping var1 and var2. Default = 5.
        method (str): Correlation method ('pearson', 'spearman', 'kendall'). Default = "pearson".
        pvalue_threshold (float): Threshold for significance in correlation. Defult = 0.05

    Returns:
        filtered_data (pd.DataFrame): Subselected data filtered based on n_minimum between vairables.
        pivot (pd.DataFrame): Pivot table of the filtered data.
        corr_masked (pd.DataFrame): Masked correlation matrix based on p-value threshold.
        correlations (pd.DataFrame): Full correlation matrix.
        pvalues (pd.DataFrame): Matrix of p-values for the correlations.
        missing_values (list): List of variables with missing data (< n_minimum).
        missing_overlap (list): List of variable pairs with insufficient data overlap.

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
    corr_masked: pd.DataFrame = field(init=False)
    correlations: pd.DataFrame = field(init=False)
    pvalues: pd.DataFrame = field(init=False)
    missing_values: list = field(init=False)
    missing_overlap: list = field(init=False)

    def __post_init__(self):
        self.filter_missing_values()
        self.pivot_data()
        self.order_columns()
        self.correlate()
        self.find_missing_overlap()
        self.process_triangle_correlogram()

    def filter_missing_values(self):
        """
        Filters out variables with occurrences less than n_minimum and updates missing_values list.
        """
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
        """
        Creates a pivot table from the filtered data.
        """
        self.pivot = self.filtered_data.pivot_table(
            values="value",
            index=self.filtered_data["mouse_id"],
            columns=[self.between, self.accross],
        )

    def order_columns(self):
        """
        Orders the columns of the pivot table based on the provided column list.
        """
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
        
    def correlate(self):
        """
        Calculates and stores correlation and p-value matrices.
        """
        self.pvalues = self.create_corr_matrix('pvalues')
        self.correlations = self.create_corr_matrix('correlations')
        self.corr_masked = self.correlations[self.pvalues < self.pvalue_threshold]
        

    def create_corr_matrix(self, result_type):
        """
        Creates a correlation matrix for either correlation values or p-values.

        Args:
            result_type (str): Type of result to return ('pvalues' or 'correlations').

        Returns:
            pd.DataFrame: A DataFrame containing the requested correlation matrix.
        """
        method = correlate(self.method, result_type)
        return self.pivot.corr(method=method, min_periods=self.n_minimum + 1).loc[
            tuple([self.var1, self.var2])
        ]

    def find_missing_overlap(self):
        """
        Identifies and reports variable pairs with insufficient data overlap.
        """
        self.missing_overlap = [
            ((self.var1, index), (self.var2, column))
            for (index, column), is_na in self.correlations.T.isna().stack().items()
            if is_na
        ]
        if self.missing_overlap:
            print(
                f"{self.grouping} insuficient overlapp for {self.missing_overlap} pairs"
            )
            print("Inspect with self.corr to adjust {columns} and redo analysis")

    def process_triangle_correlogram(self):
        """
        Masks the upper triangle of the correlogram if the matrix is not square.
        """
        if not self.is_square:
            mask = np.triu(np.ones(self.corr_masked.shape, dtype=bool), k=1)
            self.corr_masked[mask] = np.nan
            np.fill_diagonal(self.corr_masked.values, 1)

    @property
    def is_square(self):
        """
        Checks if the matrix is square (var1 and var2 are different).

        Returns:
            bool: True if the matrix is square, False otherwise.
        """
        return self.var1 != self.var2

    def get_title(self):
        """
        Generates a title for the correlogram based on the matrix configuration.

        Returns:
            str: A title string.
        """
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
        group_by (str): The column name in 'data' to group by (generally 'experiment').
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

    ## TODO: delete when new model done, should have only pertinent data
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
            *(set(matrix.corr_masked.index) for matrix in self.matrices)
        )
        conserved_cols = set.intersection(
            *(set(matrix.corr_masked.columns) for matrix in self.matrices)
        )

        for matrix in self.matrices:
            rows_to_drop = [
                row for row in matrix.corr_masked.index if row not in conserved_rows
            ]
            cols_to_drop = [
                col for col in matrix.corr_masked.columns if col not in conserved_cols
            ]
            matrix.corr_masked = matrix.corr_masked.drop(index=rows_to_drop, columns=cols_to_drop)
