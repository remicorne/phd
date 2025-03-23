import pandas as pd
import numpy as np
from scipy import stats
from functools import cached_property
from typing import List, Dict, Optional
from module.core.statistics.outliers import get_outliers
from module.core.statistics.quantitative import (
    get_quantitative_statistics_pipeline,
    TestPipeline,
)


class DataAnalyzer:

    def __init__(
        self,
        df: pd.DataFrame,
        subject_column: str,
        group_column: str,
        measurement_columns: List[str],
        independent_variables: List[str],
        value_column: str = "value",
        parametric: bool = None,
    ):
        """
        Initialize the DataAnalyzer with required parameters

        Args:
            df: Input DataFrame
            subject_column: Column name identifying subjects
            group_column: Column name for grouping similar contexts
            measurement_columns: Columns identifying different metrics
            independent_variables: Dictionary of {variable_name: column_name}
                                   for boolean flags
            value_column: Column containing measurement values. Defaults to value
        """
        # Validate columns
        required_columns = (
            [subject_column, group_column, value_column]
            + measurement_columns
            + independent_variables
        )

        if not set(required_columns).issubset(df.columns):
            missing = set(required_columns) - set(df.columns)
            raise ValueError(f"Missing required columns: {missing}")

        self.data = df
        self.working_df = self.data
        self.subject_column = subject_column
        self.group_column = group_column
        self.measurement_columns = measurement_columns
        self.independent_variables = independent_variables
        self.value_column = value_column
        self.outliers_df = pd.DataFrame()
        self.has_removed_outliers = False
        self.parametric = parametric

    @property
    def aggregate_stats(self) -> pd.DataFrame:
        """Calculate and cache aggregated statistics"""
        group_cols = [self.group_column] + self.measurement_columns
        groupby = self.working_df.groupby(group_cols)[self.value_column]
        agg_df = groupby.agg(
            mean=np.mean,
            std=np.std,
            sem=lambda x: np.std(x, ddof=1) / np.sqrt(len(x)),
            count="count",
            value=list,
        ).reset_index()

        agg_df[["is_normal", "F", "p"]] = agg_df.apply(
            lambda row: self._shapiro_test(row[self.value_column]),
            axis=1,
            result_type="expand",
        )

        return agg_df.merge(groupby.describe().reset_index())

    def remove_outliers(self, method: str = "grubbs", threshold: float = 3):
        """
        Remove outliers using specified method

        """
        group_cols = [self.group_column] + self.measurement_columns
        outliers_df = (
            self.data.groupby(group_cols)[self.value_column]
            .agg(lambda x: get_outliers(x, method, threshold))
            .reset_index()
            .rename(columns={"value": "outliers"})
        )
        self.data = self.data.merge(outliers_df, how="left")
        self.data["is_outlier"] = self.data.apply(
            lambda row: row["value"] in row["outliers"], axis=1
        )
        self.data = self.data.drop(columns=["outliers"])

        self.outliers_df = self.data[self.data.is_outlier]
        self.working_df = self.data[~self.data.is_outlier]
        return self

    def get_quantitative_stats(
        self, pipeline: Optional[List[str]] = None, p_value_threshold: float = 0.05
    ) -> pd.DataFrame:
        """
        Calculate quantitative statistics with optional processing pipeline

        Args:
            pipeline: List of processing steps to apply
            p_value_threshold: Significance threshold for statistical tests

        Returns:
            DataFrame of statistical results
        """
        if pipeline is None:
            pipeline = self.infer_quantitative_statistics_pipeline()

        test_pipeline = TestPipeline(
            pipeline,
            {
                **self.__dict__,
                "p_value_threshold": p_value_threshold,
            },
        )
        return (
            self.working_df.groupby(self.measurement_columns)
            .apply(test_pipeline.run_pipeline)
            .reset_index()
        )

    def infer_quantitative_statistics_pipeline(self):
        multiple_factors = len(self.independent_variables) >= 2
        multiple_groups = len(self.working_df[self.group_column].unique()) >= 2
        paired = len(self.working_df[self.group_column].unique()) == 1
        parametric = (
            self.aggregate_stats.is_normal.all()
            if self.parametric is None
            else self.parametric
        )

        if (
            paired
        ):  # TODO: think about pairing implemntation (add timestamp? dataset id? or timepoint col in measurment_cols)
            raise NotImplementedError("Paired design not supported")

        return get_quantitative_statistics_pipeline(
            multiple_factors,
            multiple_groups,
            paired,
            parametric,
        )

    @staticmethod
    def _shapiro_test(x: pd.Series) -> bool:
        """Helper method for normality test"""
        if len(x) < 3:
            return np.nan
        F, p = stats.shapiro(x)
        return p > 0.05, F, p

    @property
    def outliers(self) -> pd.DataFrame:
        """Get DataFrame of outliers"""
        return self.outliers_df
