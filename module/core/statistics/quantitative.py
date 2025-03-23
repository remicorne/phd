from abc import ABC, ABCMeta, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Type, Any, List, Optional, Union, Tuple
import inspect
from dataclasses import dataclass
from functools import cached_property
import pandas as pd
import numpy as np
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import pingouin as pg
from module.core.utils import parallel_process

TEST_DECISION_TREE = {
    (False, False, False, True): ["ttest"],
    (False, False, True, True): ["paired_ttest"],
    (False, True, False, True): ["one_way_anova", "tukey_hsd"],
    (False, True, True, True): ["repeated_measures_anova", "paired_ttest"],
    (True, True, False, True): ["two_way_anova", "one_way_anova", "tukey_hsd"],
}


def get_quantitative_statistics_pipeline(
    multiple_factors, multiple_groups, paired, parametric
):
    """
    Determines the statistical tests pipeline based on experiment design parameters.

    Args:
        multiple_factors (bool): Flag indicating if multiple factors are involved.
        multiple_groups (bool): Flag indicating if multiple groups are considered.
        paired (bool): Flag indicating if the design is paired.
        parametric (bool): Flag indicating if the tests should be parametric.

    Returns:
        list of str: A list of names representing the statistical tests to be applied based on the input parameters.
    """
    key = (multiple_factors, multiple_groups, paired, parametric)
    if key not in TEST_DECISION_TREE:
        raise NotImplementedError(
            f"Tests pipeline not implemented for {dict(zip(['multiple_factors', 'multiple_groups', 'paired', 'parametric'], key))}"
        )
    return TEST_DECISION_TREE[key]


class StatisticalTestMeta(ABCMeta):
    """Metaclass handling automatic test registration"""

    _registry: Dict[str, Type["AbstractStatisticalTest"]] = {}

    def __new__(cls, name, bases, attrs):
        new_class = super().__new__(cls, name, bases, attrs)
        if ABC not in bases:
            cls._registry[attrs["test_name"]] = new_class
        return new_class

    @classmethod
    def get(cls, test_name) -> Type["AbstractStatisticalTest"]:
        if test_name not in cls.list():
            raise ValueError(f"Test {test_name} not registered")
        return cls._registry[test_name]

    @classmethod
    def list(mcs):
        return mcs._registry.keys()


class BaseAbstractStatisticalTest(ABC, metaclass=StatisticalTestMeta):
    """Abstract base class for all statistical tests, defines common interface
    Used to standardize test execution and result formatting with result_table and result_string
    """

    test_name: str

    def __init__(
        self,
        group_column: str,
        value_column: str,
        p_value_threshold: float,
    ):
        self.group_column = group_column
        self.value_column = value_column
        self.p_value_threshold = p_value_threshold

    def get_results(self, data: pd.DataFrame) -> pd.Series:
        result_string, result_table = self.prepare_results(data)
        result_table["result_string"] = result_string
        return result_table

    @abstractmethod
    def prepare_results(self, data: pd.DataFrame) -> Tuple[Union[str, pd.Series]]:
        """Runs the test and handles abstract subclass specificities in result formatting"""

    @abstractmethod
    def run_test(
        self, data: pd.DataFrame
    ) -> Tuple[Union[str, pd.Series], float, int, int, float]:
        """Main method to execute the statistical test and return results"""

    def result_table(self, result, p_value, is_significant) -> pd.Series:
        """Standardized result formatting"""
        return pd.Series(
            {
                "test": self.test_name,
                "result": result,
                "p_value": p_value,
                "is_significant": is_significant,
            }
        )

    @classmethod
    def get_parameters(cls) -> List[str]:
        return inspect.signature(cls.__init__).parameters


class AbstractStatisticalGlobalTest(BaseAbstractStatisticalTest, ABC):
    """Contains specificities for all global statistical tests"""

    def prepare_results(self, data: pd.DataFrame) -> Tuple[Union[str, pd.Series]]:
        result, p_value, df1, df2, F = self.run_test(data)
        is_significant = p_value <= self.p_value_threshold
        return self.result_string(
            p_value, is_significant, df1, df2, F
        ), self.result_table(result, p_value, is_significant)

    def result_string(self, p_value, is_significant, df1, df2, F):
        return (
            f"F({df1}, {df2}) = {F:.3g}, p = {p_value:.2g} {'*' if is_significant else ''}",
        )


class OneWayANOVATest(AbstractStatisticalGlobalTest):
    test_name = "one_way_anova"

    def run_test(
        self, data: pd.DataFrame
    ) -> Tuple[Union[str, pd.Series], float, int, int, float]:
        model = ols(
            f"{self.value_column} ~ C({self.group_column})",
            data=data,
        ).fit()
        anova_table = anova_lm(model, typ=2)
        p_value = anova_table["PR(>F)"][0]
        return (
            anova_table,
            p_value,
            int(anova_table["df"][0]),
            int(anova_table["df"][1]),
            anova_table["F"][0],
        )


class TwoWayANOVATest(AbstractStatisticalGlobalTest):
    test_name = "two_way_anova"

    def __init__(
        self,
        group_column,
        value_column,
        p_value_threshold,
        independent_variables,
    ):
        super().__init__(group_column, value_column, p_value_threshold)
        self.independent_variables = independent_variables

    def run_test(
        self, data: pd.DataFrame
    ) -> Tuple[Union[str, pd.Series], float, int, int, float]:
        result = pg.anova(
            data=data,
            dv="value",
            between=self.independent_variables,
            detailed=True,
        ).round(3)

        return (
            result,
            result["p-unc"][2],
            int(result["DF"][2]),
            int(result["DF"][3]),
            result["F"][2],
        )
        # is_significant = isinstance(p_value, float) and p_value < self.p_value_threshold # Not sure if still useful


class AbstractStatisticalPostHocTest(BaseAbstractStatisticalTest, ABC):
    """Contains specificities for all post-hoc statistical tests"""

    def prepare_results(self, data: pd.DataFrame) -> Tuple[Union[str, pd.Series]]:
        result, significant_resuts = self.run_test(data)
        is_significant = bool(len(significant_resuts))
        return self.result_string(significant_resuts), self.result_table(
            result, significant_resuts, is_significant
        )

    def result_string(self, significant_pairs: pd.DataFrame):
        return "\n".join(
            significant_pairs.apply(
                lambda row: f"({row.group1}, {row.group2}) p = {row['p-adj']:.2g}",
                axis=1,
            )
        )


class TukeyHSDTest(AbstractStatisticalPostHocTest):
    test_name = "tukey_hsd"

    def run_test(
        self, data: pd.DataFrame
    ) -> Tuple[Union[str, pd.Series], float, int, int, float]:
        tukey_result = pairwise_tukeyhsd(
            endog=data[self.value_column],
            groups=data[self.group_column],
            alpha=self.p_value_threshold,
        )
        results_df = pd.DataFrame(
            tukey_result._results_table.data[1:],
            columns=tukey_result._results_table.data[0],
        )

        significant_results = results_df[results_df["reject"]][
            ["group1", "group2", "p-adj"]
        ]
        return results_df, significant_results


class TestPipeline:

    def __init__(self, tests: List[str], context: dict):
        self.tests = []
        for test in tests:
            Test: Type[BaseAbstractStatisticalTest] = StatisticalTestMeta.get(test)
            parameters = {
                key: value
                for key, value in context.items()
                if key in Test.get_parameters()
            }
            self.tests.append(Test(**parameters))

    def run_pipeline(self, data: pd.DataFrame):
        return pd.concat(
            (test.get_results(data) for test in self.tests), copy=False, axis=1
        ).T
