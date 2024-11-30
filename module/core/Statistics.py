from dataclasses import dataclass, field
from typing import ClassVar
import pandas as pd
import numpy as np
from module.core.Dataset import PickleDataset, SelectableDataFrame

# from module.core.HPLC import HPLC
from module.core.Metadata import (
    ExperimentInformation,
    ProjectInformation,
)
from tqdm import tqdm
import scipy
import pingouin as pg
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from module.core.utils import parallel_process
from IPython.display import HTML
from module.core.utils import is_array_like
import statsmodels.api as sm
from statsmodels.formula.api import ols


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
    return {
        (False, False, False, True): ["ttest"],
        (False, False, True, True): ["paired_ttest"],
        (False, True, False, True): ["one_way_anova", "tukey"],
        (False, True, True, True): ["repeated_measures_anova", "paired_ttest"],
        (True, True, False, True): ["two_way_anova", "one_way_anova", "tukey"],
    }[(multiple_factors, multiple_groups, paired, parametric)]


@dataclass
class QuantitativeStatistic:
    """
    Handles statistical analysis for quantitative data.
    Needs all the required data to perform the statistical analysis as
    well as parameters that define the statistics to be performed.
    TODO: This class is a messy and needs a lot of refactoring + finding a smarter way to do it.

    Raises:
        ValueError: _description_
        ValueError: _description_

    Returns:
        _type_: _description_
    """

    data: pd.DataFrame
    independant_variables: list[str]
    is_paired: bool
    is_parametric: bool
    p_value_threshold: float
    group_column: str = field(default="treatment", kw_only=True)
    delay_execution: bool = field(default=False, kw_only=True)
    metadata: dict = field(default_factory=dict, kw_only=True)
    pipeline: list = field(default=None, kw_only=True)

    def __post_init__(self):
        if self.delay_execution:
            self.delay_execution = False
        else:
            self.filtered_data = self.data.select(value="notna")
            # check enough data
            self.has_enough_data = all(
                [
                    group_data.value.count() >= 5
                    for _, group_data in self.filtered_data.groupby(self.group_column)
                ]
            )
            self.pipeline = self.pipeline or get_quantitative_statistics_pipeline(
                len(self.independant_variables) >= 2,
                len(self.data[self.group_column].unique()) >= 2,
                self.is_paired,
                self.is_parametric,
            )
            self.statistical_test = self.pipeline[0]
            self.post_hoc_test = self.pipeline[-1]
            self.filtered_data = self.data.select(value="notna")
            if self.has_enough_data:
                self.results = SelectableDataFrame(self.execute_stats_pipeline())
                self.significant_pairs = (
                    self.results.select(test=self.post_hoc_test).iloc[0, :].p_value
                    if self.post_hoc_test in self.results.test.to_list()
                    else None
                )
            else:
                self.results = SelectableDataFrame(
                    [
                        {
                            "test": test,
                            "is_significant": False,
                            "result": "Not enough data",
                            "result_string": "n/a",
                            "p_value_threshold": self.p_value_threshold,
                            "p_value": np.nan,
                        }
                        for test in self.pipeline + ["validation"]
                    ]
                )
            for key, val in self.metadata.items():
                self.results[key] = val
            self.is_significant = self.results.is_significant.all()

    @staticmethod
    def calculate_from_selection(
        data,
        experiments,
        p_value_threshold: float = None,
        pipeline: list = None,
    ):
        """
        Calculate statistical for data, autmaticcaly groups by experiment, compound regions.

        Args:
            data (pd.DataFrame): The name of the project.
            experiemnts (str): The experiments to calculate for.
            p_value_threshold (float, optional): The p-value threshold used for statistical analysis. Defaults to None.

        Returns:
            SelectableDataFrame: Containing the statistical results.

        """

        groupings = []

        experiments = (
            [experiments] if isinstance(experiments, pd.Series) else experiments
        )

        for experiment in experiments:
            groupings.extend(
                [
                    QuantitativeStatistic(
                        data=group_data,
                        independant_variables=experiment.independant_variables,
                        is_paired=experiment.paired,
                        is_parametric=experiment.parametric,
                        p_value_threshold=p_value_threshold,
                        pipeline=pipeline,
                        delay_execution=True,
                        metadata={
                            "experiment": experiment.label,
                            "compound": compound,
                            "region": region,
                        },
                    )
                    for (region, compound), group_data in tqdm(
                        data.select(group_id=experiment.groups).groupby(
                            ["region", "compound"]
                        ),
                        desc=f"Preparing statistical groupings for {experiment.label}",
                    )
                ]
            )

        statistics = parallel_process(groupings, description="Calculating statistics")

        results = []
        for statistic in statistics:
            result = statistic.results
            result["fully_significant"] = statistic.is_significant
            results.append(result)

        return statistics, SelectableDataFrame(pd.concat(results))

    # for parallel
    def __call__(self):
        self.__post_init__()
        return self

    def execute_stats_pipeline(self):
        results = []
        for test in self.pipeline:
            test_method_name = f"get_{test}"
            if not hasattr(self, test_method_name):
                raise ValueError(f"Unknown test: {test}")
            test_method = self.__getattribute__(test_method_name)
            test_result = test_method()
            results.append(
                {
                    **test_result,
                    "test": test,
                    "p_value_threshold": self.p_value_threshold,
                }
            )
            self.__setattr__(test, test_result["result"])
        return results

    def get_tukey(self):
        """
        Performs Tukey's HSD (Honest Significant Difference) test to identify significant differences between groups.

        Args:
            data (pd.DataFrame): The DataFrame containing the data to analyze, which must include a 'value' column and a 'treatment' column.
            p_value_threshold (float): The significance level to determine if the results are statistically significant.

        Returns:
            tuple: A tuple containing:
                - A list of tuples, where each tuple contains the pairs of treatments that have a significant difference.
                - A boolean indicating if any significant results were found.
                - A DataFrame of the complete results from the Tukey HSD test.
        """
        columns, *stats_data = pairwise_tukeyhsd(
            endog=self.filtered_data.value,
            groups=self.filtered_data[self.group_column],
            alpha=self.p_value_threshold,
        )._results_table.data
        results = pd.DataFrame(stats_data, columns=columns)
        significance_infos = pd.DataFrame(
            list(
                results[results.reject].apply(
                    lambda res: [(res.group1, res.group2), res["p-adj"]], axis=1
                )
            ),
            columns=["pairs", "p_values"],
        )
        return {
            "p_value": [
                significance_infos.pairs.tolist(),
                significance_infos.p_values.tolist(),
            ],
            "is_significant": len(significance_infos) > 0,
            "result": results,
            "result_string": significance_infos,
        }

    def get_one_way_anova(self):
        """
        Performs a one-way ANOVA test to determine if there are any statistically significant differences between the means of three or more independent (unrelated) groups.

        Args:
            data (pd.DataFrame): The DataFrame containing the data, where 'value' is the dependent variable and 'treatment' is the independent variable used to define groups.
            p_value_threshold (float): The alpha level used to determine the threshold for significance.

        Returns:
            tuple: A tuple containing:
                - The p-value from the ANOVA test.
                - A boolean indicating if the test result is significant at the given threshold.
                - A DataFrame containing the ANOVA test results (F-value and p-value).
        """
        model = ols(f"value ~ C({self.group_column})", data=self.filtered_data).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        p_value = anova_table["PR(>F)"][0]
        is_significant = p_value <= self.p_value_threshold
        F = anova_table["F"][0]
        df1, df2 = anova_table["df"][0], anova_table["df"][1]

        return {
            "p_value": p_value,
            "is_significant": p_value <= self.p_value_threshold,
            "result": anova_table,
            "result_string": f"F({int(df1)}, {int(df2)}) = {F:.3g}, p = {p_value:.2g} {'*' if is_significant else ''}",
        }

    def get_two_way_anova(self):
        """
        Performs a two-way ANOVA to evaluate the effect of two nominal predictor variables on a continuous outcome variable.

        Args:
            data (pd.DataFrame): The DataFrame containing the experimental data. 'value' is the dependent variable. The independent variables should be specified as columns following the 'experiment' column.

        Returns:
            tuple: A tuple containing:
                - The p-value of the interaction term in the ANOVA table.
                - A boolean indicating if the interaction effect is significant below the given p_value_threshold.
                - A DataFrame with detailed ANOVA results including F-values and p-values for each effect.
        """
        data = self.filtered_data.copy()
        data[self.independant_variables] = data.apply(
            lambda row: [
                variable in row.independant_variables
                for variable in self.independant_variables
            ],
            axis=1,
            result_type="expand",
        )
        results = pg.anova(
            data=data,
            dv="value",
            between=self.independant_variables,
            detailed=True,
        ).round(3)

        p_value = results["p-unc"][2]
        is_significant = isinstance(p_value, float) and p_value < self.p_value_threshold

        return {
            "p_value": p_value,
            "is_significant": is_significant,
            "result": results,
            "result_string": f"F({int(results['DF'][2])}, {int(results['DF'][3])}) = {results['F'][2]:3g}, p = {p_value:2g} {'*' if is_significant else '' }",
        }
