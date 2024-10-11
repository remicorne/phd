from dataclasses import dataclass, field
from typing import ClassVar
import pandas as pd
import numpy as np
from module.core.Dataset import PickleDataset, SelectableDataFrame
from module.core.HPLC import HPLC
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
    multiple_factors, multiple_treatments, paired, parametric
):
    """
    Determines the statistical tests pipeline based on experiment design parameters.

    Args:
        multiple_factors (bool): Flag indicating if multiple factors are involved.
        multiple_treatments (bool): Flag indicating if multiple treatments are considered.
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
    }[(multiple_factors, multiple_treatments, paired, parametric)]


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
    treatments: list[str]
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
            # check all treatments present and enough data
            self.has_enough_data = set(
                self.filtered_data[self.group_column].unique()
            ) == set(self.treatments) and all(
                [
                    group_data.value.count() >= 5
                    for _, group_data in self.filtered_data.groupby(self.group_column)
                ]
            )
            self.pipeline = self.pipeline or get_quantitative_statistics_pipeline(
                len(self.independant_variables) >= 2,
                len(self.treatments) >= 2,
                self.is_paired,
                self.is_parametric,
            )
            self.statistical_test = self.pipeline[0]
            self.post_hoc_test = self.pipeline[-1]
            self.filtered_data = self.data.select(value="notna")
            if self.has_enough_data:
                self.results = SelectableDataFrame(self.execute_stats_pipeline())
                self.significant_pairs = (
                    self.results.select(test=self.post_hoc_test).p_value
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
                        } for test in self.pipeline + ["validation"]
                    ]
                )
            for key, val in self.metadata.items():
                self.results[key] = val
            self.is_significant = self.results.is_significant.all()

    @classmethod
    def calculate(
        cls,
        project: str,
        experiment: str = None,
        compound: str = None,
        region: str = None,
        p_value_threshold: float = None,
        remove_outliers="calculated",
    ):
        """
        Calculate statistical results for a given experiment, compound, and region.
        If one ore some of the parameters are None or lists, calculate the results for all possible combinations.

        Args:
            project (str): The name of the project.
            experiment (str|list): The name of the experiment.
            compound (str|list): The name of the compound.
            region (str|list): The name of the region.
            p_value_threshold (float, optional): The p-value threshold used for statistical analysis. Defaults to None.
            remove_outliers (str, optional): Whether to remove outliers. Must be 'eliminated', 'calculated', or False. Defaults to "calculated".

        Returns:
            SelectableDataFrame: Containing the statistical results.

        Raises:
            ValueError: If remove_outliers is not 'eliminated', 'calculated', or False.
        """

        data = HPLC(project).full_df

        if compound:
            compounds = (
                compound
                if is_array_like(compound)
                else compound.replace(" ", "").split(",")
            )
        else:
            compounds = data.compound.unique()

        if region:
            regions = (
                region if is_array_like(region) else region.replace(" ", "").split(",")
            )
        else:
            regions = data.region.unique()

        if experiment:
            experiment = (
                experiment
                if is_array_like(experiment)
                else experiment.replace(", ", "").split(",")
            )
            experiments = [
                ExperimentInformation(project).select(experiment=experiment)
                for experiment in experiment
            ]
        else:
            experiments = [experiment for experiment in ExperimentInformation(project)]

        p_value_threshold = (
            p_value_threshold or ProjectInformation(project).p_value_threshold
        )

        if remove_outliers == "eliminated":
            data = data.select(outlier_status=["normal", "kept"])
        elif remove_outliers == "calculated":
            data = data.select(is_outlier=lambda val: val != True)
        elif remove_outliers is not False:
            raise ValueError(
                "remove_outliers must be 'eliminated', 'calculated', or False"
            )

        data = data.select(
            region=regions,
            compound=compounds,
        )

        groupings = []

        for experiment in experiments:
            groupings.extend(
                [
                    cls(
                        data=group_data,
                        independant_variables=experiment.independant_variables,
                        treatments=experiment.treatments,
                        is_paired=experiment.paired,
                        is_parametric=experiment.parametric,
                        p_value_threshold=p_value_threshold,
                        delay_execution=True,
                        metadata={
                            "project": project,
                            "experiment": experiment.label,
                            "compound": compound,
                            "region": region,
                        },
                    )
                    for (region, compound), group_data in tqdm(
                        data.select(treatment=experiment.treatments).groupby(
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

        return SelectableDataFrame(pd.concat(results))
    
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
        
        experiments = [experiments] if isinstance(experiments, pd.Series) else experiments

        for experiment in experiments:
            groupings.extend(
                [
                    QuantitativeStatistic(
                        data=group_data,
                        independant_variables=experiment.independant_variables,
                        treatments=experiment.treatments,
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
                        data.select(treatment=experiment.treatments).groupby(
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
        model = ols(f'value ~ C({self.group_column})', data=self.filtered_data).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        p_value = anova_table["PR(>F)"][0]
        F = anova_table['F'][0]
        df1, df2 = anova_table['df'][0], anova_table['df'][1]
    
        return {
            "p_value": p_value,
            "is_significant": p_value <= self.p_value_threshold,
            "result": anova_table,
            "result_string": f"F({df1}, {df2}) = {F}",
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
        
        
        return {
            "p_value": results["p-unc"][2],
            "is_significant": isinstance(results["p-unc"][2], float)
            and results["p-unc"][2] < self.p_value_threshold,
            "result": results,
            "result_string": f"F({results['DF'][2]}, {results['DF'][3]}) = {results['F'][2]}"
        }


@dataclass(repr=False)
class Statistics(PickleDataset):
    """
    Contains all quantitative statistical resuts for hplc

    """

    project: str
    filename: ClassVar[str] = "statistics"

    def select(self, **kwargs):
        """
        Enable selection of specific significance pairs
        Args:
            significant_pair: tuple of significant pairs

        """
        searched_pair = kwargs.pop("significant_pair", None)
        data = self.df.select(**kwargs)
        if searched_pair:
            return data[
                data.p_value.apply(
                    lambda p_value: (
                        set(searched_pair)
                        in [set(significant_pair) for significant_pair in p_value[0]]
                        if is_array_like(p_value)
                        else False
                    )
                )
            ]
        return data

    def generate(self):
        return QuantitativeStatistic.calculate(self.project).select(
            fully_significant=True
        )

    @property
    def significant_results(self):
        """Return groupings where all tests are significant.

        Returns:
            pd.DataFrame: Significant results.
        """
        return SelectableDataFrame(
            pd.concat(
                [
                    results
                    for _, results in self.df.groupby(
                        ["experiment", "compound", "region"]
                    )
                    if results.is_significant.all()
                ]
            )
        )

    @property
    def insufficent_data(self):
        """Return groupings where one group has less than 5 data points.

        Returns:
            pd.DataFrame: Insufficient data results.
        """
        try:
            invalid_groupings = self.select(test="validation")
            return invalid_groupings[invalid_groupings.result == "Not enough data"]
        except ValueError as e:
            if "EMPTY SELECTION" in str(e):
                return "No data"

        return self.df[
            (self.df.test == "validation" & self.df.result == "Not enough data")
        ]


@dataclass
class AggregateStatistics(PickleDataset):
    """
    Calculates group level descriptive statistics.
    """

    project: str
    filename: ClassVar[str] = "aggregate_statistics"

    def generate(self):
        result_ls = []
        for (treatment, region, compound), groupby_df in tqdm(
            HPLC(self.project).full_df.select(value="notna", is_outlier=False).groupby(by=["treatment", "region", "compound"]),
            desc="Calculating aggregate statistics",
        ):
            
            if len(groupby_df) >= 3:
                F, p, = scipy.stats.shapiro(groupby_df["value"])
                is_parametric = p > 0.05            
            else:
                F, p, is_parametric = np.nan, np.nan, np.nan
            
            mean, std, sem, values = [
                groupby_df.value.mean(),
                groupby_df.value.std(),
                groupby_df.value.sem(),
                groupby_df.value.values,
            ]
            result_ls.append(
                [treatment, region, compound, F, p, is_parametric, mean, std, sem, values]
            )
        return pd.DataFrame(
            result_ls,
            columns=[
                "treatment",
                "region",
                "compound",
                "shapiro_F",
                "shapiro_p",
                "is_parametric",
                "mean",
                "std",
                "sem",
                "values",
            ],
        )
        
        
