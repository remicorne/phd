from dataclasses import dataclass
from typing import ClassVar
import pandas as pd
import numpy as np
from module.core.Dataset import PickleDataset
from module.core.FullHPLC import ExperimentFullHPLC
from module.core.Metadata import ExperimentInformation, ProjectInformation
from tqdm import tqdm
import scipy
import pingouin as pg
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from module.core.utils import parallel_process


def get_tukey(data, p_value_threshold):
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
        endog=data["value"], groups=data["treatment"], alpha=p_value_threshold
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
    return (
        [significance_infos.pairs.tolist(), significance_infos.p_values.tolist()],
        len(significance_infos) > 0,
        results,
    )


def get_one_way_anova(data, p_value_threshold):
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
    F_value, p_value = scipy.stats.f_oneway(
        *[list(group_df["value"]) for _, group_df in data.groupby("treatment")]
    )
    # print(f'oneWAY_ANOVA F_value: {F_value}, p_value: {p_value}')
    return (
        p_value,
        p_value <= p_value_threshold,
        pd.DataFrame([[F_value, p_value]], columns=["F", "p_value"]),
    )


def get_two_way_anova(data, p_value_threshold):
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
    independant_variables_col_index = list(data.columns).index("experiment")
    independant_variables = list(data.columns[independant_variables_col_index + 1 :])

    results = pg.anova(
        data=data,
        dv="value",
        between=independant_variables,
        detailed=True,
    ).round(3)
    return (
        results["p-unc"][2],
        isinstance(results["p-unc"][2], float)
        and results["p-unc"][2] < p_value_threshold,
        results,
    )


QUANTITATIVE_STAT_METHODS = {
    "two_way_anova": get_two_way_anova,
    "one_way_anova": get_one_way_anova,
    "tukey": get_tukey,
}


def execute_quantitative_stats_pipeline(
    data,
    statistics_pipeline,
    p_value_threshold,  # Structured this way for parallel processing
):
    """
    Processes quantitative statistical tests based on the given parameters and dataset.

    Args:
        data__statistics_pipeline__p_value_threshold__experiment_region_compound (tuple): A tuple containing:
            - data (DataFrame): The data for a single experiment with non-NA values.
            - statistics_pipeline (list): A list of test names to be applied.
            - p_value_threshold (float): The p-value threshold for significance determination.
            - experiment_region_compound (dict): Additional experiment info mapping.

    Returns:
        list of dicts: Each dictionary contains results of statistical tests including:
            - the original experiment_region_compound info,
            - test name ('test'),
            - significance flag ('is_significant'),
            - statistical results ('result'),
            - applied p-value threshold ('p_value_threshold'),
            - computed p-value ('p_value').
    """
    test_results = []

    for test in statistics_pipeline:
        p_value, is_significant, stats_results = QUANTITATIVE_STAT_METHODS[test](
            data=data,
            p_value_threshold=p_value_threshold,
        )
        test_results.append(
            {
                "test": test,
                "is_significant": is_significant,
                "result": stats_results,
                "p_value_threshold": p_value_threshold,
                "p_value": p_value,
            }
        )

    return test_results


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
class StatisticalGrouping:

    data: pd.DataFrame
    experiment_infos: object  # Experiment
    compound: str
    region: str
    p_value_threshold: float

    def __post_init__(self):

        self.statistics_pipeline = get_quantitative_statistics_pipeline(
            len(self.experiment_infos["independant_variables"]) >= 2,
            len(self.experiment_infos["groups"]) >= 2,
            self.experiment_infos["paired"],
            self.experiment_infos["parametric"],
        )

        # Create result base
        self.grouping_information = {
            "experiment": self.experiment_infos["experiment"],
            "region": self.region,
            "compound": self.compound,
        }
        
        self.filtered_data = self.data[self.data.value.notna()].select(is_outlier=False)

        # all([]) == True if the iterable is empty
        self.has_enough_data = not self.filtered_data.empty and all(
            [
                treatment_data.value.count() >= 5
                for _, treatment_data in self.filtered_data.groupby("treatment")
            ]
        )
        self.quantitative_stats_parameters = self.filtered_data, self.statistics_pipeline, self.p_value_threshold
        # del self.experiment_infos # Avoid pickling issues with parallel execution caused by complex objects

    def calculate_results(self):
        return (
            execute_quantitative_stats_pipeline(*self.quantitative_stats_parameters)
            if self.has_enough_data
            else [{
                **self.grouping_information,
                "test": "validation",
                "is_significant": False,
                "result": "Not enough data",
                "p_value_threshold": self.p_value_threshold,
                "p_value": np.nan,
            }]
        )
        
    def get_results(self):
        return [{**self.grouping_information, **result} for result in self.calculate_results()]

    
    
def parallel_execution_wrapper(grouping: StatisticalGrouping):
    return grouping.get_results()


@dataclass(repr=False)
class Statistics(PickleDataset):
    """
    Manages the statistical analysis for a set of experiments, including outlier filtering, merging data with HPLC results,
    and generating statistical results based on defined experiments.

    Attributes:
        treatment_information (TreatmentInformation): Information about treatments in the experiments.
        experiments (list(Experiment)): A collection of `Experiment` objects detailing individual experiments.
        p_value_threshold (float): The p-value threshold used to determine the statistical significance.
        hplc (HPLC): An HPLC data object containing the HPLC results.
        outliers (Outliers): An Outliers data object for managing and filtering outlier data.
        _name (ClassVar): Static name attribute for the class, set to "statistics".

    Methods:
        generate: Processes experiments to produce a DataFrame of statistical results.
        significant_results: Returns a DataFrame of results where all tests are significant.
        insufficent_data: Returns a DataFrame of results flagged with "Not enough data".
    """

    project: str
    filename: ClassVar[str] = "statistics"

    def generate(self):
        groupings = []
        experiment_information = ExperimentInformation(self.project)
        project_information = ProjectInformation(self.project)
        for single_experiment_infos in experiment_information.list:
            experiment_df = ExperimentFullHPLC(self.project, single_experiment_infos['experiment'])
            # iterate over every data grouping to build a list of arguments for parallel processing
            for (compound, region), data in tqdm(
                experiment_df.select(experiment=single_experiment_infos['experiment']).groupby(by=["compound", "region"]),
                desc=f"Preparing stat groupings for ",
            ):
                groupings.append(StatisticalGrouping(data, single_experiment_infos, compound, region, project_information.p_value_threshold))

        # Process the statistical tests in parallel
        results = parallel_process(
            parallel_execution_wrapper,
            groupings,
            description=f"Calculating stats for each group",
        )

        # Unpack results and merge
        return pd.DataFrame(grouping_subresult for grouping_results in results for grouping_subresult in grouping_results)


    def get_quantitative_stats(self, experiment, compound, region, p_value_threshold):
        """
        Calculates the statistics for a grouping using the test pipeline
        Saves the results to the fulll statistical df
        Returns the quantitative statistical results for a specific experiment, compound, and region.

        Args:
            experiment (str): The name of the experiment.
            compound (str): The name of the compound.
            region (str): The name of the region.

        Returns:
            pd.DataFrame: A DataFrame containing the statistical results for the specified experiment, compound, and region.
        """
        data = ExperimentFullHPLC(self.project, experiment).select(compound=compound, region=region)
        experiment_infos = ExperimentInformation(self.project).get_experiment(experiment)
        return StatisticalGrouping(data, experiment_infos, compound, region, p_value_threshold).get_results()

    @property
    def significant_results(self):
        """Return groupings where all tests are significant.

        Returns:
            pd.DataFrame: Significant results.
        """
        return pd.concat(
            [
                results
                for _, results in self.df.groupby(["experiment", "compound", "region"])
                if results.is_significant.all()
            ]
        )

    @property
    def significant_tests(self):
        """Return groupings where all tests are significant.

        Returns:
            pd.DataFrame: Significant results.
        """
        return pd.concat(
            [
                results[results.is_significant]
                for _, results in self.df.groupby(["experiment", "compound", "region"])
            ]
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
