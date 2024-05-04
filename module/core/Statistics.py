from dataclasses import dataclass
from typing import ClassVar
import pandas as pd
import numpy as np
from module.core.Dataset import Dataset
from module.core.HPLC import HPLC
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


def process_quantitative_stats(
    data__test_pipeline__p_value_threshold,  # Structured this way for parallel processing
):
    """
    Processes quantitative statistical tests based on the given parameters and dataset.

    Args:
        data__test_pipeline__p_value_threshold__experiment_region_compound (tuple): A tuple containing:
            - data (DataFrame): The data for a single experiment with non-NA values.
            - test_pipeline (list): A list of test names to be applied.
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
    data, test_pipeline, p_value_threshold = data__test_pipeline__p_value_threshold
    test_results = []
    data = data[data.value.notna()]

    for test in test_pipeline:
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
class Statistics(Dataset):
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

    full_df: pd.DataFrame
    experiments: object  # list(Experiment)
    p_value_threshold: float
    _name: ClassVar = "statistics"

    def generate(self):
        results = []

        # Create df with full experiment information, duplicating shared groups
        hplc_per_experiment_list = []
        for experiment in self.experiments.values():
            experiment_hplc = self.full_df.select(
                {"is_outlier": False, "group_id": experiment.groups}
            )
            experiment_hplc.loc[:, "experiment"] = experiment.name
            hplc_per_experiment_list.append(experiment_hplc)

        hplc_per_experiment_df = pd.concat(hplc_per_experiment_list)

        # iterate over every data grouping to build a list of arguments for parallel processing
        stats_pipeline_args = []
        grouping_infos = []
        for (experiment_name, compound, region), data in tqdm(
            hplc_per_experiment_df.groupby(by=["experiment", "compound", "region"]),
            desc="Preparing experimental groups",
        ):

            # Get experiment object
            experiment = self.experiments[experiment_name]

            experiment_statistics_pipeline = get_quantitative_statistics_pipeline(
                len(experiment.independant_variables) >= 2,
                len(experiment.groups) >= 2,
                experiment.paired,
                experiment.parametric,
            )

            # Add independant variables as boolean columns
            data[experiment.independant_variables] = list(
                data.independant_variables.apply(
                    lambda group_independant_variables: [
                        experiment_variable in group_independant_variables
                        for experiment_variable in experiment.independant_variables
                    ],
                )
            )

            # Create result base
            experiment_compound_region = {
                "experiment": experiment.name,
                "region": region,
                "compound": compound,
            }

            # Groupings where there is not enough data for a group are not treated as stat tests reject them
            if any(
                [
                    treatment_data.value.count() < 5
                    for _, treatment_data in data.groupby("treatment")
                ]
            ):
                results.append(
                    {
                        **experiment_compound_region,
                        "test": "validation",
                        "is_significant": np.nan,
                        "result": "Not enough data",
                        "p_value_threshold": self.p_value_threshold,
                        "p_value": None,
                    }
                )
            # Otherwise, add to the pipeline
            else:
                # grouping info is used to merge the results back together
                grouping_infos.append(experiment_compound_region)
                stats_pipeline_args.append(
                    (
                        data,
                        experiment_statistics_pipeline,
                        self.p_value_threshold,  # Structured this way for parallel processing
                    )
                )

        # Process the statistical tests in parallel
        experiment_results = parallel_process(
            process_quantitative_stats,
            stats_pipeline_args,
            description=f"Calculating stats for each group",
        )

        # Merge results back together with grouping info and add to results
        for experiment_result, experiment_compound_region in zip(
            experiment_results, grouping_infos
        ):
            results.extend(
                [
                    {**experiment_sub_result, **experiment_compound_region}
                    for experiment_sub_result in experiment_result
                ]
            )

        return pd.DataFrame(results)

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
            invalid_groupings = self.select({"test": "validation"})
            return invalid_groupings[invalid_groupings.result == 'Not enough data']
        except ValueError as e:
            if "EMPTY SELECTION" in str(e):
                return "No data"
                
        return self.df[(self.df.test == "validation" & self.df.result == "Not enough data")]
