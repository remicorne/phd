from dataclasses import dataclass
from typing import ClassVar
import pandas as pd
from module.core.Dataset import Dataset, sub_select
from module.core.Outliers import Outliers
from module.core.HPLC import HPLC
from module.core.Metadata import TreatmentInformation
from tqdm import tqdm
import scipy
import pingouin as pg
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from module.core.utils import parallel_process


def get_tukey(data, p_value_threshold):
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
    # get independant variables
    independant_variables_col_index = list(data.columns).index('experiment')
    independant_variables = list(data.columns[independant_variables_col_index + 1:])

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


def process_quantitative_stats(data_test_pipeline_p_value_threshold):
    """_summary_

    Args:
        experiment_info (dict): experiment info mapping
        data (df): data for signle experiment
        p_value_threshold (float): threshold applied to all tests

    Returns:
        is_significant (bool), significance_infos (list, list): [treatment parings, pvalues]
    """
    data, test_pipeline, p_value_threshold = data_test_pipeline_p_value_threshold
    test_results = []
    region = data.region.iloc[0]
    compound = data.compound.iloc[0]
    data = data[data.value.notna()]

    for test in test_pipeline:
        p_value, is_significant, stats_results = QUANTITATIVE_STAT_METHODS[test](
            data=data,
            p_value_threshold=p_value_threshold,
        )
        ## TODO handle warnings:
        # RuntimeWarning: Degrees of freedom <= 0 for slice
        # RuntimeWarning: invalid value encountered in scalar divide
        # DegenerateDataWarning: all input arrays have length 1.  f_oneway requires that at least one input has length greater than 1.

        test_results.append(
            {
                "region": region,
                "compound": compound,
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
    return {
        (False, False, False, True): ["ttest"],
        (False, False, True, True): ["paired_ttest"],
        (False, True, False, True): ["one_way_anova", "tukey"],
        (False, True, True, True): ["repeated_measures_anova", "paired_ttest"],
        (True, True, False, True): ["two_way_anova", "one_way_anova", "tukey"],
    }[(multiple_factors, multiple_treatments, paired, parametric)]


@dataclass
class Statistics(Dataset):

    treatment_information: TreatmentInformation
    experiments: object  # list(Experiment)
    p_value_threshold: float
    hplc: HPLC
    outliers: Outliers
    _name: ClassVar = "statistics"

    def generate(self):
        results = []

        filtered_outliers = self.outliers.select({"is_outlier": False})
        common_columns = filtered_outliers.columns.intersection(self.hplc.df.columns)
        hplc_without_outliers = filtered_outliers.merge(self.hplc.df, on=common_columns.tolist())

        hplc_per_experiment_list = []
        # Create df with full experiment information
        for experiment in self.experiments.values():
            exepriment_hplc = sub_select(
                hplc_without_outliers, {"group_id": experiment.groups}
            ).merge(self.treatment_information.df, on="group_id")
            exepriment_hplc["experiment"] = experiment.name
            hplc_per_experiment_list.append(exepriment_hplc)
            
        hplc_per_experiment_df = pd.concat(hplc_per_experiment_list)
            
        stats_pipeline_args = []
        for (experiment, compound, region), data in tqdm(hplc_per_experiment_df.groupby(by=["experiment", "compound", "region"]), desc="Preparing experimental groups"):
            
            experiment = self.experiments[experiment]
            
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
            
            if any([treatment_data.value.count() < 5 for _, treatment_data in data.groupby("treatment")]):
                results.append(
                    {
                        "region": region,
                        "compound": compound,
                        "test": "validation",
                        "is_significant": False,
                        "result": "Not enough data",
                        "p_value_threshold": self.p_value_threshold,
                        "p_value": None,
                    }
                )

            else:
                stats_pipeline_args.append([data, experiment_statistics_pipeline, self.p_value_threshold])

        experiment_results = parallel_process(
            process_quantitative_stats,
            stats_pipeline_args,
            description=f"Calculating stats for each group",
        )

        [results.extend(experiment_result) for experiment_result in experiment_results]

        return pd.DataFrame(results)
