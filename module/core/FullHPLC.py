from module.core.HPLC import HPLC
from module.core.Outliers import Outliers
from module.core.Metadata import ExperimentInformation, TreatmentInformation
from dataclasses import dataclass


@dataclass
class FullHPLC:

    project: str

    def select(self, **selector):
        return self.df.select(**selector)

    @property
    def df(self):
        hplc = HPLC(self.project)
        outliers = Outliers(self.project)
        treatment_information = TreatmentInformation(self.project)
        common_columns = hplc.df.columns.intersection(outliers.df.columns).tolist()
        full_df = hplc.df.merge(outliers.df, on=common_columns, how="left").merge(
            treatment_information.df, on="group_id", how="left"
        )
        return full_df


@dataclass
class ExperimentFullHPLC(FullHPLC):

    experiment: str

    @property
    def df(self):
        single_experiment_info = ExperimentInformation(self.project).get_experiment(
            self.experiment
        )
        experiment_full_hplc = FullHPLC(self.project).select(
            group_id=single_experiment_info["groups"]
        )
        experiment_full_hplc.loc[:, "experiment"] = self.experiment
        experiment_full_hplc[single_experiment_info["independant_variables"]] = list(
            experiment_full_hplc.independant_variables.apply(
                lambda group_independant_variables: [
                    experiment_independant_variable in group_independant_variables
                    for experiment_independant_variable in single_experiment_info[
                        "independant_variables"
                    ]
                ],
            )
        )
        return experiment_full_hplc
