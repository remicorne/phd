from module.core.Dataset import SelectableDataFrame
from module.core.HPLC import HPLC
from module.core.Outliers import Outliers
from module.core.Metadata import ExperimentInformation, TreatmentInformation
from dataclasses import dataclass, field


@dataclass
class FullHPLC:

    project: str
    experiment: str = field(default=None)

    def select(self, **selector) -> SelectableDataFrame:
        return self.df.select(**selector)

    @property
    def df(self) -> SelectableDataFrame:
        df = HPLC(self.project).full_df.extend(Outliers(self.project))
        if self.experiment:
            
            single_experiment_info = ExperimentInformation(self.project).get_experiment(
                self.experiment
            )
            df = df.select(
                group_id=single_experiment_info["groups"]
            )
            df.loc[:, "experiment"] = self.experiment
            df[single_experiment_info["independant_variables"]] = list(
                df.independant_variables.apply(
                    lambda group_independant_variables: [
                        experiment_independant_variable in group_independant_variables
                        for experiment_independant_variable in single_experiment_info[
                            "independant_variables"
                        ]
                    ],
                )
            )
        return df
        
