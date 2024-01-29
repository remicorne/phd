import pandas as pd
from module.core.Dataset import Dataset
from dataclasses import dataclass
from typing import ClassVar


class ProjectData(Dataset):
    _type: ClassVar[str]
    _columns: ClassVar[list]

    def generate_data(self):
        return pd.DataFrame(self._temlate)

    def initialize(self):
        super().initialize()
        self.edit_excel()


class TreatmentMapping(ProjectData):
    _type: ClassVar[str] = "treatment_mapping"
    _template: ClassVar[list] = {
        "group_id": [1],
        "treatment": ["vehicles"],
        "color": ["blue"],
    }


class ExperimentInfo(ProjectData):
    _type: ClassVar[str] = "experimental_information"
    _temlate: ClassVar[list] = {
        "experiment": ["agonist antagonist"],
        "groups": [1, 2, 5, 8],
        "independant_variables": ["TCB2", "MDL"],
        "paired": [False],
        "parametric": [True],
        "outliers": ["grubbs"],
    }


@dataclass
class Metadata:
    project: str = None

    def __post_init__(self):
        self.treatment_mapping = TreatmentMapping(self.project)
        self.experiment_info = ExperimentInfo(self.project)

    def edit(self):
        self.treatment_mapping.edit_excel()
        self.experiment_info.edit_excel()

    def label_compound_data(self, compound_data):
        return compound_data.merge(self.full_mapping, on="group_id")

    @property
    def full_mapping(self):
        return (
            self.experiment_info.get()[["experiment", "groups"]]
            .explode("groups")
            .rename(columns={"groups": "group_id"})
            .merge(self.treatment_mapping, on="group_id")
        )
