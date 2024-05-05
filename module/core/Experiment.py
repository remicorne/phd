from distutils.util import strtobool
from dataclasses import dataclass


@dataclass
class Experiment:

    project: object  # Project
    experiment_information: dict

    def __post_init__(self):
        self.name = self.experiment_information["experiment"]
        self.groups = [
            int(group)
            for group in self.experiment_information["groups"]
            .replace(" ", "")
            .split(",")
        ]
        self.independant_variables = (
            self.experiment_information["independant_variables"]
            .replace(" ", "")
            .split(",")
        )
        self.paired = (
            self.experiment_information["paired"]
            if isinstance(self.experiment_information["paired"], bool)
            else strtobool(self.experiment_information["paired"])
        )
        self.parametric = (
            self.experiment_information["parametric"]
            if isinstance(self.experiment_information["parametric"], bool)
            else strtobool(self.experiment_information["parametric"])
        )
        self.location = f"{self.project.location}/{self.name}"
        del self.experiment_information

    def get_data(self, full_hplc):
        data = full_hplc.select({"group_id": self.groups})
        data.loc[:, "experiment"] = self.name
        data[self.independant_variables] = list(
            data.independant_variables.apply(
                lambda group_independant_variables: [
                    experiment_independant_variable in group_independant_variables
                    for experiment_independant_variable in self.independant_variables
                ],
            )
        )
        return data

    def __repr__(self) -> str:
        return f"""
    groups: {self.groups}
    independant variables: {self.independant_variables}
    paired: {self.paired}
    parametric: {self.parametric}
    """
