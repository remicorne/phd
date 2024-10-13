import os, re
from dataclasses import dataclass
from typing import ClassVar
import pandas as pd
import numpy as np
from module.core.Dataset import PickleDataset, SelectableDataFrame
from module.core.Metadata import (
    ProjectInformation,
    ExperimentInformation,
    TreatmentInformation,
)
from module.core.questions import yes_or_no, input_escape
from module.core.Constants import REGIONS, COMPOUNDS
from tqdm import tqdm
from outliers import smirnov_grubbs as grubbs
from module.core.utils import parallel_process


def detect_raw_data(project):
    for filename in os.listdir():
        if project.name.upper() in filename.upper():
            if os.path.splitext(filename)[1].lower() in [".csv", ".xlsx"]:
                is_correct = yes_or_no(f"Detected {filename}. Confirm?")
                if is_correct:
                    return filename
    return None


def handle_raw_col_name(column):
    match = re.match(r"(\w+)[-_ ](\w+)", column)
    while not match or not len(match.groups()) == 2:
        column = input(f"Wrong format: '{column}', input new 'compound_region': ")
        match = re.match(r"(\w+)[-_ ](\w+)", column)
    return match.groups()


@dataclass(repr=False)
class RawHPLC(PickleDataset):

    project: str
    filename: ClassVar[str] = "raw_hplc"

    def generate(self):
        project_information = ProjectInformation(self.project)
        filepath = f"{os.getcwd()}/{project_information.raw_data_filename}"
        extension = os.path.splitext(project_information.raw_data_filename)[1].lower()
        raw_data = (
            pd.read_excel(filepath) if extension == ".xlsx" else pd.read_csv(filepath)
        )
        return raw_data.rename(columns=self.get_valid_columns(raw_data.columns))

    def get_valid_columns(self, columns):
        mandatory_columns = ["mouse_id", "group_id"]
        absent_columns = [col for col in mandatory_columns if col not in columns]
        if absent_columns:
            raise ValueError(f"Absent columns: {(', ').join(absent_columns)}")
        compound_region_tuples = [
            handle_raw_col_name(col) for col in columns if col not in mandatory_columns
        ]
        invalid_compounds = {
            compound
            for compound, _ in compound_region_tuples
            if compound not in COMPOUNDS
        }
        if invalid_compounds:
            print(f"Invalid compounds: {invalid_compounds}")

        invalid_regions = {
            region for _, region in compound_region_tuples if region not in REGIONS
        }
        if invalid_regions:
            print(f"Invalid regions: {invalid_regions}")

        compound_translator = {
            invalid_compound: COMPOUNDS.get_valid_choice(invalid_compound)
            for invalid_compound in invalid_compounds
        }
        region_translator = {
            invalid_region: REGIONS.get_valid_choice(invalid_region)
            for invalid_region in invalid_regions
        }
        return {
            f"{compound}_{region}": f"{compound if compound in COMPOUNDS else compound_translator[compound]}_{region if region in REGIONS else region_translator[region]}"
            for compound, region in compound_region_tuples
        }


class ProjectSelectableDataframe(SelectableDataFrame):
    """
    Shameful hack.
    TODO: Eliminate asap
    """
    
    def __init__(self, data=None, project=None, *args, **kwargs):
        super().__init__(data, *args, **kwargs)
        self.project = project
        
    @property
    def _constructor(self):
        return lambda *args, **kwargs: ProjectSelectableDataframe(*args, project=self.project, **kwargs) if hasattr(self, "project") else SelectableDataFrame(*args, **kwargs)

    def select(self, **selector) -> SelectableDataFrame:
        experiment = selector.pop("experiment", None)
        if experiment:
            experiment = ExperimentInformation(self.project).select(label=experiment)
            selector["group_id"] = experiment.groups
        data = SelectableDataFrame(self).select(**selector)
        return ProjectSelectableDataframe(data, project=self.project) if hasattr(self, "project") else data

@dataclass(repr=False)
class HPLC(PickleDataset):

    project: str
    filename: ClassVar[str] = "hplc"

    def generate(self):
        raw_data = RawHPLC(self.project).df
        compound_data = raw_data.melt(
            id_vars=["mouse_id", "group_id"], value_vars=raw_data.columns[2:]
        )
        compound_data[["compound", "region"]] = compound_data["variable"].str.split(
            "_", expand=True
        )
        compound_data = compound_data.drop(columns=["variable"])
        ratio_data = pd.merge(
            left=compound_data,
            right=compound_data,
            on=[
                "mouse_id",
                "group_id",
                "region",
            ],
            suffixes=["_1", "_2"],
        ).reset_index(
            drop=True
        )  # merge every compound to every other for each mouse, we want to reset the index (ie make sure it has no duplicates) becaus many update operations will use it

        ratio_data = ratio_data[(ratio_data.compound_1 != ratio_data.compound_2)]

        def calculateRatio(row):
            ratio_name = f"{row.compound_1}/{row.compound_2}"
            return [
                ratio_name,
                row.value_1 / row.value_2 if row.value_1 and row.value_2 else np.nan,
            ]

        tqdm.pandas(desc="Calculating ratios", unit="ratio")

        ratio_data[["compound", "value"]] = ratio_data.progress_apply(
            calculateRatio,
            axis=1,
            result_type="expand",
        )  # calculate ratio
        # Drop duplicate columns
        compound_and_ratios_df = pd.concat(
            [
                compound_data,
                ratio_data.drop(
                    columns=["compound_1", "compound_2", "value_1", "value_2"]
                ),
            ]
        )
        return compound_and_ratios_df.replace(0, np.nan)
    
    def select(self, **selector) -> ProjectSelectableDataframe:
        return self.full_df.select(**selector)

    @property
    def full_df(self) -> ProjectSelectableDataframe:
        data = ProjectSelectableDataframe(pd.concat([self.df, TissueWeight(self.project).df]))
        data = ProjectSelectableDataframe(
            data.extend(Outliers(self.project)).extend(
                TreatmentInformation(self.project)
            ), self.project
        )
        return data[[col for col in data.columns if "Unnamed" not in col]]

    @property
    def compounds(self):
        return self.df.compound.unique()

    @property
    def regions(self):
        return self.df.region.unique()

    @property
    def compounds_and_regions(self):
        return {
            "compounds": self.df.compound.unique(),
            "regions": self.df.region.unique(),
        }


def grubbs_test(values, p_value_threshold):
    """
    Takes a list of values on which to perform the test and returns normal values
    """
    return grubbs.test(values, alpha=float(p_value_threshold))


def get_labeled_df(df__test__p_value_threshold):
    # if standar variation is 0, we can't calculate outliers
    df, test, p_value_threshold = df__test__p_value_threshold
    only_values = df[df.value != 0].dropna()
    if only_values.value.count() < 3:
        df["outlier_status"] = False
        return df
    outlier_test = OUTLIER_TESTS[test]
    normal_values = outlier_test(only_values.value.tolist(), p_value_threshold)
    df["is_outlier"] = df.value.apply(lambda value: value not in normal_values)
    df["outlier_status"] = df.is_outlier.apply(
        lambda is_outlier: "suspected" if is_outlier else "normal"
    )
    return df


OUTLIER_TESTS = {"grubbs": grubbs_test}


@dataclass(repr=False)
class Outliers(PickleDataset):

    project: str
    filename: ClassVar[str] = "outliers"

    def generate(self):
        hplc = HPLC(self.project)
        project_information = ProjectInformation(self.project)
        cases = [
            (
                subset_df,
                project_information.outlier_test,
                project_information.p_value_threshold,
            )
            for _, subset_df in hplc.df.groupby(
                ["group_id", "compound", "region", "region"]
            )
        ]
        results = parallel_process(
            cases, get_labeled_df, description="Calculating outliers"
        )
        return pd.concat(results)

    def update(self, updates):
        data = self.df.set_index(["compound", "region", "mouse_id", "group_id"])
        updates = updates.set_index(["compound", "region", "mouse_id", "group_id"])
        data.update(updates[["outlier_status"]])
        self.save(data.reset_index())

   
@dataclass 
class TissueWeight(PickleDataset):
    
    project: str
    filename: ClassVar[str] = "tissue_weight"
    
    def generate(self):
        filename = None
        while not os.path.exists(f"{os.getcwd()}/{filename}"):
            filename = input_escape("Enter tissue weight filename: ")
        filepath = f"{os.getcwd()}/{filename}"
        extension = os.path.splitext(filename)[1].lower()
        raw_data = pd.read_excel(filepath) if extension == ".xlsx" else pd.read_csv(filepath)
        raw_data.rename(columns={compound: REGIONS.get_valid_choice(compound) for compound in raw_data.columns[1:]})
        raw_data = Groups(self.project).extend(raw_data)
        raw_data = raw_data.melt(
            id_vars=["mouse_id", "group_id"], value_vars=raw_data.columns[2:]
        )
        raw_data["compound"] = "weight"
        return raw_data.rename(columns={"variable": "region"})
        

@dataclass
class Groups(PickleDataset):
    
    project: str
    filename: ClassVar[str] = "groups"
    
    def generate(self):
        hplc_data = HPLC(self.project).df[["mouse_id", "group_id"]]
        return hplc_data.drop_duplicates()