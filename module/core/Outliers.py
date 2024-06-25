from dataclasses import dataclass
from typing import ClassVar
import pandas as pd
from module.core.Dataset import PickleDataset
from module.core.Metadata import ProjectInformation, TreatmentInformation
from module.core.HPLC import HPLC
from tqdm import tqdm
from outliers import smirnov_grubbs as grubbs


def grubbs_test(values, p_value_threshold):
    """
    Takes a list of values on which to perform the test and returns normal values
    """
    return grubbs.test(values, alpha=float(p_value_threshold))


def get_labeled_df(df, test, p_value_threshold):
    # if standar variation is 0, we can't calculate outliers
    only_values = df[df.value != 0].dropna()
    if only_values.value.count() < 3:
        df["outlier_status"] = False
        return df
    outlier_test = OUTLIER_TESTS[test]
    normal_values = outlier_test(only_values.value.tolist(), p_value_threshold)
    df["is_outlier"] = df.value.apply(lambda value: value not in normal_values)
    df["outlier_status"] = df.is_outlier.apply(lambda is_outlier: "suspected" if is_outlier else "normal")
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
            (subset_df, project_information.outlier_test, project_information.p_value_threshold)
            for _, subset_df in hplc.df.groupby(["group_id", "compound", "region", "region"])
        ]
        results = [
            get_labeled_df(*case)
            for case in tqdm(cases, desc="Calculating outliers", unit="group")
        ]
        return pd.concat(results)
    
    def update(self, updates):
        columns = self.df.columns.to_list() + ["updated_outlier_status"] # Eliminate redundant columns to avoid mergs pb
        data = self.extend(updates[columns])
        data.outlier_status = data.updated_outlier_status.combine_first(data.outlier_status)
        data.drop(columns="updated_outlier_status", inplace=True)
        self.save(data)
        