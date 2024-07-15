from dataclasses import dataclass
from typing import ClassVar
import pandas as pd
from module.core.Dataset import PickleDataset
from module.core.Metadata import ProjectInformation, TreatmentInformation
from module.core.HPLC import HPLC
from tqdm import tqdm
from outliers import smirnov_grubbs as grubbs
from module.core.utils import parallel_process


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
        results = parallel_process(
            cases, get_labeled_df, description="Calculating outliers")
        return pd.concat(results)
    
    def update(self, updates):
        data = self.df.set_index(['compound', 'region', 'mouse_id', 'group_id'])
        updates = updates.set_index(['compound', 'region', 'mouse_id', 'group_id'])
        data.update(updates[["outlier_status"]])
        self.save(data.reset_index())
        