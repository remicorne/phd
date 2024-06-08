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
        df["is_outlier"] = False
        return df
    outlier_test = OUTLIER_TESTS[test]
    normal_values = outlier_test(only_values.value.tolist(), p_value_threshold)
    df["is_outlier"] = df.value.apply(lambda x: x not in normal_values)
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
        # Timin=region list comprehension with tqdm
        results = [
            get_labeled_df(*case)
            for case in tqdm(cases, desc="Calculating outliers", unit="group")
        ]
        return pd.concat(results)
    
    def select_outliers(self, treatment, compound, region):
        df = self.extend(TreatmentInformation(self.project)).select(treatment=treatment, compound=compound, region=region)
        return TreatmentOutlierFigure(df,  compound, region,treatment)
        
        