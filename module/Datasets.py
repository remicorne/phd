import pandas as pd
import numpy as np
from module.Filesystem import Filesystem
from module.Cache import Cache



class Dataset:
    """Is the parent dataset class implementing everythingthat is generic about various datasets
    NOT TO BE USED, ONLY TO BE INHERITED FROM
    """
    dataset_type = ''
    data_type = ""
    

    def __init__(self, experiment, config) -> None:
        self.experiment = experiment
        self.path = f"{Filesystem.INPUT}/{self.experiment}"
        self.cache = Cache(self.path, self.data_type)
        self.config = config
        self.data = self.cache.get_or_add(self.initialise_data)
    
    def select(self, **selectors):
        for column, value in selectors.items():
            if value not in self.data[column].unique():
                if value in self.config.region_subclassification[value]['region']:
                    value = self.config.region_subclassification[value]['region']
            df = df[
                df[column].isin(value) if isinstance(value, list) else df[column] == value
            ]
        return df

    def initialise_data(self):
        """The code is dataset specific and therefore only present in child classes
        """










class RawDataset(Dataset):
    dataset_type = ''
    raw_data_extension = '.csv'

    def initialise_data(self):
        filename =  f"{self.experiment}_{self.dataset_type}{self.raw_data_extension}"
        try:
            raw_data = pd.read_csv(f"{self.path}/{filename}", header=0).replace(0, np.nan)
            return self.format_raw_data(raw_data)
        except FileNotFoundError:
            print(f"MISSING {filename} FILE, PLEASE ADD IT TO {self.path}")
            return pd.DataFrame()
   
    def format_raw_data(self, raw_data):
        raise NotImplementedError("This method should be implemented in child classes")


 
def calculateRatio(num_ratios):
    def calculator(row):
        ratio_name = f"{row.compound_1}/{row.compound_2}"
        print("CALCULATING", row.name, "OF", num_ratios, "RATIOS")
        return [
            ratio_name,
            row.value_1 / row.value_2 if row.value_2 else np.NaN,
        ]
    return calculator
       
class HPLC(RawDataset):

    raw_dataset_type = "HPLC"

        
    def format_raw_data(self, raw_data):
        new_format_df = raw_data.melt(
        id_vars=["mouse_id", "group_id"], value_vars=raw_data.columns[2:]
        )
        new_format_df[["compound", "region"]] = new_format_df.apply(
            lambda x: x.variable.split("_"), axis=1, result_type="expand"
        )
        
        compound_df = self.config.apply_treatment_mapping(new_format_df).drop(columns=["variable"])
        ratios_df = pd.merge(
            left=compound_df,
            right=compound_df,
            on=[
                "mouse_id",
                "group_id",
                "region",
                "experiment",
                "color",
                "treatment",
            ],
            suffixes=["_1", "_2"],
        ).reset_index(
            drop=True
        )  # merge every compound to every other for each mouse, we want to reset the index (ie make sure it has no duplicates) becaus many update operations will use it

        ratios_df = ratios_df[(ratios_df.compound_1 != ratios_df.compound_2)]

        ratios_df[["compound", "value"]] = ratios_df.apply(
            calculateRatio(len(ratios_df)),
            axis=1,
            result_type="expand",
        )  # calculate ratio
        # Drop duplicate columns
        compound_and_ratios_df = pd.concat(
            [
                compound_df,
                ratios_df.drop(columns=["compound_1", "compound_2", "value_1", "value_2"]),
            ]
        )
        return compound_and_ratios_df
        
class HT(RawDataset):
    dataset_type = "HT"

    def format_raw_data(self, raw_data):
        return self.config.apply_treatment_mapping(raw_data)
        