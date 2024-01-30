import os, re
from dataclasses import dataclass
from typing import ClassVar
import pandas as pd
from module.core.Dataset import Dataset
from module.core.Questions import Questions
from module.core.Metadata import Metadata
from module.core.Constants import REGIONS, COMPOUNDS

ROOT = os.getcwd()  # This gives terminal location (terminal working dir)
CACHE = "cache"


def detect_raw_data(project):
    for filename in os.listdir(ROOT):
        if project.upper() in filename.upper():
            if os.path.splitext(filename)[1].lower() in [".csv", ".xlsx"]:
                is_correct = Questions.yes_or_no(f"Detected {filename}. Confirm?")
                if is_correct:
                    return filename
    return None


def is_valid_file(file_path):
    if not os.path.isfile(file_path):
        print("Not found", file_path)
        return False
    extension = os.path.splitext(file_path)[1].lower()
    if extension not in [".xlsx", ".csv"]:
        print("Invalid extension:", extension)
        return False

    return True

def handle_raw_col_name(column):
    match = re.match(r"(\w+)[-_ ](\w+)", column)
    while not match or not len(match.groups()) == 2:
        column = Questions.input(
            f"Wrong format: '{column}', input new 'compound_region': "
        )
        match = re.match(r"(\w+)[-_ ](\w+)", column)
    return match.groups()

@dataclass
class RawHPLC(Dataset):
    
    _type: ClassVar[str] = "raw_hplc"

    def generate_data(self):
        raw_data = self.get_raw_data()
        return raw_data.rename(columns=self.get_valid_columns(raw_data.columns))

    def get_raw_data(self):
        raw_data_filename = detect_raw_data(self.project) or Questions.input(
            "Enter HPLC excel filename"
        )
        file_path = f"{ROOT}/{raw_data_filename}"

        while not is_valid_file(file_path):
            print(raw_data_filename, "NOT FOUND")
            raw_data_filename = Questions.input("Enter excel HPLC filename")
            file_path = f"{ROOT}/{raw_data_filename}"

        extension = os.path.splitext(file_path)[1].lower()
        return pd.read_excel(file_path) if extension == ".xlsx" else pd.read_csv(file_path)

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
            invalid_region: REGIONS.get_valid_choice(invalid_region) for invalid_region in invalid_regions
        }
        return {
            f"{compound}_{region}": f"{compound if compound in COMPOUNDS else compound_translator[compound]}_{region if region in REGIONS else region_translator[region]}"
            for compound, region in compound_region_tuples
        }

@dataclass
class HPLC(Dataset):
    
    _type: ClassVar[str] = "compound_and_ratios"

    def generate_data(self):
        raw_data = RawHPLC(self.project).get()
        compound_data = raw_data.melt(
            id_vars=["mouse_id", "group_id"], value_vars=raw_data.columns[2:]
        )
        compound_data[["compound", "region"]] = compound_data.apply(
            lambda x: x.variable.split("_"), axis=1, result_type="expand"
        )
        compound_data = Metadata(self.project).label_compound_data(compound_data)
        compound_data = compound_data.drop(columns=["variable"])
        ratio_data = pd.merge(
            left=compound_data,
            right=compound_data,
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

        ratio_data = ratio_data[(ratio_data.compound_1 != ratio_data.compound_2)]

        def calculateRatio(row):
            ratio_name = f"{row.compound_1}/{row.compound_2}"
            print("CALCULATING", row.name, "OF", len(ratio_data), "RATIOS")
            return [
                ratio_name,
                row.value_1 / row.value_2 if row.value_2 else None,
            ]

        ratio_data[["compound", "value"]] = ratio_data.apply(
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
        return compound_and_ratios_df
