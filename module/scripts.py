import re
import numpy as np
from module.core.display import upload_excel
from module.core.Constants import ConstantRegistry


def raw_hplc():

    def is_valid_raw_col_name(column):
        match = re.match(r"(\w+)[-_ ](\w+)", column)
        return match and len(match.groups()) == 2

    def get_valid_columns(columns):
        column_tuples = np.transpose(column.split("_") for column in columns)
        corrected_compounds_and_regions = np.array([])
        for value_type, values in zip(
            ("compounds", "regions"), np.transpose(column_tuples)
        ):
            unique_values = set(values)
            registry = ConstantRegistry.get_registry(value_type)
            invalid_values = unique_values - set(registry)
            if invalid_values:
                print(f"Invalid {value_type}: {invalid_values}")
            correction_mapper = {
                value: registry.get_valid_choice(value) for value in unique_values
            }
            corrected_compounds_and_regions.append(
                correction_mapper[value] for value in values
            )
            return [
                "_".join(pair) for pair in corrected_compounds_and_regions.transpose()
            ]

    mandatory_columns = ["mouse_id", "group_id"]
    raw_data = upload_excel("Select raw HPLC file to import")
    absent_columns = [col for col in mandatory_columns if col not in raw_data.columns]
    if absent_columns:
        raise ValueError(f"Absent columns: {absent_columns}")
    invalid_columns = list(is_valid_raw_col_name, filter(raw_data.columns))
    if any(invalid_columns):
        raise ValueError(f"Invalid columns: {invalid_columns}, correct names and retry")
    return raw_data.rename(columns=get_valid_columns(raw_data.columns))
