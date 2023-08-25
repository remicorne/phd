import os
from module.constants import CACHE_DIR
from module.utils import getJSON

def applyTreatmentMapping(df, filename):
    filename = filename.split(".")[0]
    treatment_mapping_path = f"{CACHE_DIR}/{filename}/treatment_mapping.json"
    # Check treatment mapping is present
    if not os.path.isfile(treatment_mapping_path):
        raise Exception(
            "TREATMENT INFORMATION ABSENT, TO SAVE TREATMENT MAPPING RUN 'setTreatment(filename, treatment_mapping)'"
        )
    treatment_mapping = getJSON((treatment_mapping_path))
    # Get the future column names from one of the treatments
    new_columns = list(list(treatment_mapping.values())[0].keys())
    df[new_columns] = df.apply(
        lambda x: treatment_mapping[str(int(x["group_id"]))].values(),
        axis=1,
        result_type="expand",
    )  # Get alll the values and assign to corresponding columns
    # Duplicate rows belonging to multiple experiment so that groupby can be done later
    return df.explode("experiments").rename(columns={"experiments": "experiment"})
