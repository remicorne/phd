import os
from module.constants import CACHE_DIR
import pandas as pd
from module.utils import checkFileSystem, getJSON, saveJSON

####### FUNCTION THAT DEAL WITH SAVING INFORMATIONS ABOUT EXPERIMENTS ###########
# TODO: replace with user friendly interactive process


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

#OLD causing settingcopy warning 
#     df[new_columns] = df.apply(
#     lambda x: treatment_mapping[str(int(x["group_id"]))].values(),
#     axis=1,
#     result_type="expand",
# )
    #new
    df.loc[:, new_columns] = df.apply(
    lambda x: pd.Series(treatment_mapping[str(int(x["group_id"]))]),
    axis=1
)
    
    # Get alll the values and assign to corresponding columns
    # Duplicate rows belonging to multiple experiment so that groupby can be done later
    return df.explode("experiments").rename(columns={"experiments": "experiment"})


def saveMetadata(
    filename,
    treatment_mapping=None,
    experimental_info=None,
    region_subclassification=None,
    compound_subclassification=None,
    compound_ratio_mapping=None
):
    subcache_dir = f"{CACHE_DIR}/{filename.split('.')[0]}"
    checkFileSystem(subcache_dir)
    if treatment_mapping:
        saveJSON(f"{subcache_dir}/treatment_mapping.json", treatment_mapping)
        print(f"TREATMENT MAPPING {treatment_mapping} SAVED TO {subcache_dir} SUBCACHE")
    if experimental_info:
        saveJSON(f"{subcache_dir}/experimental_info.json", experimental_info)
        print(f"EXPERIMENTAL INFO {experimental_info} SAVED TO {subcache_dir} SUBCACHE")
    if region_subclassification:
        saveJSON(f"{subcache_dir}/region_subclassification.json", region_subclassification)
        print(f"REGION SUBCLASSIFICATION {region_subclassification} SAVED TO {subcache_dir} SUBCACHE")
    if compound_subclassification:
        saveJSON(f"{subcache_dir}/compound_subclassification.json", compound_subclassification)
        print(f"COMPOUND SUBCLASSIFICATION {compound_subclassification} SAVED TO {subcache_dir} SUBCACHE")
    if compound_ratio_mapping:
        saveJSON(f"{subcache_dir}/compound_ratio_mapping.json", compound_ratio_mapping)
        print(f"COMPOUND RATIO MAPPING {compound_ratio_mapping} SAVED TO {subcache_dir} SUBCACHE")

