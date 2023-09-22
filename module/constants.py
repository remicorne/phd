######### HERE WE DECLARE THE CONSTANTS USED BY OTHER FILES ############
# Constant are meant to be constants, the should not changed, that's what wariables or user are for

import os

### Constants that reflect the filesystem structure, used by util functions
ROOT = os.getcwd()  # This gives terminal location (terminal working dir)
INPUT_DIR = f"{ROOT}/input"
OUTPUT_DIR = f"{ROOT}/output"
CACHE_DIR = f"{INPUT_DIR}/cache"

COLUMN_ORDER = {
    "region": [
        "OF",
        "PL",
        "CC",
        "IC",
        "M",
        "SJ",
        "SL1",
        "SL6",
        "SR6",
        "SR1",
        "AC",
        "V",
        "Am",
        "dH",
        "vH",
        "NAc",
        "VM",
        "DM",
        "VL",
        "DL",
        "MD",
        "VPL",
        "VPR",
        "DG",
        "Y",
        "SC",
        "SN",
        "VTA",
        "DR",
        "MR",
        "CE",
    ],
    "compound": ["DA", "DOPAC", "HVA", "3MT", "5HT", "5HIAA", "GLU", "GLN"],
}

# Converts keyword 'region' to 'compound' and vice versa, used in conjunction with the function below
CORRELOGRAM_TYPE_CONVERTER = {"region": "compound", "compound": "region"}


def getCorrelogramColumns(correlogram_type):
    """
    Get the default correlogram colums to display based on correlogram_type in ['region', 'compound']
    """
    return COLUMN_ORDER[CORRELOGRAM_TYPE_CONVERTER[correlogram_type]]
