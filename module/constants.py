import os


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

CORRELOGRAM_TYPE_CONVERTER = {"region": "compound", "compound": "region"}


def getCorrelogramColumns(correlogram_type):
    return COLUMN_ORDER[CORRELOGRAM_TYPE_CONVERTER[correlogram_type]]
