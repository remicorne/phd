import os


ROOT : str = os.getcwd()  # This gives terminal location (terminal working dir)
INPUT_DIR : str = f"{ROOT}/input"
OUTPUT_DIR : str = f"{ROOT}/output"
CACHE_DIR : str = f"{INPUT_DIR}/cache"

COLUMN_ORDER : dict = {
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


CORRELOGRAM_COLUMN_ORDER : dict = {
    "compound": COLUMN_ORDER["region"],
    "region": COLUMN_ORDER["compound"],
}


