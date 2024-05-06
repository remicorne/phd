'''
Project level file to be filled and ran to generate required dictionaries

PROJECT : ALZHEIMERS

'''
from module.constants import CACHE_DIR
from module.utils import checkFileSystem, saveJSON

filename = "ALZHEIMERS_HPLC.csv" #file name of the excel file saved as csv
# columns : 'mouse_id'  'group_id'  'compound_regions'   ...

########## Project Prams ##########
treatment_mapping = (
    {  # TO DO : change name to treatment_info and add columns in df REMI
        1: {"treatment": "young_WT",
            "color": "white",
            "experiments": ["Alzheimers"],
            },
        2: {"treatment": "young_AD",
            "color": "bisque",
            "experiments": ["Alzheimers"],
            },
        3: {"treatment": "old_WT",
            "color": "red",
            "experiments": ["Alzheimers"],
            },
        4: {"treatment": "old_AD",
            "color": "firebrick",
            "experiments": ["Alzheimers"],
            },
    }
)

experimental_info = {
    "Alzheimers": {"groups": [1, 2, 3, 4], 
                      "independant_vars": ["old", "AD"], 
                      "paired": False,
                      "parametric": True, #idealy this would be True / False / Check : check would involve checking the data using the spearman test which should already be done then taking the majority of the data to be parametric or not chose that 
                      "outliers": 'grubbs'
                      },
}


########## KEYS ##########

compounds = {
    'NA': 'Noradrenaline',
    'DO': '3,4-dihydroxyphenylacetic acid',
    'DA': 'Dopamine',
    'HVA': 'Homovanillic acid',
    '5HIAA': '5-Hydroxyindoleacetic acid',
    '5HT': 'Serotonin'
}

regions = {
    'MO': 'Medial Orbitofrontal Cortex',
    'LO': 'Lateral Orbitofrontal Cortex',
    'M2': 'Secondary Motor Cortex',
    'PL': 'Prelimbic Cortex',
    'IL': 'Infralimbic Cortex',
    'aCg': 'Anterior Cingulate Cortex (Granular)',
    'pCg': 'Posterior Cingulate Cortex (Granular)',
    'S1': 'Primary Somatosensory Cortex',
    'ENT': 'Entorhinal Cortex',
    'A': 'Amygdala',
    'DH': 'Dorsal Hippocampus',
    'VH': 'Ventral Hippocampus',
    'SH': 'Shell of the Nucleus Accumbens',
    'CO': 'Core of the Nucleus Accumbens',
    'VM': 'Ventromedial Striatum',
    'DM': 'Dorsomedial Striatum',
    'VL': 'Ventrolateral Striatum',
    'DL': 'Dorsolateral Striatum',
    'GP': 'Globus Pallidus',
    'vlT': 'Ventrolateral Thalamus',
    'MD': 'Mediodorsal Thalamus',
    'HYP': 'Hypothalamus',
    'CB': 'Cerebellum',
    'SN': 'Substantia Nigra',
    'VTA': 'Ventral Tegmental Area',
    'DR': 'Dorsal Raphe',
    'MR': 'Medial Raphe'
}

########## SAVE ##########
project_dict = {
    'treatment_mapping': treatment_mapping,
    'experimental_info': experimental_info,
    'compounds': compounds,
    'regions': regions
}

def saveProjectDict(filename, project_dict):
    subcache_dir = f"{CACHE_DIR}/{filename.split('.')[0]}"
    checkFileSystem(subcache_dir)
    for name, data in project_dict.items():
        saveJSON(f"{subcache_dir}/{name}.json", data)
        print(f"{name} SAVED TO {subcache_dir} SUBCACHE")

#### INIT ####
saveProjectDict(filename, project_dict)
