'''
Project level file to be filled and ran to generate required dictionaries
The file must be verified to have consistent acronim use for compounds and regions

PROJECT : TCB-2

'''
from module.constants import CACHE_DIR
from module.utils import checkFileSystem, saveJSON

filename = "TCB2_data_HPLC.csv" #file name of the excel file saved as csv
# columns : 'mouse_id'  'group_id'  'compound_regions'   ...



########## Project Prams ##########

treatment_mapping = (
    {  # TO DO : change name to treatment_info and add columns in df REMI
        1: {"treatment": "vehicles",
            "color": "white",
            "experiments": ["dose_response", "agonist_antagonist"],
            },
        2: {"treatment": "0.3mg/kgTCB",
            "color": "salmon",
            "experiments": ["dose_response"],
            },
        3: {"treatment": "3mg/kgTCB",
            "color": "red",
            "experiments": ["dose_response", "agonist_antagonist"],
            },
        4: {"treatment": "10mg/kgTCB",
            "color": "firebrick",
            "experiments": ["dose_response"],
            },
        5: {"treatment": "0.2mg/kgMDL",
            "color": "grey",
            "experiments": ["agonist_antagonist"],
            },
        6: {"treatment": "TCB+MDL",
            "color": "black",
            "experiments": ["agonist_antagonist"],
            },
    }
)

experimental_info = {
    "dose_response": {"groups": [1, 2, 3, 4], 
                      "independant_vars": ["TCB2"], 
                      "paired": False,
                      "parametric": True, #idealy this would be True / False / Check : check would involve checking the data using the spearman test which should already be done then taking the majority of the data to be parametric or not chose that 
                      "outliers": 'grubbs'
                      },

    "agonist_antagonist": {"groups": [1,3,5,6],  
                         "independant_vars": ["TCB2","MDL"],
                         "paired": False,
                         "parametric": True,
                         "outliers": "grubbs"
                           },
}


########## KEYS ##########

compounds = {
    'LDOPA': 'L-3,4-dihydroxyphenylalanine',
    'NA': 'Noradrenaline',
    'A': 'Adrenaline',
    '5HTP': '5-Hydroxytryptophan',
    'DOPAC': '3,4-Dihydroxyphenylacetic acid',
    'DA': 'Dopamine',
    '5HIAA': '5-Hydroxyindoleacetic acid',
    'HVA': 'Homovanillic acid',
    '5HT': 'Serotonin',
    'VMA': 'Vanillylmandelic acid',
    '3MT': '3-Methoxytyramine',
    'ASP': 'Aspartate',
    'GLU': 'Glutamate',
    'ASPN': 'Asparagine',
    'HIS': 'Histidine',
    'LSER': 'L-Serine',
    'GLN': 'Glutamine',
    'ARG': 'Arginine',
    'GLY': 'Glycine',
    'THR': 'Threonine',
    'TAU': 'Taurine',
    'ALA': 'Alanine',
    'TYR': 'Tyrosine',
    'GABA': 'Gamma-Aminobutyric acid'
}
regions = {
  "OF": "Orbital Frontal Cortex",
    "PL": "Prelimbic Cortex",
    "CC": "Cingulate Cortex",
    "IC": "Insular Cortex",
    "M": "Motor Cortex",
    "SJ": "Primary Somatosensory Cortex",
    "SL1": "Left Somatosensory Cortex Layer 1",
    "SL6": "Left Somatosensory Cortex Layer 6",
    "SR6": "Right Somatosensory Cortex Layer 6",
    "SR1": "Right Somatosensory Cortex Layer 1",
    "AC": "Auditory Cortex",
    "V": "Visual Cortex",
    "Am": "Amygdala",
    "dH": "Dorsal Hippocampus",
    "vH": "Ventral Hippocampus",
    "NAc": "Nucleus Accumbens",
    "VM": "Ventromedial Thalamus",
    "DM": "Dorsomedial Thalamus",
    "VL": "Ventrolateral Thalamus",
    "DL": "Dorsolateral Thalamus",
    "MD": "Medio Dorsal Thalamus",
    "VPL": "Left Ventral Posterior Medial Thalamus",
    "VPR": "Right Ventral Posterior Medial Thalamus",
    "DG": "Lateral Geniculate Nucleus",
    "Y": "Lateral Hypothalamus",
    "SC": "Superior Colliculus",
    "SN": "Substantia Nigra",
    "VTA": "Ventral Tegmental Area",
    "DR": "Dorsal Raphe Nucleus",
    "MR": "Medial Raphe Nucleus",
    "CE": "Cerebellum"
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
