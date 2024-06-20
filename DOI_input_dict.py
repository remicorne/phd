'''
Project level file to be filled and ran to generate required dictionaries
The file must be verified to have consistent acronim use for compounds and regions

PROJECT : DOI

'''
from module.constants import CACHE_DIR
from module.utils import checkFileSystem, saveJSON

filename = "DOI_data_HPLC.xlsx" #file name of the excel file saved as csv
# columns : 'mouse_id'  'group_id'  'compound_regions'   ...



########## Project Prams ##########
treatment_mapping = (
    {  # TO DO : change name to treatment_info and add columns in df REMI
        1: {"treatment": "vehicle_cage",
            "color": "thistle",
            "experiments": ["environment"],
            },
        2: {"treatment": "vehicle_plastic",
            "color": "lavender",
            "experiments": ["environment", "DOI_environment"],
            },
        3: {"treatment": "vehicle_orange",
            "color": "sandybrown",
            "experiments": ["dose_response", "environment", "DOI_environment"],
            },
        4: {"treatment": "3mg/kgDOI_plastic",
            "color": "darkgreen",
            "experiments": ["DOI_environment"],
            },
        5: {"treatment": "0.03mg/kgDOI_orange",
            "color": "lightskyblue",
            "experiments": ["dose_response"],
            },
        6: {"treatment": "0.3mg/kgDOI_orange",
            "color": "cornflowerblue",
            "experiments": ["dose_response"],
            },
        7: {"treatment": "3mg/kgDOI_orange",
            "color": "darkcyan",
            "experiments": ["dose_response", "DOI_environment"],
            } } )


experimental_info = {
    "dose_response": {"groups": [3, 5, 6, 7], 
                      "independant_vars": ["DOI"], 
                      "paired": False,
                      "parametric": True, #idealy this would be True / False / Check : check would involve checking the data using the spearman test which should already be done then taking the majority of the data to be parametric or not chose that 
                      "outliers": 'grubbs'
                      },

    "environment": {"groups": [1,2,3],  
                         "independant_vars": ["environment"],
                         "paired": False,
                         "parametric": True,
                         "outliers": "grubbs"
                           },

    "DOI_environment": {"groups": [2,3,4,7],  
                         "independant_vars": ["DOI","environment"],
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
    "LC": "Locus Coeruleus",
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
