######### HERE WE DECLARE THE CONSTANTS USED BY OTHER FILES ############
# Constant are meant to be constants, the should not changed, that's what variables or user are for

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
    "compound": ['DA','NA','5HT','GLU','GABA','ASP','GLY'],#perhaps this should be a full list
}

# Converts keyword 'region' to 'compound' and vice versa, used in conjunction with the function below
CORRELOGRAM_TYPE_CONVERTER = {"region": "compound", "compound": "region"}


def getCorrelogramColumns(correlogram_type):
    """
    Get the default correlogram colums to display based on correlogram_type in ['region', 'compound']
    """
    return COLUMN_ORDER[CORRELOGRAM_TYPE_CONVERTER[correlogram_type]]


########## SUBCLASIFICATIONS ########## 

#TODO this needs to be saved and also applied but will work only after 
#the naming of regions and compounds has been set to something consistent and a verification function has been made

region_subclassification = {
    'cortex': {'regions': ['OF', 'PL', 'CC', 'IC', 'M', 'SL1', 'SR1', 'AC', 'V'], 'color': 'mediumblue'}, #SR6, SL6 ?
    'subcortical_telencephalon': {'regions': ['Am', 'dH', 'vH', 'NAc', 'VM', 'DM', 'VL', 'DL'], 'color': 'orange'},
    'diencephalon': {'regions': ['MD', 'VPL', 'VPR', 'DG', 'Y'], 'color': 'darkorchid'},
    'mesencephalon': {'regions': ['SC', 'SN', 'VTA', 'DR', 'MR'], 'color': 'forestgreen'},
    'cerebellum': {'regions': ['CE'], 'color': 'peru'}
}

compound_subclassification = {
    'monoamines':['LDOPA', 'NA', 'A', '5HTP', 'DOPAC', 'DA', '5HIAA', 'HVA', '5HT', 'VMA', '3MT'],
    'amino_acids':['ASP', 'GLU', 'ASPN', 'HIS', 'LSER', 'GLN', 'ARG', 'GLY', 'THR', 'TAU', 'ALA', 'TYR', 'GABA'],
    'neurotransmitters':['DA','NA','5HT','GLU','GABA','ASP','GLY']
}

region_circuitry = {
    'reward_system': ['PL', 'NAc', 'Am', 'VTA', 'Y'],
    'vision_equilibrium': ['V', 'DLG', 'SC', 'CE', 'SN'],
    'basal_ganglia': ['DM', 'DL', 'VM', 'VL', 'NAC', 'SN', 'VTA', 'DR', 'MR'],
    'cognitive_loops': ['DM', 'DL', 'VM', 'VL', 'NAC', 'OF', 'PL', 'CC'],
    'motor_loops': ['DM', 'DL', 'VM', 'VL', 'NAC', 'M', 'SJ'],
    'sensory_systems': ['V', 'SL1', 'SR1', 'AC', 'DLG', 'VPL', 'VPR'],
    'limbic_loops': ['Am', 'dH', 'vH', 'MD', 'OF', 'CC', 'NAC', 'VTA', 'DR', 'MR', 'Y'],
    'thalamocortical_interaction': ['OF', 'PL', 'CC', 'M', 'SJ', 'SL1', 'SR1', 'AC', 'V', 'MD', 'VPL', 'VPR', 'DLG']
}

########## STANDARD KEYS ##########
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

# TCB2 regions need to be replaced with allen or atlas regions , and aditionaly the x,y,z spacial cordinates 
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
