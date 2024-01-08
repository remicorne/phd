
#TODO remove this function and replace with a saver for this page only so it can be filled in 
#for one experiment because adjusting every json is too much 

# saveMetadata(
#     filename, treatment_mapping=treatment_mapping, experimental_info=experimental_info, region_subclassification=region_subclassification, compound_subclassification=compound_subclassification, compound_ratio_mapping=compound_ratio_mapping
# )


##### DOCS ? #####

# outlier tests  = ['grubbs', 'ROUT']  #todo add code for ROUT
# post hoc tests = [ 'pLSD' , 'tukey' ] #todo add code for
# quantitative stats = ['two_way_anova', 'one_way_anova', 'tukey'] #TODO this should not be here but with post hoc not written yet
# correlative stats = ['pearson', 'spearmam', 'kendall']


#DFs floating around
# aggergate_stats, quantative_stats, compound_and_ratio_df ++graph_stats


##### INPUT DICTIONARY #####

# treatment_mapping = (
#     {  # TO DO : change name to treatment_info and add columns in df REMI
#         1: {"treatment": "vehicles",
#             "color": "white",
#             "experiments": ["dose_response", "agonist_antagonist"],
#             },
#         2: {"treatment": "0.3mg/kgTCB",
#             "color": "salmon",
#             "experiments": ["dose_response"],
#             },
#         3: {"treatment": "3mg/kgTCB",
#             "color": "red",
#             "experiments": ["dose_response", "agonist_antagonist"],
#             },
#         4: {"treatment": "10mg/kgTCB",
#             "color": "firebrick",
#             "experiments": ["dose_response"],
#             },
#         5: {"treatment": "0.2mg/kgMDL",
#             "color": "grey",
#             "experiments": ["agonist_antagonist"],
#             },
#         6: {"treatment": "TCB+MDL",
#             "color": "black",
#             "experiments": ["agonist_antagonist"],
#             },
#     }
# )

# compound_ratio_mapping = { 
#     "DOPAC": ["DA"],
#     "5HIAA": ["5HT"],
#     "3MT": ["DA"],
#     "HVA": ["DA", "3MT", "DOPAC"],
#     "GLN": ["GLU"],
# }

# region_subclassification = {
#     'cortex': {'regions': ['OF', 'PL', 'CC', 'IC', 'M', 'SL1', 'SR1', 'AC', 'V'], 'color': 'mediumblue'}, #SR6, SL6 ?
#     'subcortical_telencephalon': {'regions': ['Am', 'dH', 'vH', 'NAc', 'VM', 'DM', 'VL', 'DL'], 'color': 'orange'},
#     'diencephalon': {'regions': ['MD', 'VPL', 'VPR', 'DG', 'Y'], 'color': 'darkorchid'},
#     'mesencephalon': {'regions': ['SC', 'SN', 'VTA', 'DR', 'MR'], 'color': 'forestgreen'},
#     'cerebellum': {'regions': ['CE'], 'color': 'peru'}
# }

# compound_subclassification = {
#     'monoamines':['LDOPA', 'NA', 'A', '5HTP', 'DOPAC', 'DA', '5HIAA', 'HVA', '5HT', 'VMA', '3MT'],
#     'amino_acids':['ASP', 'GLU', 'ASPN', 'HIS', 'LSER', 'GLN', 'ARG', 'GLY', 'THR', 'TAU', 'ALA', 'TYR', 'GABA'],
#     'neurotransmitters':['DA','NA','5HT','GLU','GABA','ASP','GLY']
# }

# experimental_info = {
#     "dose_response": {"groups": [1, 2, 3, 4], 
#                       "independant_vars": ["TCB2"], 
#                       "paired": False,
#                       "parametric": True, #idealy this would be True / False / Check : check would involve checking the data using the spearman test which should already be done then taking the majority of the data to be parametric or not chose that 
#                       "outliers": 'grubbs'
#                       },

#     "agonist_antagonist": {"groups": [1,3,5,6],  
#                          "independant_vars": ["TCB2","MDL"],
#                          "paired": False,
#                          "parametric": True,
#                          "outliers": "grubbs"
#                            },
# }

