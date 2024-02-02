from main import *

# SET FILENAME
filename = "ALZHEIMERS_HPLC.xlsx"  #in current working directory 

#LOAD METTADATA 
correlogram(filename, 
            experiment='Alzheimers', #  dose_response agonist_antagonist
            correlogram_type='compound',
            to_correlate='DA', 
            p_value_threshold=0.05, 
            n_minimum=5, 
            columns=["MO","LO","M2","IL","aCg","pCg","S1","ENT","A","DH", "VH","SH", "CO", "GP","vlT","HYP","CB"],
            from_scratch=True,
            
            )