from main import *

# SET FILENAME
filename = "ALZHEIMERS_HPLC.xlsx"  #in current working directory 

quantitativeHistogram( 
    filename,
    p_value_threshold=0.05,
    compound='5HIAA/5HT',
    region='IL',
    experiment= 'Alzheimers', 
    do_outliers=False, 
    from_scratch=True
)