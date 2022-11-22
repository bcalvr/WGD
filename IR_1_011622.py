import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from datetime import datetime as dt
from datetime import timedelta
from scipy.signal import lfilter
import pickle as pickle
import sys

# load required files

# CHANGE DATE
with open(f"../Input_files/jh_daily_ir_01_16_22.pkl",'rb') as file: #Contains only country name for each date
    jh_daily_ir = pickle.load(file)
# CHANGE DATE
with open(f"../Input_files/jh_data_01_16_22.pkl",'rb') as file:
    jh_data=pickle.load(file)
# CHANGE DATE
with open(f"../Input_files/jh_daily_cases_01_16_22.pkl",'rb') as file:
    jh_daily_cases = pickle.load(file)   
# CHANGE DATE
with open(f"../Input_files/seq_country_counts_01_16_22.pkl",'rb') as file:
    sequence_country_counts = pickle.load(file)
# CHANGE DATE
with open(f"../Input_files/basis_vector_IR_FR_PP_01_16_22.pkl",'rb') as file:
    IR = pickle.load(file)
# CHANGE DATES
with open(f"../Input_files/aggregated_mutations_window00_12_16_21_01_16_22.pkl",'rb') as file:
    aggregated_mutations_window = pickle.load(file)
    
    

#Looping over single days
list_dates=jh_daily_ir.columns[9::]

# For this step, create a few parallel scripts, each cycling through ~15 days. For example, if your update covers 30 days, you would create 2 scripts for each IR, PP and lagged_FR.
# Get the indexes corresponding to days by looking up list_dates as before.

# CHANGE INDEXES
for date in tqdm(list_dates[684:700]):
#Looping over single variants
    for var in aggregated_mutations_window[date].keys():
        #For each variant count countries, 
        #totalCounts=len(aggregated_mutations_window[date][var][0])
        for counted in Counter(aggregated_mutations_window[date][var]).items(): #counted is country and counts
            if sequence_country_counts[date][counted[0]]==0:
                IR[date][var].append((counted[1])*jh_daily_cases.loc[counted[0],date])
            else:    
                IR[date][var].append(round((counted[1]/sequence_country_counts[date][counted[0]])*jh_daily_cases.loc[counted[0],date],2))
                
# CHANGE DATES
with open(f"../Input_files/IR1_011622.pkl",'wb') as file:
    pickle.dump(IR,file,protocol=pickle.HIGHEST_PROTOCOL)

