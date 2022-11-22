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
dates = pd.read_table('dates.txt',header=None,index_col=0)
old_date = dates.loc['old_date'].values[0]
new_date = dates.loc['new_date'].values[0]

with open(f"../Input_files/jh_daily_ir_{new_date}.pkl",'rb') as file: #Contains only country name for each date
    jh_daily_ir = pickle.load(file)
with open(f"../Input_files/jh_data_{new_date}.pkl",'rb') as file:
    jh_data=pickle.load(file)
with open(f"../Input_files/jh_daily_cases_{new_date}.pkl",'rb') as file:
    jh_daily_cases = pickle.load(file)
with open(f"../Input_files/seq_country_counts_{new_date}.pkl",'rb') as file:
    sequence_country_counts = pickle.load(file)
with open(f"../Input_files/basis_vector_IR_FR_PP_{new_date}.pkl",'rb') as file:
    IR = pickle.load(file)
with open(f"../Input_files/aggregated_mutations_window00_{dt.strftime(dt.strptime(old_date,'%m_%d_%y') + timedelta(days=1),'%m_%d_%y')}_{new_date}.pkl",'rb') as file:
    aggregated_mutations_window = pickle.load(file)

# For this step, create a few parallel scripts, each cycling through ~15 days. For example, if your update covers 30 days, you would create 2 scripts for each IR, PP and lagged_FR.
# Get the indexes corresponding to days by looking up list_dates as before.

# CHANGE INDEXES

#Looping over single days
list_dates=jh_daily_ir.columns[9::]
for date in tqdm(list_dates[int(float(sys.argv[1])):int(float(sys.argv[2]))]):
#Looping over single variants
    for var in aggregated_mutations_window[date].keys():
        #For each variant count countries,
        #totalCounts=len(aggregated_mutations_window[date][var][0])
        for counted in Counter(aggregated_mutations_window[date][var]).items(): #counted is country and counts
            if sequence_country_counts[date][counted[0]]==0:
                IR[date][var].append((counted[1])*jh_daily_cases.loc[counted[0],date])
            else:
                IR[date][var].append(round((counted[1]/sequence_country_counts[date][counted[0]])*jh_daily_cases.loc[counted[0],date],2))

# CHANGE DATE
with open(f"../Input_files/IR2_{new_date.replace('_','')}.pkl",'wb') as file:
    pickle.dump(IR,file,protocol=pickle.HIGHEST_PROTOCOL)
