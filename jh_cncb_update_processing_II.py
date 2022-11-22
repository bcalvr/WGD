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
import gc

# load required files - lagged FR
dates = pd.read_table('dates.txt',header=None,index_col=0)
old_date = dates.loc['old_date'].values[0]
new_date = dates.loc['new_date'].values[0]
date_idx = pd.read_csv('dateindex_splits.csv',names=['s','e'])

# LOAD ALL OUTPUTS YOU GENERATED
lagged_FR_init = {}

for i in range(len(date_idx)):
    with open(f"../Input_files/lagged_FR{i+1}_{new_date.replace('_','')}.pkl",'rb') as file:
        lagged_FR_init[i] = pickle.load(file)
        
list_dates=list(lagged_FR_init[0].keys())

lagged_FR=lagged_FR_init[0]

for i in range(len(date_idx)):
    for date in tqdm(list_dates[int(date_idx.iloc[i]['s']):int(date_idx.iloc[i]['e'])]):

        for var in lagged_FR[date].keys():
            lagged_FR[date][var]=lagged_FR_init[i][date][var]  

# If you have more than 2 partial files, keep adding them to the aggregate list - for example, supposing you have two more partial files lagged_FR3 and lagged_FR4:

# for date in tqdm(list_dates[717:730]):   # CHANGE DATES
#    for var in lagged_FR[date].keys():#\n",
#        lagged_FR[date][var]=lagged_FR3[date][var]

# for date in tqdm(list_dates[730:745]):  # CHANGE DATES
#     for var in lagged_FR[date].keys():#\n",
#         lagged_FR[date][var]=lagged_FR4[date][var]

del lagged_FR_init # delete all partial files
gc.collect()

# IR

# LOAD ALL OUTPUTS YOU GENERATED
IR_init = {}

for i in range(len(date_idx)):
    with open(f"../Input_files/IR{i+1}_{new_date.replace('_','')}.pkl",'rb') as file:
        IR_init[i] = pickle.load(file)

IR=IR_init[0]

for i in range(len(date_idx)):
    for date in tqdm(list_dates[int(date_idx.iloc[i]['s']):int(date_idx.iloc[i]['e'])]):
        for var in IR[date].keys():
            IR[date][var]=IR_init[i][date][var]
        
# SAME AS ABOVE - EXTEND TO ALL PARTIAL FILES AS NEEDED
        
del IR_init
gc.collect()

# PP

# LOAD ALL OUTPUTS YOU GENERATED
PP_init = {}

for i in range(len(date_idx)):
    with open(f"../Input_files/PP{i+1}_{new_date.replace('_','')}.pkl",'rb') as file:
        PP_init[i] = pickle.load(file)    

PP=PP_init[0]

# SAME AS ABOVE - EXTEND TO ALL PARTIAL FILES AS NEEDED

for date in tqdm(list_dates[int(date_idx.iloc[i]['s']):int(date_idx.iloc[i]['e'])]):
    for var in PP[date].keys():
        PP[date][var]=PP_init[i][date][var]

del PP_init
gc.collect()

# save 

# CHANGE DATES
with open(f"../Input_files/lagged_FR_{old_date}_{new_date}.pkl",'wb') as file:
    pickle.dump(lagged_FR,file,protocol=pickle.HIGHEST_PROTOCOL)
# CHANGE DATES
with open(f"../Input_files/IR_{old_date}_{new_date}.pkl",'wb') as file:
    pickle.dump(IR,file,protocol=pickle.HIGHEST_PROTOCOL)
# CHANGE DATES
with open(f"../Input_files/PP_{old_date}_{new_date}.pkl",'wb') as file:
    pickle.dump(PP,file,protocol=pickle.HIGHEST_PROTOCOL)

# If you just need to load the files before for the final table writing step, uncomment this cell and execute

# with open(f"../Input_files/lagged_FR_{old_date}_{new_date}.pkl",'rb') as file:
#    lagged_FR = pickle.load(file)
    
# with open(f"../Input_files/IR_{old_date}_{new_date}.pkl",'rb') as file:
#    IR = pickle.load(file)
    
# with open(f"../Input_files/PP_{old_date}_{new_date}.pkl",'rb') as file:
#    PP = pickle.load(file)
    
# list_dates=list(lagged_FR.keys())

with open(f"../Input_files/aggregated_mutations_window00_{dt.strftime(dt.strptime(old_date,'%m_%d_%y') + timedelta(days=1),'%m_%d_%y')}_{new_date}.pkl",'rb') as file:
    aggregated_mutations_window = pickle.load(file)

# Finally, assemble and save the output tables for GPR.

# CHANGE INDEXES - these will be the same as in jh_cncb_update_processing_I.py, line 58 (AFTER YOU UPDATED THE INDEXES TO CURRENT UPDATE)
for date in tqdm(list_dates[int(date_idx.iloc[0]['s']):int(date_idx.iloc[-1]['e'])]):
    unique_muts = pd.DataFrame(list(aggregated_mutations_window[date].keys()),columns=['descriptor'])
    unique_muts = pd.DataFrame.join(unique_muts,pd.DataFrame(unique_muts['descriptor'].str.split(',').to_list())) #assign columns with parsed descriptor
    unique_muts.set_index('descriptor',inplace=True)
    unique_muts[0] = pd.to_numeric(unique_muts[0])
    unique_muts[1] = pd.to_numeric(unique_muts[1])
    unique_muts.sort_values([0,1],inplace=True)
    unique_muts.columns = ['Start','End','Ref','Alt','VEP','Variant Type','Ref_AA','Alt_AA','AA']

    counts = {}
    num_countries = {}
    infection_rate = {}
    fatality_rate = {}
    counted_countries = {}

    for desc in unique_muts.index:
        counts[desc] = len(aggregated_mutations_window[date][desc])
        num_countries[desc] = len(set(aggregated_mutations_window[date][desc]))
        ir_sum=sum(IR[date][desc])
        pp_sum=sum(PP[date][desc])
        if pp_sum==0:
            infection_rate[desc] = ir_sum
        else:
            infection_rate[desc] = ir_sum/pp_sum #/counts[desc]
        #infection_rate[desc] = np.mean(aggregated_mutations_window[date][desc][1]) #since accumulated over countries before denomenator of mean is wrong
        if counts[desc]==0:
            fatality_rate[desc] = sum(lagged_FR[date][desc])
        else:
            fatality_rate[desc] = sum(lagged_FR[date][desc])/counts[desc]#/counts[desc]
        #fatality_rate[desc] = np.mean(aggregated_mutations_window[date][desc][2]) #since accumulated over countries before denomenator of mean is wrong
        counted_countries[desc] = dict(Counter(aggregated_mutations_window[date][desc]))

    unique_muts['counts'] = unique_muts.index.to_series().map(counts)
    unique_muts['countries'] = unique_muts.index.to_series().map(num_countries)
    unique_muts['infection_rate'] = unique_muts.index.to_series().map(infection_rate)
    unique_muts['fatality_rate'] = unique_muts.index.to_series().map(fatality_rate)
    unique_muts['counted_countries'] = unique_muts.index.to_series().map(counted_countries)
    unique_muts.index = unique_muts.index.str.split(',').str[0:4].str.join('_')
    unique_muts.sort_values('counts',ascending=False)
    os.mkdir(f'../Output Files/cumulative_daily_AF_v4_lag_{new_date}')
    # CREATE A FOLDER FOR THE OUTPUT IN ../Output Files: cumulative_daily_AF_v4_lag_<your date>
    # CHANGE FOLDER NAME BELOW TO THE FOLDER YOU CREATED
    unique_muts.to_csv(f"../Output Files/cumulative_daily_AF_v4_lag_{new_date}/{date.strftime('%m_%d_%y')}.csv")