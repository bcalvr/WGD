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

with open(f"../Input_files/meta_{new_date}.pkl",'rb') as file:
    meta = pickle.load(file)

date_range = jh_daily_ir.columns[9::] #skip the first 9 (feb 1), arbitrarily since not too much data yet
sequence_country_counts = {date:dict.fromkeys(set(meta['Country'].str.lower()),0) for date in date_range} #Dictionary of dates, sub dictionary country names with counts

for date in tqdm(date_range):
    for country in set(meta['Country'].str.lower()):
        sequence_country_counts[date][country] = meta[(meta['Sample Collection Date']<=date) & (meta['Country']==country)].shape[0]


with open(f"../Input_files/seq_country_counts_{new_date}.pkl",'wb') as file:
    pickle.dump(sequence_country_counts,file,protocol=pickle.HIGHEST_PROTOCOL)
