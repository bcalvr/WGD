import pandas as pd
import numpy as np
import pickle as pickle
from datetime import datetime as dt

def index_splitter():
    # load required files
    dates = pd.read_table('dates.txt',header=None,index_col=0)
    old_date = dates.loc['old_date'].values[0]
    new_date = dates.loc['new_date'].values[0]

    with open(f"../Input_files/jh_daily_ir_{new_date}.pkl",'rb') as file: #Contains only country name for each date
        jh_daily_ir = pickle.load(file)

    # Print date start, end, and length
    list_dates=jh_daily_ir.columns[9::]
    old_idx = np.where(list_dates==dt.strptime(old_date,'%m_%d_%y'))[0][0]
    new_idx = np.where(list_dates==dt.strptime(new_date,'%m_%d_%y'))[0][0]
    idx_len = new_idx - old_idx + 1

    print(f"The old date ({old_date}) is at index {old_idx}, the new date ({new_date}) is at index {new_idx}, for a total index length of {idx_len}, requiring {int(np.floor(idx_len/15))} runs of 15 dates and one of {idx_len%15} dates.")

    # Output date range indices
    index_pairs = [tuple([old_idx + i*15,old_idx + (i+1)*15]) for i in range(int(np.floor(idx_len/15)))]
    index_pairs.append(tuple([old_idx + int(np.floor(idx_len/15))*15,new_idx + 1]))
    return index_pairs

np.savetxt('dateindex_splits.csv',index_splitter(),delimiter =',')

# print(index_splitter())