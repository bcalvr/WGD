import numpy as np
import pandas as pd
from datetime import datetime as dt

tt1 = pd.read_csv('testfull.csv')
tt1.iloc[1,1] = 45678

tt1.to_csv('testfull1.csv')

print(f'Test script 2 complete at {dt.now()}.')