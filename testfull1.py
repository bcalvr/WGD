import pandas as pd
from datetime import datetime as dt
import time

tt = pd.DataFrame([[1,2,3,7],[0,5,6,3]])

time.sleep(5)

tt.to_csv('testfull.csv')

print(f"Test script 1 complete at {dt.now()}.")