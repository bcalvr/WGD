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
import time

cell_time = time.time()

dates = pd.read_table('dates.txt',header=None,index_col=0)
old_date = dates.loc['old_date'].values[0]
new_date = dates.loc['new_date'].values[0]

date_start = '1/23/20'
date_up_to = new_date.replace('_','/').lstrip('0')
jh_path = '/gpfs/group/balch/data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/'
confirmed_global = 'time_series_covid19_confirmed_global.csv'
deaths_global = 'time_series_covid19_deaths_global.csv'

confirmed = pd.read_csv(jh_path+confirmed_global)
deaths = pd.read_csv(jh_path+deaths_global)

#Need to adjust country names to match those in CNCB
confirmed.loc[confirmed['Country/Region']=='Taiwan*','Country/Region'] = 'Taiwan'
confirmed.loc[confirmed['Country/Region']=='US','Country/Region'] = 'United States'
confirmed.loc[confirmed['Country/Region']=='Congo (Kinshasa)','Country/Region'] = 'Democratic Republic of the Congo'
confirmed.loc[confirmed['Country/Region']=='Korea, South','Country/Region'] = 'South Korea'
confirmed.loc[confirmed['Country/Region']=='Czechia','Country/Region'] = 'Czech Republic'
confirmed.loc[confirmed['Country/Region']=='Burma','Country/Region'] = 'myanmar'
confirmed.loc[confirmed['Country/Region']=='Congo (Brazzaville)','Country/Region'] = 'republic of the congo'
confirmed.loc[confirmed['Country/Region']=="Cote d'Ivoire",'Country/Region'] = 'cotedivoire'

confirmed['Country/Region'] = confirmed['Country/Region'].str.lower()
confirmed.set_index('Country/Region',inplace=True)

deaths.loc[deaths['Country/Region']=='Taiwan*','Country/Region'] = 'Taiwan'
deaths.loc[deaths['Country/Region']=='US','Country/Region'] = 'United States'
deaths.loc[deaths['Country/Region']=='Congo (Kinshasa)','Country/Region'] = 'Democratic Republic of the Congo'
deaths.loc[deaths['Country/Region']=='Korea, South','Country/Region'] = 'South Korea'
deaths.loc[deaths['Country/Region']=='Czechia','Country/Region'] = 'Czech Republic'
deaths.loc[deaths['Country/Region']=='Burma','Country/Region'] = 'myanmar'
deaths.loc[deaths['Country/Region']=='Congo (Brazzaville)','Country/Region'] = 'republic of the congo'
deaths.loc[deaths['Country/Region']=="Cote d'Ivoire",'Country/Region'] = 'cotedivoire'

deaths['Country/Region'] = deaths['Country/Region'].str.lower()
deaths.set_index('Country/Region',inplace=True)

print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())

cell_time = time.time()

with open(f"../Input_files/confirmed_{new_date}.pkl",'wb') as file:
    pickle.dump(confirmed,file,protocol=pickle.HIGHEST_PROTOCOL)
    
with open(f"../Input_files/deaths_{new_date}.pkl",'wb') as file:
    pickle.dump(deaths,file,protocol=pickle.HIGHEST_PROTOCOL)
    
print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())

cell_time = time.time()

# take latest daily report in update range - in JH data
date = f"{dt.strftime(dt.strptime(new_date,'%m_%d_%y'),'%m-%d-%Y')}.csv"
jh_path = '/gpfs/group/balch/data/COVID-19/csse_covid_19_data/csse_covid_19_daily_reports/'
jh_data = pd.read_csv(jh_path+date)

#drop rows that have nans for incidence rate
jh_data = jh_data[~jh_data['Incident_Rate'].isnull()]
jh_data = jh_data[jh_data['Incident_Rate']!=0]

#matching countries of jh data to that of cncb countries
jh_data.loc[jh_data['Country_Region']=='Taiwan*','Country_Region'] = 'Taiwan'
jh_data.loc[jh_data['Country_Region']=='US','Country_Region'] = 'United States'
jh_data.loc[jh_data['Country_Region']=='Korea, South','Country_Region'] = 'South Korea'
jh_data.loc[jh_data['Country_Region']=='Czechia','Country_Region'] = 'Czech Republic'
jh_data.loc[jh_data['Country_Region']=='Burma','Country_Region'] = 'myanmar'
jh_data.loc[jh_data['Country_Region']=='Congo (Kinshasa)','Country_Region'] = 'Democratic Republic of the Congo'
jh_data.loc[jh_data['Country_Region']=='Congo (Brazzaville)','Country_Region'] = 'republic of the congo'
jh_data.loc[jh_data['Country_Region']=="Cote d'Ivoire",'Country_Region'] = 'cotedivoire'

# Population is back calculated through the incident rate number: 
# population = confirmed_cases/incident_rate_per100k * 100,000 
jh_data['Population'] = np.ceil(jh_data['Confirmed']/jh_data['Incident_Rate']*100000)
jh_data['Country_Region'] = jh_data['Country_Region'].str.lower()

print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())

cell_time = time.time()

with open(f"../Input_files/jh_data_{new_date}.pkl",'wb') as file:
    pickle.dump(jh_data,file,protocol=pickle.HIGHEST_PROTOCOL)

print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())

cell_time = time.time()

#cncb meta, need to match countries
# previous issues seem to be fixed (previously: 1) `Host instead of Host; 2) Submission date and Submitting lab swapped)
column_names = ['Virus Strain Name','Accession ID','Related ID','Nuc.Completeness','Sequence Quality','Host','Location','Sample Collection Date','Submitting Lab']

# UPLOAD NEW META FILE FIRST (FROM CNCB) - name it accordingly and replace the filename below with the new meta filename
meta=pd.read_table(f'../Input_files/metadata_{new_date}.tsv',usecols=column_names,index_col=1) 

meta.fillna(' ',inplace=True)

#Remove low quality and partial reads, pretty sure cncb does not run variant annotation for these anyways
meta = meta[meta['Sequence Quality']!='Low']
meta = meta[meta['Nuc.Completeness']!='Partial']

#set country column and lowercase
meta['Country'] = meta['Location'].str.split('/').str[0].str.strip()
meta['Country'] = meta['Country'].str.lower()

#Adjust country typos

meta.loc[meta['Country']=='?romania','Country'] = 'romania'
meta.loc[meta['Country']=='viet nam','Country'] = 'vietnam'
meta.loc[meta['Country']=='czech repubic','Country'] = 'czech republic'
meta.loc[meta['Country']=='ivory coast','Country'] = 'cotedivoire'

meta.loc[meta['Country']=='blegium','Country'] = 'belgium'
meta.loc[meta['Country']=='frnace','Country'] = 'france'
meta.loc[meta['Country']=='morocoo','Country'] = 'morocco'
meta.loc[meta['Country']=='mosambik','Country'] = 'mozambique'
meta.loc[meta['Country']=='congo','Country'] = 'democratic republic of the congo'
meta.loc[meta['Country']=='guinea bissau','Country'] = 'guinea-bissau'
meta.loc[meta['Country']=='the republic of guinea-bissau','Country'] = 'guinea-bissau'

# #Remove Crimea and Palestine
meta = meta[meta['Country']!='crimea']
meta = meta[meta['Country']!='palestine']
meta = meta[meta['Country']!='saint martin']
meta = meta[meta['Country']!='sint maarten']
meta = meta[meta['Country']!='french guiana']
meta = meta[meta['Country']!='guadeloupe']
meta = meta[meta['Country']!='martinique']
meta = meta[meta['Country']!='gibraltar']
meta = meta[meta['Country']!='northern mariana islands']


#Filter using date as well
meta = meta[meta['Sample Collection Date']!='2020-00-00'] #bad dates in there
meta.loc[:,'Sample Collection Date'] = pd.to_datetime(meta['Sample Collection Date'], yearfirst=True)
dt_date_up_to = dt.strptime(date_up_to, "%m/%d/%y")
dt_date_start = dt.strptime(date_start, "%m/%d/%y")
meta = meta[(meta.loc[:,'Sample Collection Date']<dt_date_up_to) & (meta.loc[:,'Sample Collection Date']>dt_date_start)]

print(meta.shape)

print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())

cell_time = time.time()

# RUN THE BELOW AND MAKE SURE YOU GET set() AS RESULT
# If not, remove non-matching countries in meta until you get an empty set difference.  
print(set(meta['Country'])-set(confirmed.index.str.lower()))
print(len(set(meta['Country'])))

print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())

cell_time = time.time()

meta.to_pickle(f"../Input_files/meta_{new_date}.pkl",protocol=pickle.HIGHEST_PROTOCOL)

print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())

cell_time = time.time()

jh_daily_ir = pd.DataFrame(columns=pd.to_datetime(confirmed.columns[4::])) #start from the second date because we take difference between dates
jh_daily_fr = pd.DataFrame(columns=pd.to_datetime(confirmed.columns[4::]))
country_daily_cts = pd.DataFrame(columns=pd.to_datetime(confirmed.columns[4::]))
country_death_cts = pd.DataFrame(columns=pd.to_datetime(confirmed.columns[4::]))

jh_daily_cases = pd.DataFrame(columns=pd.to_datetime(confirmed.columns[4::])) #start from the second date because we take difference between dates
jh_daily_deaths = pd.DataFrame(columns=pd.to_datetime(confirmed.columns[4::]))

jh_country_lvl = pd.DataFrame(columns=['cases','deaths','population'])

for country in set(meta['Country'].str.lower()):
    population = jh_data.loc[jh_data['Country_Region']==country,'Population'].sum()
    summed_cases = confirmed.loc[confirmed.index==country,date_up_to].sum()
    summed_deaths = deaths.loc[deaths.index==country,date_up_to].sum()
    jh_country_lvl.loc[country] = [summed_cases,summed_deaths,population]
    
    cases_timeperiod = 1+confirmed.loc[confirmed.index==country].iloc[:,4::].sum().values
    deaths_timeperiod = deaths.loc[deaths.index==country].iloc[:,4::].sum().values
    jh_daily_ir.loc[country] = cases_timeperiod / population * 100000
    jh_daily_fr.loc[country] = np.nan_to_num(deaths_timeperiod / cases_timeperiod) * 100
    country_daily_cts.loc[country] = (confirmed.loc[[country]].iloc[:,4::].sum().values-confirmed.loc[[country]].iloc[:,3:-1].sum().values).clip(min=0)
    country_death_cts.loc[country] = (deaths.loc[[country]].iloc[:,4::].sum().values-deaths.loc[[country]].iloc[:,3:-1].sum().values).clip(min=0)

    jh_daily_cases.loc[country] = confirmed.loc[confirmed.index==country].iloc[:,4::].sum().values
    jh_daily_deaths.loc[country] = np.nan_to_num(deaths_timeperiod)
    
jh_country_lvl.sort_index(inplace=True)
jh_daily_ir.sort_index(inplace=True)
jh_daily_fr.sort_index(inplace=True)
jh_daily_cases.sort_index(inplace=True)
jh_daily_deaths.sort_index(inplace=True)

print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())

cell_time = time.time()

jh_daily_ir.to_csv(f"../Output Files/daily_counts/daily_ir_{dt.today().strftime('%m_%d_%y')}.csv")
jh_daily_fr.to_csv(f"../Output Files/daily_counts/daily_fr_{dt.today().strftime('%m_%d_%y')}.csv")
country_daily_cts.to_csv(f"../Output Files/daily_counts/daily_cases_{dt.today().strftime('%m_%d_%y')}.csv")
country_death_cts.to_csv(f"../Output Files/daily_counts/daily_deaths_{dt.today().strftime('%m_%d_%y')}.csv")

print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())

cell_time = time.time()

jh_daily_ir.to_pickle(f"../Input_files/jh_daily_ir_{new_date}.pkl",protocol=pickle.HIGHEST_PROTOCOL)

print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())

cell_time = time.time()

jh_daily_cases.to_pickle(f"../Input_files/jh_daily_cases_{new_date}.pkl",protocol=pickle.HIGHEST_PROTOCOL)

print('---Cell run time:',round(np.floor((time.time() - cell_time)/60)),'minutes,',round((time.time() - cell_time)%60),'seconds---')
print("Cell completed:",dt.now())