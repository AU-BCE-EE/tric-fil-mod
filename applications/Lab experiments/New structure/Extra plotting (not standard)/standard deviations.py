# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:54:56 2024

@author: Mortensen
"""


import pandas as pd




first = '4'
second = '1'

file_path = '../../Master_spreadsheet.xlsx'
df = pd.read_excel(file_path, sheet_name = 'Data')

params = df[df['key'] == float(first)+0.1*float(second)]

pH1 = params['pH1'].values[0]
pH2 = params['pH2'].values[0]
pH3 = params['pH3'].values[0]
cycle1 = int(params['cycle1'].values[0])
cycle2 = int(params['cycle2'].values[0])
cycle3 = int(params['cycle3'].values[0])
cycle4 = params['cycle4'].values[0]
if not pd.isna(cycle4):
    cycle4=int(cycle4)
else:
    cycle4 = 'NaN'
cycle5 = int(params['cycle5'].values[0])
length = params['length'].values[0]
vol = params['volume'].values[0]
no = params['no'].values[0]

  # Inlet concentrations load in from data treatment and definition

ex1 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')

t1norm = ex1['Time in h'] - ex1['Time in h'][cycle1]
t1 = t1norm [cycle1:cycle1+length]
C1= ex1['Concentration in g/m^3'][cycle1:cycle1+length]

df = C1[251:451]
inlet_std = df.std()

print('inlet standard deviation: '+inlet_std)


#The background could not be imported. Instead, it was done in the document found by following this path
#'..//Special experiments (1,7 and 8)/Experiment_7/Scripts for plotting/plotting background.py





