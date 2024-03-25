# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:00:48 2024

@author: Mortensen
"""

import numpy as np
import pandas as pd
import sys

#specify excel file and destination
filename = 'Background'

excel_file_path = '..//Raw_data/'+filename+'.xlsx'
csv_destination_path = '..//Processed_data/'+filename+'.csv'



#defining columns of interest 
sheetname1 = 'Time   Cycle'
columnname1 = 'Relative Time'
sheetname2 = 'Raw signal intensities'
columnname21 = 'm/z 37.00 ch8'
columnname22 = 'm/z 21.00 ch1'
sheetname3 = 'Concentration'
columnname31 = 'm/z 35.00 ch7'
columnname32 = 'm/z 49.00 ch9'
columnname33 = 'm/z 63.00 ch10'


#Loading the columns of interest
df1 = pd.read_excel(excel_file_path,sheet_name=sheetname1)
df2 = pd.read_excel(excel_file_path,sheet_name=sheetname2)
df3 = pd.read_excel(excel_file_path,sheet_name=sheetname3)


#defining the raw data
t_s = df1[columnname1] #time in s
mz37 = df2[columnname21] # m/z 37 signal (water cluster)
mz21 = df2[columnname22] # m/z 21 signal (H3O+)
mz35 = df3[columnname31]  #m/z 35 signal (H2S)
mz49 = df3[columnname32] #m/z 49 signal (methyl mercaptan)
mz63 = df3[columnname33] #m/z 63 signal (dimethyl sulfide)


#modifying and correction the raw data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Dimethyl sulfide and methyl mercaptan calibrations
#read the inlet files. Only the "concentration" sheet and m/z 49 and 63 is needed

#Methyl mercaptan

df4 = pd.read_excel('..//Raw_data/experiment_7.1.inlet1.xlsx',sheet_name=sheetname3)
df5 = pd.read_excel('..//Raw_data/experiment_7.1.inlet2.xlsx',sheet_name=sheetname3)
#load raw data from the sheet "Concentration"

inlet149 = df4.loc[[101,233],columnname32]
inlet249 = df5.loc[[57,189],columnname32]
#Choose m/z 49 (columnname32) and the right indexes (start of experiment from logbook is index 61 
#and 17 respectively + 40cycles to reach equilibrium)
average49 = (inlet149.mean()+inlet249.mean())/2 #find the average of all the chosen data points
correction49 = 100/average49 # known concentration is 100ppb (0.5L/min 5ppm + 24.5l/min air)
mm = correction49 * mz49 * 10**-9 * 1.01325 / (0.00008314*298) * 48 # concentration of methyl mercaptan in g/m^3 (started in ppb)

#dimethyl sulfide
inlet163 = df4.loc[[101,233],columnname33]
inlet263 = df5.loc[[57,189],columnname33]
#Choose m/z 49 (columnname33) and the right indexes (start of experiment from logbook is index 61 
#and 17 respectively + 40cycles to reach equilibrium)
average63 = (inlet163.mean()+inlet263.mean())/2 #finding the average
correction63 = 100/average63 # known concentration is 100ppb (0.5L/min 5ppm + 24.5l/min air)
dms = correction63 * mz63 * 10**-9 * 1.01325 / (0.00008314*298) * 62 # concentration of dimethyl sulfide in g/m^3 (started in ppb)

#H2S humidity calibration
from Calibration_func import a,b,maxhumid    #import calibration parameters from calibration script


time =  t_s / 3600 #time in h
humid = mz37 / mz21
for i in range(len(humid)):
    if humid[i] > maxhumid:
        sys.exit('Error, the humidity is too high for this calibration')
correction = a * np.log(humid) + b # correction factor, correcting for the amount of water

H2S = mz35 * correction * 10**-6 * 1.01325 / (0.00008314*298) * 34.08088 # concentration of H2S in g/m^3

results = pd.DataFrame ({'H2S concentration in g/m^3':H2S ,'DMS concentration in g/m^3':dms,'mm concentration in g/m^3':mm, 'Time in h':time})

results.to_csv(csv_destination_path, header = True, index=False)



    






