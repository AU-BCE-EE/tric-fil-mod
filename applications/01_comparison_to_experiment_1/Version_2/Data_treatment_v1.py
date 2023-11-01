# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 09:21:41 2023

@author: Mortensen
"""

##data treatment of experiment 1.1.1

#Import packages
import numpy as np
import csv
import math
import matplotlib.pyplot as plt

###input:

name='ex_1.1.1.csv' #name of file
#Calibration parameters
a=0.0082992
b=-0.0062214
#start and stop times
cycle_start=65
cycle_stop=520


#Initial arrays
time=H2S=[0.0]


#reading the file and saving time, mz35 and humid
#Remember to change the file name
with open(name) as file: 
    data=csv.reader(file,delimiter=';')
    header=next(data)
    for row in data:
        time_s= float(row[header.index('ï»¿Relative Time')])
        time = np.insert(time,0,time_s/3600)
        mz35 = float(row[header.index('m/z 35.00 ch3')])
        humid =float(row[header.index('37/21')])
        H2S = np.insert(H2S,0,float(mz35*(a*math.log(humid)+b))*10**-6*1.01325/(0.00008314*298)*34.08088) # calibration curve
#Be aware the the data is "reversed" so that time 0 is at index 860 (or similar for other data sets)

C_out=H2S[-cycle_stop:-cycle_start]

time_norm=time-time[-cycle_start]
t=time_norm[-cycle_stop:-cycle_start]

plt.plot(t,C_out)
#axes are time in s and concentration in ppm
plt.show


