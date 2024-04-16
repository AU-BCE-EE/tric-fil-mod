# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:00:48 2024

@author: Mortensen
"""

import numpy as np
import csv
import math


def Data_reciever (filename,startcycle,stopcycle,a,b):
    
    
    
    #Initial arrays
    time=H2S=[0.0]


    #reading the file and saving time, mz35 and humid
    #Remember to change the file name
    with open(filename) as file: 
        data=csv.reader(file,delimiter=';')
        header=next(data)
        for row in data:
            time_s= float(row[header.index('Relative Time')])
            time = np.insert(time,0,time_s/3600)
            mz35 = float(row[header.index('m/z 35.00 ch3')])
            humid =float(row[header.index('37/21')])
            H2S = np.insert(H2S,0,float(mz35*(a*math.log(humid)+b))*10**-6*1.01325/(0.00008314*298)*34.08088) # calibration curve and calculation from ppm to g/m3
    #Be aware the the data is "reversed" so that time 0 is at index 860 (or similar for other data sets)

    C_out=H2S[-stopcycle:-startcycle]

    time_norm=time-time[-startcycle]
    t=time_norm[-stopcycle:-startcycle]
    return{'time':t, 'concentration':C_out}





