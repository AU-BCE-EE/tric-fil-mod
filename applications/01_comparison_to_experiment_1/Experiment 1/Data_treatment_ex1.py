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
from scipy.optimize import curve_fit




##Calibration parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#from experiment in lab
#start time
startcyclecali=20
endcyclecali=700

#initial guesses
params=[0.008,-0.006]

namecali='calibration_1.csv' #name of file


#Initial arrays
humid=correction=[0.0]

#actual concentration in ppm
C_actual=40


#reading the file and saving correction factor and humid
#Remember to change the file name
with open(namecali) as filecali: 
    datacali=csv.reader(filecali,delimiter=';')
    headercali=next(datacali)
    for row in datacali:
        mz35=float(row[headercali.index('ï»¿m/z 35.00 ch3')])
        correction= np.insert(correction,0,C_actual/mz35) #define correction factor as the ratio between actual and measured concentrations
        humid =np.insert(humid,0,float(row[headercali.index('37/21')])) #Define humidity as mz37 / mz 21. Make this calculation in excel from "raw siganl intensities"
#Note that the order is reversed so that the first obtained data is last in the array
x=humid[-endcyclecali:-startcyclecali]   
y=correction[-endcyclecali:-startcyclecali]     
        
def logcurve(h,A,B):#a and b is calibration parameters, h is humidity and c is correction factor
    n=h.size
    c=np.zeros(n)
    for i in range(n):
        c[i]=A*np.log(h[i])+B
    return c

#Making the regression
copt,ccov=curve_fit(logcurve,x,y,p0=params)

a=float(copt[0])
b=float(copt[1]) 

#data set number 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

name1='ex_1.1.1.csv' #name of file

#start and stop times
cycle_start1=65
cycle_stop1=520


#Initial arrays
time1=H2S1=[0.0]


#reading the file and saving time, mz35 and humid
#Remember to change the file name
with open(name1) as file1: 
    data1=csv.reader(file1,delimiter=';')
    header1=next(data1)
    for row in data1:
        time_s1= float(row[header1.index('ï»¿Relative Time')])
        time1 = np.insert(time1,0,time_s1/3600)
        mz351 = float(row[header1.index('m/z 35.00 ch3')])
        humid1 =float(row[header1.index('37/21')])
        H2S1 = np.insert(H2S1,0,float(mz351*(a*math.log(humid1)+b))*10**-6*1.01325/(0.00008314*298)*34.08088) # calibration curve and calculation from ppm to g/m3
#Be aware the the data is "reversed" so that time 0 is at index 860 (or similar for other data sets)

C_out1=H2S1[-cycle_stop1:-cycle_start1]

time_norm1=time1-time1[-cycle_start1]
t1=time_norm1[-cycle_stop1:-cycle_start1]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#data set number 2

name2='ex_1.2.1.csv' #file name

#start and stop times
cycle_start2=70
cycle_stop2=367


#Initial arrays
time2=H2S2=[0.0]


#reading the file and saving time, mz35 and humid
#Remember to change the file name
with open(name2) as file2: 
    data2=csv.reader(file2,delimiter=';')
    header2=next(data2)
    for row in data2:
        time_s2= float(row[header2.index('ï»¿Relative Time')])
        time2 = np.insert(time2,0,time_s2/3600)
        mz352 = float(row[header2.index('m/z 35.00 ch3')])
        humid2 =float(row[header2.index('37/21')])
        H2S2 = np.insert(H2S2,0,float(mz352*(a*math.log(humid2)+b))*10**-6*1.01325/(0.00008314*298)*34.08088) # calibration curve
#Be aware the the data is "reversed" so that time 0 is at index 860 (or similar for other data sets)

C_out2=H2S2[-cycle_stop2:-cycle_start2]

time_norm2=time2-time2[-cycle_start2]
t2=time_norm2[-cycle_stop2:-cycle_start2]



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#data set number 3

name3='ex_1.3.1.csv' #file name

#start and stop times
cycle_start3=80
cycle_stop3=470


#Initial arrays
time3=H2S3=[0.0]


#reading the file and saving time, mz35 and humid
#Remember to change the file name
with open(name3) as file3: 
    data3=csv.reader(file3,delimiter=';')
    header3=next(data3)
    for row in data3:
        time_s3= float(row[header3.index('ï»¿Relative Time')])
        time3 = np.insert(time3,0,time_s3/3600)
        mz353 = float(row[header3.index('m/z 35.00 ch3')])
        humid3 =float(row[header3.index('37/21')])
        H2S3 = np.insert(H2S3,0,float(mz353*(a*math.log(humid3)+b))*10**-6*1.01325/(0.00008314*298)*34.08088) # calibration curve
#Be aware the the data is "reversed" so that time 0 is at index 860 (or similar for other data sets)

C_out3=H2S3[-cycle_stop3:-cycle_start3]

time_norm3=time3-time3[-cycle_start3]
t3=time_norm3[-cycle_stop3:-cycle_start3]


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#data set number 4

name4='ex_1.4.1.csv' #file name

#start and stop times
cycle_start4=80
cycle_stop4=450


#Initial arrays
time4=H2S4=[0.0]


#reading the file and saving time, mz35 and humid
#Remember to change the file name
with open(name4) as file4: 
    data4=csv.reader(file4,delimiter=';')
    header4=next(data4)
    for row in data4:
        time_s4= float(row[header4.index('ï»¿Relative Time')])
        time4 = np.insert(time4,0,time_s4/3600)
        mz354 = float(row[header4.index('m/z 35.00 ch3')])
        humid4 =float(row[header4.index('37/21')])
        H2S4 = np.insert(H2S4,0,float(mz354*(a*math.log(humid4)+b))*10**-6*1.01325/(0.00008314*298)*34.08088) # calibration curve
#Be aware the the data is "reversed" so that time 0 is at index 860 (or similar for other data sets)

C_out4=H2S4[-cycle_stop4:-cycle_start4]

time_norm4=time4-time4[-cycle_start4]
t4=time_norm4[-cycle_stop4:-cycle_start4]




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#data set number 5

name5='ex_1.5.1.csv' #file name

#start and stop times
cycle_start5=90
cycle_stop5=510


#Initial arrays
time5=H2S5=[0.0]


#reading the file and saving time, mz35 and humid
#Remember to change the file name
with open(name5) as file5: 
    data5=csv.reader(file5,delimiter=';')
    header5=next(data5)
    for row in data5:
        time_s5= float(row[header5.index('ï»¿Relative Time')])
        time5 = np.insert(time5,0,time_s5/3600)
        mz355 = float(row[header5.index('m/z 35.00 ch3')])
        humid5 =float(row[header5.index('37/21')])
        H2S5 = np.insert(H2S5,0,float(mz355*(a*math.log(humid5)+b))*10**-6*1.01325/(0.00008314*298)*34.08088) # calibration curve
#Be aware the the data is "reversed" so that time 0 is at index 860 (or similar for other data sets)

C_out5=H2S5[-cycle_stop5:-cycle_start5]

time_norm5=time5-time5[-cycle_start5]
t5=time_norm5[-cycle_stop5:-cycle_start5]


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#data set number 6

name6='ex_1.6.1.csv' #file name

#start and stop times
cycle_start6=100
cycle_stop6=460


#Initial arrays
time6=H2S6=[0.0]


#reading the file and saving time, mz35 and humid
#Remember to change the file name
with open(name6) as file6: 
    data6=csv.reader(file6,delimiter=';')
    header6=next(data6)
    for row in data6:
        time_s6= float(row[header6.index('ï»¿Relative Time')])
        time6 = np.insert(time6,0,time_s6/3600)
        mz356 = float(row[header6.index('m/z 35.00 ch3')])
        humid6 =float(row[header6.index('37/21')])
        H2S6 = np.insert(H2S6,0,float(mz356*(a*math.log(humid6)+b))*10**-6*1.01325/(0.00008314*298)*34.08088) # calibration curve
#Be aware the the data is "reversed" so that time 0 is at index 860 (or similar for other data sets)

C_out6=H2S6[-cycle_stop6:-cycle_start6]

time_norm6=time6-time6[-cycle_start6]
t6=time_norm6[-cycle_stop6:-cycle_start6]
