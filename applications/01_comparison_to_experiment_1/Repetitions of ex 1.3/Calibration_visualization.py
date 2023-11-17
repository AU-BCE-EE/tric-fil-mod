# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 10:19:13 2023

@author: bruger
"""

##data treatment of experiment 1.3 repetitions

#Import packages
import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#start time
startcyclecali=20
endcyclecali=700

#initial guesses
params=[0.008,-0.006]

#data set number 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
        
def logcurve(h,a,b):#a and b is calibration parameters, h is humidity and c is correction factor
    n=h.size
    c=np.zeros(n)
    for i in range(n):
        c[i]=a*np.log(h[i])+b
    return c

#Making the regression
copt,ccov=curve_fit(logcurve,x,y,p0=params)
#plotting
yfit=logcurve(x,copt[0],copt[1])
plt.plot(x,y,'o')
plt.plot(x,yfit,label=('y= %1.5f'%copt[0]+'*ln(x) %1.5f'%copt[1]))
plt.xlabel('m/z 37 / m/z 21')
plt.ylabel('correction factor []')
plt.title('Calibration Curve')
plt.legend(loc='lower right')
plt.show
A=float(copt[0])
B=float(copt[1])
