# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:46:28 2024

@author: Mortensen
"""

#load packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#moving average over 5 data points on each side
window_size = 11

# Compare different pH values_____________________________________________________


ex1 = pd.read_csv('..//Output data/experimental6.1.csv', sep = ',')
model1 = pd.read_csv('..//Output data/model_6.1.csv', sep = ',')
ex1['Moving_Average'] = ex1['experimental'].rolling(window=window_size,center = True).mean()

ex2 = pd.read_csv('..//Output data/experimental6.2.csv', sep = ',')
model2 = pd.read_csv('..//Output data/model_6.2.csv', sep = ',')
ex2['Moving_Average'] = ex2['experimental'].rolling(window=window_size,center = True).mean()

ex3 = pd.read_csv('..//Output data/experimental6.3.csv', sep = ',')
model3 = pd.read_csv('..//Output data/model_6.3.csv', sep = ',')
ex3['Moving_Average'] = ex3['experimental'].rolling(window=window_size,center = True).mean()

ex4 = pd.read_csv('..//Output data/experimental6.4.csv', sep = ',')
model4 = pd.read_csv('..//Output data/model_6.4.csv', sep = ',')
ex4['Moving_Average'] = ex4['experimental'].rolling(window=window_size,center = True).mean()




#define parameters for breakthrough curve calculation
por_l = 0.23
por_g = 0.795115372 - por_l
BT1 = 14.5*por_g/27.2991
BT2 = 14.5*por_g/54.2766
BT1label='Theoretical Breakthrough curve 6.1 and 6.4)'
BT2label='Theoretical Breakthrough curve 6.2 and 6.3)'


plt.plot(ex1['experimental time (h)'],ex1['Moving_Average'],color='c',label='experimaltal 6.1')
plt.plot(ex2['experimental time (h)'],ex2['Moving_Average'],color='k',label='experimaltal 6.2')
plt.plot(ex3['experimental time (h)'],ex3['Moving_Average'],color='y',label='experimaltal 6.3')
plt.plot(ex4['experimental time (h)'],ex4['Moving_Average'],color='m',label='experimaltal 6.4')
plt.plot(model1['model time(h)'],model1['model'],color='c',linestyle='dashed', label='model 6.1')
plt.plot(model2['model time(h)'],model2['model'],color='k',linestyle='dashed', label='model 6.2')
plt.plot(model3['model time(h)'],model3['model'],color='y',linestyle='dashed', label='model 6.3')
plt.plot(model4['model time(h)'],model4['model'],color='m',linestyle='dashed', label='model 6.4')
plt.axvline(x=BT1/60,linestyle='-',label=BT1label) #breakthrough curve
plt.axvline(x=BT2/60,linestyle='-',label=BT2label)
plt.axhline(y=0.0596,color='g',label='Expected inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0,0.35)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Velocity comparison')
plt.savefig('..//Plots/Ex 6 velocity comparison.png', bbox_inches='tight')

#Compare the two replicates (5.4 and 6.1)_________________________________________








