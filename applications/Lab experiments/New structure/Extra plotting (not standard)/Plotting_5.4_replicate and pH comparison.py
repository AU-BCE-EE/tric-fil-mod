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


ex1 = pd.read_csv('..//Output data/experimental5.1.csv', sep = ',')
model1 = pd.read_csv('..//Output data/model_5.1.csv', sep = ',')
ex1['Moving_Average'] = ex1['experimental'].rolling(window=window_size,center = True).mean()

ex2 = pd.read_csv('..//Output data/experimental5.2.csv', sep = ',')
model2 = pd.read_csv('..//Output data/model_5.2.csv', sep = ',')
ex2['Moving_Average'] = ex2['experimental'].rolling(window=window_size,center = True).mean()

ex3 = pd.read_csv('..//Output data/experimental5.3.csv', sep = ',')
model3 = pd.read_csv('..//Output data/model_5.3.csv', sep = ',')
ex3['Moving_Average'] = ex3['experimental'].rolling(window=window_size,center = True).mean()

ex4 = pd.read_csv('..//Output data/experimental5.4.csv', sep = ',')
model4 = pd.read_csv('..//Output data/model_5.4.csv', sep = ',')
ex4['Moving_Average'] = ex4['experimental'].rolling(window=window_size,center = True).mean()

ex5 = pd.read_csv('..//Output data/experimental5.5.csv', sep = ',')
model5 = pd.read_csv('..//Output data/model_5.5.csv', sep = ',')
ex5['Moving_Average'] = ex5['experimental'].rolling(window=window_size,center = True).mean()

ex6 = pd.read_csv('..//Output data/experimental5.6.csv', sep = ',')
model6 = pd.read_csv('..//Output data/model_5.6.csv', sep = ',')
ex6['Moving_Average'] = ex6['experimental'].rolling(window=window_size,center = True).mean()

ex7 = pd.read_csv('..//Output data/experimental5.7.csv', sep = ',')
model7 = pd.read_csv('..//Output data/model_5.7.csv', sep = ',')
ex7['Moving_Average'] = ex7['experimental'].rolling(window=window_size,center = True).mean()


#define parameters for breakthrough curve calculation
por_l = 0.20498
por_g = 0.795115372 - por_l
BT = 14.5*por_g/27.2991
BTlabel='Theoretical Breakthrough curve)'


plt.plot(ex7['experimental time (h)'],ex7['Moving_Average'],color='y',label='experimaltal pH=5.99')
plt.plot(ex1['experimental time (h)'],ex1['Moving_Average'],color='c',label='experimaltal pH=7.08')
plt.plot(ex2['experimental time (h)'],ex2['Moving_Average'],color='k',label='experimaltal pH=7.27')
plt.plot(ex5['experimental time (h)'],ex5['Moving_Average'],color='b',label='experimaltal pH=7.55')
plt.plot(ex3['experimental time (h)'],ex3['Moving_Average'],color='r',label='experimaltal pH=7.76')
plt.plot(ex4['experimental time (h)'],ex4['Moving_Average'],color='m',label='experimaltal pH=8.02')
plt.plot(model7['model time(h)'],model7['model'],color='y',linestyle='dashed', label='model pH=5.99')
plt.plot(model1['model time(h)'],model1['model'],color='c',linestyle='dashed', label='model pH=7.08')
plt.plot(model2['model time(h)'],model2['model'],color='k',linestyle='dashed', label='model pH=7.27')
plt.plot(model5['model time(h)'],model5['model'],color='b',linestyle='dashed', label='model pH=7.55')
plt.plot(model3['model time(h)'],model3['model'],color='r',linestyle='dashed', label='model pH=7.76')
plt.plot(model4['model time(h)'],model4['model'],color='m',linestyle='dashed', label='model pH=8.02')
plt.axvline(x=BT/60,linestyle='-',label=BTlabel) #breakthrough curve
plt.axhline(y=0.0596,color='g',label='Expected inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0,0.35)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('pH comparison')
plt.savefig('..//Plots/Ex 5 pH comparison.png', bbox_inches='tight')

#Compare the two replicates (5.4 and 6.1)_________________________________________


#Import data and make the moving average to avoid noise


#import data

ex8 = pd.read_csv('..//Output data/experimental6.1.csv', sep = ',')
model8 = pd.read_csv('..//Output data/model_6.1.csv', sep = ',')
ex8['Moving_Average'] = ex8['experimental'].rolling(window=window_size).mean()



# plot
plt.clf()
plt.plot(ex4['experimental time (h)'],ex4['Moving_Average'],color='c',label='experimaltal 5.4')
plt.plot(ex8['experimental time (h)'],ex8['Moving_Average'],color='k',label='experimaltal 6.1')
plt.plot(model4['model time(h)'],model4['model'],color='c',linestyle='dashed', label='model 5.4')
plt.plot(model8['model time(h)'],model8['model'],color='k',linestyle='dashed', label='model 6.1')
plt.axvline(x=BT/60,linestyle='-',label=BTlabel) #breakthrough curve
plt.axhline(y=0.0596,color='g',label='Expected inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0,0.35)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('replicate comparison')
plt.savefig('..//Plots/5.4 compared to 6.1.png', bbox_inches='tight')





