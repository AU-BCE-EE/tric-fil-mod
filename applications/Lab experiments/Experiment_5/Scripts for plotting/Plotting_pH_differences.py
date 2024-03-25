# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:46:28 2024

@author: Mortensen
"""

#load packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#Import data and make the moving average to avoid noise

window_size = 10

ex1 = pd.read_csv('..//Compare_pH_data/experimental5.1.csv', sep = ',')
model1 = pd.read_csv('..//Compare_pH_data/model_5.1.csv', sep = ',')
ex1['Moving_Average'] = ex1['experimental'].rolling(window=window_size).mean()

ex2 = pd.read_csv('..//Compare_pH_data/experimental5.2.csv', sep = ',')
model2 = pd.read_csv('..//Compare_pH_data/model_5.2.csv', sep = ',')
ex2['Moving_Average'] = ex2['experimental'].rolling(window=window_size).mean()

ex3 = pd.read_csv('..//Compare_pH_data/experimental5.3.csv', sep = ',')
model3 = pd.read_csv('..//Compare_pH_data/model_5.3.csv', sep = ',')
ex3['Moving_Average'] = ex3['experimental'].rolling(window=window_size).mean()

ex4 = pd.read_csv('..//Compare_pH_data/experimental5.4.csv', sep = ',')
model4 = pd.read_csv('..//Compare_pH_data/model_5.4.csv', sep = ',')
ex4['Moving_Average'] = ex4['experimental'].rolling(window=window_size).mean()

ex5 = pd.read_csv('..//Compare_pH_data/experimental5.5.csv', sep = ',')
model5 = pd.read_csv('..//Compare_pH_data/model_5.5.csv', sep = ',')
ex5['Moving_Average'] = ex5['experimental'].rolling(window=window_size).mean()

ex6 = pd.read_csv('..//Compare_pH_data/experimental5.6.csv', sep = ',')
model6 = pd.read_csv('..//Compare_pH_data/model_5.6.csv', sep = ',')
ex6['Moving_Average'] = ex6['experimental'].rolling(window=window_size).mean()

ex7 = pd.read_csv('..//Compare_pH_data/experimental5.7.csv', sep = ',')
model7 = pd.read_csv('..//Compare_pH_data/model_5.7.csv', sep = ',')
ex7['Moving_Average'] = ex7['experimental'].rolling(window=window_size).mean()





por_g = 0.77
BT = 14.5*por_g/25
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
plt.axhline(y=0.055,color='g',label='Expected inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0)
plt.ylim(0)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('pH comparison')
plt.savefig('..//Plots/pH comparison.png', bbox_inches='tight')





