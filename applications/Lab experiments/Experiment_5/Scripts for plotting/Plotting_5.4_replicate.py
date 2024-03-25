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


ex1 = pd.read_csv('..//Compare_pH_data/experimental5.4.csv', sep = ',')
model1 = pd.read_csv('..//Compare_pH_data/model_5.4.csv', sep = ',')
ex1['Moving_Average'] = ex1['experimental'].rolling(window=window_size).mean()

ex2 = pd.read_csv('..//Compare_pH_data/experimental6.1.csv', sep = ',')
model2 = pd.read_csv('..//Compare_pH_data/model_6.1.csv', sep = ',')
ex2['Moving_Average'] = ex2['experimental'].rolling(window=window_size).mean()



por_g = 0.77
BT = 14.5*por_g/25
BTlabel='Theoretical Breakthrough curve)'



plt.plot(ex1['experimental time (h)'],ex1['Moving_Average'],color='c',label='experimaltal 5.4')
plt.plot(ex2['experimental time (h)'],ex2['Moving_Average'],color='k',label='experimaltal 6.1')
plt.plot(model1['model time(h)'],model1['model'],color='c',linestyle='dashed', label='model 5.4')
plt.plot(model2['model time(h)'],model2['model'],color='k',linestyle='dashed', label='model 6.1')
plt.axvline(x=BT/60,linestyle='-',label=BTlabel) #breakthrough curve
plt.axhline(y=0.055,color='g',label='Expected inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0)
plt.ylim(0)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('replicate comparison')
plt.savefig('..//Plots/replicate comparison.png', bbox_inches='tight')





