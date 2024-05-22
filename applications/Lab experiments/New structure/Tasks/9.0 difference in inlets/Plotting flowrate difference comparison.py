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


ex1 = pd.read_csv('Output data/inlet_5.7.csv', sep = ',')


ex2 = pd.read_csv('Output data/inlet_5.4.csv', sep = ',')

diff = ex1['inlet']-ex2['inlet']
time = ex1['time(h)']

plt.plot(time*60,diff)
plt.xlabel('Time (min)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.grid(True)
plt.xlim(0,20)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Inlet difference')
plt.show()


#Compare the two replicates (5.4 and 6.1)_________________________________________








