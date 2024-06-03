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


ex1 = pd.read_csv('Output data/experimental6.1.csv', sep = ',')


ex2 = pd.read_csv('Output data/experimental5.4.csv', sep = ',')



diff_inlet = ex1['experimental inlet']-ex2['experimental inlet']
diff_outlet = ex1['experimental outlet']-ex2['experimental outlet']
time = ex1['experimental time (h)'] * 60
average = np.sum(diff_inlet[0:403])/len(diff_inlet[0:403])
print(average)
plt.plot(time,diff_outlet,label = 'outlet')
plt.plot(time,diff_inlet,label = 'inlet')
plt.xlabel('Time (min)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.grid(True)
plt.xlim(0,20)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Difference between repetitions on different days')
plt.show()


#Compare the two replicates (5.4 and 6.1)_________________________________________








