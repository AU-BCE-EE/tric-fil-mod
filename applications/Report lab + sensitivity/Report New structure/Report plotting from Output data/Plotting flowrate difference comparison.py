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
BT1label='Theoretical breakthrough curve setting 1 and 4)'
BT2label='Theoretical breakthrough curve setting 2 and 3)'

# Define y values for the expected inlet curve @40ppm
expected_inlet = np.full_like(ex1['experimental time (h)']*60, 0.0596)  

# Make the line drop to 0 at x=5, as the oulse is 5min
expected_inlet[ex1['experimental time (h)']*60 >= 5] = 0


plt.plot(ex1['experimental time (h)']*60,ex1['Moving_Average'],color='c',label='Measured outlet at setting 1')
plt.plot(ex2['experimental time (h)']*60,ex2['Moving_Average'],color='k',label='Measured outlet at setting 2')
plt.plot(ex3['experimental time (h)']*60,ex3['Moving_Average'],color='y',label='Measured outlet at setting 3')
plt.plot(ex4['experimental time (h)']*60,ex4['Moving_Average'],color='m',label='Measured outlet at setting 4')
plt.plot(ex1['experimental time (h)']*60, expected_inlet, color='g', label='Expected inlet concentration')
plt.plot(model1['model time(h)']*60,model1['model'],color='c',linestyle='dashed', label='Model outlet at setting 1')
plt.plot(model2['model time(h)']*60,model2['model'],color='k',linestyle='dashed', label='Model outlet at setting 2')
plt.plot(model3['model time(h)']*60,model3['model'],color='y',linestyle='dashed', label='Model outlet at setting 3')
plt.plot(model4['model time(h)']*60,model4['model'],color='m',linestyle='dashed', label='Model outlet at setting 4')
# plt.axvline(x=BT1,linestyle='-',label=BT1label) #breakthrough curve
# plt.axvline(x=BT2,linestyle='-',label=BT2label)
plt.xlabel('Time (min)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.2))
plt.grid(True)
plt.xlim(0,20)
plt.ylim(0,0.07)
plt.title('Velocity comparison')
#plt.show()
plt.savefig('..//Plots/Ex 6 velocity comparison.png', bbox_inches='tight')

#Compare the two replicates (5.4 and 6.1)_________________________________________








