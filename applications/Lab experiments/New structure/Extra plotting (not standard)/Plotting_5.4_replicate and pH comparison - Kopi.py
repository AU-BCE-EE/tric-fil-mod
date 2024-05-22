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
BTlabel='Theoretical Breakthrough curve'

# Define y values for the expected inlet curve @40ppm
expected_inlet = np.full_like(ex7['experimental time (h)']*60, 0.0596)  

# Make the line drop to 0 at x=5, as the oulse is 5min
expected_inlet[ex7['experimental time (h)']*60 >= 5] = 0



plt.plot(ex7['experimental time (h)']*60,ex7['Moving_Average'],color='y',label='Measured outlet at pH=5.99')
plt.plot(ex1['experimental time (h)']*60,ex1['Moving_Average'],color='c',label='Measured outlet at pH=7.08')
#plt.plot(ex2['experimental time (h)']*60,ex2['Moving_Average'],color='k',label='Measured outlet at pH=7.27')
plt.plot(ex5['experimental time (h)']*60,ex5['Moving_Average'],color='b',label='Measured outlet at pH=7.55')
plt.plot(ex3['experimental time (h)']*60,ex3['Moving_Average'],color='r',label='Measured outlet at pH=7.76')
plt.plot(ex4['experimental time (h)']*60,ex4['Moving_Average'],color='m',label='Measured outlet at pH=8.02')
plt.plot(model7['model time(h)']*60,model7['model'],color='y',linestyle='dashed', label='Model outlet at pH=5.99')
plt.plot(model1['model time(h)']*60,model1['model'],color='c',linestyle='dashed', label='Model outlet at pH=7.08')
#plt.plot(model2['model time(h)']*60,model2['model'],color='k',linestyle='dashed', label='Model outlet at pH=7.27')
plt.plot(model5['model time(h)']*60,model5['model'],color='b',linestyle='dashed', label='Model outlet at pH=7.55')
plt.plot(model3['model time(h)']*60,model3['model'],color='r',linestyle='dashed', label='Model outlet at pH=7.76')
plt.plot(model4['model time(h)']*60,model4['model'],color='m',linestyle='dashed', label='Model outlet at pH=8.02')
plt.axvline(x=BT,linestyle='-',label=BTlabel) #breakthrough curve
plt.plot(ex7['experimental time (h)']*60, expected_inlet, color='g', label='Expected inlet concentration')
plt.xlabel('Time (min)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.2))
plt.grid(True)
plt.xlim(0,20)
plt.ylim(0,0.07)
plt.title('pH comparison, baseline model')
#plt.show()
plt.savefig('..//Plots/Ex 5 pH comparison.png', bbox_inches='tight')

#Compare the two replicates (5.4 and 6.1)_________________________________________


#Import data and make the moving average to avoid noise


#import data

ex8 = pd.read_csv('..//Output data/experimental6.1.csv', sep = ',')
model8 = pd.read_csv('..//Output data/model_6.1.csv', sep = ',')
ex8['Moving_Average'] = ex8['experimental'].rolling(window=window_size).mean()



# Define y values for the expected inlet curve @40ppm
expected_inlet = np.full_like(ex4['experimental time (h)']*60, 0.0596)  

# Make the line drop to 0 at x=5, as the oulse is 5min
expected_inlet[ex4['experimental time (h)']*60 >= 5] = 0

# plot
plt.clf()
plt.plot(ex4['experimental time (h)']*60,ex4['Moving_Average'],color='c',label='Outlet gas phase concentration day 1')
plt.plot(ex8['experimental time (h)']*60,ex8['Moving_Average'],color='k',label='Outlet gas phase concentration day 2')
#plt.plot(model4['model time(h)']*60,model4['model'],color='c',linestyle='dashed', label='model 5.4')
#plt.plot(model8['model time(h)']*60,model8['model'],color='k',linestyle='dashed', label='model 6.1')
plt.axvline(x=BT/60,linestyle='-',label=BTlabel) #breakthrough curve
plt.plot(ex4['experimental time (h)']*60, expected_inlet, color='g', label='Expected inlet concentration')
plt.xlabel('Time (min)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.grid(True)
plt.xlim(0,20)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('replicate comparison')
plt.show()
#plt.savefig('..//Plots/5.4 compared to 6.1.png', bbox_inches='tight')





