# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:46:28 2024

@author: Mortensen
"""

#load packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Velocity comparison (experiment 6)

#moving average over 5 data points on each side
window_size = 11


ex1 = pd.read_csv('..//Output data/residuals_6.1.csv', sep = ',')
ex1['Moving_Average'] = ex1['residuals[g/m3]'].rolling(window=window_size,center = True).mean()

ex2 = pd.read_csv('..//Output data/residuals_6.2.csv', sep = ',')
ex2['Moving_Average'] = ex2['residuals[g/m3]'].rolling(window=window_size,center = True).mean()

ex3 = pd.read_csv('..//Output data/residuals_6.3.csv', sep = ',')
ex3['Moving_Average'] = ex3['residuals[g/m3]'].rolling(window=window_size,center = True).mean()

ex4 = pd.read_csv('..//Output data/residuals_6.4.csv', sep = ',')
ex4['Moving_Average'] = ex4['residuals[g/m3]'].rolling(window=window_size,center = True).mean()



plt.plot(ex1['time_residuals[h]']*60,ex1['Moving_Average'],color='c',label='Velocity setting 1')
plt.plot(ex2['time_residuals[h]']*60,ex2['Moving_Average'],color='k',label='Velocity setting 2')
plt.plot(ex3['time_residuals[h]']*60,ex3['Moving_Average'],color='y',label='Velocity setting 3')
plt.plot(ex4['time_residuals[h]']*60,ex4['Moving_Average'],color='m',label='Velocity setting 4')
plt.xlabel('Time (min)')
plt.ylabel('Residual (g/m3)')
plt.legend()
plt.grid(which='both', linestyle='-', linewidth=0.5)
plt.xticks(range(0,20,1))
plt.xlim(0,20)
plt.ylim(-0.015,0.015)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Velocity comparison')
plt.show()
#plt.savefig('..//Plots/Ex 6 velocity comparison residuals.png', bbox_inches='tight')
plt.close()


#pH comparison

ex5 = pd.read_csv('..//Output data/residuals_5.1.csv', sep = ',')
ex5['Moving_Average'] = ex5['residuals[g/m3]'].rolling(window=window_size,center = True).mean()

ex6 = pd.read_csv('..//Output data/residuals_5.3.csv', sep = ',')
ex6['Moving_Average'] = ex6['residuals[g/m3]'].rolling(window=window_size,center = True).mean()

ex7 = pd.read_csv('..//Output data/residuals_5.4.csv', sep = ',')
ex7['Moving_Average'] = ex7['residuals[g/m3]'].rolling(window=window_size,center = True).mean()

ex8 = pd.read_csv('..//Output data/residuals_5.5.csv', sep = ',')
ex8['Moving_Average'] = ex8['residuals[g/m3]'].rolling(window=window_size,center = True).mean()

ex9 = pd.read_csv('..//Output data/residuals_5.7.csv', sep = ',')
ex9['Moving_Average'] = ex9['residuals[g/m3]'].rolling(window=window_size,center = True).mean()

plt.clf()
plt.plot(ex9['time_residuals[h]']*60,ex9['Moving_Average'],color='y',label='pH 6')
plt.plot(ex5['time_residuals[h]']*60,ex5['Moving_Average'],color='c',label='pH 7')
plt.plot(ex8['time_residuals[h]']*60,ex8['Moving_Average'],color='b',label='pH 7.5')
plt.plot(ex6['time_residuals[h]']*60,ex6['Moving_Average'],color='r',label='pH 7.75')
plt.plot(ex7['time_residuals[h]']*60,ex7['Moving_Average'],color='m',label='pH 8')
plt.xlabel('Time (min)')
plt.ylabel('Residual (g/m3)')
plt.legend()
plt.grid(which='both', linestyle='-', linewidth=0.5)
plt.xticks(range(0,20,1))
plt.xlim(0,20)
plt.ylim(-0.03,0.03)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('pH comparison')
plt.show()
#plt.savefig('..//Plots/Ex 5 pH comparison residuals.png', bbox_inches='tight')








