# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:46:28 2024

@author: Mortensen
"""

#load packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd




first = '4'
second = '1'

file_path = '../../Master_spreadsheet.xlsx'
df = pd.read_excel(file_path, sheet_name = 'Data')

params = df[df['key'] == float(first)+0.1*float(second)]

pH1 = params['pH1'].values[0]
pH2 = params['pH2'].values[0]
pH3 = params['pH3'].values[0]
cycle1 = int(params['cycle1'].values[0])
cycle2 = int(params['cycle2'].values[0])
cycle3 = int(params['cycle3'].values[0])
cycle4 = params['cycle4'].values[0]
if not pd.isna(cycle4):
    cycle4=int(cycle4)
else:
    cycle4 = 'NaN'
cycle5 = int(params['cycle5'].values[0])
length = params['length'].values[0]
vol = params['volume'].values[0]
no = params['no'].values[0]

  # Inlet concentrations load in from data treatment and definition

ex1 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')


#second inlet profile (last data obtained)
ex5 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.inlet2.csv', sep = ',')


# The three repetitions
ex2 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.1.csv', sep = ',')
ex3 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.2.csv', sep = ',')
ex4 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.3.csv', sep = ',')

t4norm = ex4['Time in h'] - ex4['Time in h'][cycle4]
t4 = t4norm [cycle4:cycle4+length+500]
C4= ex4['Concentration in g/m^3'][cycle4:cycle4+length+500]
    



# selecting and normalising the time, because all data files are started before the H2S
# was turned on, and not after a specific time, it all depends on the individual experient. 
#However the cycle at which the H2S was turned on was written down.
t1norm = ex1['Time in h'] - ex1['Time in h'][cycle1]
t1 = t1norm [cycle1:cycle1+length]
C1= ex1['Concentration in g/m^3'][cycle1:cycle1+length]

t5norm = ex5['Time in h'] - ex5['Time in h'][cycle5]
t5 = t5norm [cycle5:cycle5+length]
C5= ex5['Concentration in g/m^3'][cycle5:cycle5+length]

t2norm = ex2['Time in h'] - ex2['Time in h'][cycle2]
t2 = t2norm [cycle2:cycle2+length+500]
C2= ex2['Concentration in g/m^3'][cycle2:cycle2+length+500]

t3norm = ex3['Time in h'] - ex3['Time in h'][cycle3]
t3 = t3norm [cycle3:cycle3+length+500]
C3= ex3['Concentration in g/m^3'][cycle3:cycle3+length+500]


if no == 1:
    por_l = 0.20498
    por_g = 0.795115372 - por_l
elif no == 2:
    por_l = 0.22494
    por_g = 0.795115372 - por_l
elif no == 3:
    por_l = 0.24056
    por_g = 0.795115372 - por_l
elif no == 4:
    por_l = 0.246304
    por_g = 0.795115372 - por_l
else: print ('Error in reading "no"') 



BT1 = 14.5*por_g/54.2766
BT2 = 14.5*por_g/27.2991
BTlabel='Theoretical Breakthrough curve'

if no == 1 or no == 4:
    BT = BT2
elif no ==2 or no == 3:
    BT= BT1
else: 
    print ('error in definition of "no". Has to be a string containing the numbers 1,2,3 or 4')

plt.clf()
plt.plot(t2*60,C2,label='Oultlet gas phase rep.1')
plt.plot(t3*60,C3,label='Oultlet gas phase rep.2')
if not cycle4 == 'NaN':
    plt.plot(t4*60,C4,label='Oultlet gas phase rep.3')
plt.plot(t1*60,C1,label='Inlet gas phase rep.1')
plt.plot(t5*60,C5,label='Inlet gas phase rep.2')
#plt.axvline(x=BT/60,linestyle='-',label=BTlabel) #breakthrough curve
plt.axhline(y=0.0596,color='g',label='Expected inlet concentration') # Average of two inlet concentrations of 43.2 and 42.3. These were for low and high gas flow rate respectively
plt.xlabel('Time (min)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.grid(True)
plt.xlim(0,20)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Initial experiment')
plt.show()
plt.savefig('..//Plots/Raw experimental data 4.1.png', bbox_inches='tight')
plt.close()




