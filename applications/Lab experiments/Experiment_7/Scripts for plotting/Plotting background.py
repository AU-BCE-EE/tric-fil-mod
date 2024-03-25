# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



# Choose experiment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# name of experiment (file named "experiment_first.second.x)
first = '7'
second = '1'

#parameters loaded. change name depending on experiment no.

from lab_parameters_71 import cycle1,cycle2,cycle3,cycle4,length 





# Inlet concentrations load in from data treatment and definition

ex1 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')


#second inlet profile (last data obtained)
ex4 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.inlet2.csv', sep = ',')


# The three repetitions
ex2 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.1.csv', sep = ',')
ex3 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.2.csv', sep = ',')

#The background
ex5 = pd.read_csv('..//Processed_data/Background.csv',sep=',')


#Methyl mercaptan___________________________________________________________________
t1norm = ex1['Time in h'] - ex1['Time in h'][cycle1]
t1 = t1norm [cycle1:cycle1+length]
C1= ex1['mm concentration in g/m^3'][cycle1:cycle1+length]

t4norm = ex4['Time in h'] - ex4['Time in h'][cycle4]
t4 = t4norm [cycle4:cycle4+length]
C4= ex4['mm concentration in g/m^3'][cycle4:cycle4+length]

t2norm = ex2['Time in h'] - ex2['Time in h'][cycle2]
t2 = t2norm [cycle2:cycle2+length+1500]
C2= ex2['mm concentration in g/m^3'][cycle2:cycle2+length+1500]

t3norm = ex3['Time in h'] - ex3['Time in h'][cycle3]
t3 = t3norm [cycle3:cycle3+length+1500]
C3= ex3['mm concentration in g/m^3'][cycle3:cycle3+length+1500]

t5 = ex5['Time in h']
C5= ex5['mm concentration in g/m^3']





# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



plt.figure()
plt.plot(t2,C2,label=first+'.'+second+'.1 experimental data')
plt.plot(t3,C3,label=first+'.'+second+'.2 experimental data')
plt.plot(t1,C1,label='inlet 1')
plt.plot(t4,C4,label='inlet 2')
plt.plot(t5,C5,label='Background')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Experiment '+first+'.'+second+' mm')
plt.savefig('..//Plots/Background_mm.png', bbox_inches='tight')






#Dimethyl sulfide____________________________________________________________________

C12= ex1['DMS concentration in g/m^3'][cycle1:cycle1+length]


C42= ex4['DMS concentration in g/m^3'][cycle4:cycle4+length]


C22= ex2['DMS concentration in g/m^3'][cycle2:cycle2+length+1500]


C32= ex3['DMS concentration in g/m^3'][cycle3:cycle3+length+1500]


C52= ex5['DMS concentration in g/m^3']


# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



plt.figure()
plt.plot(t2,C22,label=first+'.'+second+'.1 experimental data')
plt.plot(t3,C32,label=first+'.'+second+'.2 experimental data')
plt.plot(t1,C12,label='inlet 1')
plt.plot(t4,C42,label='inlet 2')
plt.plot(t5,C52,label='Background')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Experiment '+first+'.'+second+' DMS')
plt.savefig('..//Plots/Background_DMS.png', bbox_inches='tight')



#Hydrogen sulfide____________________________________________________________________

C12= ex1['H2S concentration in g/m^3'][cycle1:cycle1+length]


C42= ex4['H2S concentration in g/m^3'][cycle4:cycle4+length]


C22= ex2['H2S concentration in g/m^3'][cycle2:cycle2+length+1500]


C32= ex3['H2S concentration in g/m^3'][cycle3:cycle3+length+1500]


C52= ex5['H2S concentration in g/m^3']


# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



plt.figure()
plt.plot(t2,C22,label=first+'.'+second+'.1 experimental data')
plt.plot(t3,C32,label=first+'.'+second+'.2 experimental data')
plt.plot(t1,C12,label='inlet 1')
plt.plot(t4,C42,label='inlet 2')
plt.plot(t5,C52,label='Background')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.ylim(-0.00005)
plt.xlim(0)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Experiment '+first+'.'+second+' H2S')
plt.savefig('..//Plots/Background_H2S.png', bbox_inches='tight')

