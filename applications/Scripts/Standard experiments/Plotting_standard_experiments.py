# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shutil.copy('../../../mod_funcs.py', '.')
from mod_funcs import tfmod 

# Choose experiment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# name of experiment (file named "experiment_first.second.x)
first = '4'
second = '4'

#parameters loaded. change name depending on experiment no.

from lab_parameters_44 import pH1,pH2,pH3,cycle1,cycle2,cycle3,cycle4,cycle5,length 



# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 0.51            # Filter length/depth (m) 
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 21.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

#Setting the liquid and gas velocities using the experiment number (second number in the file name)
if second == '1' or second == '4':
    v_g = 53/3600
elif second =='2' or second == '3':
    v_g = 106/3600
else: print ('Error in reading "second"')

if second == '1' or second == '2':
    v_l = 0.4/3600
elif second =='4' or second == '3':
    v_l = 1.2/3600
else: print ('Error in reading "second"')

# Porosity and water content definition
if second == '1':
    por_l = 0.21
    por_g = 0.77
elif second == '2':
    por_l = 0.23
    por_g = 0.75
elif second == '3':
    por_l = 0.24
    por_g = 0.74
elif second == '4':
    por_l = 0.23
    por_g = 0.75
else: print ('Error in reading "second"') 


k = 0       # Reaction rate (1/s). Small because of inert carrier
                 # Reaction could be acid/base that changes the pH



# realistic pKa
pKa = 7.


## Breakthrough time calculated from volume (V=14.5L), Por_g=0.8 and volumetric velocities 25l/min (v_g=53m/h) 
# . 37.5L/min (v_g=79m/h) and 50 L/min (v_g=106m/h)

BT1 = 14.5*por_g/50
BT2 = 14.5*por_g/25
BTlabel='Theoretical Breakthrough curve)'

if second == '1' or second == '4':
    BT = BT2
elif second =='2' or second == '3':
    BT= BT1
else: 
    print ('error in definition of "second". Has to be a string containing the numbers 1,2,3 or 4')



# Inlet concentrations load in from data treatment and definition

ex1 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')


#second inlet profile (last data obtained)
ex5 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet2.csv', sep = ',')


# The three repetitions
ex2 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.1.csv', sep = ',')
ex3 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.2.csv', sep = ',')
ex4 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.3.csv', sep = ',')


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

t4norm = ex4['Time in h'] - ex4['Time in h'][cycle4]
t4 = t4norm [cycle4:cycle4+length+500]
C4= ex4['Concentration in g/m^3'][cycle4:cycle4+length+500]


    



#Average of the two inlet profiles for use in the model

c_in=np.mean([C1,C5],axis=0)

cgin = pd.DataFrame({'time': t5*3600, 
                     'cgin': c_in})

clin = 0

# Times for model output in h
tt = 0.35
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH1, temp = temp, dens_l = dens_l, times = times)
pred1label= first+'.'+second+'.1 model' #label on plots




pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH2, temp = temp, dens_l = dens_l, times = times)
pred2label=first+'.'+second+'.2 model'



pred3 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH3, temp = temp, dens_l = dens_l, times = times)
pred3label=first+'.'+second+'.3 model'




# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




plt.plot(t2,C2,label=first+'.'+second+'.1 experimental data')
plt.plot(t3,C3,label=first+'.'+second+'.2 experimental data')
plt.plot(t4,C4,label=first+'.'+second+'.3 experimental data')
plt.plot(t1,C1,label='inlet 1')
plt.plot(t5,C5,label='inlet 2')
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :],color='k',label=pred1label)
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :],color='m',label=pred2label)
plt.plot(pred3['time'] / 3600, pred3['gas_conc'][nc - 1, :],color='b',label=pred3label)
plt.axvline(x=BT/60,linestyle='-',label=BTlabel) #breakthrough curve
plt.axhline(y=0.055,color='g',label='Expected inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0)
plt.ylim(0)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Experiment '+first+'.'+second)
plt.savefig('../..//Plots/Experiment '+first+'.'+second+'.png', bbox_inches='tight')

#table with input paramters
#import module
from tabulate import tabulate
parameters = [['v_g', pred1['inputs']['v_g']], ['v_l', pred1['inputs']['v_l']], ['pH1', pH1], ['pH2', pH2], ['pH3', pH3], ['k', k], ['countercurrent', pred1['inputs']['counter']],
              ['recirculation', pred1['inputs']['recirc']], ['water content', por_l], ['porosity', por_g], ['temperature', temp]]
head = ['parameter', 'value']
table = tabulate(parameters, headers=head, tablefmt="grid")

# Create a figure and axis
fig, ax = plt.subplots()

# Hide the axes
ax.axis('off')

# Plot the table
table_ax = ax.table(cellText=parameters, colLabels=head, loc='center', cellLoc='center')

# Adjust font size
table_ax.auto_set_font_size(False)
table_ax.set_fontsize(10)

plt.title('Experiment '+first+'.'+second)

# Save the figure
plt.savefig('../..//Plots/Inputs/Input_parameters_'+first+'.'+second+'.png', bbox_inches='tight')

