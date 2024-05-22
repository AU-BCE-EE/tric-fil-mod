# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sys.path.append("../../../../../..")  # Add the directory containing mod_funcs.py to Python path
from mod_funcs import tfmod  

# Choose experiment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# name of experiment (file named "experiment_first.second.x)
first = '7'
second = '1'

#parameters loaded. change name depending on experiment no.

from lab_parameters_71 import pH1,pH2,cycle1,cycle2,cycle3,cycle4,length,vol 


# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 0.51            # Filter length/depth (m) 
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.56, 3500.)
temp = 21.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

no=1
area = (0.19/2)**2 * 3.14159265
v_res = vol * 10**(-6) / area

#Setting the liquid and gas velocities using the experiment number (second number in the file name)
if no == 1 or no == 4:
    v_g = 27.2991 / 1000 / 60 / area # gas velocity using calibration curve for mass flow controllers L/min /1000L/m3 / 60 s/min / m2 = m/s
elif no ==2 or no == 3:
    v_g = 54.2766 /1000 / 60 /area # gas velocity using calibration curve for mass flow controllers L/min /1000L/m3 / 60 s/min / m2 = m/s
else: print ('Error in reading "no"')

if no == 1 or no == 2:
    v_l = 0.390681119/3600 #in m/s, based on pump calibration number 2
elif no ==4 or no == 3:
    v_l = 1.253842217/3600 #in m/s, based on pump calibration number 2
else: print ('Error in reading "no"')

# Porosity and water content definition
#from water content experiments, as well as porosity experiments. Average of three measurements. gas phase is porosity minus water content
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


k = 0       # Reaction rate (1/s). Small because of inert carrier
                 # Reaction could be acid/base that changes the pH



# realistic pKa
pKa = 14.



## Breakthrough time calculated from volume (V=14.5L), Por_g=0.8 and volumetric velocities 25l/min (v_g=53m/h) 
# . 37.5L/min (v_g=79m/h) and 50 L/min (v_g=106m/h)

BT1 = 14.5*por_g/54.2766
BT2 = 14.5*por_g/27.2991
BTlabel='Theoretical Breakthrough curve'

if no == 1 or no == 4:
    BT = BT2
elif no ==2 or no == 3:
    BT= BT1
else: 
    print ('error in definition of "no". Has to be a string containing the numbers 1,2,3 or 4')



# Inlet concentrations load in from data treatment and definition

ex1 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')


#second inlet profile (last data obtained)
ex4 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.inlet2.csv', sep = ',')


# The three repetitions
ex2 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.1.csv', sep = ',')
ex3 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.2.csv', sep = ',')



# selecting and normalising the time, because all data files are started before the H2S
# was turned on, and not after a specific time, it all depends on the individual experient. 
#However the cycle at which the H2S was turned on was written down.
t1norm = ex1['Time in h'] - ex1['Time in h'][cycle1]
t1 = t1norm [cycle1:cycle1+length]
C1= ex1['DMS concentration in g/m^3'][cycle1:cycle1+length]

t4norm = ex4['Time in h'] - ex4['Time in h'][cycle4]
t4 = t4norm [cycle4:cycle4+length]
C4= ex4['DMS concentration in g/m^3'][cycle4:cycle4+length]

t2norm = ex2['Time in h'] - ex2['Time in h'][cycle2]
t2 = t2norm [cycle2:cycle2+length]
C2= ex2['DMS concentration in g/m^3'][cycle2:cycle2+length]

t3norm = ex3['Time in h'] - ex3['Time in h'][cycle3]
t3 = t3norm [cycle3:cycle3+length]
C3= ex3['DMS concentration in g/m^3'][cycle3:cycle3+length]




clin = 0

# Times for model output in h
tt = 0.35
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#moving average over 5 data points on each side
window_size = 11

c_out = np.mean([C2,C3],axis = 0)
# c_out= pd.DataFrame(c_out)
# c_out = c_out.rolling(window=window_size,center = True).mean()
c_in = np.mean([C1,C4],axis = 0)





cgin = pd.DataFrame({'time': t3*3600, 
                     'cgin': c_in})


pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = (pH1+pH2)/2, temp = temp, dens_l = dens_l, times = t3*3600, v_res = v_res, recirc = True, counter = True)
pred2label='Model outlet gas phase'



# Define y values for the expected inlet curve @5ppm
expected_inlet = np.full_like(ex1['Time in h']*60, 0.000254)  

# Make the line drop to 0 at x=5, as the oulse is 5min
expected_inlet[ex1['Time in h']*60 >= 5] = 0


# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


c_out_plot = pd.DataFrame(c_out).rolling(window=window_size,center = True).mean()

plt.plot(t3 * 60,c_out_plot,label='Measured outlet gas phase',color ='b')
plt.plot(pred2['time'] / 60, pred2['gas_conc'][nc - 1, :],color='b',linestyle = 'dashed',label=pred2label)
plt.axvline(x=BT,linestyle='-',label=BTlabel) #breakthrough curve
plt.plot(ex1['Time in h']*60, expected_inlet, color='g', label='Expected inlet concentration')
plt.xlabel('Time (min)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.grid(True)
plt.xlim(0,20)
plt.ylim(0)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Dimethyl sulfide')
plt.show()
#plt.savefig('..//Plots/Experiment '+first+'.'+second+' DMS.png', bbox_inches='tight')




MAE_baseline = length**-1*np.sum(abs( pred2['gas_conc'][nc-1:]-c_out ))
c_out_mean = np.mean(c_out)
ME_baseline = (np.sum((c_out-c_out_mean)**2)-np.sum((pred2['gas_conc'][nc-1:]-c_out)**2))/np.sum((c_out-c_out_mean)**2)

print('ME '+str(ME_baseline))
print('MAE '+str(MAE_baseline))




