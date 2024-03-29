# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shutil.copy('../../../mod_funcs.py', '.')
from mod_funcs import tfmod
 


# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 0.51            # Filter length/depth (m) 
por_g = 0.77      # (m3/m3) Estimated by porosity experiments
por_l = 0.22     # (m3/m3) Estimated by porosity experiments
v_g = 0.017       #not relevant as this is entered manually later
v_l = 1E-4        #not relevant as this is entered manually later
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 21.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 0        # Reaction rate (1/s). Small because of inert carrier
                 # Reaction could be acid/base that changes the pH

pH = 8.07 

# realistic pKa
pKa = 7.

#Calculating the cross sectional area, and dividing the volume by it, as required by the model
vol = 1200 #volume of reservoir in mL
area = (0.19/2)**2 * 3.14159265
v_res = vol * 10**(-6) / area

## Breakthrough time calculated from volume (V=14.5L), Por_g=0.8 and volumetric velocities 25l/min (v_g=53m/h) 
# . 37.5L/min (v_g=79m/h) and 50 L/min (v_g=106m/h)

BT1 = 14.5*por_g/50
BT1label='Theoretical breakthrough for vg=106m/h'




# Inlet concentrations load in from data treatment and definition

from DT_func import Data_reciever
from Calibration_func import constant1,constant2 
#Loading data using the DT_func. Needs a csv file, semicolon separated with the following columns: 
    #Relative Time,  m/z 35.00 ch3 and 37/21. Found in "Time", "concentration" and "raw signal intensities" (make new column that divides m/z37 with m/z21) in PTR-MS file. 
    #constant1 and constant2 are calibration constants used to correct the signal for water content
    #start and stop cycles are the lines in the excel that are imported. All other left out. 
    # note that the data is so that time goes from high to low, but the concentration follows the same pattern, so this should not matter.
    #output is "time" in h and "concentration" in g/m^3
ex1=Data_reciever(filename='experiment_3.2.inlet.csv',startcycle=120,stopcycle=16800,a=constant1,b=constant2)



#Inlet profiles for the model

cgin = pd.DataFrame({'time': ex1['time'][::-1]*3600, 
                     'cgin': ex1['concentration'][::-1]})


clin = 0

# Times for model output in h
tt = 3
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#v_g=60m/h, v_l=0.4m/h
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = 102/3600, v_l = 0.4/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times,v_res = v_res, recirc = True, counter = True)
pred1label= 'Model' #label on plots








# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Loading experimental data from experiment 2.1
ex2=Data_reciever(filename='experiment_3.2.1.csv',startcycle=75,stopcycle=16800,a=constant1,b=constant2)


#plotting


plt.plot(ex2['time'],ex2['concentration'],label='3.2.1 experimental data')
plt.plot(ex1['time'],ex1['concentration'],label='inlet')
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :],color='k',label=pred1label)
plt.axvline(x=BT1/60,linestyle='-',label=BT1label) #breakthrough curve
plt.axhline(y=0.055,color='g',label='Expected inlet concentration')
plt.xlim(0)
plt.ylim(0)
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.5)) #Moves legend out of plot
plt.title('Experiment 2.1')
plt.show()

#table with input paramters
#import module
from tabulate import tabulate
parameters=[['v_g',pred1['inputs']['v_g']],['v_l',pred1['inputs']['v_l']],['pH',pH],['k',k],['countercurrent',pred1['inputs']['counter']],
            ['recirculation',pred1['inputs']['recirc']],['water content',por_l],['porosity',por_g],['temperature',temp]]
head=['parameter','value']
print(tabulate(parameters,headers=head,tablefmt="grid"))



