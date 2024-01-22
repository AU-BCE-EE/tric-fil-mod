# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from Mod_Funcs import tfmod  


# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 0.51            # Filter length/depth (m) 
por_g = 0.80      # (m3/m3) Estimated by porosity experiments
por_l = 0.08     # (m3/m3) Estimated by porosity experiments
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

pH = 8.1 #CHANGE_ME

# realistic pKa
pKa = 7.


## Breakthrough time calculated from volume (V=14.5L), Por_g=0.8 and volumetric velocities 25l/min (v_g=53m/h) 
# . 37.5L/min (v_g=79m/h) and 50 L/min (v_g=106m/h)

BT1 = 14.5*por_g/50
BT1label='Theoretical breakthrough for vg=106m/h'
BT2 = 14.5*por_g/25
BT2label='Theoretical Breakthrough for vg=53m/h)'
BT3 = 14.5*por_g/37.5
BT3label= 'Theoretical Breakthrough for vg=79m/h)'


# Inlet concentrations load in from data treatment and definition
from Data_treatment_ex2 import t1
from Data_treatment_ex2 import C_out1
cgin = pd.DataFrame({'time': t1[::-1], 
                     'cgin': C_out1[::-1]})


clin = 0

# Times for model output in h
tt = 0.35
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#v_g=60m/h, v_l=0.4m/h
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = 53/3600, v_l = 0.4/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
Daw=pred1['pars']['Daw'] #used in equilibrium calculations
pred1label= '2.1.1 and 2.1.2 model' #label on plots

#second set of inlet profile for the third experiment
from Data_treatment_ex2 import t5
from Data_treatment_ex2 import C_out5
cgin = pd.DataFrame({'time': t5, 
                     'cgin': C_out5})

pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = 53/3600, v_l = 0.4/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
pred2label='2.1.3 model'






# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Loading experimental data from experiment 2.1
from Data_treatment_ex2 import t2
from Data_treatment_ex2 import C_out2
from Data_treatment_ex2 import t3
from Data_treatment_ex2 import C_out3
from Data_treatment_ex2 import t4
from Data_treatment_ex2 import C_out4

#plotting

plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :],label=pred1label)
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :],label=pred2label)
plt.plot(t2,C_out2,label='2.1.1 experimental data')
plt.plot(t3,C_out3,label='2.1.2 experimental data')
plt.plot(t4,C_out4,label='2.1.3 experimental data')
plt.plot(t1,C_out1,label='inlet 1')
plt.plot(t5,C_out5,label='inlet 2')
plt.axvline(x=BT1/60,linestyle='-',label=BT1label) #breakthrough curve
plt.axhline(y=0.055,color='g',label='Expected inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.5)) #Moves legend out of plot
plt.title('Experiment 2.1')
plt.show()

