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
por_g = 0.80      # (m3/m3) Estimated by volume calculations
por_l = 0.018     # (m3/m3) Estimated by volume calculations
v_g = 0.017       #not relevant as this is entered manually
v_l = 1E-4        #not relevant as this is entered manually
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 21.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 0        # Reaction rate (1/s). Small because of inert carrier
                 # Reaction could be acid/base that changes the pH

pH = 6.

# realistic pKa
pKa = 7.


# Put inlet concentrations at equilibrium
from Data_treatment_repetitions import C_inlet
cgin = float(C_inlet)  #Average measured value. Note: For future this will be a function of time. But as it is not on the same time axis as the outlet concentrations,
# an average value is just as good as anything time dependent.  
clin = 0

# Times for model output in h
tt = 0.08 
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Use only model 3!! And make it correct
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = 106/3600, v_l = 1.2/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

pred1label='Model with average measured inlet concentration' 





# Closed-form solution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# z = np.linspace(0, L, 10)
# pr = por_l / por_g
# ctin = por_g * cgin + por_l * clin
# ct = ctin * np.exp(-k  * por_l / (v_g * Kaw) * z)
# cg = Kaw * ct / (por_g * Kaw + por_l)
# cl = cg / Kaw



# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Breakthrough time calculated from volume (V=14.5L), Por_g=0.8 and volumetric velocities 25l/min (v_g=53m/h) 
# and 50 L/min (v_g=106m/h)

BT1 = 14.5*por_g/50
BT2 = 14.5*por_g/25



#All experimental data in one graph

from Data_treatment_repetitions import t1,t2,t3,t4,t5,C_out1,C_out2,C_out3,C_out4,C_out5


plt.plot(t1,C_out1,color='b', label='New setup 1')
plt.plot(t2,C_out2,color='r',label='New setup 2')
plt.plot(t3,C_out3,color='y',label='Old setup 1')
plt.plot(t4,C_out4,color='k',label='Old setup 2')
plt.plot(t5,C_out5,color='c',label='Old setup 3')
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :], color='m',label=pred1label)
plt.axvline(x=BT1/60,color='violet',linestyle='-',label='Theoretical breakthrough for vg=106m/h')
#plt.axvline(x=BT2/60,color='deeppink',linestyle='-',label='Theoretical Breakthrough for vg=53m/h)')
plt.axhline(y=cgin,color='g',label='Measured inlet concentration')
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.5))
plt.title('Experiment 3 all experimental data')
plt.xlabel('Time(h)')
plt.ylabel('concentration (g/m3)')
plt.show()
