# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from Mod_funcs import tfmod
  


# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L =0.51            # Filter length/depth (m) 
por_g = 0.73     # (m3/m3) Estimated by volume calculations
por_l = 0.08       # (m3/m3) Estimated at 0.08, changed to be able to tell the graphs apart
v_g = 106/3600       # superficial gas velocity m/s (approx. middle of ex1)
v_l = 0.8/3600        #liquid superficial velocity m/s (approx middle of ex1)
nc = 100          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 21.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 0        # Reaction rate (1/s). Small because of inert carrier
                 # Reaction could be acid/base that changes the pH

pH = 9.

# realistic pKa
pKa = 7.


# Put inlet concentrations at equilibrium
cgin = 0.05575209  #corresponding to 40ppm
clin = 0

# Times for model output in h
tt = 0.08 
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)




# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Profiles
# Gas
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, 1], color='b',label='time=start')
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 160], color='r',label='time=1/5')
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 120], color='y',label='time=2/5')
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 80], color = 'k',label='time=3/5')
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 40], color='c',label='time=4/5')
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 1], color='m',label='time=end')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.title('Gas Phase @pH=9')
plt.legend()
plt.show()

#Liquid
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, 1], color='b',label='time=start')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 160], color='r',label='time=1/5')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 120], color='y',label='time=2/5')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 80], color = 'k',label='time=3/5')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 40], color='c',label='time=4/5')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 1], color='m',label='time=end')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.title('Liquid Phase @pH=9')
plt.legend()
plt.show()

# Profiles
# Gas
plt.plot(pred1['cell_pos'], pred1['gas_mass'][:, 1], color='b',label='time=start')
plt.plot(pred1['cell_pos'], pred1['gas_mass'][:, nt - 160], color='r',label='time=1/5')
plt.plot(pred1['cell_pos'], pred1['gas_mass'][:, nt - 120], color='y',label='time=2/5')
plt.plot(pred1['cell_pos'], pred1['gas_mass'][:, nt - 80], color = 'k',label='time=3/5')
plt.plot(pred1['cell_pos'], pred1['gas_mass'][:, nt - 40], color='c',label='time=4/5')
plt.plot(pred1['cell_pos'], pred1['gas_mass'][:, nt - 1], color='m',label='time=end')
plt.xlabel('Location (m)')
plt.ylabel('Compound mass. (g)')
plt.title('Gas Phase @ pH=9')
plt.legend()
plt.show()

#Liquid

plt.plot(pred1['cell_pos'], pred1['liq_mass'][:, 1], color='b',label='time=start')
plt.plot(pred1['cell_pos'], pred1['liq_mass'][:, nt - 160], color='r',label='time=1/5')
plt.plot(pred1['cell_pos'], pred1['liq_mass'][:, nt - 120], color='y',label='time=2/5')
plt.plot(pred1['cell_pos'], pred1['liq_mass'][:, nt - 80], color = 'k',label='time=3/5')
plt.plot(pred1['cell_pos'], pred1['liq_mass'][:, nt - 40], color='c',label='time=4/5')
plt.plot(pred1['cell_pos'], pred1['liq_mass'][:, nt - 1], color='m',label='time=end')
plt.xlabel('Location (m)')
plt.ylabel('Compound mass. (g)')
plt.title('Liquid Phase @ pH=9')
plt.legend()
plt.show()


#The entire column
#plotting the total mass of sulfur and H2S in thhe entire column as a function of time    
 
plt.plot(pred1['time']/3600, pred1['H2S_mass'], color='b',label='mass of H2S')
plt.plot(pred1['time']/3600, pred1['tot_mass'], color='r',label='mass of H2S and HS')
plt.xlabel('Time (h)')
plt.ylabel('Compound mass. (g)')
plt.title('Total mass in column @ pH=9')
plt.legend()
plt.show()   






