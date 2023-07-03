# Demo of Python air filter model
# This version uses variable inlet concentrations in gas phase, read in from a file

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shutil.copy('../../mod_funcs.py', '.')
from mod_funcs import tfmod


# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 2            # Filter length/depth (m)
por_g = 0.5      # (m3/m3)
por_l = 0.25     # (m3/m3)
v_g = 0.03       
v_l = 2E-5       
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 15.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 20. / 3600  # Reaction rate (1/s)

pH = 7.
pKa = 7.

cgin = 1

# Fixed for water
clin = 0.        # Fresh water concentration (g/m3)

# Times for model output
tt = 3 
# Number of time rows
nt = 500
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reference, red line
Kga = 0.004
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Recirculation, blue line
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, recirc = True)

# Turn off reaction
# Green line in plots
k = 0. 
pred3 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, recirc = False)

# No reaction with recirculation, expect flat profiles
# Orange line in plots
pred4 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, recirc = True)


# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot outlet concentration (= 1 - removal efficiency here because cgin = 1)
# Gas concentration (outlet air) 
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :], 'r')
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :], 'b-')
plt.plot(pred3['time'] / 3600, pred3['gas_conc'][nc - 1, :], 'g')
plt.plot(pred4['time'] / 3600, pred4['gas_conc'][nc - 1, :], 'orange')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('outlet_gas_conc.png')

# Liquid concentration (in last layer)
plt.clf()
plt.plot(pred1['time'] / 3600, pred1['liq_conc'][nc - 1, :], 'r')
plt.plot(pred2['time'] / 3600, pred2['liq_conc'][nc - 1, :], 'b-')
plt.plot(pred3['time'] / 3600, pred3['liq_conc'][nc - 1, :], 'g-')
plt.plot(pred4['time'] / 3600, pred4['liq_conc'][nc - 1, :], 'orange')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('outlet_liq_conc.png')

# Profiles
# Gas
plt.clf()
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 1], 'r')
plt.plot(pred2['cell_pos'], pred2['gas_conc'][:, nt - 1], 'b-')
plt.plot(pred3['cell_pos'], pred3['gas_conc'][:, nt - 1], 'g-')
plt.plot(pred4['cell_pos'], pred4['gas_conc'][:, nt - 1], 'orange')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_gas_conc.png')

# Liquid
plt.clf()
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 1], 'r')
plt.plot(pred2['cell_pos'], pred2['liq_conc'][:, nt - 1], 'b-')
plt.plot(pred3['cell_pos'], pred3['liq_conc'][:, nt - 1], 'g-')
plt.plot(pred4['cell_pos'], pred4['liq_conc'][:, nt - 1], 'orange')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_liq_conc.png')
