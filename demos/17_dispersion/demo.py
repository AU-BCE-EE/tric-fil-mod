# Effect of dispersion (test)

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
L = 0.1            # Filter length/depth (m); zooming in on first 10 cm!
por_g = 0.5      # (m3/m3)
por_l = 0.25     # (m3/m3)
v_g = 0.03       
v_l = 2E-5       
nc = 100          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 15.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 0         # Reaction rate (1/s)

pH = 7.
pKa = 7.

cgin = 1

# Fixed for water
clin = 0.        # Fresh water concentration (g/m3)

# Times for model output
tt = 3 
# Number of time rows
nt = 100
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reference, red line, no dispersion

pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Blue
disp_gas = 1.0E-7 # m2/s
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', disp_gas = disp_gas, k = k, 
              henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Green line in plots
disp_liq = 0.0 # m2/s
pred3 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', disp_liq = disp_liq,
              k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Orange line in plots
disp_gas = 0.0
pred4 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', disp_gas = disp_gas, 
              k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
plt.figure(2)
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 1], 'r')
plt.plot(pred2['cell_pos'], pred2['liq_conc'][:, nt - 1], 'b-')
plt.plot(pred3['cell_pos'], pred3['liq_conc'][:, nt - 1], 'g-')
plt.plot(pred4['cell_pos'], pred4['liq_conc'][:, nt - 1], 'orange')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_liq_conc.png')
