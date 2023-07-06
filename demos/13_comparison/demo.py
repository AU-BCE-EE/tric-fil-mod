# Comparison of numerical Python model to closed-form solution

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
nc = 100          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 15.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 1000         # Reaction rate (1/s)

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
# Reference, red line
Kga = 0.004
k = 0.
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Blue
k = 0.01
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Green line in plots
k = 1.
pred3 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Orange line in plots
k = 1000.
pred4 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)


# Closed-form solution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z = np.linspace(0, L, 10)
cg = cgin * np.exp(-Kga / v_g * z)

# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Profiles
# Gas
plt.clf()
plt.plot(z, cg, 'bo')
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 1], 'r')
plt.plot(pred2['cell_pos'], pred2['gas_conc'][:, nt - 1], 'b-')
plt.plot(pred3['cell_pos'], pred3['gas_conc'][:, nt - 1], 'g-')
plt.plot(pred4['cell_pos'], pred4['gas_conc'][:, nt - 1], 'orange')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_gas_conc.png')

# Liquid
plt.clf()
plt.plot(z, 0 * cg, 'bo')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 1], 'r')
plt.plot(pred2['cell_pos'], pred2['liq_conc'][:, nt - 1], 'b-')
plt.plot(pred3['cell_pos'], pred3['liq_conc'][:, nt - 1], 'g-')
plt.plot(pred4['cell_pos'], pred4['liq_conc'][:, nt - 1], 'orange')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_liq_conc.png')
