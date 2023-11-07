# Comparison of numerical Python model to closed-form solution
# Idea (from Anders) is that at very high pH the volatile species concentration in the liquid phase is very low
# So reaction rate does not limit removal from the gas phase
# So we should see convergence between numerical and closed-form solution

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
# Reference, red line

pH = 7.
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Blue
pH = 10.
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Green line in plots
pH = 11.
pred3 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Orange line in plots
pH = 12.
pred4 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)


# Closed-form solution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For comparison with pH-specific Kga's: (OBS this can be done more elegantly...)
z = np.linspace(0, L, 25)
Kga = pred1['pars']['Kga']
cg_7 = cgin * np.exp(-Kga / v_g * z)
Kga = pred2['pars']['Kga']
cg_10 = cgin * np.exp(-Kga / v_g * z)
Kga = pred3['pars']['Kga']
cg_11 = cgin * np.exp(-Kga / v_g * z)
Kga = pred4['pars']['Kga']
cg_12 = cgin * np.exp(-Kga / v_g * z)

# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Profiles
# Gas
plt.clf()
plt.plot(z, cg_7, 'ro')
plt.plot(z, cg_10, 'bo')
plt.plot(z, cg_11, 'go')
plt.plot(z, cg_12, 'yo')
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 1], 'r')
plt.plot(pred2['cell_pos'], pred2['gas_conc'][:, nt - 1], 'b-')
plt.plot(pred3['cell_pos'], pred3['gas_conc'][:, nt - 1], 'g-')
plt.plot(pred4['cell_pos'], pred4['gas_conc'][:, nt - 1], 'orange')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_gas_conc.png')

# Liquid
plt.figure(2)
plt.plot(z, 0 * cg_12, 'bo')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 1], 'r')
plt.plot(pred2['cell_pos'], pred2['liq_conc'][:, nt - 1], 'b-')
plt.plot(pred3['cell_pos'], pred3['liq_conc'][:, nt - 1], 'g-')
plt.plot(pred4['cell_pos'], pred4['liq_conc'][:, nt - 1], 'orange')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_liq_conc.png')

