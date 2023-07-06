# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

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
v_l = 1E-10
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 15.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 0.1         # Reaction rate (1/s)

pH = 7.

# High pKa so we don't have to consider ionization in closed-form solution
pKa = 100.

# From output of pred1
Kaw = 0.3307 # From pred1 output

# Put inlet concentrations at equilibrium
cgin = 1
clin = cgin / Kaw

# Times for model output
tt = 3 
# Number of time rows
nt = 100
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reference, red line
Kga = 0.001
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Blue
Kga = 0.1
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Green line in plots
Kga = 10
pred3 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Orange line in plots
Kga = 100
pred4 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)


# Closed-form solution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z = np.linspace(0, L, 10)
pr = por_l / por_g
ctin = por_g * cgin + por_l * clin
ct = ctin * np.exp(-k  * por_l / (v_g * Kaw) * z)
cg = Kaw * ct / (por_g * Kaw + por_l)
cl = cg / Kaw

# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Profiles
# Gas
plt.clf()
plt.plot(z, cg, 'bo')
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 1], 'r')
plt.plot(pred2['cell_pos'], pred2['gas_conc'][:, nt - 1], 'b-')
plt.plot(pred3['cell_pos'], pred3['gas_conc'][:, nt - 1], 'g-')
plt.plot(pred4['cell_pos'], pred4['gas_conc'][:, nt - 1], color = 'orange', linestyle = 'dashed')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_gas_conc.png')

# Liquid
plt.clf()
plt.plot(z, cl, 'bo')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 1], 'r')
plt.plot(pred2['cell_pos'], pred2['liq_conc'][:, nt - 1], 'b-')
plt.plot(pred3['cell_pos'], pred3['liq_conc'][:, nt - 1], 'g-')
plt.plot(pred4['cell_pos'], pred4['liq_conc'][:, nt - 1], color = 'orange', linestyle = 'dashed')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_liq_conc.png')
