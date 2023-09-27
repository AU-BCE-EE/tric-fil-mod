# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from Annes_Playground_mod_funcs import tfmod  


# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 0.51            # Filter length/depth (m) 
por_g = 0.5      # (m3/m3) TBD
por_l = 0.25     # (m3/m3) TBD
v_g = 0.017       #not relevant
v_l = 1E-4        #not relevant
nc = 200          # Number of model cells (layers)
cg0 = 0.0086          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 21.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 0        # Reaction rate (1/s). Small because of inert carrier
                 # Reaction could be acid/base that changes the pH

pH = 7.

# realistic pKa
pKa = 7.

# From output of pred1
Kaw = 0.3659 # From pred1 output (check: should be the same as Kaw1)

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
# Reference, red line
#v_g=60m/h, v_l=0.4m/h
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = 53/3600, v_l = 0.4/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
Kaw1=pred1['Kaw']

# Blue
#v_g=60m/h, v_l=1m/h
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = 106/3600, v_l = 0.4/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Green line in plots
#v_g=100m/h, v_l=0.4m/h
pred3 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = 106/3600, v_l = 1.2/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Orange line in plots
#v_g=100m/h, v_l=1m/h
pred4 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = 53/3600, v_l = 1.2/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)


pred5 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = 53/3600, v_l = 0.8/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)



pred6 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = 106/3600, v_l = 0.8/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
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
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 1], color='b',label='v_g=53m/h, v_l=0.4m/h')
plt.plot(pred2['cell_pos'], pred2['gas_conc'][:, nt - 1], color='r',label='v_g=106m/h, v_l=0.4m/h')
plt.plot(pred3['cell_pos'], pred3['gas_conc'][:, nt - 1], color='y',label='v_g=106m/h, v_l=1.2m/h')
plt.plot(pred4['cell_pos'], pred4['gas_conc'][:, nt - 1], color = 'k',label='v_g=53m/h, v_l=1.2m/h')
plt.plot(pred5['cell_pos'], pred1['gas_conc'][:, nt - 1], color='c',label='v_g=53m/h, v_l=0.8m/h')
plt.plot(pred6['cell_pos'], pred1['gas_conc'][:, nt - 1], color='m',label='v_g=106m/h, v_l=0.8m/h')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.title('Gas Phase')
plt.legend()
plt.show()

# Liquid
plt.clf()
plt.plot(z, cl, 'bo')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 1], color='b',label='v_g=53m/h, v_l=0.4m/h')
plt.plot(pred2['cell_pos'], pred2['liq_conc'][:, nt - 1], color='r',label='v_g=106m/h, v_l=0.4m/h')
plt.plot(pred3['cell_pos'], pred3['liq_conc'][:, nt - 1], color='y',label='v_g=106m/h, v_l=1.2m/h')
plt.plot(pred4['cell_pos'], pred4['liq_conc'][:, nt - 1], color = 'k',label='v_g=53m/h, v_l=1.2m/h')
plt.plot(pred5['cell_pos'], pred1['liq_conc'][:, nt - 1], color='c',label='v_g=53m/h, v_l=0.8m/h')
plt.plot(pred6['cell_pos'], pred1['liq_conc'][:, nt - 1], color='m',label='v_g=106m/h, v_l=0.8m/h')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.title('Liquid Phase')
plt.legend()
plt.show()

#Outlet concentrations as function of time
#Gas
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :], color='b',label='v_g=53m/h, v_l=0.4m/h')
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :], color='r',label='v_g=106m/h, v_l=0.4m/h')
plt.plot(pred3['time'] / 3600, pred3['gas_conc'][nc - 1, :], color='y',label='v_g=106m/h, v_l=1.2m/h')
plt.plot(pred4['time'] / 3600, pred4['gas_conc'][nc - 1, :], color = 'k',label='v_g=53m/h, v_l=1.2m/h')
plt.plot(pred5['time'] / 3600, pred1['gas_conc'][nc - 1, :], color='c',label='v_g=53m/h, v_l=0.8m/h')
plt.plot(pred6['time'] / 3600, pred1['gas_conc'][nc - 1, :], color='m',label='v_g=106m/h, v_l=0.8m/h')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Gas Phase')
plt.show()

#Liquid
plt.plot(pred1['time'] / 3600, pred1['liq_conc'][nc - 1, :], color='b',label='v_g=53m/h, v_l=0.4m/h')
plt.plot(pred2['time'] / 3600, pred2['liq_conc'][nc - 1, :], color='r',label='v_g=106m/h, v_l=0.4m/h')
plt.plot(pred3['time'] / 3600, pred3['liq_conc'][nc - 1, :], color='y',label='v_g=106m/h, v_l=1.2m/h')
plt.plot(pred4['time'] / 3600, pred4['liq_conc'][nc - 1, :], color = 'k',label='v_g=53m/h, v_l=1.2m/h')
plt.plot(pred5['time'] / 3600, pred1['liq_conc'][nc - 1, :], color='c',label='v_g=53m/h, v_l=0.8m/h')
plt.plot(pred6['time'] / 3600, pred1['liq_conc'][nc - 1, :], color='m',label='v_g=106m/h, v_l=0.8m/h')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Liquid Phase')
plt.show()