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
L =0.51            # Filter length/depth (m) (actual 0.51, changed to be able to tell graphs apart)
por_g = 0.80      # (m3/m3) Estimated by volume calculations
por_l = 0.018       # (m3/m3) Estimated at 0.018, changed to be able to tell the graphs apart
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

pH = 6.

# realistic pKa
pKa = 7.

# From output of pred1
Kaw = 0.3732 # From pred1 output (check: should be the same as Kaw1)

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

pH=9
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
Kaw1=pred1['Kaw']


pH=8
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

pH=7.5
pred3 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

pH=7
pred4 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

pH=6.5
pred5 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)


pH=6
pred6 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
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
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 1], color='b',label='pH=7')
plt.plot(pred2['cell_pos'], pred2['gas_conc'][:, nt - 1], color='r',label='pH=6.5')
plt.plot(pred3['cell_pos'], pred3['gas_conc'][:, nt - 1], color='y',label='pH=6')
plt.plot(pred4['cell_pos'], pred4['gas_conc'][:, nt - 1], color = 'k',label='pH=5.5')
plt.plot(pred5['cell_pos'], pred1['gas_conc'][:, nt - 1], color='c',label='pH=5')
plt.plot(pred6['cell_pos'], pred1['gas_conc'][:, nt - 1], color='m',label='pH=4.5')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.title('Gas Phase')
plt.legend()
plt.show()

# Liquid
plt.clf()
plt.plot(z, cl, 'bo')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 1], color='b',label='7')
plt.plot(pred2['cell_pos'], pred2['liq_conc'][:, nt - 1], color='r',label='6.5')
plt.plot(pred3['cell_pos'], pred3['liq_conc'][:, nt - 1], color='y',label='6')
plt.plot(pred4['cell_pos'], pred4['liq_conc'][:, nt - 1], color = 'k',label='5.5')
plt.plot(pred5['cell_pos'], pred1['liq_conc'][:, nt - 1], color='c',label='5')
plt.plot(pred6['cell_pos'], pred1['liq_conc'][:, nt - 1], color='m',label='4.5')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.title('Liquid Phase')
plt.legend()
plt.show()

#Outlet concentrations as function of time
#Gas
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :], color='b',label='pH=9')
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :], color='r',label='pH=8')
plt.plot(pred3['time'] / 3600, pred3['gas_conc'][nc - 1, :], color='y',label='pH=7.5')
plt.plot(pred4['time'] / 3600, pred4['gas_conc'][nc - 1, :], color = 'k',label='pH=7')
plt.plot(pred5['time'] / 3600, pred1['gas_conc'][nc - 1, :], color='c',label='pH=6.5')
plt.plot(pred6['time'] / 3600, pred1['gas_conc'][nc - 1, :], color='m',label='pH=6')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Gas Phase')
plt.show()

#Liquid
plt.plot(pred1['time'] / 3600, pred1['liq_conc'][nc - 1, :], color='b',label='pH=9')
plt.plot(pred2['time'] / 3600, pred2['liq_conc'][nc - 1, :], color='r',label='pH=8')
plt.plot(pred3['time'] / 3600, pred3['liq_conc'][nc - 1, :], color='y',label='pH=7.5')
plt.plot(pred4['time'] / 3600, pred4['liq_conc'][nc - 1, :], color = 'k',label='pH=7')
plt.plot(pred5['time'] / 3600, pred1['liq_conc'][nc - 1, :], color='c',label='pH=6.5')
plt.plot(pred6['time'] / 3600, pred1['liq_conc'][nc - 1, :], color='m',label='pH=6')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Liquid Phase')
plt.show()