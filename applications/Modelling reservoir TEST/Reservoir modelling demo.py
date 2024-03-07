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
L = 0.48            # Filter length/depth (m) 
por_g = 0.5      # (m3/m3) TBD
por_l = 0.25     # (m3/m3) TBD
v_g = 0.017       #not relevant
v_l = 1E-4        #not relevant
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 20.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 0         # Reaction rate (1/s). Small because of inert carrier.
                    #Reaction could be acid/base that changes the pH

pH = 7.

# realistic pKa
pKa = 7.



# Put inlet concentrations at equilibrium
cgin = 0.055  #corresponding to 15ppm
clin = 0

# Times for model output in h
tt = 1 
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reference, red line
#v_g=60m/h, v_l=0.4m/h
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = 60/3600, v_l = 0.4/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times,recirc = True)
pred1label = 'v_res = 0'

# Blue
#v_g=60m/h, v_l=1m/h
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = 60/3600, v_l = 0.4/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = 0.0018, counter = False, recirc = True)
pred2label = 'v_res = 0.0018m3'

# Green line in plots
#v_g=100m/h, v_l=0.4m/h
pred3 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = 60/3600, v_l = 0.4/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = 0.01, recirc = True)
pred3label = 'v_res = 0.01m3'

# Orange line in plots
#v_g=100m/h, v_l=1m/h
pred4 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = 60/3600, v_l = 0.4/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = 0.1, recirc = True)
pred4label = 'v_res = 0.1m3'


# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Profiles
# Gas
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 1], 'r',label=pred1label)
plt.plot(pred2['cell_pos'], pred2['gas_conc'][:, nt - 1], 'b-',label=pred2label)
plt.plot(pred3['cell_pos'], pred3['gas_conc'][:, nt - 1], 'g-',label=pred3label)
plt.plot(pred4['cell_pos'], pred4['gas_conc'][:, nt - 1], color = 'orange', linestyle = 'dashed',label=pred4label)
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.title('Gas Phase')
plt.legend()
plt.savefig('gas_phase_profile.png')
plt.show()


# Liquid
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 1], 'r',label=pred1label)
plt.plot(pred2['cell_pos'], pred2['liq_conc'][:, nt - 1], 'b-',label=pred2label)
plt.plot(pred3['cell_pos'], pred3['liq_conc'][:, nt - 1], 'g-',label=pred3label)
plt.plot(pred4['cell_pos'], pred4['liq_conc'][:, nt - 1], color = 'orange', linestyle = 'dashed',label=pred4label)
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.title('Liquid Phase')
plt.legend()
plt.savefig('liq_phase_profile.png')
plt.show()

#Outlet concentrations as function of time
#Gas
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :], 'r',label=pred1label)
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :], 'r',linestyle='dashed',label=pred2label)
plt.plot(pred3['time'] / 3600, pred3['gas_conc'][nc - 1, :], 'b',label=pred3label)
plt.plot(pred4['time'] / 3600, pred4['gas_conc'][nc - 1, :], 'b',linestyle='dashed',label=pred4label)
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Gas Phase')
plt.savefig('gas_phase_time')
plt.show()

#Liquid
plt.plot(pred1['time'] / 3600, pred1['liq_conc'][nc - 1, :], 'r',label=pred1label)
plt.plot(pred2['time'] / 3600, pred2['liq_conc'][nc - 1, :], 'r',linestyle='dashed',label=pred2label)
plt.plot(pred3['time'] / 3600, pred3['liq_conc'][nc - 1, :], 'b',label=pred3label)
plt.plot(pred4['time'] / 3600, pred4['liq_conc'][nc - 1, :], 'b',linestyle='dashed',label=pred4label)
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Liquid Phase')
plt.savefig('liq_phase_time.png')
plt.show()