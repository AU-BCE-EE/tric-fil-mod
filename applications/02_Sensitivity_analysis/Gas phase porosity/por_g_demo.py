# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shutil.copy('../../../mod_funcs.py', '.')
from mod_funcs import tfmod 


# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 0.51            # Filter length/depth (m) 
por_g = 0.73      # (m3/m3) changed later
por_l = 0.08       # (m3/m3) Close to volume-based estimate
v_g = 106/3600       # superficial gas velocity m/s (chosen to represent ex1)
v_l = 0.8/3600        #liquid superficial velocity m/s (chosen to represent ex1)
nc = 100          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 21.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 0.1        # Reaction rate (1/s). Small because of inert carrier
                 # Reaction could be acid/base that changes the pH

pH = 8.

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

por_g=0.65
Kga=0.02
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
pred1label='Por_g=0.65, Kga=0.02'


por_g=0.7
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
pred2label='Por_g=0.7, Kga=onda'

por_g=0.75
pred3 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
pred3label='Por_g=0.75, Kga=0.02'

por_g=0.8
pred4 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
pred4label='Por_g=0.8, Kga=onda'

por_g=0.85
pred5 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = Kga, k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
pred5label='Por_g=0.85, Kga=0.02'


por_g=0.9
pred6 = tfmod(L = L, por_g = por_g, por_l = por_l,v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
pred6label='Por_g=0.9, Kga=onda'


# Closed-form solution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# z = np.linspace(0, L, 10)
# pr = por_l / por_g
# ctin = por_g * cgin + por_l * clin
# ct = ctin * np.exp(-k  * por_l / (v_g * Kaw) * z)
# cg = Kaw * ct / (por_g * Kaw + por_l)
# cl = cg / Kaw

# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


BT1 = 14.5*por_g/50
BT2 = 14.5*por_g/25


#Outlet concentrations as function of time
#Gas
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :], color='b',label=pred1label)
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :], color='r',label=pred2label)
plt.plot(pred3['time'] / 3600, pred3['gas_conc'][nc - 1, :], color='y',label=pred3label)
plt.plot(pred4['time'] / 3600, pred4['gas_conc'][nc - 1, :], color = 'k',label=pred4label)
plt.plot(pred5['time'] / 3600, pred5['gas_conc'][nc - 1, :], color='c',label=pred5label)
plt.plot(pred6['time'] / 3600, pred6['gas_conc'][nc - 1, :], color='m',label=pred6label)
plt.axvline(x=BT1/60,color='violet',linestyle='-',label='Theoretical breakthrough for vg=106m/h')
#plt.axvline(x=BT2/60,color='deeppink',linestyle='-',label='Theoretical Breakthrough for vg=53m/h)')
plt.axhline(y=cgin,color='g',label='inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.5))
plt.title('Gas Phase NOTE: k=0.1 and pH=8')
plt.show()

#Liquid
plt.plot(pred1['time'] / 3600, pred1['liq_conc'][nc - 1, :], color='b',label=pred1label)
plt.plot(pred2['time'] / 3600, pred2['liq_conc'][nc - 1, :], color='r',label=pred2label)
plt.plot(pred3['time'] / 3600, pred3['liq_conc'][nc - 1, :], color='y',label=pred3label)
plt.plot(pred4['time'] / 3600, pred4['liq_conc'][nc - 1, :], color = 'k',label=pred4label)
plt.plot(pred5['time'] / 3600, pred5['liq_conc'][nc - 1, :], color='c',label=pred5label)
plt.plot(pred6['time'] / 3600, pred6['liq_conc'][nc - 1, :], color='m',label=pred6label)
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Liquid Phase NOTE: k=0.1 and pH=8')
plt.show()