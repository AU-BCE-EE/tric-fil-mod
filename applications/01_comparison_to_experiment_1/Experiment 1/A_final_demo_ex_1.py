# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from Mod_Funcs import tfmod  


# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 0.51            # Filter length/depth (m) 
por_g = 0.80      # (m3/m3) Estimated by volume calculations
por_l = 0.018     # (m3/m3) Estimated by volume calculations
v_g = 0.017       #not relevant as this is entered manually
v_l = 1E-4        #not relevant as this is entered manually
nc = 200          # Number of model cells (layers)
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
# Reference, red line
#v_g=60m/h, v_l=0.4m/h
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = 53/3600, v_l = 0.4/3600, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)
Daw=pred1['pars']['Daw']

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







## Breakthrough time calculated from volume (V=14.5L), Por_g=0.8 and volumetric velocities 25l/min (v_g=53m/h) 
# and 50 L/min (v_g=106m/h)

BT1 = 14.5*por_g/50
BT2 = 14.5*por_g/25

# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Profiles
# Gas
#plt.clf()
#plt.plot(z, cg, 'bo')
#plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, nt - 1], color='b',label='v_g=53m/h, v_l=0.4m/h')
#plt.plot(pred2['cell_pos'], pred2['gas_conc'][:, nt - 1], color='r',label='v_g=106m/h, v_l=0.4m/h')
#plt.plot(pred3['cell_pos'], pred3['gas_conc'][:, nt - 1], color='y',label='v_g=106m/h, v_l=1.2m/h')
#plt.plot(pred4['cell_pos'], pred4['gas_conc'][:, nt - 1], color = 'k',label='v_g=53m/h, v_l=1.2m/h')
#plt.plot(pred5['cell_pos'], pred1['gas_conc'][:, nt - 1], color='c',label='v_g=53m/h, v_l=0.8m/h')
#plt.plot(pred6['cell_pos'], pred1['gas_conc'][:, nt - 1], color='m',label='v_g=106m/h, v_l=0.8m/h')
#plt.xlabel('Location (m)')
#plt.ylabel('Compound conc. (g/m3)')
#plt.title('Gas Phase')
#plt.legend()
#plt.show()

# Liquid
#plt.clf()
#plt.plot(z, cl, 'bo')
#plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, nt - 1], color='b',label='v_g=53m/h, v_l=0.4m/h')
#plt.plot(pred2['cell_pos'], pred2['liq_conc'][:, nt - 1], color='r',label='v_g=106m/h, v_l=0.4m/h')
#plt.plot(pred3['cell_pos'], pred3['liq_conc'][:, nt - 1], color='y',label='v_g=106m/h, v_l=1.2m/h')
#plt.plot(pred4['cell_pos'], pred4['liq_conc'][:, nt - 1], color = 'k',label='v_g=53m/h, v_l=1.2m/h')
#plt.plot(pred5['cell_pos'], pred1['liq_conc'][:, nt - 1], color='c',label='v_g=53m/h, v_l=0.8m/h')
#plt.plot(pred6['cell_pos'], pred1['liq_conc'][:, nt - 1], color='m',label='v_g=106m/h, v_l=0.8m/h')
#plt.xlabel('Location (m)')
#plt.ylabel('Compound conc. (g/m3)')
#plt.title('Liquid Phase')
#plt.legend()
#plt.show()

#Outlet concentrations as function of time
#Gas
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :], color='b',label='1.1 v_g=53m/h, v_l=0.4m/h')
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :], color='r',label='1.2 v_g=106m/h, v_l=0.4m/h')
plt.plot(pred3['time'] / 3600, pred3['gas_conc'][nc - 1, :], color='y',label='1. 3v_g=106m/h, v_l=1.2m/h')
plt.plot(pred4['time'] / 3600, pred4['gas_conc'][nc - 1, :], color = 'k',label='1.4 v_g=53m/h, v_l=1.2m/h')
plt.plot(pred5['time'] / 3600, pred5['gas_conc'][nc - 1, :], color='c',label='1.5 v_g=53m/h, v_l=0.8m/h')
plt.plot(pred6['time'] / 3600, pred6['gas_conc'][nc - 1, :], color='m',label='1.6v_g=106m/h, v_l=0.8m/h')
plt.axvline(x=BT1/60,color='violet',linestyle='-',label='Theoretical breakthrough for vg=106m/h')
plt.axvline(x=BT2/60,color='deeppink',linestyle='-',label='Theoretical Breakthrough for vg=53m/h)')
plt.axhline(y=0.055,color='g',label='inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.5))
plt.title('Gas Phase')
plt.show()

#Liquid
plt.plot(pred1['time'] / 3600, pred1['liq_conc'][nc - 1, :], color='b',label='v_g=53m/h, v_l=0.4m/h')
# plt.plot(pred2['time'] / 3600, pred2['liq_conc'][nc - 1, :], color='r',label='v_g=106m/h, v_l=0.4m/h')
# plt.plot(pred3['time'] / 3600, pred3['liq_conc'][nc - 1, :], color='y',label='v_g=106m/h, v_l=1.2m/h')
# plt.plot(pred4['time'] / 3600, pred4['liq_conc'][nc - 1, :], color = 'k',label='v_g=53m/h, v_l=1.2m/h')
# plt.plot(pred5['time'] / 3600, pred5['liq_conc'][nc - 1, :], color='c',label='v_g=53m/h, v_l=0.8m/h')
plt.plot(pred6['time'] / 3600, pred6['liq_conc'][nc - 1, :], color='m',label='v_g=106m/h, v_l=0.8m/h')
#equilibrium concentrations
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :]/Daw, color='b',linestyle='dashed')
# plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :]/Daw, color='r',linestyle='dashed')
# plt.plot(pred3['time'] / 3600, pred3['gas_conc'][nc - 1, :]/Daw, color='y',linestyle='dashed')
# plt.plot(pred4['time'] / 3600, pred4['gas_conc'][nc - 1, :]/Daw, color = 'k',linestyle='dashed')
# plt.plot(pred5['time'] / 3600, pred5['gas_conc'][nc - 1, :]/Daw, color='c',linestyle='dashed')
plt.plot(pred6['time'] / 3600, pred6['gas_conc'][nc - 1, :]/Daw, color='m',linestyle='dashed',label='All dashed are equilibrium values')
plt.axvline(x=BT1/60,color='violet',linestyle='-',label='Theoretical breakthrough for vg=106m/h')
plt.axvline(x=BT2/60,color='deeppink',linestyle='-',label='Theoretical Breakthrough for vg=53m/h)')
plt.axhline(y=0.055/Daw,color='g', linestyle='dashed',label='equilibrium with inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.5))
plt.title('Liquid Phase model and equilibrium values')
plt.show()

#Mass of H2S instead of concentration
plt.plot(pred1['time'] / 3600, pred1['gas_mass'][nc - 1, :], color='b',label='Gas Phase')
plt.plot(pred1['time'] / 3600, pred1['liq_mass'][nc - 1, :], color='r',label='Liquid Phase')
plt.legend()
plt.title('Comparison of mass flows as predicted by the model')
plt.xlabel('time(h)')
plt.ylabel('mass (g) or (g/m^2)')
plt.show()


#comparison to data experiment 1
from Data_treatment_ex1 import t1
from Data_treatment_ex1 import C_out1
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :], color='b',label='model')
plt.plot(t1,C_out1,label='experimental data')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Experiment 1.1')
plt.show()

#comparison to data experiment 2
from Data_treatment_ex1 import t2
from Data_treatment_ex1 import C_out2
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :], color='b',label='model')
plt.plot(t2,C_out2,label='experimental data')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Experiment 1.2')
plt.show()

#comparison to data experiment 3
from Data_treatment_ex1 import t3
from Data_treatment_ex1 import C_out3
plt.plot(pred3['time'] / 3600, pred3['gas_conc'][nc - 1, :], color='b',label='model')
plt.plot(t3,C_out3,label='experimental data')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Experiment 1.3')
plt.show()


#comparison to data experiment 4
from Data_treatment_ex1 import t4
from Data_treatment_ex1 import C_out4
plt.plot(pred4['time'] / 3600, pred4['gas_conc'][nc - 1, :], color='b',label='model')
plt.plot(t4,C_out4,label='experimental data')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Experiment 1.4')
plt.show()

#comparison to data experiment 5
from Data_treatment_ex1 import t5
from Data_treatment_ex1 import C_out5
plt.plot(pred5['time'] / 3600, pred5['gas_conc'][nc - 1, :], color='b',label='model')
plt.plot(t5,C_out5,label='experimental data')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Experiment 1.5')
plt.show()


#comparison to data experiment 6
from Data_treatment_ex1 import t6
from Data_treatment_ex1 import C_out6
plt.plot(pred6['time'] / 3600, pred6['gas_conc'][nc - 1, :], color='b',label='model')
plt.plot(t6,C_out6,label='experimental data')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.title('Experiment 1.6')
plt.show()


#All experimental data in one graph


plt.plot(t1,C_out1,color='b', label='1.1 v_g=53m/h, v_l=0.4m/h')
plt.plot(t2,C_out2,color='r',label='1.2 v_g=106m/h, v_l=0.4m/h')
plt.plot(t3,C_out3,color='y',label='1.3 v_g=106m/h, v_l=1.2m/h')
plt.plot(t4,C_out4,color='k',label='1.4 v_g=53m/h, v_l=1.2m/h')
plt.plot(t5,C_out5,color='c',label='1.5 v_g=53m/h, v_l=0.8m/h')
plt.plot(t6,C_out6,color='m',label='1.6 v_g=106m/h, v_l=0.8m/h')
plt.axvline(x=BT1/60,color='violet',linestyle='-',label='Theoretical breakthrough for vg=106m/h')
plt.axvline(x=BT2/60,color='deeppink',linestyle='-',label='Theoretical Breakthrough for vg=53m/h)')
plt.axhline(y=0.055,color='g',label='inlet concentration')
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.5))
plt.title('Experiment 1 all experimental data')
plt.xlabel('Time(h)')
plt.ylabel('concentration (g/m3)')
plt.show()
