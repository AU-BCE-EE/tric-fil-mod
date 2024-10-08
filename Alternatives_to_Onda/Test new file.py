# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 11:33:40 2024

@author: au611147
"""
import sys
import numpy as np
import matplotlib.pyplot as plt


# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sys.path.append("../")  # Add the directory containing mod_funcs.py to Python path
from mod_funcs import tfmod 



clin = 0 #inlet concentration for the liquid phase
cgin = 0.0596 #inlet gas phase concentration of H2S

# Times for model output in h
tt = 0.35
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600

L = 0.51            # Filter length/depth (m) 
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 21.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

#Liquid and gas phase porosities
por_l = 0.24056
por_g = 0.795115372 - por_l

#calculations regarding the reservoir used for recirculation
vol = 700 #mL of reservoir
area = (0.19/2)**2 * 3.14159265
v_res = vol * 10**(-6) / area
#Gas and liquid phase velocities, in m/s, as in ex3 in Annes Thesis
v_g = 54.2766 / 1000 / 60 / area
v_l = 1.253842217/3600


pKa = 7 #for H2S
pH = 8 #arbitrary
k = 1 # reaction rate constant in 1/s, arbitrary

#Simulations

pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = v_res, k2=k, recirc = True, counter = True)
pred1label= 'Onda prediction' #label on

pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'KD', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = v_res, k2=k, typ = 'LR', recirc = True, counter = True)
pred2label= 'Kim and Deshusses Lava rock prediction' #label on

pred3 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'KD', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = v_res, k2=k, typ = 'PUF', recirc = True, counter = True)
pred3label= 'Kim and Deshusses PUF prediction' #label on

pred4 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'KD', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = v_res, k2=k, typ = 'PR', recirc = True, counter = True)
pred4label= 'Kim and Deshusses Pall ring prediction' #label on

pred5 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'KD', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = v_res, k2=k, typ = 'PCB', recirc = True, counter = True)
pred5label= 'Kim and Deshusses porous ceramic bead prediction' #label on

pred6 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'KD', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = v_res, k2=k, typ = 'PCR', recirc = True, counter = True)
pred6label= 'Kim and Deshusses porous ceramic ring prediction' #label on


#Plotting
plt.clf()
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :],color='k',label=pred1label)
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :],color='m',label=pred2label)
plt.plot(pred3['time'] / 3600, pred3['gas_conc'][nc - 1, :],color='b',label=pred3label)
plt.plot(pred4['time'] / 3600, pred4['gas_conc'][nc - 1, :],color='c',label=pred4label)
plt.plot(pred5['time'] / 3600, pred5['gas_conc'][nc - 1, :],color='r',label=pred5label)
plt.plot(pred6['time'] / 3600, pred6['gas_conc'][nc - 1, :],color='orange',label=pred6label)
plt.axhline(y=0.0596,color='g',label='Inlet concentration') # Average of two inlet concentrations of 43.2 and 42.3. These were for low and high gas flow rate respectively
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0,0.35)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Comparison between Onda and KD')
plt.show()
#plt.savefig('Plots/Onda vs KD.png', bbox_inches='tight')
plt.close()