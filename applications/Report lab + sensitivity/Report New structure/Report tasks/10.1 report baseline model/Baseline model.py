# Making a baseline to compare changes in parameters (pH, k, kl, temp) to

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sys.path.append("../../../../..")  # Add the directory containing mod_funcs.py to Python path

from mod_funcs import tfmod

vol = 700
no = 1
pH = 8


# Times for model output in h
tt = 0.35
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600


# Define y values for the expected inlet curve @40ppm
expected_inlet = np.full_like(times / 60, 0.0596)  

# Make the line drop to 0 at x=5, as the oulse is 5min
expected_inlet[times / 60 >= 5] = 0




#Calculating the cross sectional area, and dividing the volume by it, as required by the model
area = (0.19/2)**2 * 3.14159265
v_res = vol * 10**(-6) / area

#Setting the liquid and gas velocities using the experiment number (second number in the file name)
if no == 1 or no == 4:
    v_g = 27.2991 / 1000 / 60 / area # gas velocity using calibration curve for mass flow controllers L/min /1000L/m3 / 60 s/min / m2 = m/s
elif no ==2 or no == 3:
    v_g = 54.2766 /1000 / 60 /area # gas velocity using calibration curve for mass flow controllers L/min /1000L/m3 / 60 s/min / m2 = m/s
else: print ('Error in reading "no"')

if no == 1 or no == 2:
    v_l = 0.390681119/3600 #in m/s, based on pump calibration number 2
elif no ==4 or no == 3:
    v_l = 1.253842217/3600 #in m/s, based on pump calibration number 2
else: print ('Error in reading "no"')

# Porosity and water content definition
#from water content experiments, as well as porosity experiments. Average of three measurements. gas phase is porosity minus water content
if no == 1:
    por_l = 0.20498
    por_g = 0.795115372 - por_l
elif no == 2:
    por_l = 0.22494
    por_g = 0.795115372 - por_l
elif no == 3:
    por_l = 0.24056
    por_g = 0.795115372 - por_l
elif no == 4:
    por_l = 0.246304
    por_g = 0.795115372 - por_l
else: print ('Error in reading "no"') 





## Breakthrough time calculated from volume (V=14.5L), Por_g and flow rate (L/min, from calibration of mass flow controllers)

BT1 = 14.5*por_g/54.2766
BT2 = 14.5*por_g/27.2991
BTlabel='Theoretical Breakthrough curve)'




#Defining inputs for the model______________________________________________________________________
#Time and inlet concentrations
c_in=expected_inlet

cgin = pd.DataFrame({'time': times, 
                     'cgin': c_in})

clin = 0

# Times for model output in h
tt = 0.35
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600

# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 0.51            # Filter length/depth (m) 
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 21.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)
# realistic pKa
pKa = 7.

k=0 


   
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, k2 = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = v_res, recirc = True, counter = True)
pred1label= 'Outlet concentration'

print(pred1['pars']['gas_rt'])

#Plotting__________________________________________________________________________________________________
plt.clf()
plt.plot(pred1['time'] / 60, pred1['gas_conc'][nc - 1, :],label=pred1label)
plt.xlabel('Time (min)')
plt.ylabel('Compound conc. (g/m3)')
plt.plot(times / 60, expected_inlet, color='g', label='Inlet concentration')
plt.legend()
plt.grid(True)
plt.xlim(0,20)
plt.ylim(0,0.07)
plt.title('Baseline model Gas phase')
plt.savefig('../..//Plots for report/Baseline model gas.png', bbox_inches='tight')
plt.close()


plt.clf()
plt.plot(pred1['time'] / 60, pred1['liq_conc'][1, :],label='Outlet concentration')
plt.plot(pred1['time'] / 60, pred1['liq_conc'][nc - 1, :],label='Inlet concentration',color = 'g')
plt.xlabel('Time (min)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.grid(True)
plt.xlim(0,20)
plt.ylim(0,2)
plt.title('Baseline model Liquid phase')
plt.savefig('../..//Plots for report/Baseline model liquid.png', bbox_inches='tight')
plt.close()



Daw = pred1['pars']['Daw']
times = np.linspace(0, tt, nt) * 60 #Time in minutes

plt.clf()
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, 10], label = 'c',color = 'b')
plt.plot(pred1['cell_pos'], (pred1['liq_conc'][:, 10])*Daw, label = 'c*',color='b',linestyle = 'dashed')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.grid(True)
plt.xlim(0,)
plt.ylim(0,0.07)
plt.legend()
plt.title('After %1.2f'%times[10]+' minutes')
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.savefig('../..//Plots for report/Baseline model profile at 1min.png', bbox_inches='tight')
plt.close()

plt.clf()
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, 50], label = 'c',color = 'b')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, 50]*Daw, label = 'c*',color='b',linestyle = 'dashed')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.grid(True)
plt.xlim(0,)
plt.ylim(0,0.07)
plt.legend()
plt.title('After %1.2f'%times[50]+' minutes')
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.savefig('../..//Plots for report/Baseline model profile at 5min.png', bbox_inches='tight')
plt.close()



# Gas
plt.clf()
plt.plot(pred1['cell_pos'], pred1['gas_conc'][:, 100], label = 'c',color ='b')
plt.plot(pred1['cell_pos'], pred1['liq_conc'][:, 100]*Daw, label = 'c*',color='b',linestyle = 'dashed')
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.grid(True)
plt.xlim(0,)
plt.ylim(0,0.07)
plt.legend()
plt.title('After %1.2f'%times[100]+' minutes')
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.savefig('../..//Plots for report/Baseline model profiles at 11min.png', bbox_inches='tight')
plt.close()



