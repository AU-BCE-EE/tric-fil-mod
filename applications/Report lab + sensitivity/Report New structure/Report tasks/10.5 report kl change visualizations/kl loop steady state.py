# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sys.path.append("../../../../..")  # Add the directory containing mod_funcs.py to Python path

from mod_funcs import tfmod

#select and define experimental data ______________________________________________
#selecting the experiment number and loading parameters from the master spreadsheet

vol = 700
no = 1
pH = 8


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
#c_in=np.mean([C1,C5],axis=0)

# cgin = pd.DataFrame({'time': t5*3600, 
#                      'cgin': c_in})

c_in = [0.0596 if i < 120 else 0 for i in range(220)]

clin = 0

# Times for model output in h
tt = 1
# Number of time rows
nt = 220
times = np.linspace(0, tt, nt) * 3600

cgin = pd.DataFrame({'time': times, 
                      'cgin': c_in})

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
kl_list = [0.000001, 'onda', 0.000004] 


preds = []
pred_labels = []

for i, kl in enumerate(kl_list):

    
    pred = tfmod(L=L, por_g=por_g, por_l=por_l, v_g=v_g, v_l=v_l, nc=nc, cg0=cg0, cl0=cl0,
                 cgin=cgin, clin=clin, Kga='individual', k=k, henry=henry, pKa=pKa,
                 pH=pH, temp=temp, dens_l=dens_l, times=times, v_res=v_res, kg = 'onda', kl = kl, ae=800, k2=k, 
                 recirc=True, counter=True)
    
    # Removal efficiency
    RE = (0.0596 - pred['gas_conc'][nc - 1, :][119]) / 0.0596
    formatted_RE = "{:.3g}".format(RE)

    
    label = f"kl {kl} m/s outlet gas phase model"
    
    preds.append(pred)
    pred_labels.append(label)
    


#Removal efficiency
RE = (0.0596-pred['gas_conc'][nc-1,:][119])/0.0596
formatted_RE = "{:.3g}".format(RE)

#Plotting__________________________________________________________________________________________________
plt.clf()
for pred, label in zip(preds, pred_labels):
    plt.plot(pred['time'] / 3600, pred['gas_conc'][nc - 1, :], label=label)
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.grid(True)
plt.xlim(0,1)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('pH 8, velocity setting 1')
plt.savefig('Plots/Experiment kl steady state, k = 0.png', bbox_inches='tight')
plt.close()
