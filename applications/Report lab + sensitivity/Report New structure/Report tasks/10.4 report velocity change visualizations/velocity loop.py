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


pH = 8
vol = 700
no = 1


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

k = 0
no_list = [1,2,3,4] 


preds = []
pred_labels = []

for i, no in enumerate(no_list):
    
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


    
    pred = tfmod(L=L, por_g=por_g, por_l=por_l, v_g=v_g, v_l=v_l, nc=nc, cg0=cg0, cl0=cl0,
                 cgin=cgin, clin=clin, Kga='individual', k=k, k2=k, henry=henry, pKa=pKa,
                 pH=pH, temp=temp, dens_l=dens_l, times=times, v_res=v_res, kg = 'onda', kl = 'onda', ae='onda', 
                 recirc=True, counter=True)
    
    label = f"velocity setting {no} outlet concentration"
    
    preds.append(pred)
    pred_labels.append(label)




#Plotting__________________________________________________________________________________________________
colors = ['C0','orange','r','m']
plt.clf()
for i, (pred, label) in enumerate(zip(preds, pred_labels)):
    plt.plot(pred['time'] / 60, pred['gas_conc'][nc - 1, :], label=label,color=colors[i % len(colors)])
#plt.plot(pred1['time'] / 60, pred1['gas_conc'][nc - 1, :],label=pred1label)
plt.xlabel('Time (min)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.plot(times / 60, expected_inlet, color='g', label='Inlet concentration')
plt.grid(True)
plt.xlim(0,20)
plt.ylim(0,0.07)
plt.title('Velocity variations')
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.savefig('Plots/velocity variations.png', bbox_inches='tight')
plt.close()




