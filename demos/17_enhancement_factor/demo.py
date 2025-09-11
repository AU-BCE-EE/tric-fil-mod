# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../../')
# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from mod_funcs_enhancement import tfmod

# Experimental and physical parameters ______________________________________________

vol = 1200 #Volume of reservoir in mL
pH = 8
temp = 21.       # (degrees C)
v_g = 0.01605 #Gas phase velocity in m/s
v_l = 0.0001085 # liquid phase velocity in m/s
por_l = 0.20498 # Water or liquid content in m3/m3
por_g = 0.795115372 - por_l #gas phase porosity
cgin = 0.057 #inlet gas concentration in g/m3
clin = 0    #inlet liquid concentration in g/m3

#Calculating the cross sectional area, and dividing the volume by it, as required by the model (in units pr cross sectional area)
area = (0.19/2)**2 * 3.14159265
v_res = vol * 10**(-6) / area
L = 0.51            # Filter length/depth (m) 

#Compound specific parameters___________________________________________________________

henry = (0.1, 2000.)  # compound specific Henry constant, from NIST
dens_l = 1000    # Liquid density (kg/m3)
pKa = 7.
# Reaction rate constants in 1/s for protonated (k) and deprotonated (k2) species
k= 0.0185
k2=0

#Model parameters______________________________________________________________________________

# Times for model output in h
tt = 0.35
# Number of time rows
nt = 200
times = np.linspace(0, tt, nt) * 3600
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3). Initial condition
cl0 = 0          # (g/m3). initial condition



# Running the model_______________________________________________________________________________________

pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda' , k = k, k2 = k2, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times, v_res = v_res, recirc = True, counter = True)
pred1label= 'Model 1'# label in plot


# Possibilities for Kga estimation: Onda and KD (Kim and Deshusses, 2009. Remember to specify typ and see Kim and Deshusses in the folder Kga_Models). 
#It is also possible to specify Kga in 1/s directly. 
#Also, if Kga = individual, kl, ae and kg can be specified separately (any parameter not specified will be estimated by Onda correlations)
# pH, cgin and clin can be specified as a pandas dataframe or numpy array as a function of time if needed. 
# This model assumes a fast first order reaction in the liquid phase in its calculation of the Hatta number and enhancement factor  

Name = 'Test' # Used only if saving the plots (as name of the file). A new set of plots will be saved if and only if you change this parameter

#When comparing multiple models, enable multiplot and add the models to be plotted
Multimod = True

if Multimod:
    models = []

    # First model (outside loop, your "baseline" one)
    models.append((pred1, pred1label))

    # Example: add a second model. More can be added. Remember to change the labels if needed
    pred2 = tfmod(
        L=L, por_g=por_g, por_l=por_l, v_g=v_g, v_l=v_l, nc=nc,
        cg0=cg0, cl0=cl0, cgin=cgin, clin=clin, Kga='onda', k=k,
        k2=k2, henry=henry, pKa=pKa, pH=4, temp=temp, dens_l=dens_l,
        times=times, v_res=v_res, recirc=True, counter=True
    )
    pred2label = 'Model 2'
    models.append((pred2, pred2label))
    
    # Example: add a third model. More can be added. Remember to change the labels if needed
    pred3 = tfmod(
        L=L, por_g=por_g, por_l=por_l, v_g=v_g, v_l=v_l, nc=nc,
        cg0=cg0, cl0=cl0, cgin=cgin, clin=clin, Kga='onda', k=k,
        k2=k2, henry=henry, pKa=pKa, pH=8.3, temp=temp, dens_l=dens_l,
        times=times, v_res=v_res, recirc=True, counter=True
    )
    pred3label = 'Model 3'
    models.append((pred3, pred3label))
    
else:
    models = []
    models.append((pred1, pred1label))








#Plotting__________________________________________________________________________________________________

#Concentration versus time plot
plt.clf()
for pred, label in models:
        plt.plot(
            pred['time'] / 3600,
            pred['gas_conc'][nc - 1, :],
            label=label
        )
plt.axhline(cgin,label = 'Inlet concentration',color = 'k') # Plots the inlet concewntration
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0,tt)
plt.ylim(0,)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title(Name)
plt.grid(True)
plt.show()
#plt.savefig(Name+'_concentration_vs_time.png', bbox_inches='tight')
plt.close()




#Concentration versus postition plot for pred1
plt.clf()
plot_times = [0.01,0.1,0.2] # time for the concentration profile in hours. Make sure this is smaller or equal to tt (hours you have chosen to simulate in line 30)! 
                    # It is possible to plot multiple times in the same plot by adding the, eg. [0.1,0.2,0.3]
for plot_time in plot_times:
    idx = (np.abs(times - plot_time*3600)).argmin()  # index closest to plot_time
    plt.plot(
        pred1['cell_pos'],
        pred1['gas_conc'][:, idx],
        label=f"Time = {pred1['time'][idx] / 3600:.2f}h"
    )
plt.xlabel('Height [m]')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0,L)
plt.ylim(0,)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title(pred1label)
plt.grid(True)
plt.show()
#plt.savefig(Name+'_concentration_vs_position.png', bbox_inches='tight')
plt.close()





#Plotting table with input paramters_____________________________________________________________________________
#import module
#If needed, enable plotting of table values for all models
Table = True
if Table:
    from tabulate import tabulate
    
    for pred, label in models:
        parameters =  [['v_g [m/s]', "{:.4f}".format (pred['inputs']['v_g'])], ['v_l [m/s]',"{:.6f}".format (pred['inputs']['v_l'])], ['pH', pred['inputs']['pH']], ['countercurrent', pred['inputs']['counter']],
                      ['recirculation', pred['inputs']['recirc']], ['water content [m^3/m^3]', pred['inputs']['por_l']], ['porosity [m^3/m^3]', pred['inputs']['por_g']], ['temperature [deg. C]', pred['inputs']['temp']], ['v_res [m^3/m^2]', "{:.4f}".format(pred['inputs']['v_res'])],
                      ['k [1/s]', pred['inputs']['k']], ['k2[1/s]', pred['inputs']['k2']],['Kga [1/s]',pred['pars']['Kga']]]
        head = ['Parameter', 'Value']
        table = tabulate(parameters, headers=head, tablefmt="grid")
        
        # Create a figure and axis
        fig, ax = plt.subplots()
        ax.axis('off')
        table_ax = ax.table(cellText=parameters, colLabels=head, loc='center', cellLoc='center')
        table_ax.auto_set_font_size(False)
        table_ax.set_fontsize(8)
        plt.title(label)
        plt.show()
        #plt.savefig(Name+'_'+label+'_table.png')


