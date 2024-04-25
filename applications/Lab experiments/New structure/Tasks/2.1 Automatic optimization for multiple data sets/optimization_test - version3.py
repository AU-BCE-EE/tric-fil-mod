 # Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import least_squares


# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sys.path.append("../../../../..")  # Add the directory containing mod_funcs.py to Python path
from mod_funcs import tfmod 


# Create an empty dictionary to store the DataFrames
results_dict = {}
modelresults_dict = {}

#This version only works for experiment5, which is where pH was varied. SLight modifications are needed before 
#the velocities can be changed (experimetn 6)

#Choose two experiments from experiment 5 that we want to fit. (options: 5.1 through 5.7)
exp_no1 = 4
exp_no2 = 7

IG_kl = 0.1
IG_k = 0.1
IG_k2 = 0.1
diff_step = 10
    
for i in [exp_no1,exp_no2]:
        second = str(i)
        first = '5'

        #parameters loaded. File name changes depending on experiment no.
        file_path = '../../../Master_spreadsheet.xlsx'
        df = pd.read_excel(file_path, sheet_name = 'Data')
        
        params = df[df['key'] == float(first)+0.1*float(second)]
        
        cycle1 = int(params['cycle1'].values[0])
        cycle2 = int(params['cycle2'].values[0])
        cycle3 = int(params['cycle3'].values[0])
        cycle5 = int(params['cycle5'].values[0])
        length = params['length'].values[0]
        
        
        
        # Inlet concentrations load in from data treatment and definition

        ex1 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')


        #second inlet profile (last data obtained)
        ex5 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet2.csv', sep = ',')


        # The outlet data loaded
        ex2 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.1.csv', sep = ',')
        ex3 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.2.csv', sep = ',')
        
        


        # selecting and normalising the time, because all data files are started before the H2S
        # was turned on, and not after a specific time, it all depends on the individual experient. 
        #However the cycle at which the H2S was turned on was written down.
        t1norm = ex1['Time in h'] - ex1['Time in h'][cycle1]
        t1 = t1norm [cycle1:cycle1+length]
        C1= ex1['Concentration in g/m^3'][cycle1:cycle1+length]

        t5norm = ex5['Time in h'] - ex5['Time in h'][cycle5]
        t5 = t5norm [cycle5:cycle5+length]
        C5= ex5['Concentration in g/m^3'][cycle5:cycle5+length]

        t2norm = ex2['Time in h'] - ex2['Time in h'][cycle2]
        t2 = t2norm [cycle2:cycle2+length]
        C2= ex2['Concentration in g/m^3'][cycle2:cycle2+length]

        t3norm = ex3['Time in h'] - ex3['Time in h'][cycle3]
        t3 = t3norm [cycle3:cycle3+length]
        C3= ex3['Concentration in g/m^3'][cycle3:cycle3+length]
        
        #Taking the average of both inlet and outlet measurements and exporting them to results
        c_out=np.mean([C2,C3],axis=0)
        
        c_in=np.mean([C1,C5],axis=0)

        results = pd.DataFrame ({'experimental_out':c_out , 'experimental_in':c_in, 'experimental time (h)':t3 })
        # Store the DataFrame in the dictionary with a key for each round
        results_dict[f'experimental{first}.{second}'] = results

        


#Load parameters from the Master speadsheet that the model needs
#parameters loaded. File name changes depending on experiment no.
file_path = '../../../Master_spreadsheet.xlsx'
df = pd.read_excel(file_path, sheet_name = 'Data')

first = '5'        
params = df[df['key'] == float(first)+0.1*exp_no1]
        
pH1 = (params['pH1'].values[0]+params['pH2'].values[0])/2        
vol1 = vol = params['volume'].values[0]

params = df[df['key'] == float(first)+0.1*exp_no2]
        
pH2 = (params['pH1'].values[0]+params['pH2'].values[0])/2 
vol2 = vol = params['volume'].values[0]

#Choose the velocities (800rpm and 25 L/min) Needs to be a loop when velocities start to change
no = 1

# realistic pKa
pKa = 7.


#Calculating the cross sectional area, and dividing the volume by it, as required by the model
area = (0.19/2)**2 * 3.14159265
v_res1 = vol1 * 10**(-6) / area
v_res2 = vol2 * 10**(-6) / area

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
BTlabel='Theoretical Breakthrough curve'

if no == 1 or no == 4:
    BT = BT2
elif no ==2 or no == 3:
    BT= BT1
else: 
    print ('error in definition of "no". Has to be a string containing the numbers 1,2,3 or 4')    
        
    
    
# Optimization function________________________________________________________________
# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 0.51            # Filter length/depth (m) 
nc = 200          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 21.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)
clin = 0


#Load time and inlet concentrations and outlet concentrations for both data sets
times2 = results_dict['experimental5.'+str(exp_no2)]['experimental time (h)'] *3600
c_meas2 = results_dict['experimental5.'+str(exp_no2)]['experimental_out'] 
c_in2 = results_dict['experimental5.'+str(exp_no2)]['experimental_in']
cgin2 = pd.DataFrame({'time': times2, 
                     'cgin': c_in2})

times1 = results_dict['experimental5.'+str(exp_no1)]['experimental time (h)'] *3600
c_meas1 = results_dict['experimental5.'+str(exp_no1)]['experimental_out'] 
c_in1 = results_dict['experimental5.'+str(exp_no1)]['experimental_in']
cgin1 = pd.DataFrame({'time': times1, 
                     'cgin': c_in1})


# Optimization loop res_calc________________________________________________________________________
def res_calc(x):
    # import parameters for the model into the optimization definition
    L, por_g, por_l, v_g, v_l, nc, cg0, cl0, cgin1, cgin2, clin,henry, pKa, pH1, pH2 , temp, dens_l, v_res1, v_res2, times2,times1, c_meas2,c_meas1 
    
    #Define wat is in each position in the answer array x

    kl = abs(x[0])
    k = abs(x[1])
    k2 = abs(x[2])
    
    #Model number 1 + residue calculation
    pred_opt1= tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
                  cl0 = cl0, cgin = cgin1, clin = clin, Kga = 'individual', k = k, k2 = k2, henry = henry, pKa = pKa, 
                  pH = pH1, temp = temp, dens_l = dens_l, times = times1, v_res = v_res1, kl = kl, kg = 'onda', ae = 800, recirc = True, counter = True)
    c1 = pred_opt1['gas_conc'][nc - 1,:]
    r = c1 - c_meas1
    
    #Model number 2 + residue calculation + appending the two residue series

    pred_opt2= tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
                  cl0 = cl0, cgin = cgin2, clin = clin, Kga = 'individual', k = k, k2 = k2, henry = henry, pKa = pKa, 
                  pH = pH2, temp = temp, dens_l = dens_l, times = times2, v_res = v_res2, kl = kl, kg = 'onda', ae = 800, recirc = True, counter = True)
    c2 = pred_opt2['gas_conc'][nc - 1,:]
    r = np.append(r,c2 - c_meas2)
           
    return r

#Finding the solution that gives the least combined deviation between model and experiments
sol = least_squares(res_calc, [IG_kl, IG_k, IG_k2], xtol = 3e-16, ftol = 3e-16, gtol = 3e-16, diff_step = diff_step)

copt = sol['x']




# Models for plotting_______________________________________________________________________________________
kl_fit = abs(float(copt[0]))
k_fit = abs(float(copt[1]))
k2_fit = abs(float(copt[2]))

#First experiment___________________________________________________________________________________________
second = str(exp_no1)

pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin1, clin = clin, Kga = 'individual' , k = k_fit, k2 = k2_fit, henry = henry, pKa = pKa, 
              pH = pH1, temp = temp, dens_l = dens_l, times = times1, v_res = v_res1, kl=kl_fit, kg='onda', ae=800, recirc = True, counter = True)
pred1label= first+'.'+second+'.k and Kga optimized model' #label on


k=0
times = np.linspace (0,0.35,200) *3600
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin1, clin = clin, Kga = 'onda', k = k, k2 = k, henry = henry, pKa = pKa, 
              pH = pH1, temp = temp, dens_l = dens_l, times = times, v_res = v_res1, recirc = True, counter = True)
pred2label= first+'.'+second+'k=0, onda model' #label on

#Plotting__________________________________________________________________________________________________
plt.clf()
plt.plot(times1/3600,c_meas1,label=first+'.'+second+' experimental data')
plt.plot(times1/3600,c_in1,label='inlet')
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :],color='k',label=pred1label)
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :],label=pred2label)
plt.axvline(x=BT/60,linestyle='-',label=BTlabel) #breakthrough curve
plt.axhline(y=0.0596,color='g',label='Expected inlet concentration') # Average of two inlet concentrations of 43.2ppm and 42.3ppm. These were for low and high gas flow rate respectively
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0,0.35)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Experiment '+first+'.'+second)
plt.savefig('Plots/Experiment '+first+'.'+second+'optimized [kl, k, k2], [IG] ['+str(IG_kl)+','+str(IG_k)+','+str(IG_k2)+'].png', bbox_inches='tight')
plt.close()

#Second experiment___________________________________________________________________________________________________
second = str(exp_no2)

pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin2, clin = clin, Kga = 'individual' , k = k_fit, k2 = k2_fit, henry = henry, pKa = pKa, 
              pH = pH2, temp = temp, dens_l = dens_l, times = times2, v_res = v_res2, kl=kl_fit, kg='onda', ae=800, recirc = True, counter = True)
pred1label= first+'.'+second+'.k and Kga optimized model' #label on


k=0
times = np.linspace (0,0.35,200) *3600
pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin2, clin = clin, Kga = 'onda', k = k, k2 = k, henry = henry, pKa = pKa, 
              pH = pH2, temp = temp, dens_l = dens_l, times = times, v_res = v_res2, recirc = True, counter = True)
pred2label= first+'.'+second+'k=0, onda model' #label on

#Plotting__________________________________________________________________________________________________
plt.clf()
plt.plot(times2 / 3600,c_meas2,label=first+'.'+second+'.1 experimental data')
plt.plot(times2 / 3600,c_in2,label='inlet')
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :],color='k',label=pred1label)
plt.plot(pred2['time'] / 3600, pred2['gas_conc'][nc - 1, :],label=pred2label)
plt.axvline(x=BT/60,linestyle='-',label=BTlabel) #breakthrough curve
plt.axhline(y=0.0596,color='g',label='Expected inlet concentration') # Average of two inlet concentrations of 43.2ppm and 42.3ppm. These were for low and high gas flow rate respectively
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0,0.35)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Experiment '+first+'.'+second)
plt.savefig('Plots/Experiment '+first+'.'+second+'optimized [kl, k, k2], [IG] ['+str(IG_kl)+','+str(IG_k)+','+str(IG_k2)+'].png', bbox_inches='tight')
plt.close()


#table with input paramters________________________________________________________________________________
#import module
from tabulate import tabulate
parameters =  [['v_g [m/s]', "{:.4f}".format (pred1['inputs']['v_g'])], ['v_l [m/s]',"{:.6f}".format (pred1['inputs']['v_l'])], ['pH1', "{:.2f}".format(pH1)], ['pH2', "{:.2f}".format(pH2)], ['countercurrent', pred1['inputs']['counter']],
              ['recirculation', pred1['inputs']['recirc']], ['water content [m^3/m^3]', "{:.2f}".format(por_l)], ['porosity [m^3/m^3]', "{:.2f}".format(por_g)], ['temperature [deg. C]', "{:.0f}".format(temp)], ['v_res [m^3/m^2]', "{:.4f}".format(pred1['inputs']['v_res'])],
              ['kl_fit[1/s]',(pred1['pars']['kl'])],['kg_onda[1/s]',"{:.4f}".format(pred1['pars']['kg'])],['ae[m^2/m^3]',pred1['pars']['ae']],['Kga_fit[1/s]',"{:.4f}".format(pred1['pars']['Kga'])], ['Kga_baseline[1/s]',"{:.4f}".format(pred2['pars']['Kga'])],
              ['k_fit[1/s]', k_fit], ['k & k2 baseline[1/s]', k],['k2_fit',k2_fit]]
head = ['parameter', 'value']
table = tabulate(parameters, headers=head, tablefmt="grid")


           


# Create a figure and axis
fig, ax = plt.subplots()

# Hide the axes
ax.axis('off')

# Plot the table
table_ax = ax.table(cellText=parameters, colLabels=head, loc='center', cellLoc='center')

# Adjust font size
table_ax.auto_set_font_size(False)
table_ax.set_fontsize(8)

plt.title('Experiment '+first+'.'+second)

# Save the figure
plt.savefig('Plots/Inputs/Input_parameters_optimized [kl, k, k2], [IG] ['+str(IG_kl)+','+str(IG_k)+','+str(IG_k2)+'].png', bbox_inches='tight')


