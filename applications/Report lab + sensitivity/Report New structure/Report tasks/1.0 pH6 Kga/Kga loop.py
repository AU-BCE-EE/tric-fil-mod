# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit


# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sys.path.append("../../../../..")  # Add the directory containing mod_funcs.py to Python path

from mod_funcs import tfmod

#select and define experimental data ______________________________________________
#selecting the experiment number and loading parameters from the master spreadsheet

first = '5'
second = '7'

#parameters loaded. File name changes depending on experiment no.
file_path = '../../../Master_spreadsheet.xlsx'
df = pd.read_excel(file_path, sheet_name = 'Data')

params = df[df['key'] == float(first)+0.1*float(second)]

pH1 = params['pH1'].values[0]
pH2 = params['pH2'].values[0]
pH3 = params['pH3'].values[0]
cycle1 = int(params['cycle1'].values[0])
cycle2 = int(params['cycle2'].values[0])
cycle3 = int(params['cycle3'].values[0])
cycle4 = params['cycle4'].values[0]
if not pd.isna(cycle4):
    cycle4=int(cycle4)
else:
    cycle4 = 'NaN'
cycle5 = int(params['cycle5'].values[0])
length = params['length'].values[0]
vol = params['volume'].values[0]
no = params['no'].values[0]


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



# Inlet concentrations load in from data treatment and definition

ex1 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')


#second inlet profile (last data obtained)
ex5 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet2.csv', sep = ',')


# The three repetitions of outlet data
ex2 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.1.csv', sep = ',')
ex3 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.2.csv', sep = ',')
if not cycle4 == 'NaN':
    ex4 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.3.csv', sep = ',')
    
    t4norm = ex4['Time in h'] - ex4['Time in h'][cycle4]
    t4 = t4norm [cycle4:cycle4+length+500]
    C4= ex4['Concentration in g/m^3'][cycle4:cycle4+length+500]



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
t2 = t2norm [cycle2:cycle2+length+500]
C2= ex2['Concentration in g/m^3'][cycle2:cycle2+length+500]

t3norm = ex3['Time in h'] - ex3['Time in h'][cycle3]
t3 = t3norm [cycle3:cycle3+length+500]
C3= ex3['Concentration in g/m^3'][cycle3:cycle3+length+500]

#Making an average of the outlet measurements for Output data and the optimization
if not cycle4 == 'NaN':
     minlen = min(len(C2),len(C3),len(C4))
     C2_cut = ex2['Concentration in g/m^3'].rolling(window=5,center = True).mean()[cycle2:cycle2+minlen]
     C3_cut = ex3['Concentration in g/m^3'].rolling(window=5,center = True).mean()[cycle3:cycle3+minlen]
     C4_cut = ex4['Concentration in g/m^3'].rolling(window=5,center = True).mean()[cycle4:cycle4+minlen]
     c_out = np.mean([C2_cut,C3_cut,C4_cut],axis=0)
else:
     minlen = min(len(C2),len(C3))
     C2_cut = ex2['Concentration in g/m^3'].rolling(window=5,center = True).mean()[cycle2:cycle2+minlen]
     C3_cut = ex3['Concentration in g/m^3'].rolling(window=5,center = True).mean()[cycle3:cycle3+minlen]
     c_out = np.mean([C2_cut,C3_cut],axis=0)
  



## Breakthrough time calculated from volume (V=14.5L), Por_g and flow rate (L/min, from calibration of mass flow controllers)

BT1 = 14.5*por_g/54.2766
BT2 = 14.5*por_g/27.2991
BTlabel='Theoretical Breakthrough curve)'

if no == 1 or no == 4:
    BT = BT2
elif no ==2 or no == 3:
    BT= BT1
else: 
    print ('error in definition of "no". Has to be a string containing the numbers 1,2,3 or 4')


#Defining inputs for the model______________________________________________________________________
#Time and inlet concentrations
c_in=np.mean([C1,C5],axis=0)

cgin = pd.DataFrame({'time': t5*3600, 
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
Kga_list = np.linspace(0.0025, 0.25,5) 


preds = []
pred_labels = []

for i, Kga in enumerate(Kga_list):
    pred = tfmod(L=L, por_g=por_g, por_l=por_l, v_g=v_g, v_l=v_l, nc=nc, cg0=cg0, cl0=cl0,
                 cgin=cgin, clin=clin, Kga=Kga, k=k, k2=k, henry=henry, pKa=pKa,
                 pH=pH1, temp=temp, dens_l=dens_l, times=times, v_res=v_res, 
                 recirc=True, counter=True)
    
    label = f"{first}.{second}.Kga {Kga} model"
    
    preds.append(pred)
    pred_labels.append(label)
    
    
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, k2 = k, henry = henry, pKa = pKa, 
              pH = pH1, temp = temp, dens_l = dens_l, times = times, v_res = v_res, recirc = True, counter = True)
pred1label= first+'.'+second+'. onda model' #label on



#Plotting__________________________________________________________________________________________________
plt.clf()
plt.plot(t2,C2,label=first+'.'+second+'.1 experimental data')
plt.plot(t3,C3,label=first+'.'+second+'.2 experimental data')
if not cycle4 == 'NaN':
    plt.plot(t4,C4,label=first+'.'+second+'.3 experimental data')
plt.plot(t1,C1,label='inlet 1')
plt.plot(t5,C5,label='inlet 2')
for pred, label in zip(preds, pred_labels):
    plt.plot(pred['time'] / 3600, pred['gas_conc'][nc - 1, :], label=label)
plt.plot(pred1['time'] / 3600, pred1['gas_conc'][nc - 1, :],label=pred1label)
plt.axvline(x=BT/60,linestyle='-',label=BTlabel) #breakthrough curve
plt.axhline(y=0.0596,color='g',label='Expected inlet concentration') # Average of two inlet concentrations of 43.2 and 42.3. These were for low and high gas flow rate respectively
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0,0.35)
plt.ylim(0,0.07)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('Experiment '+first+'.'+second)
plt.savefig('Plots/Experiment '+first+'.'+second+'Kga loop.png', bbox_inches='tight')
plt.close()

#table with input paramters
#import module
from tabulate import tabulate
parameters = [['v_g', pred1['inputs']['v_g']], ['v_l', pred1['inputs']['v_l']], ['pH1', pH1], ['pH2', pH2], ['pH3', pH3], ['k', k], ['countercurrent', pred1['inputs']['counter']],
              ['recirculation', pred1['inputs']['recirc']], ['water content', por_l], ['porosity', por_g], ['temperature', temp], ['v_res', pred1['inputs']['v_res']]]
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
table_ax.set_fontsize(10)

plt.title('Experiment '+first+'.'+second)

# Save the figure
plt.savefig('Plots/Inputs/Input_parameters_'+first+'.'+second+'baseline parameters.png', bbox_inches='tight')


