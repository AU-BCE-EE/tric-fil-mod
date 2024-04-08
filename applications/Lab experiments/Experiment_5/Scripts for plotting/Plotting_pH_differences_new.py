# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
import importlib
import pandas as pd
import sys
import matplotlib.pyplot as plt

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sys.path.append("../../../..")  # Add the directory containing mod_funcs.py to Python path

from mod_funcs import tfmod

# Create an empty dictionary to store the DataFrames
results_dict = {}
modelresults_dict = {}


# Choose experiment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i in range (1,8):
    first = '5'
    second = str(i)

    #parameters loaded. File name changes depending on experiment no.
    filename = f'lab_parameters_5{second}'
    
    # Dynamically import the module using importlib
    module = importlib.import_module(filename)
    
    # Access the variables from the imported module
    pH1 = module.pH1
    pH2 = module.pH2
    cycle1 = module.cycle1
    cycle2 = module.cycle2
    cycle3 = module.cycle3
    cycle4 = module.cycle4
    length = module.length
    vol = module.vol

    # Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # See notes in tfmod.py for more complete descriptions
    L = 0.51            # Filter length/depth (m) 
    nc = 200          # Number of model cells (layers)
    cg0 = 0          # (g/m3)
    cl0 = 0          # (g/m3)
    henry = (0.1, 2000.)
    temp = 21.       # (degrees C)
    dens_l = 1000    # Liquid density (kg/m3)

    #Setting the liquid and gas velocities and the porosity and water content

    v_g = 53/3600
    v_l = 0.4/3600

    por_l = 0.21
    por_g = 0.77


    #Calculating the cross sectional area, and dividing the volume by it, as required by the model
    area = (0.19/2)**2 * 3.14159265
    v_res = vol * 10**(-6) / area


    k = 0       # Reaction rate (1/s). Small because of inert carrier
    k2=k                 # Reaction could be acid/base that changes the pH



    # realistic pKa
    pKa = 7.


    ## Breakthrough time calculated from volume (V=14.5L), Por_g=0.8 and volumetric velocities 25l/min (v_g=53m/h) 
    # . 37.5L/min (v_g=79m/h) and 50 L/min (v_g=106m/h)

    BT1 = 14.5*por_g/50
    BT2 = 14.5*por_g/25
    BTlabel='Theoretical Breakthrough curve)'
    BT=BT2





    # Inlet concentrations load in from data treatment and definition

    ex1 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')


    #second inlet profile (last data obtained)
    ex5 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.inlet2.csv', sep = ',')


    # The three repetitions
    ex2 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.1.csv', sep = ',')
    ex3 = pd.read_csv('..//Processed_data/experiment_'+first+'.'+second+'.2.csv', sep = ',')


    # selecting and normalising the time, because all data files are started before the H2S
    # was turned on, and not after a specific time, it all depends on the individual experient. 
    #However the cycle at which the H2S was turned on was written down.
    t1norm = ex1['Time in h'] - ex1['Time in h'][cycle1]
    t1 = t1norm [cycle1:cycle1+length]
    C1= ex1['Concentration in g/m^3'][cycle1:cycle1+length]

    t5norm = ex5['Time in h'] - ex5['Time in h'][cycle4]
    t5 = t5norm [cycle4:cycle4+length]
    C5= ex5['Concentration in g/m^3'][cycle4:cycle4+length]

    t2norm = ex2['Time in h'] - ex2['Time in h'][cycle2]
    t2 = t2norm [cycle2:cycle2+length]
    C2= ex2['Concentration in g/m^3'][cycle2:cycle2+length]

    t3norm = ex3['Time in h'] - ex3['Time in h'][cycle3]
    t3 = t3norm [cycle3:cycle3+length]
    C3= ex3['Concentration in g/m^3'][cycle3:cycle3+length]




        



    #Average of the two inlet profiles for use in the model

    c_in=np.mean([C1,C5],axis=0)

    cgin = pd.DataFrame({'time': t5*3600, 
                         'cgin': c_in})

    clin = 0

    c_out = np.mean([C2,C3],axis=0)

    # Times for model output in h
    tt = 0.35
    # Number of time rows
    nt = 200
    times = np.linspace(0, tt, nt) * 3600
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
                  cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, k2 = k2, henry = henry, pKa = pKa, 
                  pH = pH1, temp = temp, dens_l = dens_l, times = times, v_res = v_res, counter = True, recirc = True)
    pred1label= first+'.'+second+'.1 model' #label on plots




    pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
                  cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, k2 = k2, henry = henry, pKa = pKa, 
                  pH = pH2, temp = temp, dens_l = dens_l, times = times, v_res = v_res, counter = True, recirc = True)
    pred2label=first+'.'+second+'.2 model'




    results = pd.DataFrame ({'experimental':c_out , 'experimental time (h)':t3 })

    # Store the DataFrame in the dictionary with a key for each round
    results_dict[f'experimental{first}.{second}'] = results



    modelresults = pd.DataFrame ({'model':(pred1['gas_conc'][nc - 1, :]+pred2['gas_conc'][nc - 1, :])/2 , 'model time(h)':pred1['time'] / 3600 })

    modelresults_dict[f'experimental{first}.{second}'] = modelresults

#import data and make the moving average to avoid noise

window_size = 10

ex1 = results_dict['experimental5.1']
model1 = modelresults_dict['experimental5.1']
ex1['Moving_Average'] = ex1['experimental'].rolling(window=window_size,center=True).mean()

ex2 = results_dict['experimental5.2']
model2 = modelresults_dict['experimental5.2']
ex2['Moving_Average'] = ex2['experimental'].rolling(window=window_size,center=True).mean()

ex3 = results_dict['experimental5.3']
model3 = modelresults_dict['experimental5.3']
ex3['Moving_Average'] = ex3['experimental'].rolling(window=window_size,center=True).mean()

ex4 = results_dict['experimental5.4']
model4 = modelresults_dict['experimental5.4']
ex4['Moving_Average'] = ex4['experimental'].rolling(window=window_size,center=True).mean()

ex5 = results_dict['experimental5.5']
model5 = modelresults_dict['experimental5.5']
ex5['Moving_Average'] = ex5['experimental'].rolling(window=window_size,center=True).mean()

ex6 = results_dict['experimental5.6']
model6 = modelresults_dict['experimental5.6']
ex6['Moving_Average'] = ex6['experimental'].rolling(window=window_size,center=True).mean()

ex7 = results_dict['experimental5.7']
model7 = modelresults_dict['experimental5.7']
ex7['Moving_Average'] = ex7['experimental'].rolling(window=window_size,center=True).mean()





por_g = 0.77
BT = 14.5*por_g/25
BTlabel='Theoretical Breakthrough curve)'


plt.plot(ex7['experimental time (h)'],ex7['Moving_Average'],color='y',label='experimaltal pH=5.99')
plt.plot(ex1['experimental time (h)'],ex1['Moving_Average'],color='c',label='experimaltal pH=7.08')
plt.plot(ex2['experimental time (h)'],ex2['Moving_Average'],color='k',label='experimaltal pH=7.27')
plt.plot(ex5['experimental time (h)'],ex5['Moving_Average'],color='b',label='experimaltal pH=7.55')
plt.plot(ex3['experimental time (h)'],ex3['Moving_Average'],color='r',label='experimaltal pH=7.76')
plt.plot(ex4['experimental time (h)'],ex4['Moving_Average'],color='m',label='experimaltal pH=8.02')
plt.plot(model7['model time(h)'],model7['model'],color='y',linestyle='dashed', label='model pH=5.99')
plt.plot(model1['model time(h)'],model1['model'],color='c',linestyle='dashed', label='model pH=7.08')
plt.plot(model2['model time(h)'],model2['model'],color='k',linestyle='dashed', label='model pH=7.27')
plt.plot(model5['model time(h)'],model5['model'],color='b',linestyle='dashed', label='model pH=7.55')
plt.plot(model3['model time(h)'],model3['model'],color='r',linestyle='dashed', label='model pH=7.76')
plt.plot(model4['model time(h)'],model4['model'],color='m',linestyle='dashed', label='model pH=8.02')
plt.axvline(x=BT/60,linestyle='-',label=BTlabel) #breakthrough curve
plt.axhline(y=0.055,color='g',label='Expected inlet concentration')
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.legend()
plt.xlim(0)
plt.ylim(0)
plt.subplot(111).legend(loc='upper center',bbox_to_anchor=(0.5,-0.2)) #Moves legend out of plot
plt.title('pH comparison')
plt.savefig('..//Plots/pH comparison.png', bbox_inches='tight')

    



