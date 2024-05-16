# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import importlib

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sys.path.append("../../../../..")  # Add the directory containing mod_funcs.py to Python path
from mod_funcs import tfmod 



# Create an empty dictionary to store the DataFrames
results_dict = {}
modelresults_dict = {}

for a in range (6,7):
    first = str(a)
    if first == '5':
        no_exp = 8
    else:
        no_exp = 5
    
    for i in range (1,no_exp):
        second = str(i)

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


        k = 0       # Reaction rate (1/s). Small because of inert carrier
                         # Reaction could be acid/base that changes the pH
        kl = 'onda' 
        
        # realistic pKa
        pKa = 7.


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



        # Inlet concentrations load in from data treatment and definition

        ex1 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')


        #second inlet profile (last data obtained)
        ex5 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet2.csv', sep = ',')


        # The three repetitions
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

        # Inlet concentration for steady state

        
        c_in = [0.0596 if i < 120 else 0 for i in range(220)]



        clin = 0

        # Times for model output in h
        tt = 1
        # Number of time rows
        nt = 220
        times = np.linspace(0, tt, nt) * 3600

        cgin = pd.DataFrame({'time': times, 
                             'cgin': c_in})

        clin = 0

      
        
        # Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # See notes in tfmod.py for more complete descriptions
        L = 0.51            # Filter length/depth (m) 
        nc = 200          # Number of model cells (layers)
        cg0 = 0          # (g/m3)
        cl0 = 0          # (g/m3)
        henry = (0.1, 2000.)
        temp = 21.       # (degrees C)
        dens_l = 1000    # Liquid density (kg/m3)
        
        if first == '3':
            recirc = False
        else: 
            recirc = True
        
        pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
                      cl0 = cl0, cgin = cgin, clin = clin, Kga = 'individual', k = k, henry = henry, pKa = pKa, 
                      pH = pH1, temp = temp, dens_l = dens_l, times = times, v_res = v_res, k2=k, kl = kl,recirc = recirc, counter = True)
        pred1label= first+'.'+second+'.1 model' #label on
        
        pred2 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
                      cl0 = cl0, cgin = cgin, clin = clin, Kga = 'individual', k = k, henry = henry, pKa = pKa, 
                      pH = pH2, temp = temp, dens_l = dens_l, times = times, v_res = v_res, k2 = k, kl = kl, recirc = recirc, counter = True)
        pred2label= first+'.'+second+'.2 model' #label on
        
        if not cycle4 == 'NaN':
           pred3 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
                         cl0 = cl0, cgin = cgin, clin = clin, Kga = 'individual', k = k, henry = henry, pKa = pKa, 
                         pH = pH3, temp = temp, dens_l = dens_l, times = times, v_res = v_res, k2 = k,kl = kl, recirc = recirc, counter = True)
           pred3label= first+'.'+second+'.3 model' #label on 
        
        

        
        #Taking the average of the outlet measurements and exporting to csv______________________________________
      
        
        #export data as CSV files
        if not cycle4 == 'NaN':
            minlen = min(len(C2),len(C3),len(C4))
            C2 = ex2['Concentration in g/m^3'][cycle2:cycle2+minlen]
            C3 = ex3['Concentration in g/m^3'][cycle3:cycle3+minlen]
            C4 = ex4['Concentration in g/m^3'][cycle4:cycle4+minlen]
            c_out = np.mean([C2,C3,C4],axis=0)
        else:
            minlen = min(len(C2),len(C3))
            C2 = ex2['Concentration in g/m^3'][cycle2:cycle2+minlen]
            C3 = ex3['Concentration in g/m^3'][cycle3:cycle3+minlen]
            c_out = np.mean([C2,C3],axis=0)
        
        results = pd.DataFrame ({'experimental':c_out[10:1010] , 'experimental time (h)':t3norm [cycle3:cycle3+1000] })

        results.to_csv('Output data/experimental_delay_'+first+'.'+second+'.csv', header = True, index=False)
        
        if not cycle4 == 'NaN':
            mod_out = (pred1['gas_conc'][nc - 1, :]+pred2['gas_conc'][nc - 1, :]+pred3['gas_conc'][nc - 1, :])/3
        else: 
            mod_out = (pred1['gas_conc'][nc - 1, :]+pred2['gas_conc'][nc - 1, :])/2



        results = pd.DataFrame ({'model':mod_out[0:1000] , 'model time(h)':pred1['time'][0:1000] / 3600 })

        results.to_csv('Output data/model_delay_'+first+'.'+second+'.csv', header = True, index=False)
        
        #making the residuals
        
        
        pred = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
                      cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
                      pH = (pH1+pH2)/2, temp = temp, dens_l = dens_l, times = t3norm[cycle3:cycle3+minlen]*3600, v_res = v_res, k2=k, recirc = recirc, counter = True)
        r = np.transpose(c_out[10:1010] - pred['gas_conc'][nc - 1, :][0:1000])
        
        results = pd.DataFrame({'residuals[g/m3]':r, 'time_residuals[h]':t3norm[cycle3:cycle3+1000]})
        
        results.to_csv('Output data/residuals_delay_'+first+'.'+second+'.csv', header = True, index=False)
        

        
        
        




        



   
    
    







    







    



