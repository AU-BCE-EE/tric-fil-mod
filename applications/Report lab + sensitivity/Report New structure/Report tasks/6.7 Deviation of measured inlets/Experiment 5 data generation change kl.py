# Comparison of numerical Python model to closed-form solution with instant partitioning (equilibrium everywhere)

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import pandas as pd



# Create an empty dictionary to store the DataFrames
results_dict = {}
modelresults_dict = {}

for a in range (5,7):
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

        cycle1 = int(params['cycle1'].values[0])
        cycle5 = int(params['cycle5'].values[0])
        length = params['length'].values[0]




        # Inlet concentrations load in from data treatment and definition

        ex1 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')


        #second inlet profile (last data obtained)
        ex5 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet2.csv', sep = ',')

        


        # selecting and normalising the time, because all data files are started before the H2S
        # was turned on, and not after a specific time, it all depends on the individual experient. 
        #However the cycle at which the H2S was turned on was written down.
        t1norm = ex1['Time in h'] - ex1['Time in h'][cycle1]
        t1 = t1norm [cycle1:cycle1+length]
        C1= ex1['Concentration in g/m^3'][cycle1:cycle1+length]

        t5norm = ex5['Time in h'] - ex5['Time in h'][cycle5]
        t5 = t5norm [cycle5:cycle5+length]
        C5= ex5['Concentration in g/m^3'][cycle5:cycle5+length]

    
        
        c_in=np.mean([C1,C5],axis=0)
        c_in[(t1*60 < 3.5) | (t1*60 > 4.5)] = 0
        c_in = c_in[c_in !=0]
        dev = (np.mean(c_in) - 0.0596) / 0.0596 * 100
        formatted_dev = f'{dev:.3g}'
        
        print(first+'.'+second+' '+str(formatted_dev)+'% afvigelse')
        
        

        

        
        
        




        



   
    
    







    







    



