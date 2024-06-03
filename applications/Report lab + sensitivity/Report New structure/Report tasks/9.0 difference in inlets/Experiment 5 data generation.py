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

for a in range (5,6):
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

        



        # Inlet concentrations load in from data treatment and definition

        ex1 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet1.csv', sep = ',')
        t1norm = ex1['Time in h'] - ex1['Time in h'][cycle1]
        t1 = t1norm [cycle1:cycle1+length]
        

        #second inlet profile (last data obtained)
       
        ex5 = pd.read_csv('../..//Processed_data/experiment_'+first+'.'+second+'.inlet2.csv', sep = ',')
        t5norm = ex5['Time in h'] - ex5['Time in h'][cycle5]
        t5 = t5norm [cycle5:cycle5+length]
        
        C1= ex1['Concentration in g/m^3'][cycle1:cycle1+length]
        C5= ex5['Concentration in g/m^3'][cycle5:cycle5+length]
        inlet = np.mean([C1,C5],axis=0)
        
        results = pd.DataFrame ({'inlet':inlet , 'time(h)':t5 })

        results.to_csv('Output data/inlet_'+first+'.'+second+'.csv', header = True, index=False)
        
        
        
        




        



   
    
    







    







    



