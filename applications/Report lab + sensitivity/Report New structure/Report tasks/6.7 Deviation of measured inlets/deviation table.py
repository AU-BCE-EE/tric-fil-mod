# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:21:45 2024

@author: Mortensen
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Create an empty dictionary to store the results
results = []

for a in range(5, 7):
    first = str(a)
    if first == '5':
        no_exp = 8
    else:
        no_exp = 5

    for i in range(1, no_exp):
        second = str(i)

        # parameters loaded. File name changes depending on experiment no.
        file_path = '../../../Master_spreadsheet.xlsx'
        df = pd.read_excel(file_path, sheet_name='Data')

        params = df[df['key'] == float(first) + 0.1 * float(second)]

        cycle1 = int(params['cycle1'].values[0])
        cycle5 = int(params['cycle5'].values[0])
        length = params['length'].values[0]

        # Inlet concentrations load in from data treatment and definition
        ex1 = pd.read_csv('../..//Processed_data/experiment_' + first + '.' + second + '.inlet1.csv', sep=',')
        ex5 = pd.read_csv('../..//Processed_data/experiment_' + first + '.' + second + '.inlet2.csv', sep=',')

        # Selecting and normalizing the time
        t1norm = ex1['Time in h'] - ex1['Time in h'][cycle1]
        t1 = t1norm[cycle1:cycle1 + length]
        C1 = ex1['Concentration in g/m^3'][cycle1:cycle1 + length]

        t5norm = ex5['Time in h'] - ex5['Time in h'][cycle5]
        t5 = t5norm[cycle5:cycle5 + length]
        C5 = ex5['Concentration in g/m^3'][cycle5:cycle5 + length]

        c_in = np.mean([C1, C5], axis=0)
        c_in[(t1 * 60 < 3.5) | (t1 * 60 > 4.5)] = 0
        c_in = c_in[c_in != 0]
        dev = (np.mean(c_in) - 0.0596) / 0.0596 * 100
        formatted_dev = f'{dev:.5g}'

        # Append the result to the list
        results.append([first + '.' + second, formatted_dev])

# Create a DataFrame from the results
results_df = pd.DataFrame(results, columns=['Experiment', 'Deviation (%)'])

# Plot the table
fig, ax = plt.subplots(figsize=(8, 6))  # Adjust the size as needed
ax.axis('tight')
ax.axis('off')
table = ax.table(cellText=results_df.values, colLabels=results_df.columns, cellLoc='center', loc='center')
table.auto_set_font_size(False)
table.set_fontsize(12)
table.scale(1.2, 1.2)  # Adjust scaling for better fit

# Save the table as a PNG file
#plt.savefig('deviation_table.png', bbox_inches='tight', dpi=300)

# Display the table
plt.show()
