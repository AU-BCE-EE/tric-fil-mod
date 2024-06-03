# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:27:19 2024

@author: Mortensen
"""
import math
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


#First attempt: Vant Hoffs equation

temp_list = [10, 15, 18, 21, 25] 

# Making the models

for i, temp in enumerate(temp_list):
    # realistic pKa, temperture dependent estimated by Vant Hoffs equation
    pKa298 = 7.

    Ka298 = 10**(-pKa298)
    DeltaHr = 1470000 #J/mol, from nist: 1470kJ/mol
    R = 8.314463 # gas constant in J/mol/K
    temp_K = temp + 273 
    KaT = Ka298 * math.exp((-DeltaHr / R) * (1/temp_K - 1/298))

    pKa = -math.log10(KaT)
    print('For temperature of '+str(temp)+' degrees C, Vant Hoff pKa = '+str(pKa))

#Second attempt linerar interpolation from US-EPA table values from 1978 (See issue 28 for link)


temp_list = [10, 15, 18, 21, 25] 

# Making the linear regression
#from the table: x= temp in deg. C, y = pKa
x = [10, 15, 20, 25, 30, 35, 40]
y= [7.20, 7.13, 7.06, 6.99, 6.92, 6.85, 6.78]
#initial guess for [slope,intersect]
params = [-0.14,7.48]
 
 
def linreg(h,a,b):#a and b is calibration parameters, x is temperature, y is pKa
      n=len(h)
      c=np.zeros(n)
      for i in range(n):
          c[i]=a*h[i]+b
      return c

  #Making the regression
copt,ccov=curve_fit(linreg,x,y,p0=params)

a=float(copt[0])
b=float(copt[1])    

yfit=linreg(x,copt[0],copt[1])
plt.plot(x,y,'o')
plt.plot(x,yfit,label=('y= %1.5f'%copt[0]+'*x + %1.5f'%copt[1]))
plt.xlabel('Temperature \u00b0 C')
plt.ylabel('pKa')
plt.legend(loc='upper right')
plt.show()

pKa_US_EPA = []
pKa_Blanes = []

for i, temp in enumerate(temp_list):
    
   
  pKa = a * temp + b
    
  print('For temperature of '+str(temp)+' degrees C, U.S. EPA pKa = '+str(pKa))
  
  pKa_US_EPA.append(pKa)
  


#Third attempt, from Blanes-Vidal et al (Again, see issue #28 for link and references)

for i, temp in enumerate(temp_list):
    
    T_K = temp + 273.15
    
    pKa = -(-3448.7 / T_K + 47.479 - 7.5227 * np.log(T_K) )
    
    
    print('For temperature of '+str(temp)+' degrees C, Blanes et al pKa = '+str(pKa))
    
    pKa_Blanes.append(pKa)
    

from tabulate import tabulate

# Create a table
parameters = list(zip(temp_list, pKa_US_EPA, pKa_Blanes))
headers = ['Temperature', 'US EPA', 'Blanes et al']
table = tabulate(parameters, headers=headers, tablefmt='grid')

# Plotting the table
fig, ax = plt.subplots(figsize=(8, 4))  # Adjust the figure size as needed

# Hide the axes
ax.axis('off')

# Plot the table
table_ax = ax.table(cellText=parameters, colLabels=headers, loc='center', cellLoc='center')

# Adjust font size
table_ax.auto_set_font_size(False)
table_ax.set_fontsize(10)

plt.title('pKa')
plt.show()


#Henry constants
#First one from nist

H_nist = []
H_Blanes = []

for i, temp in enumerate(temp_list):
    
    TK = temp + 273.15
    henry = (0.1, 2000.)
    R = 0.083144    # Gas constant (L bar / K-mol)
    dens_l = 1000

    kh = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15)) # mol/kg-bar as liq:gas
    kh = kh * dens_l / 1000                                # mol/L-bar
    Kaw = 1 / (kh * R * TK) 
    
    H_nist.append(Kaw)

for i, temp in enumerate(temp_list):
    TK = temp + 273.15
    H = 10**(5.703 - 884.94 / TK)
    H_Blanes.append(H)
    


# Create a table
parameters = list(zip(temp_list, H_nist, H_Blanes))
headers = ['Temperature', 'Nist', 'Blanes et al']
table = tabulate(parameters, headers=headers, tablefmt='grid')

# Plotting the table
fig, ax = plt.subplots(figsize=(8, 4))  # Adjust the figure size as needed

# Hide the axes
ax.axis('off')

# Plot the table
table_ax = ax.table(cellText=parameters, colLabels=headers, loc='center', cellLoc='center')

# Adjust font size
table_ax.auto_set_font_size(False)
table_ax.set_fontsize(10)

plt.title('Henrys constant')
plt.show()
    
    
    


