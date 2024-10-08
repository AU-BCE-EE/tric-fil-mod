# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 09:14:38 2024

@author: au611147
"""

import sys
import math
#The following script describes mass transfer from Kim and Deshusses, in their 2008 article " Determination of mass transfer coefficients for packing materials
# used in biofilters and biotrickling filters for air pollution control - 2". The calling function in the model is Kga = 'KD' when calling mod_funcs in a demo or application


#Uses different types of packing material. 
#LR = lavas rock, PUF = polyurethane foam
#PR = pall rings. PCB = porous ceramic beads
#PCR = porous ceramic rings



def Kga_KD(typ,v_l,v_g,henry,temp,dens_l,pKa,pH):
    
    R = 0.083144    # Gas constant (L bar / K-mol)
    TK = temp + 273.15
    kh = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15)) # mol/kg-bar as liq:gas
    kh = kh * dens_l / 1000                                # mol/L-bar
    Kaw = 1 / (kh * R * TK)                                # Neutral air-water distribution
    # alpha 0 (fraction as uncharged species)
    alpha0 = 1 / (1 + 10**(pH - pKa))
    Daw = alpha0 * Kaw
    
    #Selecting parameters i2 and logC2 from liquid phase velocity and filter material
    if typ == 'TBD': #TBD = to be determined
        print('typ = '+typ)
        sys.exit('Typ is not defined correctly. Needs to be LR (lava rock) PUF, PR (Pall ring), PCB (porous ceramic beads) or PCR (Porous ceramic ring). Error originaring in the Kga determination using Kim and Deshusses')
   
    if typ =='LR':
        if v_l * 3600 > 8:
            logC2 = 3.06
            i2 = 0.14
        elif v_l*3600 < 1:
            logC2 = 2.83
            i2 = 0.02
        else:
            logC2 = 2.81
            i2 = 0.24
            
    elif typ =='PUF':    
        if v_l * 3600 > 8:
            logC2 = 2.75
            i2 = 0.07
        elif v_l*3600 < 3:
            logC2 = 1.19
            i2 = 0.31
        else:
            logC2 = 1.93
            i2 = 0.29
            

    elif typ =='PR':    
        if v_l * 3600 > 8:
            logC2 = 1.86
            i2 = 0.37
        elif v_l*3600 < 3:
            logC2 = 0.83
            i2 = 0.55
        else:
            logC2 = 2.24
            i2 = 0.28

    elif typ =='PCB':    
        if v_l * 3600 > 8:
            logC2 = 2.78
            i2 = 0.35
        elif v_l*3600 < 3:
            logC2 = 2.03
            i2 = 0.46
        else:
            logC2 = 2.72
            i2 = 0.37

    elif typ =='PCR':    
        if v_l * 3600 > 8:
            logC2 = 3.08
            i2 = 0.09
        elif v_l*3600 < 3:
            logC2 = 1.47
            i2 = 0.58
        else:
            logC2 = 2.88
            i2 = 0.17
            
            
   
    #Calculating the gas side mass transfer coefficient in s^-1
    Kgaw = (10**(logC2+i2*math.log10(v_g*3600)))/3600
   #print('Kgaw = '+str(Kgaw*3600)+',filtermaterial = '+str(typ))
   
    #Selecting the parameters i3 and logC3 from filter material and gas velocity
    
    if typ =='LR':
        if v_g * 3600 < 500:
            logC3 = 1.23
            i3 = 0.84
        elif v_g * 3600 >2100:
            logC3 = 1.36
            i3 = 0.84
        else:
            logC3 = 1.29
            i3 = 0.86
            
    elif typ =='PUF':
        if v_g * 3600 < 500:
            logC3 =0.52
            i3 = 0.82
        elif v_g * 3600 >3100:
            logC3 = 0.56
            i3 = 0.87
        else:
            logC3 = 0.53
            i3 = 0.90
    
    elif typ =='PR':
        if v_g * 3600 < 500:
            logC3 =0.63
            i3 = 0.84
        elif v_g * 3600 >3100:
            logC3 = 0.76
            i3 = 0.82
        else:
            logC3 = 0.67
            i3 = 0.83
    
    elif typ =='PCB':
        if v_g * 3600 < 250:
            logC3 =1.44
            i3 = 0.94
        elif v_g * 3600 >550:
            logC3 = 1.41
            i3 = 0.95
        else:
            logC3 = 1.45
            i3 = 0.94
            
    elif typ =='PCR':
        if v_g * 3600 < 500:
            logC3 = 0.99
            i3 = 0.59
        elif v_g * 3600 >2100:
            logC3 = 1.17
            i3 = 1.03
        else:
            logC3 = 1.34
            i3 = 0.82
            
    #Calculating the liquid side mass transfer coefficient in s^-1
    Klaw = (10**(logC3+i3*math.log10(v_l*3600)))/3600
   #print('KlAw = '+str(Klaw*3600)+',filtermaterial = '+str(typ))
    #Total mass transfer coefficient calculation and export
    Rtot = 1/Kgaw + Daw / Klaw
    Kga = 1/Rtot
   #print ('KD Kga = '+str(Kga)+',filtermaterial = '+str(typ))
    return Kga
            



    