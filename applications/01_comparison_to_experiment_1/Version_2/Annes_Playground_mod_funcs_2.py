# 1D 2-phase bio/chem air filter model in Python

# Packages
import numpy as np
import math
import pandas as pd
import sys
from scipy.integrate import solve_ivp 

# Mass transfer coefficient function
# See tfmod function for units on inputs
def Kga_onda(pH, temp, henry, pKa, pres, ssa, v_g, v_l, por_g, dens_l):

    # Hard-wired constants
    g = 9.81        # m / sec^2
    Dg = 1.16E-5    # gas diffusion coefficient in m2 / sec; compound specific
    Dliq = 1.89E-9  # liquid diffusion coefficient
    sigm_c = 0.75   # critical surface tension
    sigm_l = 0.0073 # surface tension
    R = 0.083144    # Gas constant (L bar / K-mol)
    mw_g = 28.97    # Air (gas mix) molecular weight (molar mass) (g/mol)

    dp = 6 * (1 - por_g) / ssa  # characteristic packing length

    if dp < 15:
        dp_emp = 2.0
    else:
        dp_emp = 5.23

    TK = temp + 273.15

    dens_g = pres * mw_g / (R * TK) # g/L = kg/m3
    visc_g = 9.1E-8 * TK - 1.16E-5   # empirical relation for gas viscosity vs TK
    visc_l = -2.55E-5 * TK + 8.51E-3
    Re = dens_l * v_l / (ssa * visc_l)
    Fr = v_l * v_l * ssa / g
    We = v_l * v_l * dens_l / (sigm_l * ssa)
     
    ae = ssa * (1.0-2.71828**(-1.45 * (sigm_c / sigm_l)**0.75 *
                            Re**0.1 * Fr**-0.05 * We**0.2))
   
    kg = dp_emp * (v_g * dens_g / (ssa * visc_g))**0.7 \
        * (visc_g / (dens_g * Dg))**(1 / 3) * (ssa * dp)**-2 * ssa * Dg
   
    kl = 0.0051 * (v_l * dens_l / (ae * visc_l))**(2/3) * (visc_l / (dens_l * Dliq))**(-0.5) * (ssa * dp)**0.4 * (dens_l / (visc_l * g))**(-1/3)
   
    kh = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15)) # mol/kg-bar as liq:gas
    kh = kh * dens_l / 1000                                # mol/L-bar
    Kaw = 1 / (kh * R * TK)                                # Neutral air-water distribution
    # alpha 0 (fraction as uncharged species)
    alpha0 = 1 / (1 + 10**(pH - pKa))
    Daw = alpha0 * Kaw
   
    Rtot = 1 / (kg * ae) + Daw / (kl * ae)
    Kga = 1 / Rtot

    return Kga

# Rates function
# Arg order: time, state variable, then arguments
def rates(t, mc, v_g, v_l, cgin, clin, vol_gas, vol_liq, vol_tot, k, Kga, Daw, counter = True, recirc = False):

  # If time-variable concentrations coming in are given, get interpolated values
  if type(clin) is pd.core.frame.DataFrame:
    clin = np.interp(t, clin.iloc[:, 0], clin.iloc[:, 1])
  elif type(clin) is np.ndarray:
    clin = np.interp(t, clin[0], clin[1])
  elif not (type(clin) is int or type(clin) is float):
    sys.exit('Error: clin input must be float, integer, numpy array, or a pandas data frame, but is none of these')

  if type(cgin) is pd.core.frame.DataFrame:
    cgin = np.interp(t, cgin.iloc[:, 0], cgin.iloc[:, 1])
  elif type(cgin) is np.ndarray:
    cgin = np.interp(t, cgin[0], cgin[1])
  elif not (type(cgin) is int or type(cgin) is float):
    sys.exit('Error: cgin input must be float, integer, numpy array, or a pandas data frame, but is none of these')

  #breakpoint()
  
  # Number of cells (layers) (note integer division)
  nc = mc.shape[0] // 2

  # Separate gas and liquid state variables (g) (g/m2)
  mcg = mc[0:nc]
  mcl = mc[nc:(2 * nc)]

  # Concentrations (g/m3)
  ccg = mcg / vol_gas
  ccl = mcl / vol_liq

  # Get outlet liquid phase concentration to use at inlet if recirc = True
  if recirc:
      if counter:
          oi = 0
      else:
          oi = nc - 1
      clin = ccl[oi]

  # Derivatives
  # Set up empty arrays
  dmg = dml = g2l = np.zeros(nc)

  # Common term, mass transfer into liquid phase (g/s)
  #g/s  1/s    m3(t)     ----g/m3(g)-----
  g2l = Kga * vol_tot * (ccg - ccl * Daw) 

  # Gas phase derivatives (g/s)
  # No reaction in gas phase
  # cddiff = concentration double difference (g/m3)
  # cvec = array of cell concentrations with inlet air added
  # rxn = 0 for as phase
  cvec = np.insert(ccg, 0, cgin)
  advec = - v_g * np.diff(cvec)
  dmg = advec - g2l

  # Liquid phase derivatives (g/s)
  # Includes transport and reaction
  rxn = k * mcl
  if counter==True:
     cvec = np.insert(ccl, 0, clin)
  else:
     v_l = - v_l
     cvec = np.insert(ccl, nc, clin)

  advec = - v_l * np.diff(cvec)
  dml = advec + g2l - rxn

  # Combine gas and liquid
  dm = np.concatenate([dmg, dml])

  #if t / 3600 > 0.05:
  #    breakpoint()

  return dm

# Model function
def tfmod(L, por_g, por_l, v_g, v_l, nc, cg0, cl0, cgin, clin, Kga, k, henry, pKa, pH, temp, dens_l, times, pres = 1., ssa = 1100, counter = False, recirc = False):

   ## Note that units are defined per 1 m2 filter cross-sectional (total) area 
   ## Below, where 2 sets of units are given this applies to the first case
   ## For the second one, the cross-sectional area is used to normalize the unit
   ## The two are mathematically equivalent
   # L = total longitudinal length/height of reactor/filter (m)
   # por_g = gas phase porosity (m3/m3 = m3(g)/m3(t) where g = gas and t = total)
   # por_l = liquid phase content (m3/m3 = m3(l)/m3(t) where l = liquid)
   # v_g = gas flow rate (m3/s) (m3/m2-s = m3(g)/m2(t)-s = superficial velocity in m/s)
   # v_l = liquid flow rate (m3/s) (m3/m2-s = m3(g)/m2(t)-s = superficial velocity in m/s)
   # nc = number of cells (layers)
   # cg0 = initial compound concentration in gas phase (g/m3 = g(compound)/m3(g))
   # cl0 = initial compound concentration in liquid phase (g/m3)
   # cgin = compound concentration in gas inflow (g/m3)
   # clin = compound concentration in liquid inflow (g/m3) (ignored if recirc = True)
   # Kga = mass transfer coefficient for gas to liquid in gas phase units (1/s = g/s-m3(t) / g/m3(g))
   # k = first-order liquid phase reaction rate constant (1/s = g/s / g)
   # henry = Henry's law constant coefficients as [k_H at 25 C, d(ln(kH)) / d(1/T)] as in NIST Chemistry Web Book
   # temp = temperature (degrees C)
   # dens_l = solution (liquid) density (kg/m3)
   # pres = total pressure (bar?)
   # ssa = particle specific surface area (m2 surface / m3 bulk volume)
   # count = Boolean for countercurrent flow

   # Save input arguments for echoing in output
   args_in = locals()

   # Constants
   # Ideal gas constant (L bar / K-mol)
   R = 0.083144 

   # Caclulate Kga if requested
   if type(Kga) is str and Kga.lower() == 'onda':
      Kga = Kga_onda(pH = pH, temp = temp, henry = henry, pKa = pKa, pres = pres, ssa = ssa, v_g = v_g, v_l = v_l, por_g = por_g, dens_l = dens_l)

   # Retention time (s)
   rt_gas = L * por_g / v_g
   rt_liq = L * por_l / v_l
   
   # Make sure some inputs are numeric to avoid integer math bug
   cg0 = float(cg0)
   cl0 = float(cl0)
   Kga = float(Kga)
   k = float(k)
   pKa = float(pKa)
   pH = float(pH)
   temp = float(temp)
   dens_l = float(dens_l)
   times = np.array(times).astype(float)
   henry = np.array(henry).astype(float)

   # Temperature and Henry's law constant
   TK = temp + 273.15
   kh = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15)) # mol/kg-bar as liq:gas
   kh = kh * dens_l / 1000                                # mol/L-bar
   Kaw = 1 / (kh * R * TK)                                # dimensionless, gas:liq, e.g., g/L / g/L or g/m3 per g/m3
   
   # Create cells
   x = np.linspace(0, L, nc + 1)  # nc + 1 values
   dx = np.diff(x)[0]             # dx, single value, same for all cells
   x = x[1:(nc + 1)] - dx / 2     # Center position in m
   
   # Cell volumes (m3) (m3(g/l/t)/m2(t))
   # Total
   vol_tot = dx * 1
   # Gas
   vol_gas = vol_tot * por_g 
   # Liquid
   vol_liq = vol_tot * por_l 
   
   # Compound conc (g/m3) and mass (g) (g/m2)
   # cc = concentration, mc = mass [position]
   # g = gas, l = liquid
   ccg = np.full((nc), cg0)
   ccl = np.full((nc), cl0)
   mcg = ccg * vol_gas
   mcl = ccl * vol_liq

   # Initial state variable array
   y0 = np.concatenate([mcg, mcl])

   # Ionization fraction
   alpha0 = 1 / (1 + 10**(pH - pKa))
   Daw = alpha0 * Kaw
   
   # Solve/integrate
   out = solve_ivp(rates, [0, max(times)], y0 = y0, 
                   t_eval = times, 
                   args = (v_g, v_l, cgin, clin, vol_gas, vol_liq, vol_tot, k, Kga, Daw, counter),
                   method = 'LSODA')
   # Set max_step=0.1 as last argument in the solver to get smooth curves that take longer time to simulate
   
   # Extract mass of compound [position, time]
   mcgt = out.y[0:nc]
   mclt = out.y[nc:(2 * nc)]

   # Get concentrations vs. time
   # Gas
   ccgt = mcgt / np.transpose(np.tile(vol_gas, (mcgt.shape[1], 1)))
   # Liquid
   cclt = mclt / np.transpose(np.tile(vol_liq, (mclt.shape[1], 1)))
   # Total
   cctt = mclt / np.transpose(np.tile(vol_tot, (mclt.shape[1], 1)))

   mct = np.concatenate([mcgt, mclt])
   
   # Return results as a dictionary
   return {'gas_conc': ccgt, 'liq_conc': cclt, 'gas_mass': mcgt, 'liq_mass': mclt, 
           'cell_pos': x, 'time': times, 'Kaw':Kaw, 
           'inputs': args_in,
           'pars': {'gas_rt': rt_gas, 'liq_rt': rt_liq, 'Kga': Kga, 'Kaw': Kaw, 'alpha0': alpha0, 'Daw': Daw, 'y0':y0}}

