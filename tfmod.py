# 1D 2-phase bio/chem air filter model in Python

# Packages
import numpy as np
import math
from scipy.integrate import solve_ivp

# Rates function
# Arg order: time, state variable, then arguments
def rates(t, mc, Q, cgin, vg, vl, vt, k, ka, hc):

    #breakpoint()
    
    # Number of cells (layers) (note integer division)
    nc = mc.shape[0] // 2

    # Separate gas and liquid state variables (g) (g/m2)
    mcg = mc[0:nc]
    mcl = mc[nc:(2 * nc)]

    # Concentrations (g/m3)
    ccg = mcg / vg
    ccl = mcl / vl

    # Derivatives
    # Set up empty arrays
    dmg = dml = g2l = np.zeros(nc)

    # Common term, mass transfer into liquid phase (g/s)
    #g/s 1/s  m3(t) ----g/m3(g)-----
    g2l = ka * vt * (ccg - ccl / hc) 

    # Gas phase derivatives (g/s)
    # cddiff = concentration double difference (g/m3)
    # cvec = array of cell concentrations with inlet air added
    cvec = np.insert(ccg, 0, cgin)
    advec = - Q * np.diff(cvec)
    dmg = advec - g2l

    # Liquid phase derivatives (g/s)
    # NTS think about mass versus conc state variable
    rxn = k * mcl
    dml = g2l - rxn

    # Combine gas and liquid
    dm = np.concatenate([dmg, dml])

    #if t / 3600 > 0.05:
    #    breakpoint()

    return dm

# Model function
def tfmod(L, gas, liq, Q, nc, cg0, cl0, cgin, ka, k, henry, temp, dens, times):

    # L = total longitudinal length/height of reactor/filter (m)
    # gas = gas phase porosity (m3/m3 = m3(g)/m3(t))
    # liq = liquid phase content (m3/m3 = m3(l)/m3(t))
    # Q = gas flow rate (m3/s) (m3/m2-s = m3(g)/m2(t)-s)
    # nc = number of cells (layers)
    # cg0 = initial compound concentration in gas phase (g/m3)
    # cl0 = initial compound concentration in liquid phase (g/m3)
    # cgin = compound concentration in inflow (g/m3)
    # ka = mass transfer coefficient for gas to liquid in gas phase units (1/s = g/s-m3(t) / g/m3(g))
    # k = first-order liquid phase reaction rate constant (1/s = g/s / g)
    # henry = Henry's law constant coefficients as [k_H at 25 C, d(ln(kH)) / d(1/T)] as in NIST web book
    # temp = temperature (degrees C)
    # dens = solution (liquid) density (kg/m3)

    # Retention time (s)
    rt = L * gas / Q

    # Ideal gas constant (L bar / K-mol)
    R = 0.083144 
    
    # Make sure some inputs are numeric to avoid integer math bug
    cg0 = float(cg0)
    cl0 = float(cl0)
    cgin = float(cgin)
    ka = float(ka)
    k = float(k)
    temp = float(temp)
    henry = np.array(henry).astype(float)

    # Temperature and Henry's law constant
    TK = temp + 273.15
    kh = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15)) # mol/kg-bar as liq:gas
    kh = kh * dens / 1000                                  # mol/L-bar
    hc = kh * R * TK                                       # dimensionless, liq:gas, e.g., g/L / g/L
    
    # Create cells
    x = np.linspace(0, L, nc + 1)  # nc + 1 values
    dx = np.diff(x)[0]             # dx, single value, same for all cells
    x = x[1:(nc + 1)] - dx / 2     # Center position in m
    
    # Cell volumes (m3) (m3/m2)
    # Total
    vt = dx * 1
    # Gas
    vg = vt * gas 
    # Liquid
    vl = vt * liq 
    
    # Compound conc (g/m3) and mass (g) (g/m2)
    # cc = concentration, mc = mass [position]
    # g = gas, l = liquid
    ccg = np.full((nc), cg0)
    ccl = np.full((nc), cl0)
    mcg = ccg * vg
    mcl = ccl * vl

    # Initial state variable array
    y0 = np.concatenate([mcg, mcl])
    
    # Solve/integrate
    out = solve_ivp(rates, [0, max(times)], y0 = y0, 
                    t_eval = times, 
                    args = (Q, cgin, vg, vl, vt, k, ka, hc),
                    method = 'LSODA')
    
    # Extract mass of compound [position, time]
    mcgt = out.y[0:nc]
    mclt = out.y[nc:(2 * nc)]

    # Get concentrations vs. time
    # Gas
    ccgt = mcgt / np.transpose(np.tile(vg, (mcgt.shape[1], 1)))
    # Liquid
    cclt = mclt / np.transpose(np.tile(vl, (mclt.shape[1], 1)))
    # Total
    cctt = mclt / np.transpose(np.tile(vt, (mclt.shape[1], 1)))

    mct = np.concatenate([mcgt, mclt])
    
    # Return results as a tuple
    return ccgt, cclt, mcgt, mclt, x, times, rt
