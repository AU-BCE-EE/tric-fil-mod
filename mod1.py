# Simple 1D diffusion model to start to think about reaction-transport models in Python

# Packages
import numpy as np
import math
from scipy.integrate import solve_ivp

# Rates function
# Arg order: time, state variable, then arguments
def rates(t, mc, Q, cgin, vg, vl, k, ka, hc):

    #breakpoint()
    
    # Number of cells (layers) (note integer division)
    nc = mc.shape[0] // 2

    # Separate gas and liquid state variables
    # Note that moronic Python language starts at 0 and ends at position beyond end of array
    mcg = mc[0:nc]
    mcl = mc[nc:(2 * nc)]

    # Concentrations (g/L)
    ccg = mcg / vg
    ccl = mcl / vl

    # Derivatives
    # Set up empty arrays
    dmg = dml = g2l = np.zeros(nc)

    # Common term, mt into liquid phase (g/s)
    g2l = ka * (ccg - ccl / hc) 

    # Gas phase derivatives (g/s)
    # cdiff = concentration difference (g/L)
    cdiff = np.diff(np.insert(ccg, 0, cgin))
    #     gas flow in    mt
    dmg = Q * - cdiff - g2l

    # Liquid phase derivatives (g/s)
    #     mt  -reaction-
    dml = g2l - k * mcl

    # Combine gas and liquid
    dm = np.concatenate([dmg, dml])

    #breakpoint()
    return dm

# Model function
def tfmod(L, gas, liq, Q, nc, cg0, cl0, cgin, ka, k, henry, temp, dens, times):

    # L = total length/height (m)
    # gas = gas phase porosity (m3/m3)
    # liq = liquid phase content (m3/m3)
    # Q = gas flow rate (m3/s)
    # nc = number of cells (layers)
    # cg0 = initial compound concentration in gas phase (g/m3)
    # cl0 = initial compound concentration in liquid phase (g/m3)
    # cgin = compound concentration in inflow (g/m3)
    # ka = mass transfer coefficient for gas -> liquid (g/s / (g/m3) -> m3/s)
    # k = first-order liquid phase reaction rate constant (1/s)
    # henry = Henry's law constant coefficients as [k_H at 25 C, d(ln(kH)) / d(1/T)] as in NIST web book
    # temp = temperature (degrees C)
    # dens = solution (liquid) density (kg/m3)

    R = 0.083144 # L bar / K-mol
    
    # Make sure some inputs are numeric to avoid integer math bug
    cg0 = float(cg0)
    cl0 = float(cl0)
    cgin = float(cgin)
    ka = float(ka)
    k = float(k)
    temp = float(temp)

    # Temperature and Henry's law constant
    TK = temp + 273.15
    kh = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15)) # mol/kg-bar as liq:gas
    kh = kh * dens # mol/L-bar
    hc = kh * R * TK # dimensionless, liq:gas, e.g., g/L / g/L
    
    # Create cells
    x = np.linspace(0, L, nc + 1)  # Original nc + 1 values
    dx = np.diff(x)[0]             # dx, single value
    x = x[1:(nc + 1)] - dx / 2     # Center position in m
    
    # Cell gas volume (m3/m2)
    vg = dx * 1 * gas 
    # Cell liquid volume (m3/m2)
    vl = dx * 1 * liq 
    
    # Compound conc and mass (g)
    # cc = concentration, mc = mass [position]
    # g = gas, l = liquid
    ccg = np.full((nc), cg0)
    ccl = np.full((nc), cl0)
    mcg = ccg * vg
    mcl = ccl * vl

    # Initial state variable array (Python language is terrible)
    y0 = np.concatenate([mcg, mcl])
    
    # Solve/integrate
    out = solve_ivp(rates, [0, max(times)], y0 = y0, 
                    t_eval = times, 
                    args = (Q, cgin, vg, vl, k, ka, hc))
    
    # Extract mass of compound [position, time]
    mcgt = out.y[0:nc]
    mclt = out.y[nc:(2 * nc)]
    # Get concentrations
    ccgt = mcgt / np.transpose(np.tile(vg, (mcgt.shape[1], 1)))
    cclt = mclt / np.transpose(np.tile(vl, (mclt.shape[1], 1)))

    mct = np.concatenate([mcgt, mclt])
    
    # Return results as a tuple
    return ccgt, cclt, mcgt, mclt, x, times
