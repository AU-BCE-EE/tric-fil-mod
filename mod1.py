# Simple 1D diffusion model to start to think about reaction-transport models in Python

# Packages
import numpy as np
import pdb
from scipy.integrate import solve_ivp

# Rates function
# Arg order: time, state variable, then arguments
def rates(t, mc, Q, cin, gv, k):
    
    #pdb.set_trace()
    cc = mc / gv
    nc = cc.shape[0]
    dm = np.zeros(cc.shape[0])
    # First cell derivative
    dm[0] = Q * (cin - cc[0]) - k * mc[0] 
    # Other cell derivatives
    dm[1:nc] = Q * -np.diff(cc) - k * mc[1:nc] # Check nc or nc - 1
    
    return dm

# Model function
def tfmod(L, por, Q, nc, c0, cin, k, times):

    # L = total length/height (m)
    # por = gas phase porosity (m3/m3)
    # Q = gas flow rate (m3/s)
    # nc = number of cells (layers)
    # c0 = initial compound concentration (g/m3)
    # cin = compound concentration in inflow (g/m3)
    # k = first-order reaction rate (1/s)

    # Make sure cc is numeric to avoid integer math bug
    c0 = float(c0)
    
    # Create cells
    x = np.linspace(0, L, nc + 1)  # Original nc + 1 values
    dx = np.diff(x)[0]             # dx, single value
    x = x[1:(nc + 1)] - dx / 2     # Center position in m
    
    # Cell gas volume (m3/m2)
    gv = dx * 1 * por 
    
    # Compound conc and mass (g)
    # cc = concentration, mc = mass [position]
    cc = np.full((nc), c0)
    mc = cc * gv

    
    # Solve/integrate
    out = solve_ivp(rates, [0, max(times)], y0 = mc, 
                    t_eval = times, 
                    args = (Q, cin, gv, k))
    
    # Extract mass of compound [position, time]
    mct = out.y
    cct = mct / np.transpose(np.tile(gv, (mct.shape[1], 1)))
    
    # Return results as tuple
    return cct, mct, x, times
