# Simple 1D diffusion model to start to think about reaction-transport models in Python

# Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pdb

# Set inputs
L = 3     # Total length/height
nc = 10  # Number of cells
tt = 600 # Time (s)
dt = 1    # Step (s)

# Filter properties
por = 0.5 # Gas phase porosity (m3/m3)
cc = 0.0    # Initial concentration (g/m3)

# Transport parameters
#D = 1E-6  # m2/s
Q = 0.001   # m3/s

# Reaction 
k = 1.E-5  # First-order 1/s

# Dirty air chemical concentration
cin = 0.1 # g/m3

# Make sure cc is numeric to avoid integer crap
cc = float(cc)

# Number of time steps 
nts = tt // dt 

# Create cells
x = np.linspace(0, L, nc + 1)  # Position in m
dx = np.diff(x)
x = x[1:(nc + 1)] - dx / 2

# Cell gas volume (m3/m2)
gv = dx * 1 * por 

# Compound conc and mass (g)
# cc = concentration, mc = mass, dims -> rows are cells (position), columns have time [position, time]
cc = np.full((nc), cc)
mc = cc * gv

# Rates function
# Arg order: time, state variable, then arguments
def rates(t, mc, Q, cin, gv, k):
    
    cc = mc / gv
    dm = np.zeros(cc.shape[0])
    # First cell derivative
    dm[0] = Q * (cin - cc[0]) - k * cc[0] 
    # Others
    dm[1:nc] = Q * - np.diff(cc) - k * cc[1:nc]
    
    return(dm)

# Solve/integrate
out = solve_ivp(rates, [0, 86400], y0 = mc, t_eval = np.linspace(0, 24, 24) * 3600, 
        args = (Q, cin, gv, k))

# Extract mass of compound [position, time]
mct = out.y
cct = mct / np.transpose(np.tile(gv, (mct.shape[1], 1)))

plt.plot(x, cct[:, 1]);plt.show()
plt.plot(x, cct[:, 23]);plt.show()

