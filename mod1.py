# Simple 1D diffusion model to start to think about reaction-transport models in Python

# Packages
import numpy as np
#import matplotlib.pyplot as plt

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
k = 1E-5  # First-order 1/s

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
cc = np.full((nc, nts), cc)
# NTS: can I do this in 1 step?
mc = cc * np.transpose(np.tile(gv, (nts, 1)))

# Initial derivative array dm/dt [position]
dm = np.zeros(nc)

# Transport and reaction loop
for t in range(1, nts) :
    # First cell derivative
    dm[0] = Q * (cin - cc[0, t-1]) - k * cc[0, t-1] 
    # Others
    dm[1:nc] = Q * - np.diff(cc[:, t-1]) - k * cc[1:nc, t-1]
    # Update mass
    mc[:, t] = mc[:, t-1] + dm * dt
    # Concentrations
    cc[:, t] = mc[:, t] / gv
#    print(t)
#    print(dm)
#    print(cc[:, t])

print(cc[9, :])

#plt.plot(x, cc[:, 1]);plt.show()
