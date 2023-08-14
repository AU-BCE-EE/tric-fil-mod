# Check output

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shutil.copy('../../mod_funcs.py', '.')
from mod_funcs import tfmod


# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See notes in tfmod.py for more complete descriptions
L = 0.1            # Filter length/depth (m); zooming in on first 10 cm!
por_g = 0.5      # (m3/m3)
por_l = 0.25     # (m3/m3)
v_g = 0.03       
v_l = 2E-5       
nc = 5          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
henry = (0.1, 2000.)
temp = 15.       # (degrees C)
dens_l = 1000    # Liquid density (kg/m3)

k = 0         # Reaction rate (1/s)

pH = 7.
pKa = 7.

cgin = 1

# Fixed for water
clin = 0.        # Fresh water concentration (g/m3)

# Times for model output
tt = 3 
# Number of time rows
nt = 10
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pH = 7.
pred1 = tfmod(L = L, por_g = por_g, por_l = por_l, v_g = v_g, v_l = v_l, nc = nc, cg0 = cg0, 
              cl0 = cl0, cgin = cgin, clin = clin, Kga = 'onda', k = k, henry = henry, pKa = pKa, 
              pH = pH, temp = temp, dens_l = dens_l, times = times)

# Check output for inputs
pred1['inputs']

# Check por_g for example
pred1['inputs']['por_g']

print(pred1['inputs']['por_g'])




