# Demo of simple air filter model

import pdb
import numpy as np
import matplotlib.pyplot as plt

# Set model inputs
L = 2
gas = 0.5        #  (m3/m3)
liq = 0.25       #  (m3/m3)
Q = 0.01         # Air flow (m3/s)
nc = 10          # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
cgin = 2         # Dirty air compound concentration (g/L)
ka = 1E-3        # (m3/s)
k = 0.1 / 3600   # First-order degradation/removal rate (1/h -> 1/s)
henry = (1000., 8000.)
temp = 15
dens = 1000

# Times for model run
tt = 200 # hours
times = np.linspace(0, tt, 500) * 3600

from mod1 import tfmod
pred = tfmod(L, gas, liq, Q, nc, cg0, cl0, cgin, ka, k, henry, temp, dens, times)

## Plot concentration within reactor
#plt.plot(pred[2], pred[0][:, 0])
#plt.plot(pred[2], pred[0][:, 1])

# Plot outlet concentration (= 1 - removal efficiency here)
plt.plot(pred[5] / 3600, pred[0][nc - 1, :])
plt.savefig('outlet_gas_conc.png')

plt.plot(pred[5] / 3600, pred[1][nc - 1, :])
plt.savefig('outlet_liq_conc.png')


