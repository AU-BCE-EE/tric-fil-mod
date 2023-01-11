# Demo of Python air filter model

# Import necessary packages
import shutil
import numpy as np
import matplotlib.pyplot as plt

# Import our model
# Easiest if we have the module in this directory, so first the current version is copied in
# We will probably come up with something more sophisticated eventually
shutil.copy('../../tfmod.py', '.')
from tfmod import tfmod

# Set model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L = 2            # Filter length/depth (m)
gas = 0.5        # (m3/m3)
liq = 0.25       # (m3/m3)
Q = 0.001        # Air flow (m3/s)
nc = 8           # Number of model cells (layers)
cg0 = 0          # (g/m3)
cl0 = 0          # (g/m3)
cgin = 1         # Dirty air compound concentration (g/m3)
henry = (0.1, 2000.)
temp = 15        # (degrees C)
dens = 1000      # Liquid density (kg/m3)

Kga = 1E-3       # Mass transfer coefficient (1/s)
k = 0. / 3600    # First-order degradation/removal rate (1/h -> 1/s)

# Times for model output, calculated from tt and nt here but could be set directly
# Total duration (hours)
tt = 2 
# Number of time rows
nt = 500
times = np.linspace(0, tt, nt) * 3600
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sim 1 no reaction ~~~~~
Kga = 1E-3 
k = 0. / 3600
pred1 = tfmod(L, gas, liq, Q, nc, cg0, cl0, cgin, Kga, k, henry, temp, dens, times)

# Sim 2 ~~~~~~~~~~~~~~~~~
Kga = 1E-3
k = 1. / 3600
pred2 = tfmod(L, gas, liq, Q, nc, cg0, cl0, cgin, Kga, k, henry, temp, dens, times)

# Sim 3 ~~~~~~~~~~~~~~~~~
Kga = 1E-3
k = 10. / 3600
pred3 = tfmod(L, gas, liq, Q, nc, cg0, cl0, cgin, Kga, k, henry, temp, dens, times)

# Sim 4 ~~~~~~~~~~~~~~~~~
Kga = 1E-5       # (1/s)
k = 10. / 3600   # First-order degradation/removal rate (1/h -> 1/s)
pred4 = tfmod(L, gas, liq, Q, nc, cg0, cl0, cgin, Kga, k, henry, temp, dens, times)

# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot outlet concentration (= 1 - removal efficiency here because cgin = 1)
# Gas concentration (outlet air) 
plt.plot(pred1[5] / 3600, pred1[0][nc - 1, :])
plt.plot(pred2[5] / 3600, pred2[0][nc - 1, :])
plt.plot(pred3[5] / 3600, pred3[0][nc - 1, :])
plt.plot(pred4[5] / 3600, pred4[0][nc - 1, :])
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('outlet_gas_conc.png')

# Liquid concentration (in last layer)
plt.clf()
plt.plot(pred1[5] / 3600, pred1[1][nc - 1, :])
plt.plot(pred2[5] / 3600, pred2[1][nc - 1, :])
plt.plot(pred3[5] / 3600, pred3[1][nc - 1, :])
plt.plot(pred4[5] / 3600, pred4[1][nc - 1, :])
plt.xlabel('Time (h)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('outlet_liq_conc.png')

# Profiles
# Gas
plt.clf()
plt.plot(pred1[4], pred1[0][:, nt - 1])
plt.plot(pred2[4], pred2[0][:, nt - 1])
plt.plot(pred3[4], pred3[0][:, nt - 1])
plt.plot(pred4[4], pred4[0][:, nt - 1])
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_gas_conc.png')

# Liquid
plt.clf()
plt.plot(pred1[4], pred1[1][:, nt - 1])
plt.plot(pred2[4], pred2[1][:, nt - 1])
plt.plot(pred3[4], pred3[1][:, nt - 1])
plt.plot(pred4[4], pred4[1][:, nt - 1])
plt.xlabel('Location (m)')
plt.ylabel('Compound conc. (g/m3)')
plt.savefig('profile_liq_conc.png')
