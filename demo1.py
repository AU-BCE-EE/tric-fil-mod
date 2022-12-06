#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Call of air filter model
# def tfmod(L, por, Q, nc, c0, cin, k, times):

import numpy as np
import matplotlib.pyplot as plt

from mod1 import tfmod
    
# Model inputs
L = 2
por = 0.5        # Gas-phase porosity (m3/m3)
Q = 0.01         # Air flow (m3/s)
nc = 10          # Number of model cells (layers)
c0 = 0           # Initial compound concentration (g/m3)
cin = 1          # Dirty air compound concentration (g/m3)
k = 0.5 / 3600   # First-order degradation/removal rate (1/h -> 1/s)

# Times for model run
tt = 0.1 # hours
times = np.linspace(0, tt, 500) * 3600

pred = tfmod(L, por, Q, nc, c0, cin, k, times)

## Plot concentration within reactor
#plt.plot(pred[2], pred[0][:, 0])
#plt.plot(pred[2], pred[0][:, 1])

# Plot outlet concentration (= 1 - removal efficiency here)
plt.plot(pred[3] / 3600, pred[0][nc - 1, :])
plt.savefig('outlet_conc.png')


