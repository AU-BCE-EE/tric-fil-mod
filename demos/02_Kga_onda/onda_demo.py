# Demo on Kga_onda function

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil

# Import our model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Easiest if we have the module in this directory, so first the current version is copied in
# We will probably come up with something more sophisticated eventually
shutil.copy('../../Kga_onda.py', '.')
from Kga_onda import Kga_onda

# Set input values
pH = 7.0
TK = 298.0
KH = 10.0
pKa = 7.0
P = 1.0   # bar
ssa = 1084.0    # m2/m3
Qstd = 2400.0   # m3/hour
Qliq = 20.0     # L/min
por = 0.28      # porosity
Ax = 21.0       # cross sectional area in m2
ssa = 1084.0    # m2/m3

Kga = Kga_onda(pH, TK, KH, pKa, P, ssa, Qstd, Qliq, por, Ax)
print('Kga(onda) =', "{:.2e}".format(Kga), ' sec-1')

