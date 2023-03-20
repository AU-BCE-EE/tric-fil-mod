# Demo on Kga_onda function

# Import necessary packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import shutil

# Import our model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Easiest if we have the module in this directory, so first the current version is copied in
# We will probably come up with something more sophisticated eventually
shutil.copy('../../Kga_onda.py', '.')
from Kga_onda import Kga_onda

# Set input values in function call
Kga = Kga_onda(pH = 7, temp = 20, henry = (0.1, 2000.), pKa = 7., P = 1., ssa = 1084., v_g = 0.03, 
               v_l = 0.003, por_g = 0.28)

print('Kga(onda) =', "{:.2e}".format(Kga), ' sec-1')

file1 = open('result.txt', 'w')
file1.write(str(Kga))
file1.close()
