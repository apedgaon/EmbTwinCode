# Finds max twin nucleation probability at the latest increment
# Helpful for quickly querrying current max probability
# Usage: python max_prob.py
import numpy as np

t = np.loadtxt('kinc_time.dat')
kinc = int(t[-1, 0])

fnuc = 'kinc_nuc_prob_' + str(kinc) + '.dat'
p = np.loadtxt(fnuc)
print('kinc = ', kinc)
print('Max probability = ', np.max(p))
print('Element number = ', np.argmax(p) + 1)
