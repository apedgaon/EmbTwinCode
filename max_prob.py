# Finds max twin nucleation probability at a given increment number
# Usage: python max_prob.py 142000
import numpy as np
import sys

fnum = sys.argv[1]
fpath = 'kinc_nuc_prob_' + fnum + '.dat'

p = np.loadtxt(fpath)
print('Max nucleation probability = ', np.max(p))
print('Element number = ', np.argmax(p) + 1)