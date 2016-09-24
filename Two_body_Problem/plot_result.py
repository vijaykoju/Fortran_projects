from matplotlib import pyplot as plt
from pylab import *
import pylab
import numpy as np

try:
	fname1 = 'positionData.dat'
	fname2 = 'energyData.dat'
	a, b, c, d = np.loadtxt(fname1, usecols = range(4), unpack=True)
	n, k, p, t = np.loadtxt(fname2, unpack=True)
except:
	err_msg = "Could not load data from file %s." % (fname1) or (fname2) or (fname1 and fname2)
	raise Exception(err_msg)

#plt.figure(1)
#plt.clf()
plt.plot(n,k)
plt.draw()
plt.show()
