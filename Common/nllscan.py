import ROOT
import os
import sys
import h5py
import math
import numpy as np
# import jax
# import jax.numpy as jnp
import matplotlib.pyplot as plt
import mplhep as hep
import scipy.optimize as opt

def func(x, a):
     return a * x**2

scans = np.empty((11,))
for i in range(5,16):
    fIn = ROOT.TFile.Open('../Fit/fit_WPlus_mass{}.root'.format(i))
    fitresults = fIn.Get('fitresults')

    for ev in fitresults:
        nllfull = eval('ev.nllvalfull')
        scans[i-5]=nllfull
# bins = np.concatenate((np.linspace(-100,0,11),np.linspace(10,100,10)))
bins = np.concatenate((np.linspace(-50,0,6),np.linspace(10,50,5)))
print(bins)
scans-=scans[5]
scans*=2.

fig, ax1 = plt.subplots()
optimizedParameters, pcov = opt.curve_fit(func, bins, scans)
plt.scatter(bins,scans)
plt.plot(bins, func(bins, *optimizedParameters), label="fit")
plt.plot(bins,np.ones((11,)))
# ax1.set_ylim(0.,10.)
# ax1.set_xlim(-10.,10.)
print(*optimizedParameters)
plt.show()