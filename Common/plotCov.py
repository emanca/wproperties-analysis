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
sys.path.append('data/')
from binning import yBins, qtBins, ptBins, etaBins, mTBins, isoBins, chargeBins
from root_numpy import hist2array
import argparse

# matplotlib stuff
plt.style.use([hep.style.CMS])

parser = argparse.ArgumentParser('')
parser.add_argument('-asimov', '--asimov', default=False, action='store_true', help='plot asimov result')

args = parser.parse_args()
asimov = args.asimov
# file with fit results
f = ROOT.TFile.Open('../Fit/FitRes/fit_WPlus_{}.root'.format('asimov' if asimov else 'data'))

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
yBins_red = np.array(yBins[:threshold_y+1])
qtBins_red = np.array(qtBins[:threshold_qt+1])
yBinsC = 0.5*(yBins_red[1:]+yBins_red[:-1])
qtBinsC = 0.5*(qtBins_red[1:]+qtBins_red[:-1])

nBins=len(yBinsC)*len(qtBinsC)
bins=np.linspace(0,nBins+1,nBins+1)
yticks = np.round(qtBins[:threshold_qt+1],0)

# helXsecs = ['helicity L', 'helicity I', 'helicity T', 'helicity A', 'helicity P','helicity UL']
helXsecs = ['unpolarised cross section','A0','A1','A2','A3','A4']

hcov = f.Get('correlation_matrix_channelhelpois')
cov = hist2array(hcov)[:,:]

# retrieve impacts per group of POI
for i,hel in enumerate(helXsecs):
    vcovreduced = cov[i*nBins:(i+1)*nBins,i*nBins:(i+1)*nBins]

    fig, ax1 = plt.subplots()
    plt. text(0.1, 0.9,"{}".format(hel), ha='center', va='center', transform=ax1.transAxes, color="k", weight="bold")
    hep.hist2dplot(vcovreduced,bins,bins, cmap ='jet', vmin =-1, vmax=1)
    plt.tight_layout()
    plt.savefig("cov{}.png".format(hel))
    plt.savefig("cov{}.pdf".format(hel))
    plt.clf()

# hcov = f.Get('correlation_matrix_channelhelmetapois')
# cov = hist2array(hcov)[:,:]

# vcovreduced = cov[0:6*6,0:6*6]
# bins=np.linspace(0,36+1,36+1)

# fig, ax1 = plt.subplots()
# hep.hist2dplot(vcovreduced,bins,bins, cmap ='jet', vmin =-1, vmax=1)
# plt.tight_layout()
# plt.savefig("covy.png")
# plt.savefig("covy.pdf")
# plt.clf()

# vcovreduced = cov[6*6:84,6*6:84]
# bins=np.linspace(0,49,49)

# fig, ax1 = plt.subplots()
# hep.hist2dplot(vcovreduced,bins,bins, cmap ='jet', vmin =-1, vmax=1)
# plt.tight_layout()
# plt.savefig("covqt.png")
# plt.savefig("covqt.pdf")
# plt.clf()

vcovreduced = cov[0:48*6,0:48*6]
bins=np.linspace(0,48*6+1,48*6+1)

fig, ax1 = plt.subplots()
hep.hist2dplot(vcovreduced,bins,bins, cmap ='jet', vmin =-1, vmax=1)
plt.xticks(np.arange(min(bins), max(bins)+1, 48))
plt.yticks(np.arange(min(bins), max(bins)+1, 48))
# Turn off tick labels
ax1.set_yticklabels([])
ax1.set_xticklabels([])
plt.tight_layout()
plt.savefig("covtot.png")
plt.savefig("covtot.pdf")
plt.clf()