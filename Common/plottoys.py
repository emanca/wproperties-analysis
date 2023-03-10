import ROOT
import os
import sys
import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from root_numpy import array2hist
sys.path.append('data/')
from binning import yBins, qtBins, ptBins, etaBins, mTBins, isoBins, chargeBins

ROOT.gROOT.SetBatch(True)

plt.style.use([hep.style.ROOT])
# hep.cms.label(loc=0, year=2016, lumi=35.9, data=True)
helicities = ['UL','I', 'T', 'A', 'P']
coefficients = ["A0", "A1", "A2", "A3", "A4",'unpolarizedxsec']

threshold_y = np.digitize(2.4,yBins)
threshold_qt = np.digitize(60.,qtBins)
yBins = np.array(yBins[:threshold_y],dtype='float64')
qtBins = np.array(qtBins[:threshold_qt],dtype='float64')
print(yBins,qtBins)
yBinsC = 0.5*(yBins[1:]+yBins[:-1])
qtBinsC = 0.5*(qtBins[1:]+qtBins[:-1])

# file with reweighted coefficients
file_preVFP = '/scratchnvme/emanca/wproperties-analysis/config/powheg_acc_rew_preVFP/WPlusJetsToMuNu_helweights.hdf5'
file_postVFP = '/scratchnvme/emanca/wproperties-analysis/config/powheg_acc_rew_postVFP/WPlusJetsToMuNu_helweights.hdf5'

f_preVFP = h5py.File(file_preVFP, mode='r+')
f_postVFP = h5py.File(file_postVFP, mode='r+')
htot_preVFP = f_preVFP['totxsecs'][:]
htot_postVFP = f_postVFP['totxsecs'][:]
h_preVFP = f_preVFP['xsecs'][:]
h_postVFP = f_postVFP['xsecs'][:]

htot = htot_preVFP+htot_postVFP
h = h_preVFP+h_postVFP
# shape h: y, qt, weights, pdf
# shape tot: y, qt, pdf
factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
factors = factors[np.newaxis,np.newaxis,...]
factors_hel = np.array([2.,2*math.sqrt(2),4.,4.*math.sqrt(2),2.,2.,2.*math.sqrt(2),4.*math.sqrt(2),1.])
factors_hel = factors_hel[np.newaxis,np.newaxis,...]

h = (h/htot[...,np.newaxis]+factors[...,1])*factors[...,0]

# h = h/factors_hel
h = h[:threshold_y,:threshold_qt,:]
h[...,-1] = 3./(16.*math.pi)*htot[:threshold_y,:threshold_qt]


hists = {}
canvs = {}
hpull = ROOT.TH1D("mass_pull","mass_pull", 100, -5.,5.)
canv = ROOT.TCanvas('mass','mass')

hists['mass']=hpull
canvs['mass']=canv
for k,c in enumerate(coefficients):
    for i in range(len(yBinsC)):
        for j in range(len(qtBinsC)):
            hpull = ROOT.TH1D("pull_y_{i}_qt_{j}_{c}".format(c=c, j=j, i=i),"pull_y_{i}_qt_{j}_{c}".format(c=c, j=j, i=i), 100, -5.,5.)
            canv = ROOT.TCanvas('pull_y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i))
            hists['y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i)]=hpull
            canvs['y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i)]=canv
for ifile in range(40):
    fIntoy = ROOT.TFile.Open('../Fit/FitRes/fit_WPlus_data_{}.root'.format(ifile))
    fitresultstoy = fIntoy.Get('fitresults')
    for ev in fitresultstoy:
        for key, hist in hists.items():
            if 'mass' in key: 
                par=(eval('ev.{}'.format(key)))
            else:
                i = int(key.split('_')[1])
                j = int(key.split('_')[3])
                coeff = coefficients.index(key.split('_')[4])
                print(key, i,j,coeff)
                if not 'unpol' in key:
                    par=(eval('ev.{}'.format(key)))-h[i,j,coeff]
                else:
                    par=eval('ev.{}'.format(key))-h[i,j,-1]
                    print(par/par_err)
            par_err=(eval('ev.{}_err'.format(key)))
            hist.Fill(par/par_err)

ROOT.gStyle.SetOptFit(1)
for i,c in canvs.items():
    c.cd()
    hists[i].Fit('gaus')
    hists[i].Draw()
    c.SaveAs('toys/{}.png'.format(hists[i].GetName()))
