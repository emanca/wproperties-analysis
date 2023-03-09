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
from root_numpy import hist2array, array2hist
import argparse
import time


# matplotlib stuff
plt.style.use([hep.style.CMS])

coefficients = ['unpolarizedxsec']

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
yBins = np.array(yBins[:threshold_y+1],dtype='float64')
qtBins = np.array(qtBins[:threshold_qt+1],dtype='float64')
yBinsC = 0.5*(yBins[1:]+yBins[:-1])
qtBinsC = 0.5*(qtBins[1:]+qtBins[:-1])
yBinsS = yBins[1:]-yBins[:-1]
qtBinsS = qtBins[1:]-qtBins[:-1]

# file with fit results
fM = ROOT.TFile.Open('../Fit/FitRes/fit_WMinus_data.root')
fitresultsM = fM.Get('fitresults')

fP = ROOT.TFile.Open('../Fit/FitRes/fit_WPlus_data.root')
fitresultsP = fP.Get('fitresults')

hcovM = fM.Get('covariance_matrix_channelhelmetapois')
covM = hist2array(hcovM)[:,:]

hcovP = fP.Get('covariance_matrix_channelhelmetapois')
covP = hist2array(hcovP)[:,:]

print(hcovP.GetXaxis().GetBinLabel(7),hcovP.GetXaxis().GetBinLabel(8+6+1))

# get value of fitted parameters
hcoeffM = np.zeros([len(qtBinsC)])
hcoeff_genM = np.zeros([len(qtBinsC)])
hcoeff_errM = np.zeros([len(qtBinsC)])
for ev in fitresultsM: #dummy because there's one event only
    # for i in range(len(yBinsC)):
        for j in range(len(qtBinsC)):
            try:
                coeff = eval('ev.qt_{j}_helmeta_{c}'.format(c=coefficients[0], j=j))
                coeff_gen = eval('ev.qt_{j}_helmeta_{c}_gen'.format(c=coefficients[0], j=j))
                coeff_err = eval('ev.qt_{j}_helmeta_{c}_err'.format(c=coefficients[0], j=j))
                # coeff = eval('ev.y_{i}_qt_{j}_{c}'.format(c=coefficients[0], j=j, i=i))
                # coeff_gen = eval('ev.y_{i}_qt_{j}_{c}_gen'.format(c=coefficients[0], j=j, i=i))
                # coeff_err = eval('ev.y_{i}_qt_{j}_{c}_err'.format(c=coefficients[0], j=j, i=i))
                # coeff = coeff/(3./16./math.pi)/35.9/qtBinsS[j]
                # coeff_gen = coeff_gen/(3./16./math.pi)/35.9/qtBinsS[j]
                # coeff_err = coeff_err/(3./16./math.pi)/35.9/qtBinsS[j]
                hcoeffM[j]=coeff
                hcoeff_genM[j]=coeff_gen
                hcoeff_errM[j]=coeff_err
            except AttributeError:
                pass

hcoeffP = np.zeros([len(qtBinsC)])
hcoeff_genP = np.zeros([len(qtBinsC)])
hcoeff_errP = np.zeros([len(qtBinsC)])
for ev in fitresultsP: #dummy because there's one event only
    # for i in range(len(yBinsC)):
        for j in range(len(qtBinsC)):
            try:
                coeff = eval('ev.qt_{j}_helmeta_{c}'.format(c=coefficients[0], j=j))
                coeff_gen = eval('ev.qt_{j}_helmeta_{c}_gen'.format(c=coefficients[0], j=j))
                coeff_err = eval('ev.qt_{j}_helmeta_{c}_err'.format(c=coefficients[0], j=j))
                # coeff = eval('ev.y_{i}_qt_{j}_{c}'.format(c=coefficients[0], j=j, i=i))
                # coeff_gen = eval('ev.y_{i}_qt_{j}_{c}_gen'.format(c=coefficients[0], j=j, i=i))
                # coeff_err = eval('ev.y_{i}_qt_{j}_{c}_err'.format(c=coefficients[0], j=j, i=i))
                # coeff = coeff/(3./16./math.pi)/35.9/qtBinsS[j]
                # coeff_gen = coeff_gen/(3./16./math.pi)/35.9/qtBinsS[j]
                # coeff_err = coeff_err/(3./16./math.pi)/35.9/qtBinsS[j]
                hcoeffP[j]=coeff
                hcoeff_genP[j]=coeff_gen
                hcoeff_errP[j]=coeff_err
            except AttributeError:
                pass

fig, (ax1,ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
hep.cms.text('work in progress', loc=0, ax=ax1)
ax2.set_xlabel('$q_T$ (GeV)')
ax1.set_ylabel(r'$\frac{d\sigma^{U+L}}{dq_T} (fb/GeV)$')
ax1.errorbar(qtBinsC,hcoeffM.ravel(), xerr=qtBinsS/2,marker = 'o', label="$W^-$ data", fmt='o')
ax1.errorbar(qtBinsC,hcoeffP.ravel(), xerr=qtBinsS/2,marker = 'o', label="$W^+$ data", fmt='o')
bins=np.append(np.tile(qtBins[:-1],len(yBinsC)),60.)
x=np.array(range(len(bins)))
# hep.histplot(hcoeff_errM.ravel()/hcoeffM.ravel(),bins = x, label="$W^-$", ax=ax1)
# hep.histplot(hcoeff_errP.ravel()/hcoeffP.ravel(),bins = x, label="$W^+$", ax=ax1)
# print('cov Wm:',covM[7-1:8+6,7-1:8+6])
# print('cov Wp:',covP[7-1:8+6,7-1:8+6])
# print('qt Wm:', hcoeffM)
# print('qt Wp:', hcoeffP)

tot_M = np.sum(hcoeffM)
tot_M_gen = np.sum(hcoeff_genM)
err_M = np.sqrt(np.matmul(np.ones(hcoeffM.shape[0]), np.matmul(covM[7-1:8+6,7-1:8+6],np.ones(hcoeffM.shape[0]))))

tot_P = np.sum(hcoeffP)
tot_P_gen = np.sum(hcoeff_genP)
err_P = np.sqrt(np.matmul(np.ones(hcoeffP.shape[0]), np.matmul(covP[7-1:8+6,7-1:8+6],np.ones(hcoeffP.shape[0]))))

print("tot W-:",tot_M/1000, "+/-", err_M/1000, "prediction:", tot_M_gen/1000)
print("tot W+:",tot_P/1000, "+/-", err_P/1000, "prediction:", tot_P_gen/1000)



ratio = hcoeffM/hcoeffP
err_ratio = np.hypot(hcoeff_errM/hcoeffM, hcoeff_errP/hcoeffP)*ratio
ax2.errorbar(qtBinsC,ratio,marker = 'o', color = "k", fmt='o', xerr=qtBinsS/2,yerr=err_ratio)
ax2.plot(qtBinsC,hcoeff_genM/hcoeff_genP) 
ax1.legend(loc='upper right',frameon=False)
plt.tight_layout()
plt.savefig("ratio_qt.png")

# hratio = ROOT.TH1D("ratio", "ratio", len(qtBins)-1, qtBins)
# array2hist(ratio,hratio, errors=err_ratio)
# print(hcoeff_genM/hcoeff_genP)
# fit=hratio.Fit('pol0', 'S')
# hratio.Draw()
# time.sleep(100)