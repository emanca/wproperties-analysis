import ROOT
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import sys
import h5py
import math
import numpy as np
from scipy import stats
# import jax
# import jax.numpy as jnp
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import mplhep as hep
sys.path.append('data/')
from binning import yBins, qtBins, ptBins, etaBins, mTBins, isoBins, chargeBins, yBins_val, qtBins_val
import argparse

parser = argparse.ArgumentParser('')
parser.add_argument('-asimov', '--asimov', default=False, action='store_true', help='plot asimov result')

args = parser.parse_args()
asimov = args.asimov

label_data = "asimov dataset" if asimov else "data"
# file with fit results
fIn = ROOT.TFile.Open('../Fit/FitRes/fit_Z_asimov_.root')
fitresults = fIn.Get('fitresults')


# matplotlib stuff
plt.style.use([hep.style.CMS])
helicities = ['UL','L', 'I', 'T', 'A', 'P']
coefficients = ['unpolarizedxsec', "A0", "A1", "A2", "A3", "A4"]

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
yBins = np.array(yBins)
qtBins = np.array(qtBins)
yBinsC = 0.5*(yBins[1:]+yBins[:-1])
qtBinsC = 0.5*(qtBins[1:]+qtBins[:-1])
yBinsS = yBins[1:]-yBins[:-1]
qtBinsS = qtBins[1:]-qtBins[:-1]


for k,c in enumerate(coefficients):
    hcoeff = np.zeros([len(yBinsC),len(qtBinsC)])
    hcoeff_gen = np.zeros([len(yBinsC),len(qtBinsC)])
    hcoeff_err = np.zeros([len(yBinsC),len(qtBinsC)])
    hcoeff_hel = np.zeros([len(yBinsC),len(qtBinsC)])
    hcoeff_hel_gen = np.zeros([len(yBinsC),len(qtBinsC)])
    hcoeff_hel_err = np.zeros([len(yBinsC),len(qtBinsC)])
    bins=np.append(np.tile(qtBins[:-1],len(yBinsC)),60.)
    x=np.array(range(len(bins)))
    binsC = 0.5*(x[1:]+x[:-1])
    print('analysing', c, helicities[k])
    for ev in fitresults: #dummy because there's one event only
        for i in range(len(yBinsC)):
            for j in range(len(qtBinsC)):
                try:
                    coeff = eval('ev.y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i))
                    coeff_gen = eval('ev.y_{i}_qt_{j}_{c}_gen'.format(c=c, j=j, i=i))
                    coeff_err = eval('ev.y_{i}_qt_{j}_{c}_err'.format(c=c, j=j, i=i))
                    coeff_hel = eval('ev.helXsecs{c}_y_{i}_qt_{j}_pmaskedexp'.format(c=helicities[k], j=j, i=i))
                    coeff_hel_err = eval('ev.helXsecs{c}_y_{i}_qt_{j}_pmaskedexp_err'.format(c=helicities[k], j=j, i=i))
                    coeff_hel_gen = eval('ev.helXsecs{c}_y_{i}_qt_{j}_pmaskedexp_gen'.format(c=helicities[k], j=j, i=i))
                    if 'unpolarizedxsec' in c:
                        coeff = coeff/(3./16./math.pi)/16.8/yBinsS[i]/qtBinsS[j]
                        coeff_gen = coeff_gen/(3./16./math.pi)/16.8/yBinsS[i]/qtBinsS[j]
                        coeff_err = coeff_err/(3./16./math.pi)/16.8/yBinsS[i]/qtBinsS[j]
                    hcoeff[i,j]=coeff
                    hcoeff_gen[i,j]=coeff_gen
                    hcoeff_err[i,j]=coeff_err
                    hcoeff_hel[i,j]=coeff_hel
                    hcoeff_hel_gen[i,j]=coeff_hel_gen
                    hcoeff_hel_err[i,j]=coeff_hel_err
                except AttributeError:
                    pass
    fig, (ax1, ax2) = plt.subplots(nrows=2,figsize=(12, 10),gridspec_kw={'height_ratios': [3, 1]})
    # fig, ax1 = plt.subplots(figsize=(48, 10))
    hep.cms.text('work in progress', loc=1, ax=ax1)
    if 'unpol' in c:
        ax1.set_ylabel(r'$\frac{d\sigma^{U+L}}{dq_Td|y|} (fb/GeV)$', fontsize=30)
    else:
        ax1.set_ylabel(r'$A_{}$'.format(k-1))
    # ax1.set_xlabel('$q_T$ (GeV)')
    # ax1.set_xticks(x) # set tick positions
    # ticks = ["{:d}".format(int(v)) for v in bins]
    # for itick,tick in enumerate(ticks):
    #     if not itick==1:
    #         if tick=="0": ticks[itick] = "60_0"
    # ax1.set_xticklabels(ticks)
    # print(hcoeff.ravel(), "+/-", hcoeff_err.ravel())
    hep.histplot(hcoeff.ravel(),bins = x, yerr = hcoeff_err.ravel(),histtype = 'errorbar', color = "k", stack = False, ax=ax1, label=label_data)
    hep.histplot(hcoeff_gen.ravel(),bins =x, color = "r", stack = False, ax=ax1, label="prediction")

    if 'unpol' in c:
        ratio = hcoeff.ravel()/hcoeff_gen.ravel()
        ratio_err= hcoeff_err.ravel()/hcoeff_gen.ravel()
        ax2.set_ylabel('data/prediction')
        ax2.set_ylim([0, 2])
    else:
        ratio = hcoeff.ravel()-(hcoeff_gen.ravel())
        ratio_err = hcoeff_err.ravel()
        ax1.set_ylim(-5,5)
        ax2.set_ylabel('data-prediction')
        ax2.set_ylim([-2, 2])
    hep.histplot(ratio,bins = x, yerr=ratio_err.ravel(), histtype = 'errorbar', color = "k", stack = False, ax=ax2)
    ax2.set_xlabel('unrolled $q_T$-y bins')
    ax1.legend(loc='upper right', frameon=False)
    plt.tight_layout()
    plt.savefig('POIplots/fit{}.png'.format(c),dpi=300)
    plt.savefig('POIplots/fit{}.pdf'.format(c),dpi=300)
    plt.clf()
    # plot helicity cross sections
    fig, (ax1, ax2) = plt.subplots(nrows=2,figsize=(48, 10),gridspec_kw={'height_ratios': [3, 1]})
    ax1.set_title("fitted {}".format(helicities[k]), fontsize=18)
    ax1.set_ylabel('')
    ax1.set_xlabel('a.u.')
    ax1.set_xticks(x) # set tick positions
    # ax1.set_xticklabels(ticks)
    # print(hcoeff.ravel(), "+/-", hcoeff_err.ravel())
    # print(hcoeff_hel.ravel(), "+/-", hcoeff_hel_err.ravel())
    
    hep.histplot(hcoeff_hel.ravel(),bins = x, yerr = hcoeff_hel_err.ravel(),histtype = 'errorbar', color = "k", stack = False, ax=ax1, label=label_data)
    hep.histplot(hcoeff_hel_gen.ravel(),bins =x, color = "r", stack = False, ax=ax1, label="prediction")

    ax1.legend(loc='upper right', frameon=False)
    plt.tight_layout()
    plt.savefig('POIplots/fit{}.png'.format(helicities[k]),dpi=300)
    plt.savefig('POIplots/fit{}.pdf'.format(helicities[k]),dpi=300)
    plt.clf()

# integrated coefficients
for k,c in enumerate(coefficients):
    print('analysing', c, helicities[k])
    hcoeff = np.zeros([len(qtBinsC)])
    hcoeff_gen = np.zeros([len(qtBinsC)])
    hcoeff_err = np.zeros([len(qtBinsC)])
    hcoeff_y = np.zeros([len(yBinsC)])
    hcoeff_gen_y = np.zeros([len(yBinsC)])
    hcoeff_err_y = np.zeros([len(yBinsC)])
    hcoeff_hel = np.zeros([len(qtBinsC)])
    hcoeff_hel_gen = np.zeros([len(qtBinsC)])
    hcoeff_hel_err = np.zeros([len(qtBinsC)])
    for ev in fitresults: #dummy because there's one event only
        for j in range(len(qtBinsC)):
            try:
                coeff = eval('ev.qt_{j}_helmeta_{c}'.format(c=c, j=j))
                coeff_gen = eval('ev.qt_{j}_helmeta_{c}_gen'.format(c=c, j=j))
                coeff_err = eval('ev.qt_{j}_helmeta_{c}_err'.format(c=c, j=j))
                coeff_hel = eval('ev.helXsecs{c}_qt_{j}_sumxsec'.format(c=helicities[k], j=j))
                coeff_hel_gen = eval('ev.helXsecs{c}_qt_{j}_sumxsec_gen'.format(c=helicities[k], j=j))
                coeff_hel_err = eval('ev.helXsecs{c}_qt_{j}_sumxsec_err'.format(c=helicities[k], j=j))
                if 'unpol' in c:
                    coeff = coeff/(3./16./math.pi)/16.8/qtBinsS[j]
                    coeff_gen = coeff_gen/(3./16./math.pi)/16.8/qtBinsS[j]
                    coeff_err = coeff_err/(3./16./math.pi)/16.8/qtBinsS[j]
                hcoeff[j]=coeff
                hcoeff_gen[j]=coeff_gen
                hcoeff_err[j]=coeff_err
                hcoeff_hel[j]=coeff_hel
                hcoeff_hel_gen[j]=coeff_hel_gen
                hcoeff_hel_err[j]=coeff_hel_err
            except AttributeError:
                pass
        for i in range(len(yBinsC)):
            try:
                coeff_y = eval('ev.y_{i}_helmeta_{c}'.format(c=c, j=j, i=i))
                coeff_gen_y = eval('ev.y_{i}_helmeta_{c}_gen'.format(c=c, j=j, i=i))
                coeff_err_y = eval('ev.y_{i}_helmeta_{c}_err'.format(c=c, j=j, i=i))
                if 'unpol' in c:
                    coeff_y = coeff_y/(3./16./math.pi)/16.8/yBinsS[i]
                    coeff_gen_y = coeff_gen_y/(3./16./math.pi)/16.8/yBinsS[i]
                    coeff_err_y = coeff_err_y/(3./16./math.pi)/16.8/yBinsS[i]
                hcoeff_y[i]=coeff_y
                hcoeff_gen_y[i]=coeff_gen_y
                hcoeff_err_y[i]=coeff_err_y
            except AttributeError:
                pass
    fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
    hep.cms.text('work in progress', loc=1, ax=ax1)
    if 'unpol' in c:
        ax1.set_ylabel(r'$\frac{d\sigma^{U+L}}{dq_T} (fb/GeV)$', fontsize=30)
    else:
        ax1.set_ylabel(r'$A_{}$'.format(k-1))
    ax2.set_xlabel('$q_T$ (GeV)')
    ax1.errorbar(qtBinsC,hcoeff.ravel(), xerr=qtBinsS/2, yerr = hcoeff_err.ravel(),marker = 'o',color = "k", label=label_data, fmt='o')
    ax1.errorbar(qtBinsC,hcoeff_gen.ravel(), xerr=qtBinsS/2)

    plt.tight_layout()
    plt.savefig('POIplots/fitintegratedZ{}.png'.format(c),dpi=300)
    plt.savefig('POIplots/fitintegratedZ{}.pdf'.format(c),dpi=300)
    plt.clf()

    # plots y-integrated
    fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
    hep.cms.text('work in progress', loc=1, ax=ax1)
    if 'unpol' in c:
        ax1.set_ylabel(r'$\frac{d\sigma^{U+L}}{dy} (fb)$', fontsize=30)
    else:
        ax1.set_ylabel(r'$A_{}$'.format(k-1))
    ax2.set_xlabel('$y$')
    print(yBinsC,'bin centers')
    ax1.errorbar(yBinsC,hcoeff_y.ravel(), xerr=yBinsS/2, yerr = hcoeff_err_y.ravel(),marker = 'o',color = "k", label=label_data, fmt='o')
    ax1.errorbar(yBinsC,hcoeff_gen_y.ravel(), xerr=yBinsS/2)
    plt.tight_layout()
    plt.savefig('POIplots/fityintegratedZ{}.png'.format(c),dpi=300)
    plt.savefig('POIplots/fityintegratedZ{}.pdf'.format(c),dpi=300)
    plt.clf()
