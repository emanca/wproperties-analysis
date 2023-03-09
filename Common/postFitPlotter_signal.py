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
from root_numpy import hist2array
from binning import yBins, qtBins, ptBins, etaBins, mTBins, isoBins, chargeBins, yBins_val, qtBins_val
import argparse

def make_error_boxes(ax, xdata, ydata, xerror, yerror, facecolor='green',
                     edgecolor='none', alpha=0.5, label=""):
                     
    # Loop over data points; create box from errors at each point
    errorboxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
                  for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T)]

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor)

    # Add collection to axes
    ax.add_collection(pc)

    # Plot errorbars
    artists = ax.plot(xdata, ydata,
                        linestyle = 'None')
    return artists

parser = argparse.ArgumentParser('')
parser.add_argument('-asimov', '--asimov', default=False, action='store_true', help='plot asimov result')

args = parser.parse_args()
asimov = args.asimov

label_data = "asimov dataset" if asimov else "data"
# file with fit results
fIn = ROOT.TFile.Open('../Fit/FitRes/fit_WPlus_{}_blockAi.root'.format('asimov' if asimov else 'data'))
fitresults = fIn.Get('fitresults')

hcov = fIn.Get('covariance_matrix_channelhelmetapois')
cov = hist2array(hcov)[:,:]

# matplotlib stuff
plt.style.use([hep.style.CMS])
# helicities = ['UL','L', 'I', 'T', 'A', 'P']
# coefficients = ['unpolarizedxsec', "A0", "A1", "A2", "A3", "A4"]
helicities = ['UL']
coefficients = ['unpolarizedxsec']


threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1

# fIntoy = ROOT.TFile.Open('../Fit/fit_Wplus.root')
# fitresultstoy = fIntoy.Get('fitresults')

factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
factors = factors[np.newaxis,np.newaxis,...]
factors_hel = np.array([2.,2*math.sqrt(2),4.,4.*math.sqrt(2),2.,2.,2.*math.sqrt(2),4.*math.sqrt(2),1.])
factors_hel = factors_hel[np.newaxis,np.newaxis,...]

charge = "WPlus"
f_aMC = ROOT.TFile.Open('/scratchnvme/wmass/REWEIGHT/genInfo_syst.root')
qt_aMC = np.sum(hist2array(f_aMC.Get('angularCoefficients_{}/YqTcT'.format("WPlus" if charge=="WMinus" else "Wplus"))),axis=-1)
fcoeffs = ROOT.TFile.Open('/scratchnvme/wmass/REWEIGHT/genInput_v7_syst_{}.root'.format("Wminus" if charge=="WMinus" else "Wplus"))

hists_aMC = []
for i in range(5):
    tmp = hist2array(fcoeffs.Get('angularCoefficients/harmonicsA{}_nom_nom'.format(i)))
    hists_aMC.append(tmp[:threshold_y,:threshold_qt])
h = np.stack(hists_aMC,axis=-1)

htot = hist2array(f_aMC.Get('angularCoefficients_{}/mapTot'.format("Wminus" if charge=="WMinus" else "Wplus")))
htot_qt = np.sum(htot[:threshold_y,:threshold_qt],axis=0)
htot_y = np.sum(htot[:threshold_y,:threshold_qt],axis=1)

# integrated coeffs
hists_aMC = []
for i in range(5):
    tmp = hist2array(fcoeffs.Get('angularCoefficients/harmonicsPtA{}_nom_nom'.format(i)))
    hists_aMC.append(tmp[:threshold_qt])
h_qt = np.stack(hists_aMC,axis=-1)

# integrated coeffs
hists_aMC = []
for i in range(5):
    tmp = hist2array(fcoeffs.Get('angularCoefficients/harmonicsYA{}_nom_nom'.format(i)))
    hists_aMC.append(tmp[:threshold_y])
h_y = np.stack(hists_aMC,axis=-1)


# file with nominal coefficients
file_preVFP_powheg = '/scratchnvme/wmass/REWEIGHT/outputW_sroychow_preVFP/WPlusJetsToMuNu_helweights.hdf5'
file_postVFP_powheg = '/scratchnvme/wmass/REWEIGHT/outputW_sroychow_postVFP/WPlusJetsToMuNu_helweights.hdf5'

# pdfs
f_preVFP_powheg = h5py.File(file_preVFP_powheg, mode='r+')
f_postVFP_powheg = h5py.File(file_postVFP_powheg, mode='r+')
htot_preVFP_powheg = f_preVFP_powheg['totxsecs_LHEPdfWeight'][:]
htot_postVFP_powheg = f_postVFP_powheg['totxsecs_LHEPdfWeight'][:]
h_preVFP_powheg = f_preVFP_powheg['xsecs_LHEPdfWeight'][:]
h_postVFP_powheg = f_postVFP_powheg['xsecs_LHEPdfWeight'][:]

htot_powheg = htot_preVFP_powheg+htot_postVFP_powheg
h_powheg = h_preVFP_powheg+h_postVFP_powheg

factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
factors = factors[np.newaxis,np.newaxis,...]
factors = factors[...,np.newaxis]
h_powheg = h_powheg.reshape(len(yBins)-1, len(qtBins)-1, 9, 103)
h_powheg = (h_powheg/htot_powheg[:,:,np.newaxis,...]+factors[:,:,:,1,...])*factors[:,:,:,0,...]
# shape h: y, qt, weights, pdf
# shape tot: y, qt, pdf

h_powheg = h_powheg[:threshold_y,:threshold_qt,:,:]
h_powheg[...,-1,:] = htot_powheg[:threshold_y,:threshold_qt,:]

print('MC@NLO', h[2,5,3])
print('powheg', h_powheg[...,3,0])

# integrated coeffs qt
htot_qt_powheg = np.sum(htot_powheg[:threshold_y,:threshold_qt,:],axis=0)
h_qt_powheg = np.sum((h_preVFP_powheg[:threshold_y,:threshold_qt,:]+h_postVFP_powheg[:threshold_y,:threshold_qt,:]).reshape(len(yBins[:threshold_y+1])-1, len(qtBins[:threshold_qt+1])-1, 9, 103),axis=0)
factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
factors = factors[np.newaxis,...]
factors = factors[...,np.newaxis]

h_qt_powheg = (h_qt_powheg/htot_qt_powheg[...,np.newaxis,:]+factors[...,1,:])*factors[...,0,:]

# integrated coeffs y
htot_y_powheg = np.sum(htot_powheg[:threshold_y,:threshold_qt,:],axis=1)
h_y_powheg = np.sum((h_preVFP_powheg[:threshold_y,:threshold_qt,:]+h_postVFP_powheg[:threshold_y,:threshold_qt,:]).reshape(len(yBins[:threshold_y+1])-1, len(qtBins[:threshold_qt+1])-1, 9, 103),axis=1)
factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
factors = factors[np.newaxis,...]
factors = factors[...,np.newaxis]

h_y_powheg = (h_y_powheg/htot_y_powheg[...,np.newaxis,:]+factors[...,1,:])*factors[...,0,:]

# scales
f_preVFP_powheg = h5py.File(file_preVFP_powheg, mode='r+')
f_postVFP_powheg = h5py.File(file_postVFP_powheg, mode='r+')
htot_preVFP_powheg_scales = f_preVFP_powheg['totxsecs_LHEScaleWeight'][:]
htot_postVFP_powheg_scales = f_postVFP_powheg['totxsecs_LHEScaleWeight'][:]
h_preVFP_powheg_scales = f_preVFP_powheg['xsecs_LHEScaleWeight'][:]
h_postVFP_powheg_scales = f_postVFP_powheg['xsecs_LHEScaleWeight'][:]

htot_powheg_scales = htot_preVFP_powheg_scales+htot_postVFP_powheg_scales
h_powheg_scales = h_preVFP_powheg_scales+h_postVFP_powheg_scales

factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
factors = factors[np.newaxis,np.newaxis,...]
factors = factors[...,np.newaxis]
h_powheg_scales = h_powheg_scales.reshape(len(yBins)-1, len(qtBins)-1, 9, 9)
h_powheg_scales = (h_powheg_scales/htot_powheg_scales[:,:,np.newaxis,...]+factors[:,:,:,1,...])*factors[:,:,:,0,...]
# shape h: y, qt, weights, pdf
# shape tot: y, qt, pdf

h_powheg_scales = h_powheg_scales[:threshold_y,:threshold_qt,:,:]
h_powheg_scales[...,-1,:] = htot_powheg_scales[:threshold_y,:threshold_qt,:]

qcdsyst = [0, 1, 3, 5, 7, 8]
# qcdsyst = [1, 5, 7, 8]
h_powheg_scales = h_powheg_scales[...,qcdsyst]
htot_powheg_scales = htot_powheg_scales[:threshold_y,:threshold_qt,qcdsyst]
h_powheg_scales_vars = np.zeros((htot_powheg_scales.shape[0],htot_powheg_scales.shape[1],len(helicities)-1,36))
print(h_powheg_scales_vars.shape,htot_powheg_scales.shape,h_powheg_scales.shape)
# for i in range(len(helicities)-1):
#     for k in range(htot_powheg_scales.shape[0]):
#         for j in range(htot_powheg_scales.shape[1]):
#             h_powheg_scales_vars[k,j,i,:] = ((np.outer(h_powheg_scales[k,j,i,:],np.reciprocal(htot_powheg_scales[k,j,:])).ravel()+factors[...,i,1,:])*factors[...,i,0,:])
factors_hel = factors_hel[...,np.newaxis]
hhel_powheg_scales = 3./(16.*math.pi)*h_powheg_scales*htot_powheg_scales[...,np.newaxis,:]/factors_hel

# # integrated coeffs
# htot_qt_powheg_scales = np.sum(htot_powheg_scales[:threshold_y,:threshold_qt,:],axis=0)
# h_qt_powheg_scales = np.sum((h_preVFP_powheg_scales[:threshold_y,:threshold_qt,:]+h_postVFP_powheg_scales[:threshold_y,:threshold_qt,:]).reshape(len(yBins[:threshold_y+1])-1, len(qtBins[:threshold_qt+1])-1, 9, 9),axis=0)
# factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])

# print(htot_qt_powheg_scales.shape)

# h_qt_powheg_scales_vars = np.zeros((htot_qt_powheg_scales.shape[0],len(helicities)-1,36))
# for i in range(len(helicities)-1):
#     for j in range(htot_qt_powheg_scales.shape[0]):
#         h_qt_powheg_scales_vars[j,i] = ((np.outer(h_qt_powheg_scales[j,i,:],np.reciprocal(htot_qt_powheg_scales[j,:])).ravel()+factors[i,1])*factors[i,0])

# # y-integrated coeffs
# htot_y_powheg_scales = np.sum(htot_powheg_scales[:threshold_y,:threshold_qt,:],axis=1)
# h_y_powheg_scales = np.sum((h_preVFP_powheg_scales[:threshold_y,:threshold_qt,:]+h_postVFP_powheg_scales[:threshold_y,:threshold_qt,:]).reshape(len(yBins[:threshold_y+1])-1, len(qtBins[:threshold_qt+1])-1, 9, 9),axis=1)
# factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])

# print(htot_y_powheg_scales.shape)

# h_y_powheg_scales = h_y_powheg_scales[...,qcdsyst]
# h_y_powheg_scales_vars = np.zeros((htot_y_powheg_scales.shape[0],len(helicities)-1,36))
# for i in range(len(helicities)-1):
#     for j in range(htot_y_powheg_scales.shape[0]):
#         h_y_powheg_scales_vars[j,i] = ((np.outer(h_y_powheg_scales[j,i,:],np.reciprocal(htot_y_powheg_scales[j,:])).ravel()+factors[i,1])*factors[i,0])

# integrated coeffs qt
htot_qt_powheg_scales = np.sum(htot_powheg_scales[:threshold_y,:threshold_qt,:],axis=0)
h_qt_powheg_scales = np.sum((h_preVFP_powheg_scales[:threshold_y,:threshold_qt,:]+h_postVFP_powheg_scales[:threshold_y,:threshold_qt,:]).reshape(len(yBins[:threshold_y+1])-1, len(qtBins[:threshold_qt+1])-1, 9, 9),axis=0)[...,qcdsyst]
factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
factors = factors[np.newaxis,...]
factors = factors[...,np.newaxis]

h_qt_powheg_scales = (h_qt_powheg_scales/htot_qt_powheg_scales[...,np.newaxis,:]+factors[...,1,:])*factors[...,0,:]

# integrated coeffs y
htot_y_powheg_scales = np.sum(htot_powheg_scales[:threshold_y,:threshold_qt,:],axis=1)
h_y_powheg_scales = np.sum((h_preVFP_powheg_scales[:threshold_y,:threshold_qt,:]+h_postVFP_powheg_scales[:threshold_y,:threshold_qt,:]).reshape(len(yBins[:threshold_y+1])-1, len(qtBins[:threshold_qt+1])-1, 9, 9),axis=1)[...,qcdsyst]
factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
factors = factors[np.newaxis,...]
factors = factors[...,np.newaxis]

h_y_powheg_scales = (h_y_powheg_scales/htot_y_powheg_scales[...,np.newaxis,:]+factors[...,1,:])*factors[...,0,:]

yBins = np.array(yBins[:threshold_y+1],dtype='float64')
qtBins = np.array(qtBins[:threshold_qt+1],dtype='float64')
print(yBins,qtBins)
yBinsC = 0.5*(yBins[1:]+yBins[:-1])
qtBinsC = 0.5*(qtBins[1:]+qtBins[:-1])
yBinsS = yBins[1:]-yBins[:-1]
qtBinsS = qtBins[1:]-qtBins[:-1]

# for k,c in enumerate(coefficients):
#     hcoeff = np.zeros([len(yBinsC),len(qtBinsC)])
#     hcoeff_gen = np.zeros([len(yBinsC),len(qtBinsC)])
#     hcoeff_err = np.zeros([len(yBinsC),len(qtBinsC)])
#     hcoeff_hel = np.zeros([len(yBinsC),len(qtBinsC)])
#     hcoeff_hel_gen = np.zeros([len(yBinsC),len(qtBinsC)])
#     hcoeff_hel_err = np.zeros([len(yBinsC),len(qtBinsC)])
#     bins=np.append(np.tile(qtBins[:-1],len(yBinsC)),60.)
#     x=np.array(range(len(bins)))
#     binsC = 0.5*(x[1:]+x[:-1])
#     print('analysing', c, helicities[k])
#     for ev in fitresults: #dummy because there's one event only
#         for i in range(len(yBinsC)):
#             for j in range(len(qtBinsC)):
#                 try:
#                     coeff = eval('ev.y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i))
#                     coeff_gen = eval('ev.y_{i}_qt_{j}_{c}_gen'.format(c=c, j=j, i=i))
#                     coeff_err = eval('ev.y_{i}_qt_{j}_{c}_err'.format(c=c, j=j, i=i))
#                     coeff_hel = eval('ev.helXsecs{c}_y_{i}_qt_{j}_pmaskedexp'.format(c=helicities[k], j=j, i=i))
#                     coeff_hel_err = eval('ev.helXsecs{c}_y_{i}_qt_{j}_pmaskedexp_err'.format(c=helicities[k], j=j, i=i))
#                     coeff_hel_gen = eval('ev.helXsecs{c}_y_{i}_qt_{j}_pmaskedexp_gen'.format(c=helicities[k], j=j, i=i))
#                     if 'unpolarizedxsec' in c:
#                         coeff = coeff/(3./16./math.pi)/35.9/yBinsS[i]/qtBinsS[j]
#                         coeff_gen = coeff_gen/(3./16./math.pi)/35.9/yBinsS[i]/qtBinsS[j]
#                         coeff_err = coeff_err/(3./16./math.pi)/35.9/yBinsS[i]/qtBinsS[j]
#                     hcoeff[i,j]=coeff
#                     hcoeff_gen[i,j]=coeff_gen
#                     hcoeff_err[i,j]=coeff_err
#                     hcoeff_hel[i,j]=coeff_hel
#                     hcoeff_hel_gen[i,j]=coeff_hel_gen
#                     hcoeff_hel_err[i,j]=coeff_hel_err
#                 except AttributeError:
#                     pass
#     fig, (ax1, ax2) = plt.subplots(nrows=2,figsize=(12, 10),gridspec_kw={'height_ratios': [3, 1]})
#     # fig, ax1 = plt.subplots(figsize=(48, 10))
#     hep.cms.text('work in progress', loc=1, ax=ax1)
#     if 'unpol' in c:
#         ax1.set_ylabel(r'$\frac{d\sigma^{U+L}}{dq_Td|y|} (fb/GeV)$', fontsize=30)
#     else:
#         ax1.set_ylabel(r'$A_{}$'.format(k-1))
#     # ax1.set_xlabel('$q_T$ (GeV)')
#     # ax1.set_xticks(x) # set tick positions
#     # ticks = ["{:d}".format(int(v)) for v in bins]
#     # for itick,tick in enumerate(ticks):
#     #     if not itick==1:
#     #         if tick=="0": ticks[itick] = "60_0"
#     # ax1.set_xticklabels(ticks)
#     # print(hcoeff.ravel(), "+/-", hcoeff_err.ravel())
#     hep.histplot(hcoeff.ravel(),bins = x, yerr = hcoeff_err.ravel(),histtype = 'errorbar', color = "k", stack = False, ax=ax1, label=label_data)
#     hep.histplot(hcoeff_gen.ravel(),bins =x, color = "r", stack = False, ax=ax1, label="prediction")

#     if 'unpol' in c:
#         # hep.histplot(htot_reshaped/35.9/yBinsS[:,np.newaxis,np.newaxis]/qtBinsS[np.newaxis,:,np.newaxis].ravel(),bins = x, color = "b", stack = False, ax=ax1)
#         if not asimov:
#             hep.histplot((htot[...,:threshold_y,:threshold_qt]/35.9/yBinsS[:,np.newaxis]/qtBinsS[np.newaxis,:]).ravel(),bins = x, color = "blue", stack = False, ax=ax1, label="aMC@NLO")

#         # ax2.fill_between(binsC, ratio-(hcoeff_err.ravel()/htot_reshaped_new.ravel()), ratio+(hcoeff_err.ravel()/htot_reshaped_new.ravel()), alpha=0.3, color="orange")

#         err_pdf = np.sqrt(np.sum(np.square(h_powheg[...,-1,0][...,np.newaxis]-h_powheg[...,-1,1:]),axis=-1))/35.9/yBinsS[:,np.newaxis]/qtBinsS[np.newaxis,:]
#         err_scales_up = np.abs(h_powheg_scales[...,-1,:].max(axis=-1))/35.9/yBinsS[:,np.newaxis]/qtBinsS[np.newaxis,:]
#         err_scales_down = np.abs(h_powheg_scales[...,-1,:].min(axis=-1))/35.9/yBinsS[:,np.newaxis]/qtBinsS[np.newaxis,:]
#         # ax1.fill_between(binsC, (err_scales_down).ravel(), (err_scales_up).ravel(), alpha=0.2, color="green",label="qcd scale uncertainty")
#         # ax1.fill_between(binsC, (hcoeff_gen-err_pdf).ravel(), (hcoeff_gen+err_pdf).ravel(), alpha=1, color="orange",label="pdf uncertainty")
#         ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#         # ratio = hcoeff.ravel()/hcoeff_gen.ravel()
#         # ratio_err= hcoeff_err.ravel()/hcoeff_gen.ravel()
#         ratio = hcoeff.ravel()/(htot[...,:threshold_y,:threshold_qt]/35.9/yBinsS[:,np.newaxis]/qtBinsS[np.newaxis,:]).ravel()
#         ratio_err= hcoeff_err.ravel()/(htot[...,:threshold_y,:threshold_qt]/35.9/yBinsS[:,np.newaxis]/qtBinsS[np.newaxis,:]).ravel()
#         ax2.set_ylabel('data/aMC@NLO')
#         ax2.set_ylim([0, 2])
#     else:
#         if not asimov:
#             hep.histplot(h[...,:threshold_y,:threshold_qt,k-1].ravel(),bins = x, color = "blue", stack = False, ax=ax1,label="aMC@NLO")
#         # hep.histplot(h_powheg[...,k-1,0].ravel(),bins =x, stack = False, ax=ax1, label="prediction 2")
#         err_pdf = np.sqrt(np.sum(np.square(h_powheg[...,k-1,0][...,np.newaxis]-h_powheg[...,k-1,1:]),axis=-1))
#         err_scales_up = np.abs(h_powheg[...,k-1,0]-h_powheg_scales_vars[...,k-1,:].max(axis=-1))
#         err_scales_down = np.abs(h_powheg[...,k-1,0]-h_powheg_scales_vars[...,k-1,:].min(axis=-1))
#         # ax1.fill_between(binsC, (hcoeff_gen-err_scales_down).ravel(), (hcoeff_gen+err_scales_up).ravel(), alpha=0.2, color="green",label="qcd scale uncertainty")
#         # ax1.fill_between(binsC, (hcoeff_gen-err_pdf).ravel(), (hcoeff_gen+err_pdf).ravel(), alpha=1, color="orange",label="pdf uncertainty")
#         ratio = hcoeff.ravel()-(hcoeff_gen.ravel())
#         ratio_err = hcoeff_err.ravel()
#         ax1.set_ylim(-5,5)
#         ax2.set_ylabel('data-prediction')
#         ax2.set_ylim([-2, 2])
#     hep.histplot(ratio,bins = x, yerr=ratio_err.ravel(), histtype = 'errorbar', color = "k", stack = False, ax=ax2)
#     ax2.set_xlabel('unrolled $q_T$-y bins')
#     ax1.legend(loc='upper right', frameon=False)
#     plt.tight_layout()
#     plt.savefig('POIplots/fit{}.png'.format(c),dpi=300)
#     plt.savefig('POIplots/fit{}.pdf'.format(c),dpi=300)
#     plt.clf()
#     # plot helicity cross sections
#     fig, (ax1, ax2) = plt.subplots(nrows=2,figsize=(48, 10),gridspec_kw={'height_ratios': [3, 1]})
#     ax1.set_title("fitted {}".format(helicities[k]), fontsize=18)
#     ax1.set_ylabel('')
#     ax1.set_xlabel('a.u.')
#     ax1.set_xticks(x) # set tick positions
#     # ax1.set_xticklabels(ticks)
#     # print(hcoeff.ravel(), "+/-", hcoeff_err.ravel())
#     # print(hcoeff_hel.ravel(), "+/-", hcoeff_hel_err.ravel())
    
#     hep.histplot(hcoeff_hel.ravel(),bins = x, yerr = hcoeff_hel_err.ravel(),histtype = 'errorbar', color = "k", stack = False, ax=ax1, label=label_data)
#     hep.histplot(hcoeff_hel_gen.ravel(),bins =x, color = "r", stack = False, ax=ax1, label="prediction")
#     # if not asimov:
#     #     hep.histplot(hhel[...,:threshold_y,:threshold_qt,k-1].ravel(),bins = x, color = "blue", stack = False, ax=ax1,label="aMC@NLO")
    
#     err_scales_up = np.abs(hcoeff_hel_gen-hhel_powheg_scales[...,k-1,:].max(axis=-1))
#     err_scales_down = np.abs(hcoeff_hel_gen-hhel_powheg_scales[...,k-1,:].min(axis=-1))
#     # ax1.fill_between(binsC, hhel_powheg_scales[...,k-1,:].max(axis=-1).ravel(), hhel_powheg_scales[...,k-1,:].min(axis=-1).ravel(), alpha=0.2, color="green",label="qcd scale uncertainty")
#     # ax1.fill_between(binsC, 10*hhel_powheg_scales[...,k-1,:].max(axis=-1).ravel(), 1./10*hhel_powheg_scales[...,k-1,:].min(axis=-1).ravel(), alpha=0.2, color="magenta",label="qcd scale uncertainty x10")

#     ax1.legend(loc='upper right', frameon=False)
#     plt.tight_layout()
#     plt.savefig('POIplots/fit{}.png'.format(helicities[k]),dpi=300)
#     plt.savefig('POIplots/fit{}.pdf'.format(helicities[k]),dpi=300)
#     plt.clf()

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
                coeff = eval('ev.qt_{j}_helmeta_{c}'.format(c=c, j=j, i=i))
                coeff_gen = eval('ev.qt_{j}_helmeta_{c}_gen'.format(c=c, j=j, i=i))
                coeff_err = eval('ev.qt_{j}_helmeta_{c}_err'.format(c=c, j=j, i=i))
                coeff_hel = eval('ev.helXsecs{c}_qt_{j}_sumxsec'.format(c=helicities[k], j=j, i=i))
                coeff_hel_gen = eval('ev.helXsecs{c}_qt_{j}_sumxsec_gen'.format(c=helicities[k], j=j, i=i))
                coeff_hel_err = eval('ev.helXsecs{c}_qt_{j}_sumxsec_err'.format(c=helicities[k], j=j, i=i))
                if 'unpol' in c:
                    coeff = coeff/(3./16./math.pi)/35.9/qtBinsS[j]
                    coeff_gen = coeff_gen/(3./16./math.pi)/35.9/qtBinsS[j]
                    coeff_err = coeff_err/(3./16./math.pi)/35.9/qtBinsS[j]
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
                    coeff_y = coeff_y/(3./16./math.pi)/35.9/yBinsS[i]
                    coeff_gen_y = coeff_gen_y/(3./16./math.pi)/35.9/yBinsS[i]
                    coeff_err_y = coeff_err_y/(3./16./math.pi)/35.9/yBinsS[i]
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
    mainplot=ax1.errorbar(qtBinsC,hcoeff.ravel(), xerr=qtBinsS/2, yerr = hcoeff_err.ravel(),marker = 'o',color = "k", label=label_data, fmt='o')

    if 'unpol' in c:
        # if not asimov:
        #     toy,=ax1.plot(qtBinsC[:threshold_qt+1], htot_qt/35.9/qtBinsS.ravel(), color = "blue", label="aMC@NLO")

        err_scales_up = np.abs(htot_qt_powheg_scales[:,:].max(axis=-1)/35.9/qtBinsS-hcoeff_gen.ravel())
        err_scales_down = np.abs(htot_qt_powheg_scales[:,:].min(axis=-1)/35.9/qtBinsS-hcoeff_gen.ravel())

        err_pdf = np.sqrt(np.sum(np.square(htot_qt_powheg[...,0][:,np.newaxis]-htot_qt_powheg[...,1:]),axis=-1))/35.9/qtBinsS
        errUp = np.hypot(err_pdf, err_scales_up)
        errDown = np.hypot(err_pdf, err_scales_down)

        _ = make_error_boxes(ax=ax1, xdata=qtBinsC, ydata=hcoeff_gen.ravel(), xerror=np.array([qtBinsS/2,qtBinsS/2]), yerror=np.array([errDown,errUp]))
        patch = mpatches.Patch(color='green', alpha=0.5,label="powheg MiNNLO")
        ax1.legend(handles=[patch,mainplot],loc='upper right',frameon=False)
        ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        ax1.set_ylim(0,5e5)
        ratio = hcoeff.ravel()/(hcoeff_gen.ravel())
        err_ratio = hcoeff_err.ravel()/hcoeff_gen.ravel()
        # ratio = hcoeff.ravel()/(htot_qt/35.9/qtBinsS).ravel()
        # err_ratio= hcoeff_err.ravel()/(htot_qt/35.9/qtBinsS).ravel()
        # ax2.set_ylabel('data/aMC@NLO')
        ax2.set_ylabel('data/prediction')
    else:
        if not asimov:
            toy, =ax1.plot(qtBinsC[:threshold_qt+1],h_qt[...,k-1].ravel(), color = "blue",label="aMC@NLO")
        
        err_scales_up = np.abs(h_qt_powheg[...,k-1,0]-h_qt_powheg_scales[:,k-1,:].max(axis=-1))
        err_scales_down = np.abs(h_qt_powheg[...,k-1,0]-h_qt_powheg_scales[:,k-1,:].min(axis=-1))

        err_pdf = np.sqrt(np.sum(np.square(h_qt_powheg[...,k-1,0][:,np.newaxis]-h_qt_powheg[...,k-1,1:]),axis=-1))
        errUp = np.hypot(err_pdf, err_scales_up)
        errDown = np.hypot(err_pdf, err_scales_down)

        _ = make_error_boxes(ax=ax1, xdata=qtBinsC, ydata=hcoeff_gen.ravel(), xerror=np.array([qtBinsS/2,qtBinsS/2]), yerror=np.array([errDown,errUp]))
        patch = mpatches.Patch(color='green', alpha=0.5,label="powheg MiNNLO")
        ax1.legend(handles=[patch,mainplot],loc='upper right',frameon=False)
        ratio = hcoeff.ravel()-(hcoeff_gen.ravel())
        err_ratio = hcoeff_err.ravel()
        ax2.set_ylabel('data-prediction')

        invcov = np.linalg.inv(cov[k*8+(k+1)*6:(k+1)*8+(k+1)*6,k*8+(k+1)*6:(k+1)*8+(k+1)*6])
        print(hcov.GetXaxis().GetBinLabel(k*8+(k+1)*6+1),hcov.GetXaxis().GetBinLabel((k+1)*8+(k+1)*6))
        chi2 = np.dot(ratio.T, np.dot(invcov,ratio))
        p_value = 1 - stats.chi2.cdf(chi2, 8)
        print(chi2,p_value)
        # ax2.text(0.95, 0.95, "chi2/ndof: {:.2f}".format(chi2/8),verticalalignment='top', horizontalalignment='right',transform=ax2.transAxes,color='black')

    # print(qtBinsC.shape, (ratio-hcoeff_err.ravel()).shape)
    if not 'unpol' in c:
        ax1.set_ylim([-2, 2])
        ax2.set_ylim([-2, 2])
        _ = make_error_boxes(ax=ax2, xdata=qtBinsC, ydata=np.zeros_like(ratio), xerror=np.array([qtBinsS/2,qtBinsS/2]), yerror=np.array([errDown,errUp]), label="powheg MiNNLO")
    else:
        _ = make_error_boxes(ax=ax2, xdata=qtBinsC, ydata=np.ones_like(ratio), xerror=np.array([qtBinsS/2,qtBinsS/2]), yerror=np.array([errDown/hcoeff_gen.ravel(),errUp/hcoeff_gen.ravel()]), label="powheg MiNNLO")
    ax2.errorbar(qtBinsC,ratio,marker = 'o', color = "k", fmt='o', xerr=qtBinsS/2,yerr=err_ratio)
    plt.tight_layout()
    plt.savefig('POIplots/fitintegrated{}.png'.format(c),dpi=300)
    plt.savefig('POIplots/fitintegrated{}.pdf'.format(c),dpi=300)
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
    mainplot=ax1.errorbar(yBinsC,hcoeff_y.ravel(), xerr=yBinsS/2, yerr = hcoeff_err_y.ravel(),marker = 'o',color = "k", label=label_data, fmt='o')

    if 'unpol' in c:
        # if not asimov:
        #     toy,=ax1.plot(yBinsC[:threshold_y+1],htot_y/35.9/yBinsS.ravel(), color = "blue", label="aMC@NLO")

        err_pdf = np.sqrt(np.sum(np.square(htot_y_powheg[...,0][:,np.newaxis]-htot_y_powheg[...,1:]),axis=-1))/35.9/yBinsS
        err_scales_up = np.abs(htot_y_powheg_scales[:,:].max(axis=-1)/35.9/yBinsS-hcoeff_gen_y.ravel())
        err_scales_down = np.abs(htot_y_powheg_scales[:,:].min(axis=-1)/35.9/yBinsS-hcoeff_gen_y.ravel())

        err_pdf = np.sqrt(np.sum(np.square(htot_y_powheg[...,0][:,np.newaxis]-htot_y_powheg[...,1:]),axis=-1))/35.9/yBinsS
        errUp = np.hypot(err_pdf, err_scales_up)
        errDown = np.hypot(err_pdf, err_scales_down)

        _ = make_error_boxes(ax=ax1, xdata=yBinsC, ydata=hcoeff_gen_y.ravel(), xerror=np.array([yBinsS/2,yBinsS/2]), yerror=np.array([errDown,errUp]))
        patch = mpatches.Patch(color='green', alpha=0.5,label="powheg MiNNLO")
        ax1.legend(handles=[patch,mainplot],loc='upper right',frameon=False)
        # for isyst in range(htot_y_powheg_scales.shape[-1]):
        #     hep.histplot(htot_y_powheg_scales[...,isyst]/35.9/yBinsS.ravel(),bins = yBins[:threshold_y+1],ax=ax1)

        ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        ratio = hcoeff_y.ravel()/(hcoeff_gen_y.ravel())
        err_ratio = hcoeff_err_y.ravel()/hcoeff_gen_y.ravel()
        # ratio = hcoeff_y.ravel()/(htot_y/35.9/yBinsS).ravel()
        # err_ratio= hcoeff_err_y.ravel()/(htot_y/35.9/yBinsS).ravel()
        # ax2.set_ylabel('data/aMC@NLO')

    else:
        if not asimov:
            hep.histplot(h_y[...,k-1].ravel(),bins = yBins[:threshold_y+1], color = "blue", stack = False, ax=ax1,label="aMC@NLO")
        err_scales_up = np.abs(h_y_powheg[...,k-1,0]-h_y_powheg_scales[:,k-1,:].max(axis=-1))
        err_scales_down = np.abs(h_y_powheg[...,k-1,0]-h_y_powheg_scales[:,k-1,:].min(axis=-1))

        err_pdf = np.sqrt(np.sum(np.square(h_y_powheg[...,k-1,0][:,np.newaxis]-h_y_powheg[...,k-1,1:]),axis=-1))
        errUp = np.hypot(err_pdf, err_scales_up)
        errDown = np.hypot(err_pdf, err_scales_down)
        # for isyst in range(h_y_powheg_scales.shape[-1]):
        #     hep.histplot(h_y_powheg_scales[...,k-1,isyst].ravel(),bins = yBins[:threshold_y+1],ax=ax1)
        _ = make_error_boxes(ax=ax1, xdata=yBinsC, ydata=hcoeff_gen_y.ravel(), xerror=np.array([yBinsS/2,yBinsS/2]), yerror=np.array([errDown,errUp]))
        patch = mpatches.Patch(color='green', alpha=0.5,label="powheg MiNNLO")
        ax1.legend(handles=[patch,mainplot],loc='upper right',frameon=False)
        ratio = hcoeff_y.ravel()-(hcoeff_gen_y.ravel())
        err_ratio = hcoeff_err_y.ravel()
        ax2.set_ylabel('data-prediction')

        invcov = np.linalg.inv(cov[k*6+k*8:(k+1)*6+k*8,k*6+k*8:(k+1)*6+k*8])
        print(hcov.GetXaxis().GetBinLabel(k*6+k*8+1),hcov.GetXaxis().GetBinLabel((k+1)*6+k*8))
        chi2 = np.dot(ratio.T, np.dot(invcov,ratio))
        p_value = 1 - stats.chi2.cdf(chi2, 6)
        print(chi2,p_value)
        # ax2.text(0.95, 0.95, "chi2/ndof: {:.2f}".format(chi2/6),verticalalignment='top', horizontalalignment='right',transform=ax2.transAxes,color='black')
    ax2.errorbar(yBinsC,ratio,marker = 'o', color = "k", fmt='o', xerr=yBinsS/2,yerr=err_ratio)
    if not 'unpol' in c:
        ax1.set_ylim([-2, 2])
        ax2.set_ylim([-2, 2])
        _ = make_error_boxes(ax=ax2, xdata=yBinsC, ydata=np.zeros_like(ratio), xerror=np.array([yBinsS/2,yBinsS/2]), yerror=np.array([errDown,errUp]), label="powheg MiNNLO")
    else:
        _ = make_error_boxes(ax=ax2, xdata=yBinsC, ydata=np.ones_like(ratio), xerror=np.array([yBinsS/2,yBinsS/2]), yerror=np.array([errDown/hcoeff_gen_y.ravel(),errUp/hcoeff_gen_y.ravel()]), label="powheg MiNNLO")
        ax1.set_ylim([0, 4e6])
        pass
    plt.tight_layout()
    plt.savefig('POIplots/fityintegrated{}.png'.format(c),dpi=300)
    plt.savefig('POIplots/fityintegrated{}.pdf'.format(c),dpi=300)
    plt.clf()
