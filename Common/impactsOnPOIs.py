import ROOT
import os
import sys
import h5py
import math
import numpy as np
# import jax
# import jax.numpy as jnp
import matplotlib
# Set the backend to Agg
matplotlib.use('Agg')

import matplotlib.pyplot as plt
# import mplhep as hep
sys.path.append('data/')
from binning import yBins, qtBins, ptBins, etaBins, mTBins, isoBins, chargeBins
# import narf
import argparse

# matplotlib stuff
# plt.style.use([hep.style.CMS])

parser = argparse.ArgumentParser('')
parser.add_argument('-asimov', '--asimov', default=False, action='store_true', help='plot asimov result')

args = parser.parse_args()
asimov = args.asimov

helicities = ['UL','L', 'I', 'T', 'A', 'P']
coefficients = ['unpolarizedxsec', "A0", "A1", "A2", "A3", "A4"]

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
yBins = np.array(yBins[:threshold_y+1],dtype='float64')
qtBins = np.array(qtBins[:threshold_qt+1],dtype='float64')
yBinsC = 0.5*(yBins[1:]+yBins[:-1])
qtBinsC = 0.5*(qtBins[1:]+qtBins[:-1])
yBinsS = yBins[1:]-yBins[:-1]
qtBinsS = qtBins[1:]-qtBins[:-1]

# file with fit results
f = ROOT.TFile.Open('../Fit/FitRes/fit_Wplus_asimov_BBB.root'.format('asimov' if asimov else 'data'))
fitresults = f.Get('fitresults')

# get value of fitted parameters
for k,c in enumerate(coefficients):
    if not 'unpolarizedxsec' in c: continue
    hcoeff = np.zeros([len(yBinsC),len(qtBinsC)])
    hcoeff_gen = np.zeros([len(yBinsC),len(qtBinsC)])
    hcoeff_err = np.zeros([len(yBinsC),len(qtBinsC)])
    print('analysing', c, helicities[k])
    for ev in fitresults: #dummy because there's one event only
        for i in range(len(yBinsC)):
            for j in range(len(qtBinsC)):
                try:
                    coeff = getattr(ev,'y_{:.1f}_qt_{:.1f}_WplusmunuPostVFP_{}'.format(yBinsC[i] , qtBinsC[j] , c))
                    coeff_gen = getattr(ev,'y_{:.1f}_qt_{:.1f}_WplusmunuPostVFP_{}_gen'.format(yBinsC[i] , qtBinsC[j] , c))
                    coeff_err = getattr(ev,'y_{:.1f}_qt_{:.1f}_WplusmunuPostVFP_{}_err'.format(yBinsC[i] , qtBinsC[j] , c))
                    if 'unpolarizedxsec' in c:
                        coeff = coeff
                        coeff_gen = coeff_gen
                        coeff_err = coeff_err
                    hcoeff[i,j]=coeff
                    hcoeff_gen[i,j]=coeff_gen
                    hcoeff_err[i,j]=coeff_err
                except AttributeError:
                    pass

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
yBins_red = np.array(yBins[:threshold_y+1])
qtBins_red = np.array(qtBins[:threshold_qt+1])
yBinsC = 0.5*(yBins_red[1:]+yBins_red[:-1])
qtBinsC = 0.5*(qtBins_red[1:]+qtBins_red[:-1])

nBins=len(yBinsC)*len(qtBinsC)
bins = np.linspace(0,nBins+1,nBins+1)
binsC = 0.5*(bins[1:]+bins[:-1])
helXsecs = ['helXsec_L', 'helXsec_I', 'helXsec_T', 'helXsec_A', 'helXsec_P','helXsec_UL']

processes =[]
for hel in helXsecs:
    for i in range(len(yBinsC)):
        for j in range(len(qtBinsC)):
            proc = 'helXsecs' + hel + '_y_{}'.format(i)+'_qt_{}'.format(j)
            processes.append(proc)

# hcov = f.Get('covariance_matrix_channelhelpois')
# cov = hist2array(hcov)[:,:]
# print(cov.shape)

impacts = []
# retrieve impacts per group of POI

# for i in range(1,6):
#     for j in range(nBins):
#         #impact is generalization of per-nuisance impacts above v^T C^-1 v
#         #where v is the matrix of poi x nuisance correlations within the group
#         #and C is is the subset of the covariance matrix corresponding to the nuisances in the group
#         vcovreduced = cov[j,i*nBins:(i+1)*nBins]
#         groupmcov = np.linalg.inv(cov[i*nBins:(i+1)*nBins,i*nBins:(i+1)*nBins])
#         # print(hcov.GetXaxis().GetBinLabel(i*nBins+1),hcov.GetXaxis().GetBinLabel((i+1)*nBins+1))
#         vimpact = np.sqrt(np.matmul(np.matmul(vcovreduced.T,groupmcov),vcovreduced))
#         impacts.append(vimpact)
# impacts =np.stack(impacts).reshape(nBins,5)/hcoeff.ravel()[:,np.newaxis]
# print(impacts.shape)
# # impacts per POI
# plt.rcdefaults()
# fig, ax = plt.subplots(figsize=(5, 5))
# hep.cms.text('work in progress', loc=0)
# for i in range(impacts.shape[1]):
#     ax.plot(binsC, impacts[:,i], label=coefficients[i+1])
# for j in range(0, len(yBinsC)+1):
#     plt.axvline(x=j*len(qtBinsC)+binsC[0],linewidth=1,linestyle='dashed')
# ax.legend(bbox_to_anchor=(1.04,1), borderaxespad=0, frameon=False)
# ax.set_yscale('log')
# ax.set_ylabel('relative impact on $\sigma^{U+L}$')
# ax.set_xlabel('unrolled $q_T$-y bins')
# plt.tight_layout()
# plt.savefig('impactsOnPOIs_other_{}.pdf'.format("asimov" if asimov else "data"),dpi=300)

hcov = f.Get('correlation_matrix_channelmu')
# cov = narf.root_to_hist(hcov).values()
cov = np.asarray(hcov).reshape(289+2,289+2)[1:-1,1:-1]
hel_idx = {}
for hel in helXsecs:
    hel_idx[hel] = []
    for i in range(1,hcov.GetNbinsX()+1):
        if hel in hcov.GetXaxis().GetBinLabel(i):
            print(i,hel,hcov.GetXaxis().GetBinLabel(i))
            hel_idx[hel].append(i-1)
        print(i, hcov.GetXaxis().GetBinLabel(i))


impact_names = []
# get stat impact
# impact_names.extend(helXsecs)
# retrieve impacts per group of nuisance
hgroupimp = f.Get('nuisance_group_impact_helpois')
for ix in range(1,hgroupimp.GetNbinsX()+1):
    for iy in range(1,hgroupimp.GetNbinsY()+1):
        print(hgroupimp.GetBinContent(ix,iy))
# groupimp = narf.root_to_hist(hgroupimp).values()
print(hgroupimp.GetNbinsX(),hgroupimp.GetNbinsY())
groupimp = np.array(hgroupimp).reshape(288+2,3+2)[1:-1,1:-1]
print(groupimp.shape,hcoeff.shape)
print(len(hel_idx['helXsec_UL']))
print(groupimp[hel_idx['helXsec_UL'],:].shape)
impacts = groupimp[hel_idx['helXsec_UL'],:]/hcoeff.T.ravel()[:,np.newaxis]
print(groupimp,hcoeff.ravel(),hcoeff.T.ravel())
for j in range(nBins): #loop over UL POIs
    for i in range(groupimp.shape[1]): #loop over groups
        impact_names.append(hgroupimp.GetYaxis().GetBinLabel(i+1))


# # print(np.sqrt(np.sum(np.dot(impacts,impacts))), 'total error summing impacts in quadrature')
# print(qtBins)
# print(qtBinsS)
# plt.rcdefaults()
fig, ax = plt.subplots(figsize=(7, 5))
# hep.cms.text('work in progress', loc=0)

for i in range(impacts.shape[1]):
    ax.plot(binsC, impacts[:,i], label=impact_names[i])
for j in range(0, len(yBinsC)+1):
    plt.axvline(x=j*len(qtBinsC)+binsC[0],linewidth=1,linestyle='dashed')
ax.legend(bbox_to_anchor=(1.04,1), borderaxespad=0, frameon=False)
ax.set_yscale('log')
ax.set_ylabel('relative impact on $\sigma^{U+L}$')
ax.set_xlabel('unrolled $q_T$-y bins')
plt.tight_layout()
plt.savefig('impactsOnPOIs_syst_{}.png'.format("asimov" if asimov else "data"),dpi=300)

# helicities = ['UL']
# coefficients = ['unpolarizedxsec']
# # integrated coefficients
# for k,c in enumerate(coefficients):
#     print('analysing', c, helicities[k])
#     hcoeff = np.zeros([len(qtBinsC)])
#     hcoeff_gen = np.zeros([len(qtBinsC)])
#     hcoeff_err = np.zeros([len(qtBinsC)])
#     hcoeff_y = np.zeros([len(yBinsC)])
#     hcoeff_gen_y = np.zeros([len(yBinsC)])
#     hcoeff_err_y = np.zeros([len(yBinsC)])
#     hcoeff_hel = np.zeros([len(qtBinsC)])
#     hcoeff_hel_gen = np.zeros([len(qtBinsC)])
#     hcoeff_hel_err = np.zeros([len(qtBinsC)])
#     for ev in fitresults: #dummy because there's one event only
#         for j in range(len(qtBinsC)):
#             try:
#                 coeff = eval('ev.qt_{j}_helmeta_{c}'.format(c=c, j=j, i=i))
#                 coeff_gen = eval('ev.qt_{j}_helmeta_{c}_gen'.format(c=c, j=j, i=i))
#                 coeff_err = eval('ev.qt_{j}_helmeta_{c}_err'.format(c=c, j=j, i=i))
#                 coeff_hel = eval('ev.helXsecs{c}_qt_{j}_sumxsec'.format(c=helicities[k], j=j, i=i))
#                 coeff_hel_gen = eval('ev.helXsecs{c}_qt_{j}_sumxsec_gen'.format(c=helicities[k], j=j, i=i))
#                 coeff_hel_err = eval('ev.helXsecs{c}_qt_{j}_sumxsec_err'.format(c=helicities[k], j=j, i=i))
#                 hcoeff[j]=coeff
#                 hcoeff_gen[j]=coeff_gen
#                 hcoeff_err[j]=coeff_err
#                 hcoeff_hel[j]=coeff_hel
#                 hcoeff_hel_gen[j]=coeff_hel_gen
#                 hcoeff_hel_err[j]=coeff_hel_err
#             except AttributeError:
#                 pass
#         for i in range(len(yBinsC)):
#             try:
#                 coeff_y = eval('ev.y_{i}_helmeta_{c}'.format(c=c, j=j, i=i))
#                 coeff_gen_y = eval('ev.y_{i}_helmeta_{c}_gen'.format(c=c, j=j, i=i))
#                 coeff_err_y = eval('ev.y_{i}_helmeta_{c}_err'.format(c=c, j=j, i=i))
#                 hcoeff_y[i]=coeff_y
#                 hcoeff_gen_y[i]=coeff_gen_y
#                 hcoeff_err_y[i]=coeff_err_y
#             except AttributeError:
#                 pass

# impact_names = []
# # get stat impact
# # impact_names.extend(helXsecs)
# # retrieve impacts per group of nuisance
# hgroupimp = f.Get('nuisance_group_impact_helmetapois')
# # for i in range(1,14+1):
# #     print(hgroupimp.GetXaxis().GetBinLabel(i))

# groupimp = hist2array(hgroupimp)[:,:]
# impacts = groupimp[:6,:]/hcoeff_y.ravel()[:,np.newaxis]
# print(impacts.shape)
# for j in range(6): #loop over UL POIs
#     for i in range(groupimp.shape[1]): #loop over groups
#         impact_names.append(hgroupimp.GetYaxis().GetBinLabel(i+1))

# # print(np.sqrt(np.sum(np.dot(impacts,impacts))), 'total error summing impacts in quadrature')
# print(qtBins)
# print(qtBinsS)
# plt.rcdefaults()
# fig, ax = plt.subplots(figsize=(7, 5))
# hep.cms.text('work in progress', loc=0)

# for i in range(impacts.shape[1]):
#     ax.plot(yBinsC, impacts[:,i], label=impact_names[i])
# # for j in range(0, len(yBinsC)+1):
# #     plt.axvline(x=j*len(qtBins)+binsC[0],linewidth=1,linestyle='dashed')
# ax.legend(bbox_to_anchor=(1.04,1), borderaxespad=0, frameon=False)
# ax.set_yscale('log')
# ax.set_ylabel('relative impact on $W y$')
# ax.set_xlabel('$y$')
# plt.tight_layout()
# plt.savefig('impactsOnPOIs_syst_yint_{}.pdf'.format("asimov" if asimov else "data"),dpi=300)