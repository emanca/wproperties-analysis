import ROOT
ROOT.gROOT.SetBatch(True)
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
f = ROOT.TFile.Open('../Fit/FitRes/fit_WMinus_{}.root'.format('asimov' if asimov else 'data'))
fitresults = f.Get('fitresults')
hcov = f.Get('covariance_matrix_channelmu')
hcorr = f.Get('correlation_matrix_channelmu')

cov = hist2array(hcov)[:,:]

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
yBins_red = np.array(yBins[:threshold_y+1])
qtBins_red = np.array(qtBins[:threshold_qt+1])
yBinsC = 0.5*(yBins_red[1:]+yBins_red[:-1])
qtBinsC = 0.5*(qtBins_red[1:]+qtBins_red[:-1])

processes = []
processes_names = []
# for i in range(1,5):
#     processes.extend(["LHEScaleWeight{}_muRmuF".format(i), "LHEScaleWeight{}_muR".format(i),"LHEScaleWeight{}_muF".format(i)])
#     processes_names.extend([r'$q_T^W~bin~{}~\mu_R\mu_F$'.format(i), r'$q_T^W~bin~{}~\mu_R$'.format(i),r'$q_T^W~bin~{}~\mu_F$'.format(i)])
# processes.extend(["pdf{}".format(i) for i in range(1,103)])
processes.extend(["SFall{}".format(i) for i in range(624)]+["SFiso{}".format(i) for i in range(624)]+["SFSyst"])
# processes.extend(["SFall","SFiso"])
# processes.extend(['fakeShapeBin{}'.format(i) for i in range(48*30)]+["fakesNormLowMt","fakesNormHighMt"])

# processes.extend(["jes", "uncl", "prefire"])
# processes.extend(["muRmuF", "muR","muF"])

nuis = np.zeros((len(processes)))
nuis_err = np.zeros((len(processes)))

for iproc,proc in enumerate(processes):
    for ev in fitresults:
        nuis[iproc]=eval('ev.{}'.format(proc))
    ibin = hcov.GetXaxis().FindBin(proc)
    nuis_err[iproc]=math.sqrt(hcov.GetBinContent(ibin,ibin))

ccorr=ROOT.TCanvas("corr","",2000,2000)
mat = ROOT.TH2D("corr","", len(processes)-1, 0, len(processes), hcorr.GetNbinsX(), 0, hcorr.GetNbinsX()+1)
for i in range(1,hcorr.GetNbinsX()+1):
    for j in range(1,hcorr.GetNbinsY()+1):
        if "SFiso" in hcorr.GetXaxis().GetBinLabel(i):
            print(hcorr.GetXaxis().GetBinLabel(i))
            print(hcorr.GetBinContent(i,j))
            mat.SetBinContent(i,j,hcorr.GetBinContent(i,j))
mat.Draw("colz")
ccorr.SaveAs("corr.png")
assert(0)
# fig, ax = plt.subplots(figsize=(10, 10))
# hep.cms.label(loc=0,year=2016, data=False)

# y_pos = np.arange(len(processes))

# ax.hlines(y_pos, 0., nuis_err, color='green')
# ax.plot(nuis_err, y_pos, '|', color='green')  # Stem ends
# ax.set_yticks(y_pos)
# ax.set_yticklabels(processes_names, fontweight='medium')
# ax.invert_yaxis()  # labels read top-to-bottom
# ax.set_xlabel('constraint on nuisance', fontweight='medium')
# ax.set_xlim(0,1)

# # Add annotation to bars
# for i in range(len(nuis_err)):
#     plt.text(nuis_err[i],y_pos[i]-0.065, 
#              str(round((nuis_err[i]), 2)),
#              fontsize=13,
#              color='black', alpha=0.7)

# plt.tight_layout()
# plt.savefig('nuis_{}.png'.format("asimov" if asimov else "data"),dpi=300)
# # plt.show()
bins=np.linspace(0, nuis.shape[0]+1,nuis.shape[0]+1)
fig, ax1 = plt.subplots()
hep.cms.text('work in progress', loc=0, ax=ax1)
hep.histplot(nuis/nuis_err,bins,yerr=nuis_err, histtype = 'errorbar', ax=ax1)
ax1.set_ylabel("nuisance value")
ax1.set_xlabel("nuisance")
plt.tight_layout()
plt.savefig('nuis_WMinus_{}.png'.format("asimov" if asimov else "data"),dpi=300)

ptBins_SF = np.array([25.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 47.0, 50.0, 55.0])
etaBins = np.array([-2.4+i*0.1 for i in range(49)])

# nuis_map = np.zeros((etaBins.shape[0],ptBins_SF.shape[0]))
fig, ax1 = plt.subplots()
# for j in range(48): # for each eta bin
#     for i in range(len(ptBins_SF)-1): #create one histogram per variation in each macro bin of SF
#         nuis_map[j,i]=

hep.hist2dplot(nuis[:624].reshape((etaBins.shape[0]-1,ptBins_SF.shape[0]-1)),etaBins,ptBins_SF, cmap = 'jet')
ax1.set_xlabel("$\eta$")
ax1.set_ylabel("$p_T$")
plt.tight_layout()
plt.savefig('nuismap_WMinus_{}.png'.format("asimov" if asimov else "data"))