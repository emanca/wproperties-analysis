import ROOT
import os
import sys
import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
sys.path.append('data/')
from binning import yBins, qtBins, ptBins, etaBins, mTBins, isoBins, chargeBins

ROOT.gROOT.SetBatch(True)

plt.style.use([hep.style.ROOT])
# hep.cms.label(loc=0, year=2016, lumi=35.9, data=True)
helicities = ['UL','I', 'T', 'A', 'P']
# coefficients = ["A0", "A1", "A2", "A3", "A4",'unpolarizedxsec']
coefficients = ['unpolarizedxsec']

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
yBins_red = np.array(yBins[:threshold_y+1])
qtBins_red = np.array(qtBins[:threshold_qt+1])
yBinsC = 0.5*(yBins_red[1:]+yBins_red[:-1])
qtBinsC = 0.5*(qtBins_red[1:]+qtBins_red[:-1])


hists = {}
canvs = {}
hpull = ROOT.TH1D("mass","mass", 100, -10.,10.)
canv = ROOT.TCanvas('mass','mass')

hists['mass_var_mass_var_0.5_noi']=hpull
canvs['mass_var_mass_var_0.5_noi']=canv
for k,c in enumerate(coefficients):
    for i in range(len(yBinsC)):
        for j in range(len(qtBinsC)):
            hpull = ROOT.TH1D("pull_y_{i}_qt_{j}_{c}".format(c=c, j=j, i=i),"pull_y_{i}_qt_{j}_{c}".format(c=c, j=j, i=i), 100, -10,10)
            canv = ROOT.TCanvas('pull_y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i))
            hists['y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i)]=hpull
            canvs['y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i)]=canv
for ifile in range(1):
    # fIntoy = ROOT.TFile.Open('../Fit/FitRes/fit_WPlus_data_{}.root'.format(ifile))
    fIntoy = ROOT.TFile.Open('../Fit/FitRes/fit_Wplus_asimovtoys.root')
    fitresultstoy = fIntoy.Get('fitresults')
    for ev in fitresultstoy:
        for key, hist in hists.items():
            if 'mass' in key: 
                continue
                par=(getattr(ev,'{}'.format(key)))
                par_err=(getattr(ev,'{}_err'.format(key)))
                hist.Fill(par/par_err)
                # hist.Fill(par)
            else:
                i = int(key.split('_')[1])
                j = int(key.split('_')[3])
                coeff = coefficients.index(key.split('_')[4])
                print(key, i,j,coeff)
                coeff = getattr(ev,'y_{:.1f}_qt_{:.1f}_WplusmunuPostVFP_{}'.format(yBinsC[i] , qtBinsC[j] , c))
                coeff_gen = getattr(ev,'y_{:.1f}_qt_{:.1f}_WplusmunuPostVFP_{}_gen'.format(yBinsC[i] , qtBinsC[j] , c))
                coeff_err = getattr(ev,'y_{:.1f}_qt_{:.1f}_WplusmunuPostVFP_{}_err'.format(yBinsC[i] , qtBinsC[j] , c))
                hist.Fill((coeff-coeff_gen)/coeff_err)
                # hist.Fill(coeff-coeff_gen)
                print(coeff-coeff_gen)

ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetOptStat(111111)

means=[]
sigmas=[]
sigmasErr =[]
for i,c in canvs.items():
    c.cd()
    fit=hists[i].Fit('gaus','LS')
    hists[i].Draw()
    c.SaveAs('toys/{}.png'.format(hists[i].GetName()))
    means.append(fit.Get().Parameter(1))
    sigmas.append(fit.Get().Parameter(2))
    sigmasErr.append(fit.Get().ParError(2))

# means=np.array(means)
# sigmas=np.array(sigmas)
# sigmasErr=np.array(sigmasErr)
# bins=np.linspace(0,means.shape[0]+1,means.shape[0]+1)
# fig, ax1 = plt.subplots()
# hep.cms.text('work in progress', loc=0, ax=ax1)
# # hep.histplot(means,bins,yerr=sigmas,histtype = 'errorbar', ax=ax1)
# hep.histplot(sigmas,bins,yerr=sigmasErr,histtype = 'errorbar', color='red', ax=ax1)
# ax1.set_ylabel("$\mu_{pull}\pm\sigma_{pull}$")
# ax1.set_xlabel("unrolled $y-q_T$ bin")
# plt.tight_layout()
# plt.savefig("pars_halfStat.pdf")
# plt.savefig("pars_halfStat.png")