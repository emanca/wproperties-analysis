import ROOT
import os
import sys
import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import scipy.stats as stats
sys.path.append('data/')
from root_numpy import hist2array
from binning import yBins, qtBins, ptBins, etaBins, mTBins, isoBins, chargeBins
from math import pi
import math 

def normal(mu,sigma):
    def f(x):
        z = 1.0*(x-mu)/sigma
        e = math.e**(-0.5*z**2)
        C = math.sqrt(2*math.pi)*sigma
        return 1.0*e
    return f


ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()

plt.style.use([hep.style.ROOT])
# hep.cms.label(loc=0, year=2016, lumi=35.9, data=True)
helicities = ['UL','I', 'T', 'A', 'P']
# coefficients = ['unpolarizedxsec']
coefficients = ['A0']


fIntoy = ROOT.TFile.Open('../Fit/FitRes/fit_WPlus_toy.root')
fitresultstoy = fIntoy.Get('fitresults')

charge = "WPlus"
f_aMC = ROOT.TFile.Open('/scratchnvme/wmass/REWEIGHT/genInfo_syst.root')
qt_aMC = np.sum(hist2array(f_aMC.Get('angularCoefficients_{}/YqTcT'.format("Wminus" if charge=="WMinus" else "Wplus"))),axis=-1)
fcoeffs = ROOT.TFile.Open('/scratchnvme/wmass/REWEIGHT/genInput_v7_syst_{}.root'.format("Wminus" if charge=="WMinus" else "Wplus"))

hists_aMC = []
for i in range(5):
    tmp = hist2array(fcoeffs.Get('angularCoefficients/harmonicsA{}_nom_nom'.format(i)))
    hists_aMC.append(tmp)
coeffs_aMC = np.stack(hists_aMC,axis=-1)

htot = hist2array(f_aMC.Get('angularCoefficients_{}/mapTot'.format("Wminus" if charge=="WMinus" else "Wplus")))

d = ROOT.RDataFrame('fitresults', '../Fit/FitRes/fit_WPlus_toy.root')

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
yBins = np.array(yBins[:threshold_y+1])
qtBins = np.array(qtBins[:threshold_qt+1])
yBinsC = 0.5*(yBins[1:]+yBins[:-1])
qtBinsC = 0.5*(qtBins[1:]+qtBins[:-1])
print(yBins,qtBins)
print(yBinsC,qtBinsC)

hists = []
canvs = []
bins = np.linspace(-5.,5.,101)
for k,c in enumerate(coefficients):
    for i in range(len(yBinsC)):
        for j in range(len(qtBinsC)):
            try:
                # htmp = d.Define('pull_y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i),'(y_{i}_qt_{j}_{c}/{fact}-{htot})/(y_{i}_qt_{j}_{c}_err/{fact})'.format(c=c, j=j, i=i, fact=((3./16./math.pi)*35.9),htot=htot[i,j]/35.9)).Histo1D('pull_y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i))
                htmp = d.Filter('status==0').Define('pull_y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i),'(y_{i}_qt_{j}_{c}+(-1)*{htot})/(y_{i}_qt_{j}_{c}_err)'.format(c=c, j=j, i=i, fact=((3./16./math.pi)*35.9),htot=coeffs_aMC[i,j,0])).Histo1D('pull_y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i))
                # htmp = d.Filter('status==0').Define('pull_y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i),'(y_{i}_qt_{j}_{c}-y_{i}_qt_{j}_{c}_gen)/(y_{i}_qt_{j}_{c}_err)'.format(c=c, j=j, i=i, fact=((3./16./math.pi)*35.9),htot=coeffs_aMC[i,j,0])).Histo1D('pull_y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i))
                # hists.append(htmp)
                canv = ROOT.TCanvas('y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i))
                # canvs.append(canv)
            except AttributeError:
                pass
htmp = d.Define('pull_mass'.format(c=c, j=j, i=i),'mass/mass_err/{}'.format(0.4158529)).Histo1D(('pull_mass', 'pull_mass',100, -5,5),'pull_mass')
hists.append(htmp)
canv = ROOT.TCanvas('mass')
canvs.append(canv)


means=[]
sigmas=[]
for i,c in enumerate(canvs):
    c.cd()
    fit=hists[i].Fit('gaus', 'S')
    means.append(fit.Get().Parameter(1))
    sigmas.append(fit.Get().Parameter(2))
    if means[i]>10: print(hists[i].GetName())
    # hists[i].Draw()
    # c.SaveAs('FIT/{}.png'.format(hists[i].GetName()))
    binspull = np.linspace(-5,5,101)
    fig, ax1 = plt.subplots()
    hep.cms.text('work in progress', loc=0, ax=ax1)
    hep.histplot(np.array(hists[i].GetValue())[1:-1],binspull, ax=ax1)
    f = normal(fit.Get().Parameter(1), fit.Get().Parameter(2))
    plt.plot(binspull, fit.Get().Parameter(0)*f(binspull))
    ax1.text(0.95, 0.95, "mean: {:.2f}+/-{:.2f} \n $\sigma$: {:.2f}+/-{:.2f}\n".format(means[i],fit.Get().ParError(1),sigmas[i],fit.Get().ParError(2)),verticalalignment='top', horizontalalignment='right',
        transform=ax1.transAxes,
        color='black')
    ax1.set_ylabel("counts")
    ax1.set_xlabel("pull mass")
    plt.tight_layout()
    plt.savefig("pullsmass.pdf")
    plt.savefig("pullsmass.png")
means=np.array(means)
sigmas=np.array(sigmas)
bins=np.linspace(0,means.shape[0]+1,means.shape[0]+1)
fig, ax1 = plt.subplots()
hep.cms.text('work in progress', loc=0, ax=ax1)
hep.histplot(means,bins,yerr=sigmas,histtype = 'errorbar', ax=ax1)
ax1.set_ylabel("$\mu_{pull}\pm\sigma_{pull}$")
ax1.set_xlabel("unrolled $y-q_T$ bin")
plt.tight_layout()
plt.savefig("pullsA0.pdf")
plt.savefig("pullsA0.png")

# #nuisance

# processes = []
# processes_names = []
# for i in range(1,5):
#     processes.extend(["LHEScaleWeight{}_muRmuF".format(i), "LHEScaleWeight{}_muR".format(i),"LHEScaleWeight{}_muF".format(i)])
#     processes_names.extend([r'$q_T^W~bin~{}~\mu_R\mu_F$'.format(i), r'$q_T^W~bin~{}~\mu_R$'.format(i),r'$q_T^W~bin~{}~\mu_F$'.format(i)])
# # processes.extend(["pdf{}".format(i) for i in range(1,103)])

# nuis = np.zeros((len(processes)))
# nuis_err = np.zeros((len(processes)))

# hists = []
# canvs = []
# for iproc,proc in enumerate(processes):
#     print(proc)
#     htmp = d.Filter('status==0').Define('pull_{c}'.format(c=proc),'({c})/({c}_err)'.format(c=proc)).Histo1D('pull_{c}'.format(c=proc))
#     hists.append(htmp)
#     canv = ROOT.TCanvas('{c}'.format(c=proc))
#     canvs.append(canv)

# means=[]
# sigmas=[]
# for i,c in enumerate(canvs):
#     c.cd()
#     fit=hists[i].Fit('gaus', 'S')
#     means.append(fit.Get().Parameter(1))
#     sigmas.append(fit.Get().Parameter(2))
#     hists[i].Draw()
#     c.SaveAs('FIT/{}.png'.format(hists[i].GetName()))
# means=np.array(means)
# sigmas=np.array(sigmas)
# bins=np.linspace(0,means.shape[0]+1,means.shape[0]+1)
# fig, ax1 = plt.subplots()
# hep.cms.text('work in progress', loc=0, ax=ax1)
# hep.histplot(means,bins,yerr=sigmas,histtype = 'errorbar', ax=ax1)
# ax1.set_ylabel("$\mu_{pull}\pm\sigma_{pull}$")
# ax1.set_xlabel("nuisance")
# plt.tight_layout()
# plt.show()