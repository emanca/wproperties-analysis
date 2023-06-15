import ROOT
from IPython.display import Image
import pandas as pd
import matplotlib.pyplot as plt
import argparse

'''This script makes pull distributions for each of the angular coefficients for toy datasets. One plot is made per angular coefficient,
plus the mass_noi distribution. The individual process distributions are extracted and fit to a gaussian. The fit gaussian parameters 
are gathered and plotted for each angular coefficient'''




fIntoy = ROOT.TFile.Open('../Fit/FitRes/fit_Wlike_1ktoys_newbins.root')
fitresultstoy = fIntoy.Get('fitresults')

branch_names = []
for branch in fitresultstoy.GetListOfBranches():
    branch_names.append(branch.GetName())

#angular coefficients of interest
coeffs = ['mass_var_mass_var_0.5_noi','unpolarizedxsec','A0','A1','A2','A3','A4']


#names of processes in branch. contains unwanted information
names = pd.DataFrame({'names':branch_names})

#filtering by coefficients 
filtered_coeffs = [names.query("names.str.contains(@coeffs[@icoeff]) & \
        ~names.str.contains('err|gen|minos|hel')")['names'].values for icoeff in range(len(coeffs))]


for icoeff,coeff in enumerate(filtered_coeffs):
    hists  = {}
    for iproc, proc in enumerate(coeff):
        hists[proc]=ROOT.TH1D(proc,proc, 100, -5.,5.)
        ROOT.gStyle.SetOptStat(111100) #this line adds underflow/overflow bins


    for ev in fitresultstoy:
        for key , hist in hists.items():
            par = getattr(ev,key)
            par_err = getattr(ev,key+'_err')
            par_gen = getattr(ev,key+'_gen')
            hist.Fill((par-par_gen)/par_err)

    means = []
    sigmas = []
    
    ROOT.gROOT.SetBatch(True) #disables display of plots as code runs
    ROOT.gStyle.SetOptFit(1)  #enables the display of fit information
    for name,hist in hists.items():
        fit = hists[name].Fit('gaus','LS')
        means.append(fit.Get().Parameter(1))
        sigmas.append(fit.Get().Parameter(2))
    plt.rcParams['font.size']=30
    plt.figure(figsize = (15,7))
    plt.errorbar(x=range(len(means)) , y=means , yerr=sigmas, fmt='o')
    plt.ylim(-2,2)
    plt.ylabel(r'$\mu_{pull}\pm{}\sigma_{pull}$')
    plt.xlabel(r'unrolled $y-q_T$ bin')
    plt.title('{} pull for 1000 toys'.format(coeffs[icoeff]))
    plt.xticks(ticks=range(len(means)),labels=coeff,rotation=90,fontsize=10)
    plt.tight_layout()
    plt.savefig('{}_pull_1000toys.png'.format(coeffs[icoeff]))
    plt.show()
    plt.close()