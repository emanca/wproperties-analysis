import ROOT
import pandas as pd
from IPython.display import Image
import numpy as np



'''plots every single process distribution from toy dataset.'''


fIntoy = ROOT.TFile.Open('../Fit/FitRes/fit_Wlike_1ktoys_newbins.root')
fitresultstoy = fIntoy.Get('fitresults')
suffix = 'Bnom'

branch_names = []
for branch in fitresultstoy.GetListOfBranches():
    branch_names.append(branch.GetName())

#angular coefficients of interest
#distribution plots will be made for all processes
coeffs = ['mass_var_mass_var_0.5_noi','unpolarizedxsec','A0','A1','A2','A3','A4']

#names of processes in pranch. contains unwanted information
names = pd.DataFrame({'names':branch_names})

#filtering by coefficients 
filtered_coeffs = [names.query("names.str.contains(@coeffs[@icoeff]) & \
        ~names.str.contains('err|gen|minos|hel')")['names'].values for icoeff in range(len(coeffs))]

print(filtered_coeffs)

ROOT.gROOT.SetBatch(True)
for icoeff,coeff in enumerate(filtered_coeffs):
    print(icoeff,coeff)
    for iproc, proc in enumerate(coeff):
        print(iproc,proc)
        hists = {}
        canvs = {}
        hpull = ROOT.TH1D(proc+'_pull',proc+"_pull", 100, -5.,5.)
        canv = ROOT.TCanvas(proc,proc)
        ROOT.gStyle.SetOptStat(111100)
        hists[proc]=hpull
        canvs[proc]=canv
        for ifile in range(100):
            fIntoy = ROOT.TFile.Open('../Fit/FitRes/fit_Wlike_iteration_{}_{}.root'.format(ifile,suffix))
            fitresultstoy = fIntoy.Get('fitresults')
            for ev in fitresultstoy:
                for key , hist in hists.items():
                    par = getattr(ev,key)
                    par_err = getattr(ev,key+'_err')
                    par_gen = getattr(ev,key+'_gen')
                    hist.Fill((par-par_gen)/par_err)
            ROOT.gStyle.SetOptFit(1)
            for i,c in canvs.items():
                c.cd()
                hists[i].Fit('gaus','L')
                #hists[i].Draw()
                c.SaveAs('{}_1ktoys.png'.format(hists[i].GetName()))