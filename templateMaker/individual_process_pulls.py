import ROOT
import pandas as pd
from IPython.display import Image
import numpy as np

'''plots every single process distribution from toy dataset.'''


fIntoy = ROOT.TFile.Open('../Fit/FitRes/fit_Wlike_1ktoys_newbins.root')
fitresultstoy = fIntoy.Get('fitresults')

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

#filtered_coeffs = [names.query("names.str.contains('_L_|_T_|_I_y_1.8_qt_44.8')\
#                & names.str.contains('L_y_1.8|qt_44.8')\
#                &~names.str.contains('sum|mask|minos|err|gen')")['names'].values]

ROOT.gROOT.SetBatch(True)
for icoeff,coeff in enumerate(filtered_coeffs):
    for iproc, proc in enumerate(coeff):
        hists = {}
        canvs = {}
        hpull = ROOT.TH1D(proc+'_pull',proc+"_pull", 100, -5.,5.)
        canv = ROOT.TCanvas(proc,proc)
        ROOT.gStyle.SetOptStat(111100)


        hists[proc]=hpull
        canvs[proc]=canv

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