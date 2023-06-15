import ROOT
import pandas as pd
from IPython.display import Image
import numpy as np
import time


fIntoy = ROOT.TFile.Open('../Fit/FitRes/fit_Wlike_1ktoys_newbins.root')
fitresultstoy = fIntoy.Get('fitresults')

proc = 'mass_var_mass_var_0.5_noi'


#ROOT.gROOT.SetBatch(True)


hists = {}
canvs = {}

hpull    = ROOT.TH1D(proc+'_pull',proc+"_pull", 100, -5.,5.)
pullcanv = ROOT.TCanvas(proc+'_pull',proc+'_pull')

hmass     = ROOT.TH1D(proc,proc, 100, -5.,5.)
masscanv = ROOT.TCanvas(proc,proc)



ROOT.gStyle.SetOptStat(111100)


hists[proc+'_pull']=hpull
canvs[proc+'_pull']=pullcanv
hists[proc]=hmass
canvs[proc]=masscanv


for ev in fitresultstoy:
    for key , hist in hists.items():
        par = getattr(ev,proc)
        par_err = getattr(ev,proc+'_err')
        par_gen = getattr(ev,proc+'_gen')
        if key == proc:
            hist.Fill(par)
        else:
            hist.Fill((par-par_gen)/par_err)

ROOT.gStyle.SetOptFit(1)
for i,c in canvs.items():
    c.cd()
    hists[i].Fit('gaus','L')
    c.Draw()
    c.SaveAs('{}_1ktoys.png'.format(hists[i].GetName()))