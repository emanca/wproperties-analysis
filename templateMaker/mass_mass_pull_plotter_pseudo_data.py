import ROOT
import pandas as pd
from IPython.display import Image
import numpy as np
import time
import argparse

parser = argparse.ArgumentParser("")
parser.add_argument('-n',type=int, default=100, help='iteration of boostrap to run')
args = parser.parse_args()
n = args.n


proc = 'mass_var_mass_var_0.5_noi'
suffix = 'Bnom'
#suffix = 'Bnom'


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(111100)
ROOT.gStyle.SetOptFit(1)

hists = {}
canvs = {}

hpull    = ROOT.TH1D(proc+'_pull',proc+"_pull", 100, -5.,5.)
pullcanv = ROOT.TCanvas(proc+'_pull',proc+'_pull')

hmass     = ROOT.TH1D(proc,proc, 100, -5.,5.)
masscanv = ROOT.TCanvas(proc,proc)

hists[proc+'_pull']=hpull
canvs[proc+'_pull']=pullcanv
hists[proc]=hmass
canvs[proc]=masscanv

for ifile in range(n):
    fIntoy = ROOT.TFile.Open('../Fit/FitRes/Helicity_tests/fit_Wlike_iteration_{}_twocharges_nominal.root'.format(ifile))
    #fIntoy = ROOT.TFile.Open('../Fit/FitRes/fit_Wlike_iteration_{}_Bnom.root'.format(ifile))
    fitresultstoy = fIntoy.Get('fitresults')
    for ev in fitresultstoy:
        for key , hist in hists.items():
            par = getattr(ev,proc)
            par_err = getattr(ev,proc+'_err')
            par_gen = getattr(ev,proc+'_gen')
            if key == proc:
                hist.Fill(par)
            else:
                hist.Fill((par-par_gen)/par_err)


for i,c in canvs.items():
    c.cd()
    hists[i].Fit('gaus','L')
    hists[i].Draw()
    #c.SaveAs('{}_{}.png'.format(hists[i].GetName(),suffix))
    c.SaveAs('Helicity_tests/{}_{}.png'.format(hists[i].GetName(),suffix))

