import ROOT
from array import array
import math
import copy
import argparse
import os
import numpy as np 
import root_numpy as rootnp
import ctypes
import scipy
ROOT.gROOT.SetBatch(True)
# ROOT.TH1.AddDirectory(False)
# ROOT.TH2.AddDirectory(False)

def cloneEmpty2D(histo,name) :
    hout = histo.Clone(name)
    for i in range(0, histo.GetNbinsX()+2) : #also overflow and underflow
        for j in range(0, histo.GetNbinsY()+2) :
            hout.SetBinContent(i,j,0)
            hout.SetBinError(i,j, 0)
    return hout 

Nrand = 1000 # number of template-toys produced
AC_list = ['L',"I","T","A","P","7","8","9","UL"]
Nqt = 11
Ny=6
# BKG_list = ["DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"] 
BKG_list=[]

outFile = ROOT.TFile('study_sumOfTemplate.root',"recreate")


chargeList = ['plus','minus']
for s in chargeList :
    inFileName = '/scratch/bertacch/wmass/wproperties-analysis/Fit/OUTPUT_25Apr_noCF_qtMax60_fixLowAcc/W'+s+'_reco.root'
    inFile = ROOT.TFile.Open(inFileName)
    
    data_obs = cloneEmpty2D(inFile.Get("data_obs"), 'data_obs')
    for q in range(1,Nqt+1) :
        for y in range (1,Ny+1) :
            for c in AC_list :
                h = inFile.Get('helXsecs'+c+'_y_'+str(y)+'_qt_'+str(q))
                data_obs.Add(h)
    for c in BKG_list :
        h = inFile.Get(c)
        data_obs.Add(h)
                                
    outFile.cd()
    data_obs.Write()
    
    # assert(0)#no minus
        
        
    
            