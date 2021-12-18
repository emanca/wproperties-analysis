import os
import ROOT
import copy
import sys
import argparse
import math
from ROOT import gStyle
ROOT.gROOT.SetBatch(True)
# ROOT.TH1.AddDirectory(False)
# ROOT.TH2.AddDirectory(False)
gStyle.SetOptStat('mr')

parser = argparse.ArgumentParser("")
parser.add_argument('-i','--input', type=str, default='/scratchssd/sroychow/NanoAOD2016-V2/',help="name of the input direcory")
parser.add_argument('-p', '--project',type=int, default=True, help="projection from input")

args = parser.parse_args()
inputDir = args.input
PROJECT = args.project

sampleDict = {
    # 'Data' : ['SingleMuon_Run2016',1, 'B_ver2']
    'Data' : ['SingleMuon_Run2016',1, 'B_ver2','C','D','E','F','G','H'],
    'MC_WToMu' : ['WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',2,'','_ext2']
    # 'MC_WToMu' : ['WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',2,'_ext2']
}



# preselection = "(Vtype==0 or Vtype=1) && HLT_SingleMu24 && Mu1_pt>25.0 && Mu1_pt<55.0 && MET_filters==1 && nVetoElectrons==0 && 1"
preselection = "(Vtype==0 || Vtype==1) && HLT_SingleMu24 && Muon_corrected_pt[Idx_mu1]>25.0 && Muon_corrected_pt[Idx_mu1]<55.0 && MET_filters==1 && nVetoElectrons==0 && 1"
mtString = 'sqrt(2*Muon_corrected_pt[Idx_mu1]*MET_pt_nom*(1.0-cos(Muon_phi[Idx_mu1]-MET_phi_nom)))'
signDict = {
    'Plus' : '&& Muon_charge[Idx_mu1]>0',
    'Minus' : ' && Muon_charge[Idx_mu1]<0'
}


if PROJECT :
    NbinMt=300
    NbinRelIso=500

    histoDict = {}
    output = ROOT.TFile("bkg_ABCD_plotter.root", "recreate")

    for sample, info in sampleDict.items() :
        for part in info[2:] : 
            inFile = ROOT.TFile.Open(inputDir+'/'+info[0]+part+'/tree.root')
            tree = inFile.Get('Events')
            for s,sstring in signDict.items() :
            
                histoDict[sample+part+s] = ROOT.TH2F('h_RelIso_vs_Mt'+sample+part+s,'h_RelIso_vs_Mt'+sample+part+s,NbinMt,0,150,NbinRelIso,0,0.6)
                tree.Project('h_RelIso_vs_Mt'+sample+part+s,'Muon_pfRelIso04_all[Idx_mu1]:'+mtString,preselection+sstring,"")#1000
                
                if sample!='Data' :
                    xsec = 61526.7/ 0.001
                    lumi = 35.9
                    lumiWeight = (xsec*lumi)/histoDict[sample+part+s].Integral()
                    histoDict[sample+part+s].Scale(lumiWeight)
                print(info[0], sample, part, histoDict[sample+part+s].GetEntries())

                output.cd()
                histoDict[sample+part+s].Write()
    output.Close()

step2input = ROOT.TFile.Open('bkg_ABCD_plotter.root')
step2output = ROOT.TFile("bkg_ABCD_plotter_added.root", "recreate")

cmslab = "#bf{CMS} #scale[0.7]{#it{Work in progress}}"
lumilab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
cmsLatex = ROOT.TLatex()

addedDict = {}  
for sample, info in sampleDict.items() :
    for s,sstring in signDict.items() :
        addedDict[sample+s] = step2input.Get('h_RelIso_vs_Mt'+sample+info[3]+s).Clone('h_RelIso_vs_Mt'+sample+s)
        for part in info[3:] : 
            temph = step2input.Get('h_RelIso_vs_Mt'+sample+part+s).Clone('h_RelIso_vs_Mt'+sample+s+part)
            addedDict[sample+s].Add(temph)
        addedDict[sample+s].GetXaxis().SetTitle("m_{T} [GeV]")
        addedDict[sample+s].GetYaxis().SetTitle("RelIso")
        addedDict[sample+s].SetTitle("Base-selection ("+sample+")")
        addedDict[sample+s].GetXaxis().SetTitleSize(0.05)
        addedDict[sample+s].GetYaxis().SetTitleSize(0.05)
        addedDict[sample+s].GetXaxis().SetTitleOffset(0.9)
        addedDict[sample+s].GetYaxis().SetTitleOffset(0.9)
        step2output.cd()
        addedDict[sample+s].Write()
    
        if sample=='Data' and s=='Plus' :
            can = ROOT.TCanvas('c_RelIso_vs_Mt'+sample+s,'c_RelIso_vs_Mt'+sample+s,1600,1200)
            can.cd()
            addedDict[sample+s].SetMarkerColorAlpha(1,0.2)
            addedDict[sample+s].GetYaxis().SetRangeUser(0.0012,0.4)
            addedDict[sample+s].Draw()
            addedDict[sample+s].SetTitle('')
            addedDict[sample+s].SetStats(0)
            
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-can.GetRightMargin(),1-0.8*can.GetTopMargin(),lumilab)

            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(can.GetLeftMargin(),1-0.8*can.GetTopMargin(),cmslab)
            
            cmsLatex.SetTextAlign(31) 
            namelab = " #scale[0.99]{Data (single #mu^{+})}"
            # namelab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
            cmsLatex.DrawLatex(1-1.05*can.GetRightMargin(),1-1.6*can.GetTopMargin(), namelab)
            
            can.SaveAs('bkg_ABCD_dataPlus.pdf')
            can.SaveAs('bkg_ABCD_dataPlus.png')
            
        
step2output.Close()
            
