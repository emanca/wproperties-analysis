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


sampleList =  ['WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8','','_ext2']

scaleDict = {   
    "Nominal"                      : '1',
    "LHEScaleWeight_muR0p5_muF0p5" : "LHEScaleWeight[0]",
    "LHEScaleWeight_muR0p5_muF1p0" : "LHEScaleWeight[1]",
    "LHEScaleWeight_muR1p0_muF0p5" : "LHEScaleWeight[3]",
    "LHEScaleWeight_muR1p0_muF2p0" : "LHEScaleWeight[5]",
    "LHEScaleWeight_muR2p0_muF1p0" : "LHEScaleWeight[7]",
    "LHEScaleWeight_muR2p0_muF2p0" : "LHEScaleWeight[8]"    
}


mtString = 'sqrt(2*Muon_corrected_pt[Idx_mu1]*MET_pt_nom*(1.0-cos(Muon_phi[Idx_mu1]-MET_phi_nom)))'
preselection = "Vtype==0 && HLT_SingleMu24 && Muon_corrected_pt[Idx_mu1]>25.0 && Muon_corrected_pt[Idx_mu1]<55.0 && MET_filters==1 && nVetoElectrons==0 && (abs(genVtype) == 14) && 1 "+" && "+ mtString+">40"
signDict = {
    'Plus' : '&& Muon_charge[Idx_mu1]>0',
    'Minus' : ' && Muon_charge[Idx_mu1]<0'
}
weiString = 'float(puWeight*lumiweight*WHSF*weightPt*weightY)'


if PROJECT :
    histoDict = {}
    output = ROOT.TFile("bkg_wpt_MCscale.root", "recreate")

    for part in sampleList[1:] : 
        inFile = ROOT.TFile.Open(inputDir+'/'+sampleList[0]+part+'/tree.root')
        tree = inFile.Get('Events')
        for s,sstring in signDict.iteritems() :
            for scaleName, scaleNum in scaleDict.iteritems() :
                histoDict[part+s+scaleName] = ROOT.TH1F('h_Wpt'+part+s+scaleName,'h_Wpt'+part+s+scaleName,600,0,200)
                tree.Project('h_Wpt'+part+s+scaleName,'TMath::Abs(GenV_preFSR[0])*'+scaleNum,preselection+sstring,"")#,1000000)
                
                xsec = 61526.7/ 0.001
                lumi = 35.9
                lumiWeight = (xsec*lumi)/histoDict[part+s+scaleName].Integral()
                histoDict[part+s+scaleName].Scale(lumiWeight)
                
                output.cd()
                histoDict[part+s+scaleName].Write()
                print "done ", part, s, scaleName
    output.Close()

step2input = ROOT.TFile.Open('bkg_wpt_MCscale.root')
step2output = ROOT.TFile("bkg_wpt_MCscale_added.root", "recreate")

addedDict = {}  
for s,sstring in signDict.iteritems() :
    colNum = 1
    for scaleName, scaleNum in scaleDict.iteritems() :
        addedDict[s+scaleName] = step2input.Get('h_Wpt'+sampleList[1]+s+scaleName).Clone('h_Wpt'+s+scaleName)
        for part in sampleList[2:] : 
            print "done ", part, s, scaleName, 'h_Wpt'+part+s+scaleName
            temph = step2input.Get('h_Wpt'+part+s+scaleName).Clone('h_Wpt'+s+part+scaleName)
            addedDict[s+scaleName].Add(temph)
        addedDict[s+scaleName].GetXaxis().SetTitle("p_{T} [GeV]")
        addedDict[s+scaleName].GetYaxis().SetTitle("dN/dp_{T} [1/"+str(addedDict[s+scaleName].GetBinWidth(1))+" GeV^{-1}]")
        addedDict[s+scaleName].SetTitle("W "+s+ " transverse momentum")
        addedDict[s+scaleName].SetLineWidth(2)
        addedDict[s+scaleName].SetLineColor(colNum)
        colNum+=1
        step2output.cd()
        addedDict[s+scaleName].Write()
    
    #ratios
    for scaleName, scaleNum in scaleDict.iteritems() :    
        addedDict[s+scaleName+'ratio'] = addedDict[s+scaleName].Clone(addedDict[s+scaleName].GetName()+'ratio')
        addedDict[s+scaleName+'ratio'].Divide(addedDict[s+'Nominal'])
        
        
    
    


        
    addedDict['can'+s] = ROOT.TCanvas('h_Wpt'+s,'h_Wpt'+s,1600,1400)
    addedDict['can'+s].cd()
    addedDict['padH'+s] = ROOT.TPad("pad_histo_"+s, "pad_histo_"+s,0,0.2,1,1)
    addedDict['padR'+s] = ROOT.TPad("pad_ratio_"+s, "pad_ratio_"+s,0,0,1,0.2)
    addedDict['padH'+s].SetBottomMargin(0.02)
    addedDict['padH'+s].Draw()
    addedDict['padR'+s].SetTopMargin(0)
    addedDict['padR'+s].SetBottomMargin(0.25)
    addedDict['padR'+s].Draw()
    
    addedDict['leg'+s] = ROOT.TLegend(0.5, 0.5, 0.90, 0.85)
    
    addedDict['padH'+s].cd() 
    addedDict['padH'+s].SetGridx()
    addedDict['padH'+s].SetGridy()
    addedDict[s+'Nominal'].Draw('hist')
    addedDict['leg'+s].AddEntry(addedDict[s+'Nominal'],'Nominal')
    addedDict['leg'+s].Draw("SAME")
    for scaleName, scaleNum in scaleDict.iteritems() :
        if scaleName == 'Nominal' : continue
        addedDict[s+scaleName].Draw("hist same")
        addedDict['leg'+s].AddEntry(addedDict[s+scaleName],scaleName.replace('LHEScaleWeight_','').replace('p','.'))
        # addedDict[sample+s].GetYaxis().SetRangeUser(0.0012,0.4)
        # addedDict['can'+s].SaveAs('bkg_ABCD_dataPlus.pdf')
        # addedDict['can'+s].SaveAs('bkg_ABCD_dataPlus.png')
        
    addedDict['padR'+s].cd()
    addedDict['padR'+s].SetGridx()
    addedDict['padR'+s].SetGridy()
    addedDict[s+'Nominal'+'ratio'].Draw('hist')
    for scaleName, scaleNum in scaleDict.iteritems() :
        if scaleName == 'Nominal' : continue
        addedDict[s+scaleName+'ratio'].Draw("hist same")
        
    addedDict['can'+s].Write()        
        
step2output.Close()
            
