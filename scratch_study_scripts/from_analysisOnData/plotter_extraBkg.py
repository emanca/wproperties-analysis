import os
import ROOT
import copy
import sys
import argparse
import math
ROOT.gROOT.SetBatch()

ROOT.gStyle.SetOptStat(0)


def CloneEmpty(histo,name) :
    hout = histo.Clone(name)
    for i in range(0, histo.GetNbinsX()+2) :
        hout.SetBinContent(i,0)
        hout.SetBinError(i,0)
    return hout 

parser = argparse.ArgumentParser("")

parser.add_argument('-o','--output', type=str, default='plot_extraBkg',help="name of the output file")
parser.add_argument('-i','--input', type=str, default='stackPlots.root',help="name input prefit plot file")
parser.add_argument('-s','--save', type=int, default=False,help="save .png and .pdf canvas")
parser.add_argument('-n','--noratio', type=int, default=False,help="no ratio W-/W+ but direct subtraction")

args = parser.parse_args()
OUTPUT = args.output
INPUT = args.input
SAVE= args.save
NORATIO= args.noratio

inFile = ROOT.TFile.Open(INPUT)


#build the histograms
varList = ['Mu1_pt','MT','Mu1_eta']
signList = ['plus','minus']
obj = {}
for v in varList :
    for s in signList : 
        obj['can'+v+s] = inFile.Get('prefit_'+v+'_'+s)
        obj['pad'+v+s] = obj['can'+v+s].GetPrimitive('pad_histo_'+v+'_'+s)
        obj['st'+v+s] = obj['pad'+v+s].GetPrimitive('stack_'+v+'_'+s)
        obj['data'+v+s] = obj['pad'+v+s].GetPrimitive('hData_'+v+'_'+s)
        obj['hs'+v+s] = obj['st'+v+s].GetHists()
        obj['W'+v+s] = obj['hs'+v+s].At(5)
        # print (v,s,obj['W'+v+s].GetMean())
    
    #build ratio
    if not NORATIO : 
        obj['Rmp'+v] = obj['W'+v+'minus'].Clone('ratio'+v)
        obj['Rmp'+v].Divide(obj['W'+v+'plus'])
        # obj['Rmp'+v] = obj['W'+v+'minus'].GetEntries()/obj['W'+v+'plus'].GetEntries()
        # if v!='Mu1_eta' : 
        #     nbin = obj['W'+v+'minus'].FindBin(42)
        #     print("dentro v", v)
        # else :
        #     nbin=1 
        # obj['Rmp'+v] = obj['W'+v+'minus'].GetBinContent(nbin)/obj['W'+v+'plus'].GetBinContent(nbin)

        
        
        #mutiply
        obj['st'+v+'plus'+'mult'] = obj['st'+v+'plus'].Clone('stack'+v+'plus'+'mult')
        for hh in range(0,6) :  
            obj['st'+v+'plus'+'mult'].GetHists().At(hh).Multiply(obj[' Rmp'+v])
            # obj['st'+v+'plus'+'mult'].GetHists().At(hh).Scale(obj['Rmp'+v])
        obj['data'+v+'plus'+'mult'] = obj['data'+v+'plus'].Clone('data'+v+'plus'+'mult')
        obj['data'+v+'plus'+'mult'].Multiply(obj['Rmp'+v])
        # obj['data'+v+'plus'+'mult'].Scale(obj['Rmp'+v])
        
        #subtract
        obj['st'+v+'minus'+'sub'] = obj['st'+v+'minus'].Clone('stack'+v+'minus'+'sub')
        for hh in range(0,6) :  
            obj['st'+v+'minus'+'sub'].GetHists().At(hh).Add(obj['st'+v+'plus'+'mult'].GetHists().At(hh),-1)
        obj['data'+v+'minus'+'sub'] = obj['data'+v+'minus'].Clone('data'+v+'minus'+'sub')
        obj['data'+v+'minus'+'sub'].Add(obj['data'+v+'plus'+'mult'],-1)
    else :
        obj['st'+v+'minus'+'sub'] = obj['st'+v+'minus'].Clone('stack'+v+'minus'+'sub')
        obj['st'+v+'minus'+'sub'].GetHists().At(5).Add(obj['st'+v+'minus'].GetHists().At(5),-1)
        obj['st'+v+'minus'+'sub'] = obj['st'+v+'minus'].Clone('stack'+v+'minus'+'sub')
        obj['data'+v+'minus'+'sub'] = obj['data'+v+'minus'].Clone('data'+v+'minus'+'sub')
        obj['data'+v+'minus'+'sub'].Add(obj['st'+v+'minus'].GetHists().At(5),-1)





#------------------- plotting part -------------------------------

sampleOrder = ['DiBoson','WtoTau','Top','DYJets','SIGNAL_Fake','WToMu']
colorList = [ROOT.kViolet+2,ROOT.kSpring+9,ROOT.kGreen+3,ROOT.kAzure+2,ROOT.kGray,ROOT.kRed+2]

outFile =  ROOT.TFile(OUTPUT+'.root', "RECREATE")

for v in varList :
    
    legend = ROOT.TLegend(0.5, 0.5, 0.90, 0.85)            
    hStack = ROOT.THStack('stack_'+v,'stack_'+v)
    hSum = CloneEmpty(obj['data'+v+'minus'+'sub'],'hsum_'+v) #sum of MC for ratio evaluation
    hData = obj['data'+v+'minus'+'sub'].Clone('hData_'+v)

    #build stacked plot
    numLim = 6
    if NORATIO :
        numLim =5 
    for hh in range(0,numLim) :  
        obj['st'+v+'minus'+'sub'].GetHists().At(hh).SetFillStyle(1001)
        obj['st'+v+'minus'+'sub'].GetHists().At(hh).SetLineWidth(1)
        obj['st'+v+'minus'+'sub'].GetHists().At(hh).SetFillColor(colorList[hh])
        hStack.Add(obj['st'+v+'minus'+'sub'].GetHists().At(hh))
        hSum.Add(obj['st'+v+'minus'+'sub'].GetHists().At(hh))
        legend.AddEntry(obj['st'+v+'minus'+'sub'].GetHists().At(hh), sampleOrder[hh],  "f")                    
    legend.AddEntry(hData, 'Data', "PE1")
    
    #build ratio plot Data/pred
    hRatio = hData.Clone('hRatio_'+v)
    hRatio.Divide(hSum) #nominal ratio
    hRatioDict = {} #syst ratio

    #build the canvas
    can = ROOT.TCanvas('prefit_'+v,'prefit_'+v,800,700)
    can.cd()
    pad_histo = ROOT.TPad("pad_histo_"+v, "pad_histo_"+v,0,0.2,1,1)
    pad_ratio = ROOT.TPad("pad_ratio_"+v, "pad_ratio_"+v,0,0,1,0.2)
    pad_histo.SetBottomMargin(0.02)
    pad_histo.Draw()
    pad_ratio.SetTopMargin(0)
    pad_ratio.SetBottomMargin(0.25)
    pad_ratio.Draw()
    
    pad_histo.cd()
    pad_histo.SetGridx()
    pad_histo.SetGridy()
    hStack.Draw("HIST")
    hData.Draw("PE1SAME")
    legend.Draw("SAME")
    
    pad_ratio.cd()
    pad_ratio.SetGridx()
    pad_ratio.SetGridy()
    hRatio.Draw("PE1")
    
    #aesthetic features
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetNColumns(2)
    
    hStack.SetMaximum(1.7*max(hData.GetMaximum(),hStack.GetMaximum())) 
    hStack.GetYaxis().SetTitle('N events')
    hStack.GetYaxis().SetTitleOffset(1.3)
    hStack.GetXaxis().SetTitle(v)
    hStack.GetXaxis().SetTitleOffset(3)
    hStack.GetXaxis().SetLabelOffset(3)
    
    hData.SetMarkerStyle(20)
    hData.SetMarkerColor(1)
    
    hRatio.SetMarkerStyle(1)#20
    hRatio.SetLineColor(1)
    hRatio.SetLineWidth(2)
    
    hRatio.SetTitle("")
    hRatio.GetYaxis().SetTitle('Data/Pred')
    hRatio.GetYaxis().SetTitleOffset(0.25)
    hRatio.GetYaxis().SetNdivisions(506)
    hRatio.SetTitleSize(0.15,'y')
    hRatio.SetLabelSize(0.12,'y')
    hRatio.GetXaxis().SetTitle(v)
    hRatio.GetXaxis().SetTitleOffset(0.6)
    hRatio.SetTitleSize(0.18,'x')
    hRatio.SetLabelSize(0.15,'x')

    outFile.cd()
    can.Write()
    
outFile.Close()