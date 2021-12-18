import ROOT
import math
from array import array
from ROOT import gStyle, gPad
ROOT.gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("2.3f")

cmslab = "#bf{CMS} #scale[0.7]{#it{Simulation Work in progress}}"
lumilab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
cmsLatex = ROOT.TLatex()


inpath = '/scratch/bertacch/wmass/wproperties-analysis/Fit/OUTPUT_25Apr_noCF_qtMax60_fixLowAcc/'
inFile = ROOT.TFile.Open(inpath+'accMap.root')


hDict = {}
cDict = {}
signList = ['Plus','Minus']
signLabel = ["W^{+}#rightarrow #mu^{+}#nu", "W^{-}#rightarrow #mu^{-}#bar{#nu}"]
for s in signList :
    hDict['acc'+s] = inFile.Get('clos'+s)
    cDict['acc'+s] = ROOT.TCanvas("acceptance"+s,"acceptance"+s,1600,1200)
    hDict['acc'+s].Draw("colz text")
    
    hDict['acc'+s].SetTitle('')
    hDict['acc'+s].GetXaxis().SetRangeUser(0,2.4)
    hDict['acc'+s].GetYaxis().SetRangeUser(0,60)
    hDict['acc'+s].GetXaxis().SetTitle('|Y_{W}|')
    hDict['acc'+s].GetYaxis().SetTitle('q_{T}^{W} [GeV]')
    hDict['acc'+s].GetZaxis().SetTitle(signLabel[signList.index(s)]+" acceptance")
    hDict['acc'+s].GetZaxis().SetTitleOffset(1.0)
    hDict['acc'+s].GetZaxis().SetTitleSize(0.05)
    hDict['acc'+s].SetContour(40)
    hDict['acc'+s].SetStats(0)
    
    cDict['acc'+s].SetRightMargin(0.15)
    gPad.Update()
    palette =    hDict['acc'+s].GetListOfFunctions().FindObject("palette")
    palette.SetX1NDC(0.875)
    # cDict['acc'+s].SetGridx()
    cDict['acc'+s].SetGridx()
    hDict['acc'+s].GetXaxis().SetNdivisions(506,ROOT.kFALSE)
    hDict['acc'+s].GetXaxis().SetTickLength(0)
    hDict['acc'+s].GetYaxis().SetTickLength(0)
    
    cmsLatex.SetNDC()
    cmsLatex.SetTextFont(42)
    cmsLatex.SetTextColor(ROOT.kBlack)
    cmsLatex.SetTextAlign(31) 
    cmsLatex.DrawLatex(1-cDict['acc'+s].GetRightMargin(),1-0.8*cDict['acc'+s].GetTopMargin(),lumilab)
    cmsLatex.SetTextAlign(11) 
    cmsLatex.DrawLatex(cDict['acc'+s].GetLeftMargin(),1-0.8*cDict['acc'+s].GetTopMargin(),cmslab)
    
    
    #horizontal lines
    hLines = []
    qtArr =[2.,  4., 6.,  8., 10., 12., 16., 20., 26., 36.]
    for qt in qtArr :
        hLines.append(ROOT.TLine(0,qt,2.4,qt))
        hLines[-1].SetLineStyle(ROOT.kDotted)
        hLines[-1].DrawLine(0,qt,2.4,qt)
        
    

outfile = ROOT.TFile('Fit_acceptance.root',"recreate")
outfile.cd()

for s in signList :
    cDict['acc'+s].Write()
    cDict['acc'+s].SaveAs("Fit_acceptance_"+s+".png")
    cDict['acc'+s].SaveAs("Fit_acceptance_"+s+".pdf")
    

