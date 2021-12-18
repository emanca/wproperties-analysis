
#info:  ######################################################################## 
# - python ../bkg_QCD_hadder.py
# - cd ../../analysisOnData
# - python runAnalysisOnQCD.py --inputDir /scratch/bertacch/wmass/wproperties-analysis/thesis_plot/bkg_QCD_iso/ --outputDir test_bkg_QCD_iso --ncores 128
# - cd test_bkg_QCD_iso
# - hadd -f QCD.root QCD_Pt*.root
# - cd ../../thesis_plot/bkg_QCD_iso
# - python ../bkg_QCD_iso.py
####### ######################################################################## 

import ROOT
import math
from ROOT import gStyle
ROOT.gROOT.SetBatch(True)
gStyle.SetOptStat(0)

cmslab = "#bf{CMS} #scale[0.7]{#it{Simulation Work in progress}}"
lumilab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
cmsLatex = ROOT.TLatex()


inFile = ROOT.TFile.Open('../../analysisOnData/test_bkg_QCD_iso/QCD.root')

legDict = {
    'int': ['25 GeV<p_{T}^{#mu}<55 GeV, |#eta^{#mu}|<2.4',ROOT.kRed+2],
    'pt25' : ['25 GeV<p_{T}^{#mu}<35 GeV, |#eta^{#mu}|<2.4',ROOT.kBlue-4],
    'pt50' : ['45 GeV<p_{T}^{#mu}<55 GeV, |#eta^{#mu}|<2.4',ROOT.kOrange+1],
    'eta2' : ['25 GeV<p_{T}^{#mu}<55 GeV, 1.6<|#eta^{#mu}|<2.4',ROOT.kGreen+2],
    'eta0' : ['25 GeV<p_{T}^{#mu}<55 GeV, |#eta^{#mu}|<0.4',ROOT.kAzure+9]
}

fakeDict ={}

for k,info in legDict.items() :
    hiso = inFile.Get('QCDhisto_iso/Nominal/MT_'+k)
    haiso = inFile.Get('QCDhisto_aiso/Nominal/MT_'+k)
    hden = hiso.Clone('den')
    hden.Add(haiso)
    fakeDict[k] = hiso.Clone('fake_'+k)
    fakeDict[k].Divide(hden)
    fakeDict[k].GetXaxis().SetTitle("m_{T} [GeV]")
    fakeDict[k].GetYaxis().SetTitle("Isolation Cut Efficiency")
    fakeDict[k].GetYaxis().SetRangeUser(0,1.05)
    fakeDict[k].GetXaxis().SetRangeUser(0,120)
    # fakeDict[k].SetTitle("Fake Rate from QCD MC, Base-selection W^{+}")
    fakeDict[k].SetTitle("")
    fakeDict[k].SetLineWidth(2)
    fakeDict[k].SetLineColor(info[1])
    fakeDict[k].SetFillStyle(0)
    fakeDict[k].SetMarkerStyle(20)
    fakeDict[k].SetMarkerColor(info[1])
    fakeDict[k].SetMarkerSize(1.8)

#fit 
linFunc = ROOT.TF1("linFunc", '[0]+[1]*x',0,120,2)
linFunc.SetParameters(0.5,0)
linFunc.SetParNames("offset","slope")     
constFunc = ROOT.TF1("constFunc", '[0]',0,120,2)
constFunc.SetParameter(0,0.5)
constFunc.SetParNames("offset")     
linFunc.SetLineColor(ROOT.kBlack)
linFunc.SetLineWidth(2)
constFunc.SetLineWidth(2)

result_linFit = fakeDict['int'].Fit(linFunc,"S","",0,120)
result_constFit = fakeDict['int'].Fit(constFunc,"S+","",0,120)

const_chi2 = str(round(constFunc.GetChisquare()/constFunc.GetNDF(),1))
lin_chi2 = str(round(linFunc.GetChisquare()/linFunc.GetNDF(),1))
lin_slope = str(round(linFunc.GetParameter(1),4))
lin_slopeErr = str(round(linFunc.GetParError(1),4))

    
#plot
can = ROOT.TCanvas('c_fakerate','c_fakerate',1600,1200)
# leg = ROOT.TLegend(0.12,0.7,0.5,0.9)
leg = ROOT.TLegend(0.12,0.72,0.43,0.9)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
# leg.SetNColumns(2)

can.cd()
can.SetGridx()
can.SetGridy()
fakeDict['int'].Draw()
fakeDict['int'].SetLineWidth(3)
fakeDict['int'].SetMarkerColor(ROOT.kRed+3)
for k,info in legDict.items() :
    if 'int' in k : continue
    fakeDict[k].Draw("same")
    leg.AddEntry(fakeDict[k],info[0])
fakeDict['int'].Draw("same")
   
legFit = ROOT.TLegend(0.43,0.72,0.83,0.9)
legFit.SetFillStyle(0)
legFit.SetBorderSize(0)
legFit.AddEntry(fakeDict['int'],legDict['int'][0])
legFit.AddEntry(constFunc, 'constant Fit, #chi^{2}/Ndf='+const_chi2)
legFit.AddEntry(linFunc, 'linear Fit, #chi^{2}/Ndf='+lin_chi2+', slope='+lin_slope+'#pm'+lin_slopeErr)
leg.Draw("same")  
legFit.Draw("same")  

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-can.GetRightMargin(),1-0.8*can.GetTopMargin(),lumilab)

cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(can.GetLeftMargin(),1-0.8*can.GetTopMargin(),cmslab)

#output
output = ROOT.TFile("bkg_QCD_iso.root", "recreate")
can.Write()
can.SaveAs("bkg_QCD_isoMT.png")
can.SaveAs("bkg_QCD_isoMT.pdf")
