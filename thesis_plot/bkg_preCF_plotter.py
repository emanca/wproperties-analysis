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
gStyle.SetOptStat(0)

parser = argparse.ArgumentParser("")
parser.add_argument('-i','--input', type=str, default='../bkgAnalysis/test',help="name of the input direcory")

args = parser.parse_args()
inputDir = args.input

inFile = ROOT.TFile.Open(inputDir+'/final_plots.root')

canvasDict = {}
legDict = {}

cmslab = "#bf{CMS} #scale[0.7]{#it{Work in progress}}"
lumilab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
cmsLatex = ROOT.TLatex()

# -------------------------------- plot fake prompt ---------------------------

can_fr = inFile.Get('Plus_eta0/c_comparison_Plus_0')
fr = can_fr.GetPrimitive('hFakes_pt_fake_Plus_0')
frFit = can_fr.GetPrimitive('hFakes_pt_fakeFit_Plus_0')
pr = can_fr.GetPrimitive('hFakes_pt_prompt_Plus_0')
prFit = can_fr.GetPrimitive('hFakes_pt_promptFit_Plus_0')

canvasDict['fakerate'] = ROOT.TCanvas('c_fakerate','c_fakerate',1600,1200)
legDict['fakerate'] = ROOT.TLegend(0.55,0.55,0.9,0.85)
canvasDict['fakerate'].cd()
canvasDict['fakerate'].SetGridx()
canvasDict['fakerate'].SetGridy()

frFit.SetLineWidth(1)
frFit.SetLineStyle(1)
frFit.SetFillStyle(3001)
frFit.SetFillColor(632-7)
prFit.SetLineWidth(1)
prFit.SetLineStyle(1)
prFit.SetFillStyle(3001)
prFit.SetFillColor(600-9)
fr.GetFunction('fitFake').SetLineColor(632)
pr.GetFunction('fitFake').SetLineColor(600+3)
fr.GetYaxis().SetRangeUser(0,1.03)
# fr.SetTitle('Fake and Prompt Rate, 0<#eta<0.1, W^{+}')
fr.SetTitle('')
fr.SetFillStyle(0)
pr.SetFillStyle(0)

fr.Draw()
fr.SetStats(0)
frFit.Draw("same E3")
pr.Draw("same")
prFit.Draw("same E3")

legDict['fakerate'].SetFillStyle(0)
legDict['fakerate'].SetBorderSize(0)
legDict['fakerate'].AddEntry(fr,'Fake Rate')
legDict['fakerate'].AddEntry(frFit,'Fake Rate fit')
legDict['fakerate'].AddEntry(pr,'Prompt Rate')
legDict['fakerate'].AddEntry(prFit,'Prompt Rate fit')
legDict['fakerate'].Draw("same")
legDict['fakerate'].SetHeader("W^{+}#rightarrow #mu^{+}#nu_{#mu}, 0<#eta^{#mu}<0.1")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['fakerate'].GetRightMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['fakerate'].GetLeftMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),cmslab)

canvasDict['fakerate'].SaveAs("bkg_preCF_fakerate.png")
canvasDict['fakerate'].SaveAs("bkg_preCF_fakerate.pdf")



# -------------------------------plot template --------------------------------

can_template = inFile.Get('Plus_eta0/c_template_Plus_0')
qcd = can_template.GetPrimitive('htempl_pt_fake_Plus_0')
ewk = can_template.GetPrimitive('htempl_pt_prompt_Plus_0')

canvasDict['template'] = ROOT.TCanvas('c_template','c_template',1600,1200)
legDict['template'] = ROOT.TLegend(0.6,0.68,0.95,0.9)
canvasDict['template'].cd()
canvasDict['template'].SetGridx()
canvasDict['template'].SetGridy()
canvasDict['template'].SetLogy()

ewk.SetFillStyle(0)
qcd.SetFillStyle(0)
# ewk.SetTitle('Transverse momentum distribution, 0<#eta<0.1, W^{+}')
ewk.SetTitle('')
ewk.GetYaxis().SetTitle('Events / '+str(ewk.GetBinWidth(1))+' GeV')
# ewk.GetYaxis().SetTitle('Events/1 [GeV^{-1}]')


ewk.Draw()
ewk.SetStats(0)
qcd.Draw("same")
ewk.GetYaxis().SetRangeUser(100,500000)

legDict['template'].SetFillStyle(0)
legDict['template'].SetBorderSize(0)
legDict['template'].AddEntry(qcd,'QCD estimation')
legDict['template'].AddEntry(ewk,'EWK MC')
legDict['template'].Draw("same")
legDict['template'].SetHeader("W^{+}#rightarrow #mu^{+}#nu_{#mu}, 0<#eta^{#mu}<0.1")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['fakerate'].GetRightMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['fakerate'].GetLeftMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),cmslab)


canvasDict['template'].SaveAs("bkg_preCF_template.png")
canvasDict['template'].SaveAs("bkg_preCF_template.pdf")


# -------------------------- plot chi2 ----------------------------------------

chi_fr_Plus = inFile.Get('Plus_Fit_parameters/c_parameters_Plus_fake_chi2red').GetPrimitive('parameters_fake_chi2red_Plus')
chi_fr_Minus = inFile.Get('Minus_Fit_parameters/c_parameters_Minus_fake_chi2red').GetPrimitive('parameters_fake_chi2red_Minus')
chi_pr_Plus = inFile.Get('Plus_Fit_parameters/c_parameters_Plus_prompt_chi2red').GetPrimitive('parameters_prompt_chi2red_Plus')
chi_pr_Minus = inFile.Get('Minus_Fit_parameters/c_parameters_Minus_prompt_chi2red').GetPrimitive('parameters_prompt_chi2red_Minus')

canvasDict['chi2'] = ROOT.TCanvas('c_chi2','c_chi2',1600,1200)
legDict['chi2'] = ROOT.TLegend(0.6,0.7,0.95,0.88)
canvasDict['chi2'].cd()
canvasDict['chi2'].SetGridx()
canvasDict['chi2'].SetGridy()

chi_pr_Plus.SetLineColor(600-4) #blue
chi_pr_Minus.SetLineColor(860+10) #azure
chi_fr_Plus.SetLineColor(632+2) #red
chi_fr_Minus.SetLineColor(800+7) #orange
chi_pr_Plus.SetStats(0)
chi_pr_Minus.SetStats(0)
chi_fr_Plus.SetStats(0)
chi_fr_Minus.SetStats(0)
chi_pr_Plus.GetYaxis().SetTitle("#chi^{2}/Ndf")
# chi_pr_Plus.SetTitle('Reduced #chi^{2}, W^{+}')
chi_pr_Plus.SetTitle('')
chi_pr_Plus.GetXaxis().SetTitle("#eta^{#mu}")
chi_pr_Plus.SetFillStyle(0)
chi_pr_Minus.SetFillStyle(0)
chi_fr_Plus.SetFillStyle(0)
chi_fr_Minus.SetFillStyle(0)
chi_pr_Plus.GetYaxis().SetRangeUser(0,3.3)

chi_pr_Plus.Draw()
chi_pr_Plus.SetStats(0)
# chi_pr_Minus.Draw("same")
chi_fr_Plus.Draw("same")
# chi_fr_Minus.Draw("same")

legDict['chi2'].SetFillStyle(0)
legDict['chi2'].SetBorderSize(0)
legDict['chi2'].AddEntry(chi_fr_Plus,'Fake Rate')
legDict['chi2'].AddEntry(chi_pr_Plus,'Prompt Rate')
# legDict['chi2'].AddEntry(chi_fr_Minus,'Fake Rate, W^{-}')
# legDict['chi2'].AddEntry(chi_pr_Minus,'Prompt Rate, W^{-}')
legDict['chi2'].Draw("same")
legDict['chi2'].SetHeader("W^{+}#rightarrow #mu^{+}#nu_{#mu}")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['fakerate'].GetRightMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['fakerate'].GetLeftMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),cmslab)

canvasDict['chi2'].SaveAs("bkg_preCF_chi2.png")
canvasDict['chi2'].SaveAs("bkg_preCF_chi2.pdf")


#output
output = ROOT.TFile("bkg_preCF.root", "recreate")
for i, c in canvasDict.items() :
    c.Write()
 