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

inFile = ROOT.TFile.Open(inputDir+'/final_plots_CFstatAna.root')

canvasDict = {}
legDict = {}

# -------------------------------- plot fake prompt ---------------------------

can_fr = inFile.Get('Plus_eta0/c_comparison_Plus_0')
fr = can_fr.GetPrimitive('hFakes_pt_fake_Plus_0')
frFit = can_fr.GetPrimitive('hFakes_pt_fakeFit_Plus_0')
pr = can_fr.GetPrimitive('hFakes_pt_prompt_Plus_0')
prFit = can_fr.GetPrimitive('hFakes_pt_promptFit_Plus_0')
fr_err = can_fr.GetPrimitive('hFakes_pt_fake_Plus_0_error')
frFit_err = can_fr.GetPrimitive('hFakes_pt_fakeFit_Plus_0_error')
pr_err = can_fr.GetPrimitive('hFakes_pt_prompt_Plus_0_error')
prFit_err = can_fr.GetPrimitive('hFakes_pt_promptFit_Plus_0_error')

canvasDict['fakerate'] = ROOT.TCanvas('c_fakerate','c_fakerate',1600,1200)
legDict['fakerate'] = ROOT.TLegend(0.4,0.6,0.9,0.85)
canvasDict['fakerate'].cd()
canvasDict['fakerate'].SetGridx()
canvasDict['fakerate'].SetGridy()

frFit.SetLineWidth(1)
frFit.SetLineStyle(1)
frFit.SetFillStyle(3001)
frFit.SetFillColor(ROOT.kRed-4)
prFit.SetLineWidth(1)
prFit.SetLineStyle(1)
prFit.SetFillStyle(3001)
prFit.SetFillColor(ROOT.kBlue-4)
fr.GetFunction('fitFake').SetLineWidth(0)
pr.GetFunction('fitFake').SetLineWidth(0)
fr.GetYaxis().SetRangeUser(0,1.03)
fr.SetTitle('Fake and Prompt Rate, 0<#eta<0.1, W^{+}')
fr.SetFillStyle(0)
pr.SetFillStyle(0)
fr_err.SetFillStyle(3005)
fr_err.SetFillColor(fr.GetLineColor())
fr_err.SetLineColor(fr.GetLineColor())
pr_err.SetFillStyle(3005)
pr_err.SetFillColor(pr.GetLineColor())
pr_err.SetLineColor(pr.GetLineColor())
frFit_err.SetFillStyle(3002)
frFit_err.SetFillColor(ROOT.kRed)
frFit_err.SetLineColor(ROOT.kRed)
prFit_err.SetFillStyle(3002)
prFit_err.SetFillColor(ROOT.kBlue-9)
prFit_err.SetLineColor(ROOT.kBlue-9)

fr.Draw()
fr.SetStats(0)
fr_err.Draw("same OP5")
frFit.Draw("same E3")
frFit_err.Draw("same E3L")
pr.Draw("same")
pr_err.Draw("same OP5")
prFit.Draw("same E3")
prFit_err.Draw("same E3L")


legDict['fakerate'].SetFillStyle(0)
legDict['fakerate'].SetBorderSize(0)
legDict['fakerate'].SetNColumns(2)
legDict['fakerate'].AddEntry(fr,'Fake Rate')
legDict['fakerate'].AddEntry(fr_err,'Fake Rate (syst. unc.)')
legDict['fakerate'].AddEntry(frFit,'Fake Rate fitted')
legDict['fakerate'].AddEntry(frFit_err,'Fake Rate fitted (syst. unc.)')
legDict['fakerate'].AddEntry(pr,'Prompt Rate')
legDict['fakerate'].AddEntry(pr_err,'Prompt Rate (syst. unc.)')
legDict['fakerate'].AddEntry(prFit,'Prompt Rate fitted')
legDict['fakerate'].AddEntry(prFit_err,'Prompt Rate fitted (syst. unc.)')
legDict['fakerate'].Draw("same")

canvasDict['fakerate'].SaveAs("bkg_CF_fakerate.png")
canvasDict['fakerate'].SaveAs("bkg_CF_fakerate.pdf")



# -------------------------------plot template --------------------------------

can_template = inFile.Get('Plus_eta0/c_template_Plus_0')
qcd = can_template.GetPrimitive('htempl_pt_fake_Plus_0')
ewk = can_template.GetPrimitive('htempl_pt_prompt_Plus_0')
qcd_err = can_template.GetPrimitive('htempl_pt_fake_Plus_0_error')
ewk_err = can_template.GetPrimitive('htempl_pt_prompt_Plus_0_error')

canvasDict['template'] = ROOT.TCanvas('c_template','c_template',1600,1200)
legDict['template'] = ROOT.TLegend(0.55,0.75,0.9,0.9)
canvasDict['template'].cd()
canvasDict['template'].SetGridx()
canvasDict['template'].SetGridy()
canvasDict['template'].SetLogy()

ewk.SetFillStyle(0)
qcd.SetFillStyle(0)
ewk.SetTitle('Transverse momentum distribution, 0<#eta<0.1, W^{+}')
ewk.GetYaxis().SetTitle('Events')
# ewk.GetYaxis().SetTitle('Events/1 [GeV^{-1}]')
qcd_err.SetFillStyle(3005)
ewk_err.SetFillStyle(3005)

ewk.Draw()
ewk.SetStats(0)
qcd.Draw("same")
qcd_err.Draw("same 0P5")
ewk_err.Draw("same 0P5")
ewk.GetYaxis().SetRangeUser(100,500000)

legDict['template'].SetFillStyle(0)
legDict['template'].SetBorderSize(0)
legDict['template'].AddEntry(qcd,'QCD estimation')
legDict['template'].AddEntry(qcd_err,'QCD estimation (syst. unc.)')
legDict['template'].AddEntry(ewk,'EWK MC')
legDict['template'].AddEntry(ewk_err,'EWK MC (syst. unc.)')
legDict['template'].Draw("same")

canvasDict['template'].SaveAs("bkg_CF_template.png")
canvasDict['template'].SaveAs("bkg_CF_template.pdf")


# -------------------------- plot chi2 ----------------------------------------

chi_fr_Plus = inFile.Get('Plus_Fit_parameters/c_parameters_Plus_fake_chi2red').GetPrimitive('parameters_fake_chi2red_Plus')
chi_fr_Minus = inFile.Get('Minus_Fit_parameters/c_parameters_Minus_fake_chi2red').GetPrimitive('parameters_fake_chi2red_Minus')
chi_pr_Plus = inFile.Get('Plus_Fit_parameters/c_parameters_Plus_prompt_chi2red').GetPrimitive('parameters_prompt_chi2red_Plus')
chi_pr_Minus = inFile.Get('Minus_Fit_parameters/c_parameters_Minus_prompt_chi2red').GetPrimitive('parameters_prompt_chi2red_Minus')
chi_fr_Plus_err = inFile.Get('Plus_Fit_parameters/c_parameters_Plus_fake_chi2red').GetPrimitive('parameters_fake_chi2red_Plus_error')
chi_fr_Minus_err = inFile.Get('Minus_Fit_parameters/c_parameters_Minus_fake_chi2red').GetPrimitive('parameters_fake_chi2red_Minus_error')
chi_pr_Plus_err = inFile.Get('Plus_Fit_parameters/c_parameters_Plus_prompt_chi2red').GetPrimitive('parameters_prompt_chi2red_Plus_error')
chi_pr_Minus_err = inFile.Get('Minus_Fit_parameters/c_parameters_Minus_prompt_chi2red').GetPrimitive('parameters_prompt_chi2red_Minus_error')

canvasDict['chi2'] = ROOT.TCanvas('c_chi2','c_chi2',1600,1200)
legDict['chi2'] = ROOT.TLegend(0.6,0.75,0.9,0.9)
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
chi_pr_Plus.SetTitle('Reduced #chi^{2}, W^{+}')
chi_pr_Plus.SetFillStyle(0)
chi_pr_Minus.SetFillStyle(0)
chi_fr_Plus.SetFillStyle(0)
chi_fr_Minus.SetFillStyle(0)
chi_pr_Plus_err.SetFillStyle(3005)
chi_pr_Plus_err.SetFillColor(chi_pr_Plus.GetLineColor())
chi_pr_Plus_err.SetLineColor(chi_pr_Plus.GetLineColor())
chi_fr_Plus_err.SetFillStyle(3005)
chi_fr_Plus_err.SetFillColor(chi_fr_Plus.GetLineColor())
chi_fr_Plus_err.SetLineColor(chi_fr_Plus.GetLineColor())
chi_pr_Plus.GetYaxis().SetRangeUser(0,2.7)

chi_pr_Plus.Draw()
chi_pr_Plus.SetStats(0)
chi_pr_Plus_err.Draw("same 0P5")
# chi_pr_Minus.Draw("same")
chi_fr_Plus.Draw("same")
chi_fr_Plus_err.Draw("same 0P5")
# chi_fr_Minus.Draw("same")

legDict['chi2'].SetFillStyle(0)
legDict['chi2'].SetBorderSize(0)
legDict['chi2'].AddEntry(chi_fr_Plus,'Fake Rate')
legDict['chi2'].AddEntry(chi_fr_Plus_err,'Fake Rate (syst. unc.)')
legDict['chi2'].AddEntry(chi_pr_Plus,'Prompt Rate')
legDict['chi2'].AddEntry(chi_pr_Plus_err,'Prompt Rate (syst. unc.)')
# legDict['chi2'].AddEntry(chi_fr_Minus,'Fake Rate, W^{-}')
# legDict['chi2'].AddEntry(chi_pr_Minus,'Prompt Rate, W^{-}')
legDict['chi2'].Draw("same")

canvasDict['chi2'].SaveAs("bkg_CF_chi2.png")
canvasDict['chi2'].SaveAs("bkg_CF_chi2.pdf")


mainhisto= {}
# -------------------------- grouped fake ----------------------------------------
canvasDict['group_fr'] = inFile.Get('Plus_eta0/c_groupSyst_Plus_0_comparison_fakeFit')
mainhisto['group_fr'] = canvasDict['group_fr'].GetPrimitive('comparison_hFakes_pt_fakeFit_Plus_0_sum_ratio')
mainhisto['group_fr'].SetStats(0)
mainhisto['group_fr'].SetTitle('Systematic uncertainties shift (Fake Rate), 0<#eta<0.1, W^{+}')
mainhisto['group_fr'].GetYaxis().SetTitle('|Var-Nom| / Nom')
canvasDict['group_fr'].SaveAs("bkg_CF_groupSyst_fr.png")
canvasDict['group_fr'].SaveAs("bkg_CF_groupSyst_fr.pdf")

# -------------------------- grouped prompt ----------------------------------------
canvasDict['group_pr'] = inFile.Get('Plus_eta0/c_groupSyst_Plus_0_comparison_promptFit')
mainhisto['group_pr'] = canvasDict['group_pr'].GetPrimitive('comparison_hFakes_pt_promptFit_Plus_0_sum_ratio')
mainhisto['group_pr'].SetStats(0)
mainhisto['group_pr'].SetTitle('Systematic uncertainties shift (Prompt Rate), 0<#eta<0.1, W^{+}')
mainhisto['group_pr'].GetYaxis().SetTitle('|Var-Nom| / Nom')
canvasDict['group_pr'].SaveAs("bkg_CF_groupSyst_pr.png")
canvasDict['group_pr'].SaveAs("bkg_CF_groupSyst_pr.pdf")

# -------------------------- grouped template ----------------------------------------
canvasDict['group_template'] = inFile.Get('Plus_eta0/c_groupSyst_Plus_0_template_fake')
mainhisto['group_template'] = canvasDict['group_template'].GetPrimitive('template_htempl_pt_fake_Plus_0_sum_ratio')
mainhisto['group_template'].SetStats(0)
mainhisto['group_template'].SetTitle('Systematic uncertainties shift (QCD p_{T} spectrum), 0<#eta<0.1, W^{+}')
mainhisto['group_template'].GetYaxis().SetTitle('|Var-Nom| / Nom')
canvasDict['group_template'].SaveAs("bkg_CF_groupSyst_template.png")
canvasDict['group_template'].SaveAs("bkg_CF_groupSyst_template.pdf")

# -------------------------- grouped fake unrolled ----------------------------------------
canvasDict['group_fr_unrolled'] = inFile.Get('Plus_unrolled/c_groupSyst_unrolled_Plus_comparison_fakeFit')
mainhisto['group_fr_unrolled'] = canvasDict['group_fr_unrolled'].GetPrimitive('comparison_hFakes_pt_fakeFit_Plus_0_sum_ratio_group_unrolled')
mainhisto['group_fr_unrolled'].SetStats(0)
mainhisto['group_fr_unrolled'].SetTitle('Systematic uncertainties shift (Fake Rate), unrolled #eta,p_{T} , W^{+}')
mainhisto['group_fr_unrolled'].GetYaxis().SetTitle('|Var-Nom| / Nom')
mainhisto['group_fr_unrolled'].GetYaxis().SetRangeUser(0,0.22)
canvasDict['group_fr_unrolled'].SaveAs("bkg_CF_groupSyst_fr_unrolled.png")
canvasDict['group_fr_unrolled'].SaveAs("bkg_CF_groupSyst_fr_unrolled.pdf")

# -------------------------- grouped prompt unrolled ----------------------------------------
canvasDict['group_pr_unrolled'] = inFile.Get('Plus_unrolled/c_groupSyst_unrolled_Plus_comparison_promptFit')
mainhisto['group_pr_unrolled'] = canvasDict['group_pr_unrolled'].GetPrimitive('comparison_hFakes_pt_promptFit_Plus_0_sum_ratio_group_unrolled')
mainhisto['group_pr_unrolled'].SetStats(0)
mainhisto['group_pr_unrolled'].SetTitle('Systematic uncertainties shift (Prompt Rate), unrolled #eta,p_{T} , W^{+}')
mainhisto['group_pr_unrolled'].GetYaxis().SetTitle('|Var-Nom| / Nom')
mainhisto['group_pr_unrolled'].GetYaxis().SetRangeUser(0,0.0025)
canvasDict['group_pr_unrolled'].SaveAs("bkg_CF_groupSyst_pr_unrolled.png")
canvasDict['group_pr_unrolled'].SaveAs("bkg_CF_groupSyst_pr_unrolled.pdf")

# -------------------------- grouped template unrolled ----------------------------------------
canvasDict['group_template_unrolled'] = inFile.Get('Plus_unrolled/c_groupSyst_unrolled_Plus_template_fake')
mainhisto['group_template_unrolled'] = canvasDict['group_template_unrolled'].GetPrimitive('template_htempl_pt_fake_Plus_0_sum_ratio_group_unrolled')
mainhisto['group_template_unrolled'].SetStats(0)
mainhisto['group_template_unrolled'].SetTitle('Systematic uncertainties shift (QCD spectrum), unrolled #eta,p_{T} , W^{+}')
mainhisto['group_template_unrolled'].GetYaxis().SetTitle('|Var-Nom| / Nom')
mainhisto['group_template_unrolled'].GetYaxis().SetRangeUser(0,0.4)
canvasDict['group_template_unrolled'].SaveAs("bkg_CF_groupSyst_template_unrolled.png")
canvasDict['group_template_unrolled'].SaveAs("bkg_CF_groupSyst_template_unrolled.pdf")



#output
output = ROOT.TFile("bkg_postCF.root", "recreate")
for i, c in canvasDict.items() :
    c.Write()
 