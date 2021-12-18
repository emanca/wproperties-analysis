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
parser.add_argument('-i','--input', type=str, default='../../bkgAnalysis/OUTPUT_1Apr_noCF',help="name of the input direcory")

args = parser.parse_args()
inputDir = args.input

# inFile = ROOT.TFile.Open(inputDir+'/final_plots_CFstatAna.root')
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
fr_err = can_fr.GetPrimitive('hFakes_pt_fake_Plus_0_error')
frFit_err = can_fr.GetPrimitive('hFakes_pt_fakeFit_Plus_0_error')
pr_err = can_fr.GetPrimitive('hFakes_pt_prompt_Plus_0_error')
prFit_err = can_fr.GetPrimitive('hFakes_pt_promptFit_Plus_0_error')

canvasDict['fakerate'] = ROOT.TCanvas('c_fakerate','c_fakerate',1600,1200)
legDict['fakerate'] = ROOT.TLegend(0.35,0.55,0.83,0.85)
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
fr.SetTitle('')
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
legDict['fakerate'].AddEntry(frFit,'Fake Rate fit')
legDict['fakerate'].AddEntry(frFit_err,'Fake Rate fit (syst. unc.)')
legDict['fakerate'].AddEntry(pr,'Prompt Rate')
legDict['fakerate'].AddEntry(pr_err,'Prompt Rate (syst. unc.)')
legDict['fakerate'].AddEntry(prFit,'Prompt Rate fit')
legDict['fakerate'].AddEntry(prFit_err,'Prompt Rate fit (syst. unc.)')
legDict['fakerate'].Draw("same")
legDict['fakerate'].SetHeader("W^{+}#rightarrow #mu^{+}#nu_{#mu}, 0<#eta^{#mu}<0.1")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['fakerate'].GetRightMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['fakerate'].GetLeftMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),cmslab)

canvasDict['fakerate'].SaveAs("bkg_CF_fakerate.png")
canvasDict['fakerate'].SaveAs("bkg_CF_fakerate.pdf")



# -------------------------------plot template --------------------------------

can_template = inFile.Get('Plus_eta0/c_template_Plus_0')
qcd = can_template.GetPrimitive('htempl_pt_fake_Plus_0')
ewk = can_template.GetPrimitive('htempl_pt_prompt_Plus_0')
qcd_err = can_template.GetPrimitive('htempl_pt_fake_Plus_0_error')
ewk_err = can_template.GetPrimitive('htempl_pt_prompt_Plus_0_error')

canvasDict['template'] = ROOT.TCanvas('c_template','c_template',1600,1200)
legDict['template'] = ROOT.TLegend(0.6,0.65,0.89,0.9)
canvasDict['template'].cd()
canvasDict['template'].SetGridx()
canvasDict['template'].SetGridy()
canvasDict['template'].SetLogy()

ewk.SetFillStyle(0)
qcd.SetFillStyle(0)
ewk.SetTitle('Transverse momentum distribution, 0<#eta<0.1, W^{+}')
ewk.SetTitle('')
ewk.GetYaxis().SetTitle('Events / '+str(ewk.GetBinWidth(1))+' GeV')
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
legDict['template'].SetHeader("W^{+}#rightarrow #mu^{+}#nu_{#mu}, 0<#eta^{#mu}<0.1")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['fakerate'].GetRightMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['fakerate'].GetLeftMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),cmslab)

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
# chi_pr_Plus.SetTitle('Reduced #chi^{2}, W^{+}')
chi_pr_Plus.SetTitle('')
chi_pr_Plus.GetXaxis().SetTitle("#eta^{#mu}")
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
legDict['chi2'].SetHeader("W^{+}#rightarrow #mu^{+}#nu_{#mu}")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['fakerate'].GetRightMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['fakerate'].GetLeftMargin(),1-0.8*canvasDict['fakerate'].GetTopMargin(),cmslab)


canvasDict['chi2'].SaveAs("bkg_CF_chi2.png")
canvasDict['chi2'].SaveAs("bkg_CF_chi2.pdf")


mainhisto= {}
# -------------------------- grouped fake ----------------------------------------
canvasDict['group_fr'] = inFile.Get('Plus_eta0/c_groupSyst_Plus_0_comparison_fakeFit')
mainhisto['group_fr'] = canvasDict['group_fr'].GetPrimitive('comparison_hFakes_pt_fakeFit_Plus_0_sum_ratio')
mainhisto['group_fr'].SetStats(0)
# mainhisto['group_fr'].SetTitle('Systematic uncertainties shift (Fake Rate), 0<#eta<0.1, W^{+}')
mainhisto['group_fr'].SetTitle('')
mainhisto['group_fr'].GetXaxis().SetTitle("p_{T}^{#mu} [GeV]")
mainhisto['group_fr'].GetYaxis().SetTitle('|Var-Nom| / Nom (Fake Rate)')
mainhisto['group_fr'].GetYaxis().SetRangeUser(0.000002,20)

# mainhisto['group_fr'+'leg'] = canvasDict['group_fr'].GetPrimitive('TPave')
# mainhisto['group_fr'+'leg'].SetHeader("W^{+}#rightarrow #mu^{+}#nu_{#mu}, 0<#eta^{#mu}<0.1")
cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['group_fr'].GetRightMargin(),1-0.8*canvasDict['group_fr'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['group_fr'].GetLeftMargin(),1-0.8*canvasDict['group_fr'].GetTopMargin(),cmslab)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['group_fr'].GetRightMargin()-0.03,canvasDict['group_fr'].GetBottomMargin()+0.03,"W^{+}#rightarrow #mu^{+}#nu_{#mu}, 0<#eta^{#mu}<0.1")


canvasDict['group_fr'].SaveAs("bkg_CF_groupSyst_fr.png")
canvasDict['group_fr'].SaveAs("bkg_CF_groupSyst_fr.pdf")

# -------------------------- grouped prompt ----------------------------------------
canvasDict['group_pr'] = inFile.Get('Plus_eta0/c_groupSyst_Plus_0_comparison_promptFit')
mainhisto['group_pr'] = canvasDict['group_pr'].GetPrimitive('comparison_hFakes_pt_promptFit_Plus_0_sum_ratio')
mainhisto['group_pr'].SetStats(0)
# mainhisto['group_pr'].SetTitle('Systematic uncertainties shift (Prompt Rate), 0<#eta<0.1, W^{+}')
mainhisto['group_pr'].SetTitle('')
mainhisto['group_pr'].GetXaxis().SetTitle("p_{T}^{#mu} [GeV]")
mainhisto['group_pr'].GetYaxis().SetTitle('|Var-Nom| / Nom (Prompt Rate)')
mainhisto['group_pr'].GetYaxis().SetRangeUser(0.000002,20)

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['group_pr'].GetRightMargin(),1-0.8*canvasDict['group_pr'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['group_pr'].GetLeftMargin(),1-0.8*canvasDict['group_pr'].GetTopMargin(),cmslab)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['group_pr'].GetRightMargin()-0.03,canvasDict['group_pr'].GetBottomMargin()+0.03,"W^{+}#rightarrow #mu^{+}#nu_{#mu}, 0<#eta^{#mu}<0.1")

canvasDict['group_pr'].SaveAs("bkg_CF_groupSyst_pr.png")
canvasDict['group_pr'].SaveAs("bkg_CF_groupSyst_pr.pdf")

# -------------------------- grouped template ----------------------------------------
canvasDict['group_template'] = inFile.Get('Plus_eta0/c_groupSyst_Plus_0_template_fake')
mainhisto['group_template'] = canvasDict['group_template'].GetPrimitive('template_htempl_pt_fake_Plus_0_sum_ratio')
mainhisto['group_template'].SetStats(0)
# mainhisto['group_template'].SetTitle('Systematic uncertainties shift (QCD p_{T} spectrum), 0<#eta<0.1, W^{+}')
mainhisto['group_template'].SetTitle('')
mainhisto['group_template'].GetXaxis().SetTitle("p_{T}^{#mu} [GeV]")
mainhisto['group_template'].GetYaxis().SetTitle('|Var-Nom| / Nom (QCD yield)')
mainhisto['group_template'].GetYaxis().SetRangeUser(0.000002,20)

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['group_template'].GetRightMargin(),1-0.8*canvasDict['group_template'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['group_template'].GetLeftMargin(),1-0.8*canvasDict['group_template'].GetTopMargin(),cmslab)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['group_template'].GetRightMargin()-0.03,canvasDict['group_template'].GetBottomMargin()+0.04,"W^{+}#rightarrow #mu^{+}#nu_{#mu}, 0<#eta^{#mu}<0.1")

canvasDict['group_template'].SaveAs("bkg_CF_groupSyst_template.png")
canvasDict['group_template'].SaveAs("bkg_CF_groupSyst_template.pdf")

# -------------------------- grouped fake unrolled ----------------------------------------
canvasDict['group_fr_unrolled'] = inFile.Get('Plus_unrolled/c_groupSyst_unrolled_Plus_comparison_fakeFit')
mainhisto['group_fr_unrolled'] = canvasDict['group_fr_unrolled'].GetPrimitive('comparison_hFakes_pt_fakeFit_Plus_0_sum_ratio_group_unrolled')
mainhisto['group_fr_unrolled'].SetStats(0)
# mainhisto['group_fr_unrolled'].SetTitle('Systematic uncertainties shift (Fake Rate), unrolled #eta,p_{T} , W^{+}')
mainhisto['group_fr_unrolled'].SetTitle('')
mainhisto['group_fr_unrolled'].GetYaxis().SetTitle('|Var-Nom| / Nom (Fake Rate)' )
mainhisto['group_fr_unrolled'].GetYaxis().SetRangeUser(0,0.30)
mainhisto['group_fr_unrolled'].GetXaxis().SetTitle(mainhisto['group_fr_unrolled'].GetXaxis().GetTitle().replace("p_{T}","p_{T}^{#mu}"))


cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['group_fr_unrolled'].GetRightMargin(),1-0.8*canvasDict['group_fr_unrolled'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['group_fr_unrolled'].GetLeftMargin(),1-0.8*canvasDict['group_fr_unrolled'].GetTopMargin(),cmslab)
cmsLatex.DrawLatex(canvasDict['group_fr_unrolled'].GetLeftMargin()+0.03,1-canvasDict['group_fr_unrolled'].GetTopMargin()-0.04,"W^{+}#rightarrow #mu^{+}#nu_{#mu}")


canvasDict['group_fr_unrolled'].SaveAs("bkg_CF_groupSyst_fr_unrolled.png")
canvasDict['group_fr_unrolled'].SaveAs("bkg_CF_groupSyst_fr_unrolled.pdf")

# -------------------------- grouped prompt unrolled ----------------------------------------
canvasDict['group_pr_unrolled'] = inFile.Get('Plus_unrolled/c_groupSyst_unrolled_Plus_comparison_promptFit')
mainhisto['group_pr_unrolled'] = canvasDict['group_pr_unrolled'].GetPrimitive('comparison_hFakes_pt_promptFit_Plus_0_sum_ratio_group_unrolled')
mainhisto['group_pr_unrolled'].SetStats(0)
# mainhisto['group_pr_unrolled'].SetTitle('Systematic uncertainties shift (Prompt Rate), unrolled #eta,p_{T} , W^{+}')
mainhisto['group_pr_unrolled'].SetTitle('')
mainhisto['group_pr_unrolled'].GetYaxis().SetTitle('|Var-Nom| / Nom (Prompt Rate)')
mainhisto['group_pr_unrolled'].GetYaxis().SetRangeUser(0,0.0025)
mainhisto['group_pr_unrolled'].GetXaxis().SetTitle(mainhisto['group_pr_unrolled'].GetXaxis().GetTitle().replace("p_{T}","p_{T}^{#mu}"))


cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['group_pr_unrolled'].GetRightMargin(),1-0.8*canvasDict['group_pr_unrolled'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['group_pr_unrolled'].GetLeftMargin(),1-0.8*canvasDict['group_pr_unrolled'].GetTopMargin(),cmslab)
cmsLatex.DrawLatex(canvasDict['group_pr_unrolled'].GetLeftMargin()+0.03,1-canvasDict['group_pr_unrolled'].GetTopMargin()-0.04,"W^{+}#rightarrow #mu^{+}#nu_{#mu}")

canvasDict['group_pr_unrolled'].SaveAs("bkg_CF_groupSyst_pr_unrolled.png")
canvasDict['group_pr_unrolled'].SaveAs("bkg_CF_groupSyst_pr_unrolled.pdf")

# -------------------------- grouped template unrolled ----------------------------------------
canvasDict['group_template_unrolled'] = inFile.Get('Plus_unrolled/c_groupSyst_unrolled_Plus_template_fake')
mainhisto['group_template_unrolled'] = canvasDict['group_template_unrolled'].GetPrimitive('template_htempl_pt_fake_Plus_0_sum_ratio_group_unrolled')
mainhisto['group_template_unrolled'].SetStats(0)
# mainhisto['group_template_unrolled'].SetTitle('Systematic uncertainties shift (QCD spectrum), unrolled #eta,p_{T} , W^{+}')
mainhisto['group_template_unrolled'].SetTitle('')
mainhisto['group_template_unrolled'].GetYaxis().SetTitle('|Var-Nom| / Nom (QCD yield)')
mainhisto['group_template_unrolled'].GetYaxis().SetRangeUser(0,0.55)
mainhisto['group_template_unrolled'].GetXaxis().SetTitle(mainhisto['group_template_unrolled'].GetXaxis().GetTitle().replace("p_{T}","p_{T}^{#mu}"))

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['group_template_unrolled'].GetRightMargin(),1-0.8*canvasDict['group_template_unrolled'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['group_template_unrolled'].GetLeftMargin(),1-0.8*canvasDict['group_template_unrolled'].GetTopMargin(),cmslab)
cmsLatex.DrawLatex(canvasDict['group_template_unrolled'].GetLeftMargin()+0.03,1-canvasDict['group_template_unrolled'].GetTopMargin()-0.04,"W^{+}#rightarrow #mu^{+}#nu_{#mu}")

canvasDict['group_template_unrolled'].SaveAs("bkg_CF_groupSyst_template_unrolled.png")
canvasDict['group_template_unrolled'].SaveAs("bkg_CF_groupSyst_template_unrolled.pdf")



# -------------------------------- fake unrolled ---------------------------

can_fr = inFile.Get('Plus_unrolled/c_unrolled_comparison_Plus')
fr = can_fr.GetPrimitive('comparison_hFakes_pt_fake_Plus_0_unrolled')
frFit = can_fr.GetPrimitive('comparison_hFakes_pt_fakeFit_Plus_0_unrolled')
pr = can_fr.GetPrimitive('comparison_hFakes_pt_prompt_Plus_0_unrolled')
prFit = can_fr.GetPrimitive('comparison_hFakes_pt_promptFit_Plus_0_unrolled')
fr_err = can_fr.GetPrimitive('comparison_hFakes_pt_fake_Plus_0_unrolled_error')
frFit_err = can_fr.GetPrimitive('comparison_hFakes_pt_fakeFit_Plus_0_unrolled_error')
pr_err = can_fr.GetPrimitive('comparison_hFakes_pt_prompt_Plus_0_unrolled_error')
prFit_err = can_fr.GetPrimitive('comparison_hFakes_pt_promptFit_Plus_0_unrolled_error')

canvasDict['fakerate_unr'] = ROOT.TCanvas('c_fakerate_unr','c_fakerate_unr',4800,2400)
legDict['fakerate_unr'] = ROOT.TLegend(0.4,0.75,0.9,0.9)
canvasDict['fakerate_unr'].cd()
canvasDict['fakerate_unr'].SetGridx()
canvasDict['fakerate_unr'].SetGridy()
canvasDict['fakerate_unr'].SetRightMargin(0.02)
canvasDict['fakerate_unr'].SetLeftMargin(0.05)
canvasDict['fakerate_unr'].SetBottomMargin(0.25)


frFit.SetLineWidth(2)
frFit.SetLineStyle(1)
frFit.SetFillStyle(1)
frFit.SetFillColor(ROOT.kRed+4)
prFit.SetLineWidth(2)
prFit.SetLineStyle(1)
prFit.SetFillStyle(1)
prFit.SetFillColor(ROOT.kBlue-4)
# fr.GetFunction('fitFake').SetLineWidth(0)
# pr.GetFunction('fitFake').SetLineWidth(0)
fr.GetYaxis().SetRangeUser(0,1.4)
# fr.SetTitle('Fake and Prompt Rate, unrolled #eta^{#mu}(p_{T}^{#mu}), W^{+}')
# fr.SetTitle('Fake and Prompt Rate, unrolled #eta^{#mu}(p_T^{#mu}), W^{+}')
fr.SetTitle('')
fr.GetXaxis().SetTitle(fr.GetXaxis().GetTitle().replace("p_{T}","p_{T}^{#mu}"))
fr.SetFillStyle(0)
pr.SetFillStyle(0)
fr.SetLineWidth(1)
fr.SetLineColor(ROOT.kRed-6)
pr.SetLineWidth(1)
fr_err.SetFillStyle(3244)
fr_err.SetFillColor(fr.GetLineColor())
fr_err.SetLineColor(fr.GetLineColor())
pr_err.SetFillStyle(3244)
pr_err.SetFillColor(pr.GetLineColor())
pr_err.SetLineColor(pr.GetLineColor())
frFit_err.SetFillStyle(3001)
frFit_err.SetFillColor(ROOT.kRed)
frFit_err.SetLineColor(ROOT.kRed)
prFit_err.SetFillStyle(3001)
prFit_err.SetFillColor(ROOT.kBlue-7)
prFit_err.SetLineColor(ROOT.kBlue-7)

fr.GetYaxis().SetTitleOffset(0.7)

fr.Draw()
fr.SetStats(0)
fr_err.Draw("same OP5")
frFit_err.Draw("same E3L")
frFit.Draw("same E3")
pr.Draw("same")
pr_err.Draw("same OP5")
prFit.Draw("same E3")
prFit_err.Draw("same E3L")


legDict['fakerate_unr'].SetFillStyle(0)
legDict['fakerate_unr'].SetBorderSize(0)
legDict['fakerate_unr'].SetNColumns(2)
legDict['fakerate_unr'].AddEntry(fr,'Fake Rate')
legDict['fakerate_unr'].AddEntry(fr_err,'Fake Rate (syst. unc.)')
legDict['fakerate_unr'].AddEntry(frFit,'Fake Rate fitted')
legDict['fakerate_unr'].AddEntry(frFit_err,'Fake Rate fitted (syst. unc.)')
legDict['fakerate_unr'].AddEntry(pr,'Prompt Rate')
legDict['fakerate_unr'].AddEntry(pr_err,'Prompt Rate (syst. unc.)')
legDict['fakerate_unr'].AddEntry(prFit,'Prompt Rate fitted')
legDict['fakerate_unr'].AddEntry(prFit_err,'Prompt Rate fitted (syst. unc.)')
legDict['fakerate_unr'].Draw("same")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['fakerate_unr'].GetRightMargin(),1-0.8*canvasDict['fakerate_unr'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['fakerate_unr'].GetLeftMargin(),1-0.8*canvasDict['fakerate_unr'].GetTopMargin(),cmslab)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['fakerate_unr'].GetRightMargin()-0.02,1-canvasDict['fakerate_unr'].GetTopMargin()-0.04,"W^{+}#rightarrow #mu^{+}#nu_{#mu}")

canvasDict['fakerate_unr'].SaveAs("bkg_CF_fakerate_unr.png")
canvasDict['fakerate_unr'].SaveAs("bkg_CF_fakerate_unr.pdf")






# -------------------------- plot slope ----------------------------------------

slope_fr_Plus = inFile.Get('Plus_Fit_parameters/c_parameters_Plus_fake_slope').GetPrimitive('parameters_fake_slope_Plus')
slope_fr_Minus = inFile.Get('Minus_Fit_parameters/c_parameters_Minus_fake_slope').GetPrimitive('parameters_fake_slope_Minus')
slope_fr_Plus_err = inFile.Get('Plus_Fit_parameters/c_parameters_Plus_fake_slope').GetPrimitive('parameters_fake_slope_Plus_error')
slope_fr_Minus_err = inFile.Get('Minus_Fit_parameters/c_parameters_Minus_fake_slope').GetPrimitive('parameters_fake_slope_Minus_error')

canvasDict['slope'] = ROOT.TCanvas('c_slope','c_slope',1600,1200)
legDict['slope'] = ROOT.TLegend(0.22,0.72,0.8,0.9)
canvasDict['slope'].cd()
canvasDict['slope'].SetGridx()
canvasDict['slope'].SetGridy()

slope_fr_Plus.SetLineColor(ROOT.kRed+2) #red
slope_fr_Minus.SetLineColor(ROOT.kBlue-4) #orange
slope_fr_Plus.SetStats(0)
slope_fr_Minus.SetStats(0)
slope_fr_Plus.SetFillStyle(0)
slope_fr_Minus.SetFillStyle(0)
slope_fr_Plus_err.SetFillStyle(3005)
slope_fr_Plus.GetYaxis().SetTitle("Fake Rate slope")
# slope_fr_Plus.SetTitle('Slope parameter of the fake rate')
slope_fr_Plus.SetTitle('')
slope_fr_Plus.GetXaxis().SetTitle("#eta^{#mu}")
slope_fr_Plus_err.SetFillColor(slope_fr_Plus.GetLineColor())
slope_fr_Minus_err.SetFillColor(slope_fr_Minus.GetLineColor())
slope_fr_Plus_err.SetLineColor(slope_fr_Plus.GetLineColor())
slope_fr_Minus_err.SetLineColor(slope_fr_Minus.GetLineColor())
slope_fr_Plus.GetYaxis().SetRangeUser(-0.004,0.008)
slope_fr_Plus.GetYaxis().SetTitleOffset(1.4)

slope_fr_Plus.Draw("same")
slope_fr_Plus_err.Draw("same 0P5")
slope_fr_Minus.Draw("same")
slope_fr_Minus_err.Draw("same 0P5")


legDict['slope'].SetFillStyle(0)
legDict['slope'].SetBorderSize(0)
legDict['slope'].AddEntry(slope_fr_Plus,'W^{+}#rightarrow #mu^{+}#nu_{#mu}')
legDict['slope'].AddEntry(slope_fr_Plus_err,'W^{+}#rightarrow #mu^{+}#nu_{#mu} (syst. unc.)')
legDict['slope'].AddEntry(slope_fr_Minus,'W^{-}#rightarrow #mu^{-}#bar#nu_{#mu}')
legDict['slope'].AddEntry(slope_fr_Minus_err,'W^{-}#rightarrow #mu^{-}#bar#nu_{#mu} (syst. unc.)')
legDict['slope'].Draw("same")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['slope'].GetRightMargin(),1-0.8*canvasDict['slope'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['slope'].GetLeftMargin(),1-0.8*canvasDict['slope'].GetTopMargin(),cmslab)


canvasDict['slope'].SaveAs("bkg_CF_slope.png")
canvasDict['slope'].SaveAs("bkg_CF_slope.pdf")


# -------------------------- plot offset ----------------------------------------

offset_fr_Plus = inFile.Get('Plus_Fit_parameters/c_parameters_Plus_fake_offset').GetPrimitive('parameters_fake_offset_Plus')
offset_fr_Minus = inFile.Get('Minus_Fit_parameters/c_parameters_Minus_fake_offset').GetPrimitive('parameters_fake_offset_Minus')
offset_fr_Plus_err = inFile.Get('Plus_Fit_parameters/c_parameters_Plus_fake_offset').GetPrimitive('parameters_fake_offset_Plus_error')
offset_fr_Minus_err = inFile.Get('Minus_Fit_parameters/c_parameters_Minus_fake_offset').GetPrimitive('parameters_fake_offset_Minus_error')

canvasDict['offset'] = ROOT.TCanvas('c_offset','c_offset',1600,1200)
legDict['offset'] = ROOT.TLegend(0.22,0.72,0.8,0.9)
canvasDict['offset'].cd()
canvasDict['offset'].SetGridx()
canvasDict['offset'].SetGridy()

offset_fr_Plus.SetLineColor(ROOT.kRed+2) #red
offset_fr_Minus.SetLineColor(ROOT.kBlue-4) #orange
offset_fr_Plus.SetStats(0)
offset_fr_Minus.SetStats(0)
offset_fr_Plus.SetFillStyle(0)
offset_fr_Minus.SetFillStyle(0)
offset_fr_Plus_err.SetFillStyle(3005)
offset_fr_Plus.GetYaxis().SetTitle("Fake Rate offset")
# offset_fr_Plus.SetTitle('Offset parameter of the fake rate')
offset_fr_Plus.SetTitle('')
offset_fr_Plus.GetXaxis().SetTitle("#eta^{#mu}")
offset_fr_Plus_err.SetFillColor(offset_fr_Plus.GetLineColor())
offset_fr_Minus_err.SetFillColor(offset_fr_Minus.GetLineColor())
offset_fr_Plus_err.SetLineColor(offset_fr_Plus.GetLineColor())
offset_fr_Minus_err.SetLineColor(offset_fr_Minus.GetLineColor())
# offset_fr_Plus.GetYaxis().SetRangeUser(-0.004,0.008)

offset_fr_Plus.Draw("same")
offset_fr_Plus_err.Draw("same 0P5")
offset_fr_Minus.Draw("same")
offset_fr_Minus_err.Draw("same 0P5")


legDict['offset'].SetFillStyle(0)
legDict['offset'].SetBorderSize(0)
legDict['offset'].AddEntry(offset_fr_Plus,'W^{+}#rightarrow #mu^{+}#nu_{#mu}')
legDict['offset'].AddEntry(offset_fr_Plus_err,'W^{+}#rightarrow #mu^{+}#nu_{#mu} (syst. unc.)')
legDict['offset'].AddEntry(offset_fr_Minus,'W^{-}#rightarrow #mu^{-}#bar#nu_{#mu}')
legDict['offset'].AddEntry(offset_fr_Minus_err,'W^{-}#rightarrow #mu^{-}#bar#nu_{#mu} (syst. unc.)')
legDict['offset'].Draw("same")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-canvasDict['slope'].GetRightMargin(),1-0.8*canvasDict['slope'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(canvasDict['slope'].GetLeftMargin(),1-0.8*canvasDict['slope'].GetTopMargin(),cmslab)

canvasDict['offset'].SaveAs("bkg_CF_offset.png")
canvasDict['offset'].SaveAs("bkg_CF_offset.pdf")



















#output
output = ROOT.TFile("bkg_postCF.root", "recreate")
for i, c in canvasDict.items() :
    c.Write()
 