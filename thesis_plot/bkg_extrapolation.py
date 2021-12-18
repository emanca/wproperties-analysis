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

inFileConst = ROOT.TFile.Open(inputDir+'/extrapolation_syst/extrapPlots_SigRegMean67GeV_constFit_CFstatAna.root')
inFileLinear = ROOT.TFile.Open(inputDir+'/extrapolation_syst/extrapPlots_SigRegMean67GeV_linearFit_CFstatAna.root')

canvasDict = {}
legDict = {}

#----------------- chi2 comparison -----------------------

chi2_const = inFileConst.Get('chi2_extrap').ProjectionX("chi2_plus_const",1,1)
chi2_linear = inFileLinear.Get('chi2_extrap').ProjectionX("chi2_plus_linear",1,1)

canvasDict['chi2'] = ROOT.TCanvas('c_chi2','c_chi2',1600,1200)
legDict['chi2'] = ROOT.TLegend(0.13,0.7,0.4,0.9)
canvasDict['chi2'].cd()
canvasDict['chi2'].SetGridx()
canvasDict['chi2'].SetGridy()

chi2_const.SetLineWidth(3)
chi2_linear.SetLineWidth(3)
chi2_const.GetXaxis().SetTitle("#eta")
chi2_const.GetYaxis().SetTitle("#chi^{2}/ndf")
chi2_const.SetTitle("Reduced #chi^{2} of f(p_{T}#times m_{T}) fit, W^{+}")
chi2_const.SetLineColor(ROOT.kBlue-7)
chi2_linear.SetLineColor(ROOT.kRed-7)
chi2_const.SetFillStyle(0)
chi2_linear.SetFillStyle(0)
chi2_const.SetStats(0)
chi2_const.Draw()
chi2_linear.Draw("SAME")

legDict['chi2'].AddEntry(chi2_const,"m_{T}-constant Fit")
legDict['chi2'].AddEntry(chi2_linear,"m_{T}-linear Fit")
legDict['chi2'].SetFillStyle(0)
legDict['chi2'].SetBorderSize(0)
legDict['chi2'].Draw("SAME")

canvasDict['chi2'].SaveAs("bkg_extrapolation_chi2.png")
canvasDict['chi2'].SaveAs("bkg_extrapolation_chi2.pdf")

#--------------- discrepnacy map -----------------

discrepancy = inFileConst.Get('discrepancyMap_Plus')

canvasDict['discrepancy'] = ROOT.TCanvas('c_discrepancy','c_discrepancy',1600,1200)
canvasDict['discrepancy'].cd()
canvasDict['discrepancy'].SetGridx()
canvasDict['discrepancy'].SetGridy()

discrepancy.SetTitle("Fake rate discrepancy extrapolated/nominal, W^{+}")
discrepancy.SetStats(0)
discrepancy.Draw("colz")

canvasDict['discrepancy'].Update() 
palette =discrepancy.GetListOfFunctions().FindObject("palette")
# palette.SetX1NDC(0.9)
palette.SetX2NDC(0.93)
# palette.SetY1NDC(0.2)
# palette.SetY2NDC(0.8)
canvasDict['discrepancy'].Modified()
canvasDict['discrepancy'].Update()

canvasDict['discrepancy'].SaveAs("bkg_extrapolation_discrepancy.png")
canvasDict['discrepancy'].SaveAs("bkg_extrapolation_discrepancy.pdf")



#output
output = ROOT.TFile("bkg_extrapolationPlots.root", "recreate")
for i, c in canvasDict.iteritems() :
    c.Write()
 