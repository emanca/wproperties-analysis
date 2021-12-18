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
# gStyle.SetOptStat('mr')
gStyle.SetOptStat(0)

mass = 80.419


#definition of the function
def pt_vs_eta(x,par) :
    m = par[0]
    qt = par[1]
    y = par[2]
    phi = par[3]
    # return m+qt+y+phi*x[0]
    return 0.5*m/( (math.sqrt(m**2+qt**2)/m)*math.cosh(y-x[0])-(qt/m)*math.cos(phi)  )
    


funcDict = {}
#define the functions

funcDict['y0_qt6_phi0'] = ROOT.TF1('y0_qt6_phi0',pt_vs_eta,-3,3,4)
funcDict['y0_qt6_phi0'].SetParameters(mass,6.5,0,0)

funcDict['y0_qt6_phi3'] = ROOT.TF1('y0_qt6_phi3',pt_vs_eta,-3,3,4)
funcDict['y0_qt6_phi3'].SetParameters(mass,6.5,0,math.pi)


funcDict['y0_qt1_phi0'] = ROOT.TF1('y0_qt1_phi0',pt_vs_eta,-3,3,4)
funcDict['y0_qt1_phi0'].SetParameters(mass,0.5,0,0)

funcDict['y0_qt1_phi3'] = ROOT.TF1('y0_qt1_phi3',pt_vs_eta,-3,3,4)
funcDict['y0_qt1_phi3'].SetParameters(mass,0.5,0,math.pi)


funcDict['y1_qt6_phi0'] = ROOT.TF1('y1_qt6_phi0',pt_vs_eta,-3,3,4)
funcDict['y1_qt6_phi0'].SetParameters(mass,6.5,1.1,0)

funcDict['y1_qt6_phi3'] = ROOT.TF1('y1_qt6_phi3',pt_vs_eta,-3,3,4)
funcDict['y1_qt6_phi3'].SetParameters(mass,6.5,1.1,math.pi)


funcDict['ym1p7_qt19_phi0'] = ROOT.TF1('ym1p7_qt19_phi0',pt_vs_eta,-3,3,4)
funcDict['ym1p7_qt19_phi0'].SetParameters(mass,19.5,-1.6,0)

funcDict['ym1p7_qt19_phi3'] = ROOT.TF1('ym1p7_qt19_phi3',pt_vs_eta,-3,3,4)
funcDict['ym1p7_qt19_phi3'].SetParameters(mass,19.5,-1.6,math.pi)

# c = ROOT.TCanvas("c","c",800,600)
# c.cd()
# hFunc.Draw()
# c.SaveAs('test.png')
# print(hFunc.Eval(2))





can = ROOT.TCanvas('can'+'template','can'+'template',1200,900)
# leg = ROOT.TLegend(0.43,0.73,0.89,0.89)
leg = ROOT.TLegend(0.37,0.70,0.89,0.89)
can.cd()
can.SetGridx()
can.SetGridy()

funcDict['y0_qt6_phi0'].SetLineColor(ROOT.kRed-7)
funcDict['y1_qt6_phi0'].SetLineColor(ROOT.kGreen-4)
funcDict['y0_qt1_phi0'].SetLineColor(ROOT.kBlue-7)
funcDict['ym1p7_qt19_phi0'].SetLineColor(ROOT.kAzure+10)
# funcDict['y0_qt6_phi0'].SetFillColor(ROOT.kRed-4)
# funcDict['y1_qt6_phi0'].SetFillColor(ROOT.kGreen+2)
# funcDict['y0_qt1_phi0'].SetFillColor(ROOT.kBlue+1)
# funcDict['ym1p7_qt19_phi0'].SetFillColor(ROOT.kAzure+10)
# funcDict['y0_qt6_phi0'].SetFillStyle(3001)#3018
# funcDict['y1_qt6_phi0'].SetFillStyle(3001)#3017
# funcDict['y0_qt1_phi0'].SetFillStyle(1001)
# funcDict['ym1p7_qt19_phi0'].SetFillStyle(3001)

funcDict['y0_qt6_phi3'].SetLineColor(ROOT.kRed+2)
funcDict['y1_qt6_phi3'].SetLineColor(ROOT.kGreen+3)
funcDict['y0_qt1_phi3'].SetLineColor(ROOT.kBlue+2)
funcDict['ym1p7_qt19_phi3'].SetLineColor(ROOT.kAzure-8)
# funcDict['y0_qt6_phi3'].SetFillColor(ROOT.kRed-4)
# funcDict['y1_qt6_phi3'].SetFillColor(ROOT.kGreen+2)
# funcDict['y0_qt1_phi3'].SetFillColor(ROOT.kBlue+1)
# funcDict['ym1p7_qt19_phi3'].SetFillColor(ROOT.kAzure+10)
# funcDict['y0_qt6_phi3'].SetFillStyle(3001)#3018
# funcDict['y1_qt6_phi3'].SetFillStyle(3001)#3017
# funcDict['y0_qt1_phi3'].SetFillStyle(1001)
# funcDict['ym1p7_qt19_phi3'].SetFillStyle(3001)

funcDict['y0_qt6_phi0'].GetXaxis().SetTitle('Gen #eta^{#mu}')
funcDict['y0_qt6_phi0'].GetYaxis().SetTitle('Gen p_{T}^{#mu} [GeV]')
# funcDict['y0_qt6_phi0'].SetTitle('Templates ditributions at Generator level')
funcDict['y0_qt6_phi0'].SetTitle('')

# maxVal = 1
# for i,h in hFunc.items() :
#     if 'y' not in i: continue
#     maxTemp = h.GetMaximum()
#     if maxTemp>maxVal : maxVal = maxTemp 

funcDict['y0_qt6_phi0'].SetLineWidth(3)
funcDict['y0_qt6_phi3'].SetLineWidth(3)
funcDict['y0_qt6_phi0'].SetLineWidth(3)
funcDict['y0_qt6_phi3'].SetLineWidth(3)
funcDict['y1_qt6_phi0'].SetLineWidth(3)
funcDict['y1_qt6_phi3'].SetLineWidth(3)
funcDict['y1_qt6_phi0'].SetLineWidth(3)
funcDict['y1_qt6_phi3'].SetLineWidth(3)
funcDict['ym1p7_qt19_phi0'].SetLineWidth(3)
funcDict['ym1p7_qt19_phi3'].SetLineWidth(3)
funcDict['ym1p7_qt19_phi0'].SetLineWidth(3)
funcDict['ym1p7_qt19_phi3'].SetLineWidth(3)
funcDict['y0_qt1_phi0'].SetLineWidth(3)
funcDict['y0_qt1_phi3'].SetLineWidth(3)
funcDict['y0_qt1_phi0'].SetLineWidth(3)
funcDict['y0_qt1_phi3'].SetLineWidth(3)

# funcDict['y0_qt6_phi0'].SetFillSyle(0)
# funcDict['y0_qt6_phi3'].SetFillSyle(0)
# funcDict['y0_qt6_phi0'].SetFillSyle(0)
# funcDict['y0_qt6_phi3'].SetFillSyle(0)
# funcDict['y1_qt6_phi0'].SetFillSyle(0)
# funcDict['y1_qt6_phi3'].SetFillSyle(0)
# funcDict['y1_qt6_phi0'].SetFillSyle(0)
# funcDict['y1_qt6_phi3'].SetFillSyle(0)
# funcDict['ym1p7_qt19_phi0'].SetFillSyle(0)
# funcDict['ym1p7_qt19_phi3'].SetFillSyle(0)
# funcDict['ym1p7_qt19_phi0'].SetFillSyle(0)
# funcDict['ym1p7_qt19_phi3'].SetFillSyle(0)
# funcDict['y0_qt1_phi0'].SetFillSyle(0)
# funcDict['y0_qt1_phi3'].SetFillSyle(0)
# funcDict['y0_qt1_phi0'].SetFillSyle(0)
# funcDict['y0_qt1_phi3'].SetFillSyle(0)

funcDict['y0_qt6_phi0'].GetXaxis().SetRangeUser(-3,3)
funcDict['y0_qt6_phi0'].GetYaxis().SetRangeUser(0,60)
        
funcDict['y0_qt6_phi0'].Draw()
funcDict['y0_qt6_phi3'].Draw("same")
funcDict['y0_qt6_phi0'].Draw("same")
funcDict['y0_qt6_phi3'].Draw("same")
funcDict['y1_qt6_phi0'].Draw("same")
funcDict['y1_qt6_phi3'].Draw("same")
funcDict['y1_qt6_phi0'].Draw("same")
funcDict['y1_qt6_phi3'].Draw("same")
funcDict['ym1p7_qt19_phi0'].Draw("same")
funcDict['ym1p7_qt19_phi3'].Draw("same")
funcDict['ym1p7_qt19_phi0'].Draw("same")
funcDict['ym1p7_qt19_phi3'].Draw("same")
funcDict['y0_qt1_phi0'].Draw("same")
funcDict['y0_qt1_phi3'].Draw("same")
funcDict['y0_qt1_phi0'].Draw("same")
funcDict['y0_qt1_phi3'].Draw("same")

boxDict = {}
boxDict['backward'] =  ROOT.TBox(-3,0,-2.4,60)
boxDict['low'] =  ROOT.TBox(-2.4,0,2.4,25)
boxDict['forward'] =  ROOT.TBox(2.4,0,3,60)
for i, h in boxDict.items() :
    h.SetFillColorAlpha(ROOT.kGray+1,0.2)
    h.Draw()




leg.AddEntry(funcDict['y0_qt1_phi0'],     'Y_{W}=0,    q_{T}^{W}=1       GeV, #phi=0      ')
leg.AddEntry(funcDict['y0_qt1_phi3'],     '#phi=#pi ')
# leg.AddEntry(funcDict['y0_qt1_phi3'],     'Y_{W}=0,    q_{T}^{W}=1    GeV, #phi=#pi')
leg.AddEntry(funcDict['y0_qt6_phi0'],     'Y_{W}=0,    q_{T}^{W}=6.5    GeV, #phi=0      ')
leg.AddEntry(funcDict['y0_qt6_phi3'],     '#phi=#pi ')
# leg.AddEntry(funcDict['y0_qt6_phi3'],     'Y_{W}=0,    q_{T}^{W}=6.5  GeV, #phi=#pi')
leg.AddEntry(funcDict['y1_qt6_phi0'],     'Y_{W}=1.1,  q_{T}^{W}=6.5   GeV, #phi=0      ')
leg.AddEntry(funcDict['y1_qt6_phi3'],     '#phi=#pi ')
# leg.AddEntry(funcDict['y1_qt6_phi3'],     'Y_{W}=1.1,  q_{T}^{W}=6.5  GeV, #phi=#pi')
leg.AddEntry(funcDict['ym1p7_qt19_phi0'], 'Y_{W}=-1.6, q_{T}^{W}=19.5 GeV, #phi=0      ')
leg.AddEntry(funcDict['ym1p7_qt19_phi3'], '#phi=#pi ')
# leg.AddEntry(funcDict['ym1p7_qt19_phi3'], 'Y_{W}=-1.6, q_{T}^{W}=19.5 GeV, #phi=#pi')
leg.SetNColumns(2)
leg.Draw("same")



can.SaveAs('Analysis_templateFunc.pdf')
can.SaveAs('Analysis_templateFunc.png')

