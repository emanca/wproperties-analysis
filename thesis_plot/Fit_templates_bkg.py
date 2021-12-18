import ROOT
import math
from ROOT import gStyle
ROOT.gROOT.SetBatch(True)
gStyle.SetOptStat(0)

signList = ['plus', 'minus']
output = ROOT.TFile("Fit_templates_bkg.root", "recreate")


for s in signList :
    inFile = ROOT.TFile.Open('../../Fit/OUTPUT_1Apr_noCF_qtMax60/W'+s+'_reco.root')
    QCDtemplate = inFile.Get('Fake') 
    can = ROOT.TCanvas('c_QCDtemplate_'+s,'c_QCDtemplate_'+s,1600,1200)
    can.cd()
    can.SetGridx()
    can.SetGridy()
    can.SetRightMargin(0.15)
    QCDtemplate.GetXaxis().SetTitle('#eta^{#mu}')
    QCDtemplate.GetYaxis().SetTitle('p_{T}^{#mu}')
    QCDtemplate.GetZaxis().SetTitle('Events')
    QCDtemplate.GetZaxis().SetTitleOffset(1.4)
    if s=='plus' :
            QCDtemplate.SetTitle('QCD background template,W^{+}')
    if s=='minus' :
            QCDtemplate.SetTitle('QCD background template,W^{-}')
    QCDtemplate.Draw("colz")
    can.Write()
    can.SaveAs("Fit_QCD_template_"+s+".png")
    can.SaveAs("Fit_QCD_template_"+s+".pdf")
