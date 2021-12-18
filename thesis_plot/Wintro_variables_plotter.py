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

parser = argparse.ArgumentParser("")
parser.add_argument('-i','--input', type=str, default='/scratchnvme/wmass/WJetsNoCUT_v2/',help="name of the input direcory")
parser.add_argument('-r','--inputReco', type=str, default='/scratchnvme/wmass/NanoAOD2016-V2/',help="name of the input direcory")
parser.add_argument('-p', '--project',type=int, default=True, help="projection from input")

args = parser.parse_args()
inputDir = args.input
PROJECT = args.project

Wcharge = {
    "plus"     :"GenPart_preFSRMuonIdx>-99 && GenPart_pdgId[GenPart_preFSRMuonIdx]<0",
    "minus"    :"GenPart_preFSRMuonIdx>-99 && GenPart_pdgId[GenPart_preFSRMuonIdx]>0"}

if PROJECT :
    # Nev = 100000000
    # Nev = 100000000

    inFile = ROOT.TFile.Open(inputDir+'/tree_0_2.root')
    tree = inFile.Get('Events')
    output = ROOT.TFile("Wintro_variables_plotter.root", "recreate")
    hDict = {}
        
    for s, sCut in Wcharge.items() :  
        print(s)
        
        #Wqt
        hDict[s+'Wqt'] = ROOT.TH1F('h_Wqt_'+s,'h_Wqt_'+s, 250,0,50)
        # tree.Project('h_Wqt_'+s,'Wpt_preFSR',sCut,"",Nev)
        tree.Project('h_Wqt_'+s,'Wpt_preFSR',sCut)
        output.cd()
        hDict[s+'Wqt'].Write()
        print("Wqt done")
        
        #template
        if s=='plus' :
            hDict[s+'template_y0_qt6'] = ROOT.TH2F('h_template_y0_qt6'+s,'h_template_y0_qt6'+s, 100,-3,3,100,0,60)
            hDict[s+'template_y1_qt6'] = ROOT.TH2F('h_template_y1_qt6'+s,'h_template_y1_qt6'+s, 100,-3,3,100,0,60)
            hDict[s+'template_y0_qt1'] = ROOT.TH2F('h_template_y0_qt1'+s,'h_template_y0_qt1'+s, 100,-3,3,100,0,60)
            hDict[s+'template_ym1p7_qt19'] = ROOT.TH2F('h_template_ym1p7_qt19'+s,'h_template_ym1p7_qt19'+s, 100,-3,3,100,0,60)
            tree.Project('h_template_y0_qt6'+s,'GenPart_pt[GenPart_preFSRMuonIdx]:GenPart_eta[GenPart_preFSRMuonIdx]',sCut+' && abs(Wrap_preFSR)<0.1 && Wpt_preFSR>6 && Wpt_preFSR<7',"")
            print("template_y0_qt6 done")
            tree.Project('h_template_y1_qt6'+s,'GenPart_pt[GenPart_preFSRMuonIdx]:GenPart_eta[GenPart_preFSRMuonIdx]',sCut+' && Wrap_preFSR>1 && Wrap_preFSR<1.2 && Wpt_preFSR>6 && Wpt_preFSR<7',"")
            print("template_y1_qt6 done")
            tree.Project('h_template_y0_qt1'+s,'GenPart_pt[GenPart_preFSRMuonIdx]:GenPart_eta[GenPart_preFSRMuonIdx]',sCut+' && abs(Wrap_preFSR)<0.1 && abs(Wpt_preFSR)<1', "")
            print("template_y0_qt1 done")
            tree.Project('h_template_ym1p7_qt19'+s,'GenPart_pt[GenPart_preFSRMuonIdx]:GenPart_eta[GenPart_preFSRMuonIdx]',sCut+' && Wrap_preFSR>-1.7 && Wrap_preFSR<-1.5 && Wpt_preFSR>19 && Wpt_preFSR<20',"")
            print("h_template_ym1p7_qt19 done")
            hDict[s+'template_y0_qt6'].Write()
            hDict[s+'template_y1_qt6'].Write()
            hDict[s+'template_y0_qt1'].Write()
            hDict[s+'template_ym1p7_qt19'].Write()
    
    
    

    output.Close()


#plotting only

cmslab = "#bf{CMS} #scale[0.7]{#it{Simulation Work in progress}}"
# lumilab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
lumilab = " #scale[0.7]{13 TeV}"
cmsLatex = ROOT.TLatex()


step2input = ROOT.TFile.Open('Wintro_variables_plotter.root')
step2output = ROOT.TFile("Wintro_variables_plotter_canvas.root", "recreate")

hDict2 = {}
for s, sCut in Wcharge.items() :  
    hDict2[s+'Wqt'] = step2input.Get('h_Wqt_'+s)
    if s=='plus' :
         hDict2[s+'template_y0_qt6'] = step2input.Get('h_template_y0_qt6'+s)
         hDict2[s+'template_y1_qt6'] = step2input.Get('h_template_y1_qt6'+s)
         hDict2[s+'template_y0_qt1'] = step2input.Get('h_template_y0_qt1'+s)
         hDict2[s+'template_ym1p7_qt19'] = step2input.Get('h_template_ym1p7_qt19'+s)


#w qt canvas
hDict2['can'+'Wqt'] = ROOT.TCanvas('can'+'Wqt','can'+'Wqt',1200,900)
hDict2['leg'+'Wqt'] = ROOT.TLegend(0.65,0.65,1,0.85)
hDict2['can'+'Wqt'].cd()
hDict2['can'+'Wqt'].SetGridx()
hDict2['can'+'Wqt'].SetGridy()
# hDict2['plus'+'Wqt'].SetMarkerColor(ROOT.kBlue+1)
# hDict2['minus'+'Wqt'].SetMarkerColor(ROOT.kRed+1)
# hDict2['plus'+'Wqt'].SetMarkerStyle(22)
# hDict2['minus'+'Wqt'].SetMarkerStyle(20)
hDict2['plus'+'Wqt'].SetLineColor(ROOT.kBlue-4)
hDict2['minus'+'Wqt'].SetLineColor(ROOT.kRed-4)
hDict2['plus'+'Wqt'].SetLineWidth(3)
hDict2['minus'+'Wqt'].SetLineWidth(3)
hDict2['plus'+'Wqt'].SetFillColor(ROOT.kBlue-4)
hDict2['minus'+'Wqt'].SetFillColor(ROOT.kRed-4)
hDict2['plus'+'Wqt'].SetFillStyle(3004)
hDict2['minus'+'Wqt'].SetFillStyle(3005)
# hDict2['plus'+'Wqt'].Sumw2()
# hDict2['minus'+'Wqt'].Sumw2()
hDict2['plus'+'Wqt'].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
# hDict2['plus'+'Wqt'].GetYaxis().SetTitle('Events/{w} GeV'.format(w=hDict2['plus'+'Wqt'].GetBinWidth(1)))
hDict2['plus'+'Wqt'].GetYaxis().SetTitle('arbitrary units')
hDict2['plus'+'Wqt'].GetYaxis().SetTitleOffset(1.5)
# hDict2['plus'+'Wqt'].Scale(1/hDict2['plus'+'Wqt'].Integral())
# hDict2['minus'+'Wqt'].Scale(1/hDict2['minus'+'Wqt'].Integral())


# hDict2['plus'+'Wqt'].SetTitle('W boson transverse momentum distribution (MC)')
hDict2['plus'+'Wqt'].SetTitle('')
hDict2['plus'+'Wqt'].Draw()
hDict2['minus'+'Wqt'].Draw("same")
hDict2['leg'+'Wqt'].AddEntry(hDict2['plus'+'Wqt'], 'W^{+}')
hDict2['leg'+'Wqt'].AddEntry(hDict2['minus'+'Wqt'], 'W^{-}')
hDict2['leg'+'Wqt'].SetLineWidth(0)
hDict2['leg'+'Wqt'].SetFillStyle(0)
hDict2['leg'+'Wqt'].Draw("same")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-hDict2['can'+'Wqt'].GetRightMargin(),1-0.8*hDict2['can'+'Wqt'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(0.07+hDict2['can'+'Wqt'].GetLeftMargin(),1-0.8*hDict2['can'+'Wqt'].GetTopMargin(),cmslab)


hDict2['can'+'Wqt'].SaveAs('W_qtMC.pdf')
hDict2['can'+'Wqt'].SaveAs('W_qtMC.png')



#template canvas
hDict2['can'+'template'] = ROOT.TCanvas('can'+'template','can'+'template',1200,900)
hDict2['leg'+'template'] = ROOT.TLegend(0.43,0.73,0.89,0.89)
hDict2['can'+'template'].cd()
hDict2['can'+'template'].SetGridx()
hDict2['can'+'template'].SetGridy()

hDict2['plus'+'template_y0_qt6'].SetLineColor(ROOT.kRed-4)
hDict2['plus'+'template_y1_qt6'].SetLineColor(ROOT.kGreen+2)
hDict2['plus'+'template_y0_qt1'].SetLineColor(ROOT.kBlue+1)
hDict2['plus'+'template_ym1p7_qt19'].SetLineColor(ROOT.kAzure+10)
hDict2['plus'+'template_y0_qt6'].SetFillColor(ROOT.kRed-4)
hDict2['plus'+'template_y1_qt6'].SetFillColor(ROOT.kGreen+2)
hDict2['plus'+'template_y0_qt1'].SetFillColor(ROOT.kBlue+1)
hDict2['plus'+'template_ym1p7_qt19'].SetFillColor(ROOT.kAzure+10)
hDict2['plus'+'template_y0_qt6'].SetFillStyle(3001)#3018
hDict2['plus'+'template_y1_qt6'].SetFillStyle(3001)#3017
hDict2['plus'+'template_y0_qt1'].SetFillStyle(1001)
hDict2['plus'+'template_ym1p7_qt19'].SetFillStyle(3001)

hDict2['plus'+'template_y0_qt6'].GetXaxis().SetTitle('Gen #eta^{#mu}')
hDict2['plus'+'template_y0_qt6'].GetYaxis().SetTitle('Gen p_{T}^{#mu} [GeV]')
# hDict2['plus'+'template_y0_qt6'].SetTitle('Templates ditributions at Generator level')
hDict2['plus'+'template_y0_qt6'].SetTitle('')
maxVal = 1
for i,h in hDict2.items() :
    if 'template_y' not in i: continue
    maxTemp = h.GetMaximum()
    if maxTemp>maxVal : maxVal = maxTemp 
        
hDict2['plus'+'template_y0_qt6'].GetZaxis().SetRangeUser(2,maxVal+1)
hDict2['plus'+'template_y0_qt6'].DrawCopy("BOX")
hDict2['plus'+'template_y0_qt6'].SetFillStyle(0)
hDict2['plus'+'template_y0_qt6'].DrawCopy("BOX SAME")
hDict2['plus'+'template_y1_qt6'].DrawCopy("BOX SAME")
hDict2['plus'+'template_y1_qt6'].SetFillStyle(0)
hDict2['plus'+'template_y1_qt6'].DrawCopy("BOX SAME")
hDict2['plus'+'template_ym1p7_qt19'].DrawCopy("BOX SAME")
hDict2['plus'+'template_ym1p7_qt19'].SetFillStyle(0)
hDict2['plus'+'template_ym1p7_qt19'].DrawCopy("BOX SAME")
hDict2['plus'+'template_y0_qt1'].DrawCopy("BOX SAME")
hDict2['plus'+'template_y0_qt1'].SetFillStyle(0)
hDict2['plus'+'template_y0_qt1'].DrawCopy("BOX SAME")
hDict2['plus'+'template_y0_qt6'].SetFillStyle(3001)
hDict2['plus'+'template_y1_qt6'].SetFillStyle(3001)
hDict2['plus'+'template_y0_qt1'].SetFillStyle(3001)
hDict2['plus'+'template_ym1p7_qt19'].SetFillStyle(3001)

boxDict = {}
boxDict['backward'] =  ROOT.TBox(-3,0,-2.4,60)
boxDict['low'] =  ROOT.TBox(-2.4,0,2.4,25)
boxDict['forward'] =  ROOT.TBox(2.4,0,3,60)
for i, h in boxDict.items() :
    h.SetFillColorAlpha(ROOT.kGray+1,0.2)
    h.Draw()

hDict2['leg'+'template'].AddEntry(hDict2['plus'+'template_y0_qt1'],     '-0.1<Y_{W}<0.1,    0 GeV <q_{T}^{W}<1   GeV')
hDict2['leg'+'template'].AddEntry(hDict2['plus'+'template_y0_qt6'],     '-0.1<Y_{W}<0.1,    6 GeV <q_{T}^{W}<7   GeV')
hDict2['leg'+'template'].AddEntry(hDict2['plus'+'template_y1_qt6'],     ' 1.0<Y_{W}<1.2,    6 GeV <q_{T}^{W}<7   GeV')
hDict2['leg'+'template'].AddEntry(hDict2['plus'+'template_ym1p7_qt19'], '-1.7<Y_{W}<-1.5, 19 GeV <q_{T}^{W}<20 GeV')
hDict2['leg'+'template'].Draw("same")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-hDict2['can'+'template'].GetRightMargin(),1-0.8*hDict2['can'+'template'].GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(hDict2['can'+'template'].GetLeftMargin(),1-0.8*hDict2['can'+'template'].GetTopMargin(),cmslab)



hDict2['can'+'template'].SaveAs('Analysis_templateGen.pdf')
hDict2['can'+'template'].SaveAs('Analysis_templateGen.png')



step2output.cd()    
for i,h in hDict2.items() :
    if 'can' in i :
        h.Write() 
step2output.Close()    
    