import ROOT
import math
from ROOT import gStyle, gPad
ROOT.gROOT.SetBatch(True)
# gStyle.SetOptStat(0)
gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("2.2f")



cmslab = "#bf{CMS} #scale[0.7]{#it{Simulation Work in progress}}"
lumilab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
cmsLatex = ROOT.TLatex()

maxQt = 6
maxY=4

inFileName = '../../regularization/OUTPUT_1Apr/regularizationFit_range11____nom_nom_qt7_y5_limqt60_noC1.root'
inFile = ROOT.TFile.Open(inFileName)

optDict = {
    'A0' : [2,3],
    'A1' : [2,5],
    'A2' : [1,4],
    'A3' : [2,4],
    'A4' : [3,5] 
}

subDict = {
    'A0' : [2,3],
    'A1' : [2,3],
    'A2' : [2,3],
    'A3' : [2,3],
    'A4' : [3,3] 
}

nameDict ={
    'A0' : 'A_{0}',
    'A1' : 'A_{1}',
    'A2' : 'A_{2}',
    'A3' : 'A_{3}',
    'A4' : 'A_{4}'  
}

hin = {}
histo = {}
canv = {}

########################### chi2 plots
for c,deg in optDict.items() :
   hin['chi2'+c] = inFile.Get('h2_chi2_Plus_'+c)
   histo['chi2'+c] = ROOT.TH2F('chi2'+c,'chi2'+c,maxY,0.5,maxY+0.5,maxQt,0.5,maxQt+0.5) #needed to shift to actual maximum deg, and not "first excluded"
   for y in range(1, histo['chi2'+c].GetNbinsX()+1) :
        for q in range(1, histo['chi2'+c].GetNbinsY()+1) :
             histo['chi2'+c].SetBinContent(y,q,hin['chi2'+c].GetBinContent(y,q))
             histo['chi2'+c].SetBinError(y,q,hin['chi2'+c].GetBinError(y,q))
             
   canv['chi2'+c] = ROOT.TCanvas('regul_chi2_'+c,'regul_chi2_'+c, 1200,800)
   histo['chi2'+c].SetTitle('')
   histo['chi2'+c].GetXaxis().SetTitle('|Y_{W}| polynomial degree')
   histo['chi2'+c].GetYaxis().SetTitle('q_{T}^{W} polynomial degree')
   histo['chi2'+c].GetZaxis().SetTitle("W^{+} "+nameDict[c]+' polynomial fit #chi^{2}/Ndf' )
   histo['chi2'+c].GetXaxis().SetNdivisions(507)
   histo['chi2'+c].SetMarkerColor(ROOT.kRed-7)
   histo['chi2'+c].SetMarkerSize(1.8)
   histo['chi2'+c].Draw("colztexte")
   canv['chi2'+c].SetRightMargin(0.15)
   gPad.Update()
   palette =  histo['chi2'+c].GetListOfFunctions().FindObject("palette")
   palette.SetX1NDC(0.875)
   
   histo['box'+'opt'+c] = ROOT.TBox(deg[0]-0.5,deg[1]-0.5,deg[0]+0.5,deg[1]+0.5)
   histo['box'+'opt'+c].SetLineColor(ROOT.kRed)
   histo['box'+'opt'+c].SetLineWidth(7)
   histo['box'+'opt'+c].SetFillStyle(0)
   histo['box'+'opt'+c].Draw("same")
   
   histo['box'+'sub'+c] = ROOT.TBox(subDict[c][0]-0.5,subDict[c][1]-0.5,subDict[c][0]+0.5,subDict[c][1]+0.5)
   histo['box'+'sub'+c].SetLineColor(ROOT.kGreen)
   histo['box'+'sub'+c].SetLineWidth(4)
   histo['box'+'sub'+c].SetFillStyle(0)
   histo['box'+'sub'+c].Draw("same")
   

   cmsLatex.SetNDC()
   cmsLatex.SetTextFont(42)
   cmsLatex.SetTextColor(ROOT.kBlack)
   cmsLatex.SetTextAlign(31) 
   cmsLatex.DrawLatex(1-canv['chi2'+c].GetRightMargin(),1-0.8*canv['chi2'+c].GetTopMargin(),lumilab)
   cmsLatex.SetTextAlign(11) 
   cmsLatex.DrawLatex(canv['chi2'+c].GetLeftMargin(),1-0.8*canv['chi2'+c].GetTopMargin(),cmslab)

outfile = ROOT.TFile('Fit_regularization_plotter.root',"recreate")
outfile.cd()
for c,deg in optDict.items() :
    canv['chi2'+c].Write()
    canv['chi2'+c].SaveAs("Fit_regularization_chi2_"+c+"_plus.png")
    canv['chi2'+c].SaveAs("Fit_regularization_chi2_"+c+"_plus.pdf")


############################ pulls ###################################

for dicc, name in zip([optDict,subDict],['opt','sub']) : 

    for c,deg in dicc.items() :
        hin[name+'pull2D'+c] = inFile.Get('Y'+str(dicc[c][0]+1)+'_Qt'+str(dicc[c][1]+1)+'/hPulls2D_WtoMuP_'+c+'_'+str(dicc[c][0]+1)+str(dicc[c][1]+1))
        hin[name+'pull1D'+c] = inFile.Get('Y'+str(dicc[c][0]+1)+'_Qt'+str(dicc[c][1]+1)+'/hPulls1D_WtoMuP_'+c+'_'+str(dicc[c][0]+1)+str(dicc[c][1]+1))
        
        canv[name+'pull'+c] = ROOT.TCanvas('regul_'+name+'Pull_'+c,'regul_'+name+'Pull_'+c, 1200,800)
        canv[name+'pull'+c].SetGridx()
        canv[name+'pull'+c].SetGridy()
        hin[name+'pull2D'+c].SetTitle('')
        hin[name+'pull2D'+c].GetXaxis().SetNdivisions(6,False)
        hin[name+'pull2D'+c].GetXaxis().SetTitle('|Y_{W}| ')
        hin[name+'pull2D'+c].GetYaxis().SetTitle('q_{T}^{W} [GeV]')
        hin[name+'pull2D'+c].GetZaxis().SetTitle("W^{+}, "+nameDict[c]+' (fit-MC)/#sigma_{fit}')
        hin[name+'pull2D'+c].GetZaxis().SetRangeUser(-5,5)
        gStyle.SetPaintTextFormat("2.1f")
        canv[name+'pull'+c].SetRightMargin(0.15)
        # hin[name+'pull2D'+c].Draw("colztext")
        hin[name+'pull2D'+c].Draw("colz")
        gPad.Update()
        palette =   hin[name+'pull2D'+c].GetListOfFunctions().FindObject("palette")
        palette.SetX1NDC(0.875)
        
        cmsLatex.SetNDC()
        cmsLatex.SetTextFont(42)
        cmsLatex.SetTextColor(ROOT.kBlack)
        cmsLatex.SetTextAlign(31) 
        cmsLatex.DrawLatex(1-canv[name+'pull'+c].GetRightMargin(),1-0.8*canv[name+'pull'+c].GetTopMargin(),lumilab)
        cmsLatex.SetTextAlign(11) 
        cmsLatex.DrawLatex(canv[name+'pull'+c].GetLeftMargin(),1-0.8*canv[name+'pull'+c].GetTopMargin(),cmslab)
        
        mean = hin[name+'pull1D'+c].GetMean()
        var = hin[name+'pull1D'+c].GetStdDev()
        gussfit = hin[name+'pull1D'+c].GetFunction('gaus')
        meanfit = gussfit.GetParameter(1)
        varfit = gussfit.GetParameter(2)
        
        labFitGauss = '#scale[1.5]{#mu='+str(round(meanfit,2))+', #sigma='+str(round(varfit,2))+'}'
        labDegOrder = "#scale[1.5]{~ Y^{%s} #times q_{T}^{%s}}"%(str(dicc[c][0]),str(dicc[c][1]))
        # cmsLatex.SetTextAlign(31) 
        # cmsLatex.DrawLatex(1-canv[name+'pull'+c].GetRightMargin()-0.02,1-canv[name+'pull'+c].GetTopMargin()-0.2,labFitGauss)
        # cmsLatex.DrawLatex(1-canv[name+'pull'+c].GetRightMargin()-0.02,1-canv[name+'pull'+c].GetTopMargin()-0.1,labDegOrder)
        
        boxLabel = ROOT.TPaveText(1-canv[name+'pull'+c].GetRightMargin()-0.01, 1-canv[name+'pull'+c].GetTopMargin()-0.01,1-canv[name+'pull'+c].GetRightMargin()-0.2,1-canv[name+'pull'+c].GetTopMargin()-0.15,"NDC")
        boxLabel.AddText(labDegOrder)
        boxLabel.AddText(labFitGauss)
        boxLabel.SetFillColorAlpha(ROOT.kWhite,0.6)
        boxLabel.Draw("same")

        canv[name+'pull'+c].Write()
        canv[name+'pull'+c].SaveAs("Fit_regularization_pull_"+name+'_'+c+"_plus.png")
        canv[name+'pull'+c].SaveAs("Fit_regularization_pull_"+name+'_'+c+"_plus.pdf")
   
   
   




   
   
      



