import math
import ROOT
import sys
import argparse
ROOT.gROOT.SetBatch(True)

        
parser = argparse.ArgumentParser("")
parser.add_argument('-o','--output_name', type=str, default='SF_comparison',help="name of the output root file")

args = parser.parse_args()
OUTPUT = args.output_name

#input files
fileDict = {}
fileDict['Trigger'+'Plus'] = ROOT.TFile.Open('data/SF_WHelicity_Trigger_Plus.root')
fileDict['Trigger'+'Minus'] = ROOT.TFile.Open('data/SF_WHelicity_Trigger_Minus.root')
fileDict['Reco'] = ROOT.TFile.Open('data/SF_WHelicity_Reco.root')
fileDict['Trigger'+'BCDEF'] = ROOT.TFile.Open('data/SF_POG_Trigger_BCDEF.root')
fileDict['Trigger'+'GH'] = ROOT.TFile.Open('data/SF_POG_Trigger_GH.root')
fileDict['Iso'+'BCDEF'] = ROOT.TFile.Open('data/SF_POG_Iso_BCDEF.root')
fileDict['Iso'+'GH'] = ROOT.TFile.Open('data/SF_POG_Iso_GH.root')
fileDict['Id'+'BCDEF'] = ROOT.TFile.Open('data/SF_POG_Id_BCDEF.root')
fileDict['Id'+'GH'] = ROOT.TFile.Open('data/SF_POG_Id_GH.root')


outDict = {}

#----------------- GET WHelicity SF ---------------#

outDict['Trigger'+'Plus'] = fileDict['Trigger'+'Plus'].Get('scaleFactor')
outDict['Trigger'+'Minus'] = fileDict['Trigger'+'Minus'].Get('scaleFactor')

#----------------- rebin WHelicity SF ---------------#
outDict['Reco'] = fileDict['Reco'].Get('scaleFactor_etaInterpolated')
for ind, h in outDict.iteritems() :
    if 'Reco' in ind :
        h.Rebin2D(3,5)
        rebFactor=15
    else :
        h.RebinY(5)
        rebFactor=5
    for xx in range(1,h.GetXaxis().GetNbins()+1) :
        for yy in range(1,h.GetYaxis().GetNbins()+1) :
            h.SetBinContent(xx,yy,h.GetBinContent(xx,yy)/rebFactor)

#remove non useful errorbars
nEtaBins = outDict['Reco'].GetXaxis().GetNbins()
nPtBins = outDict['Reco'].GetYaxis().GetNbins()
for eta in range(1,nEtaBins+1) :
        for pt in range(1,nPtBins+1) :
            outDict['Trigger'+'Minus'].SetBinError(eta,pt,0)
            outDict['Trigger'+'Plus'].SetBinError(eta,pt,0)
            outDict['Reco'].SetBinError(eta,pt,0)

#----------------- get POG SF ---------------#
SFlist = ['Trigger','Id','Iso']
era_ratios = [0.548,0.452]
outDict['Iso'+'BCDEF']=fileDict['Iso'+'BCDEF'].Get('NUM_TightRelIso_DEN_MediumID_eta_pt')
outDict['Iso'+'GH']=fileDict['Iso'+'GH'].Get('NUM_TightRelIso_DEN_MediumID_eta_pt')
outDict['Id'+'BCDEF']=fileDict['Id'+'BCDEF'].Get('NUM_MediumID_DEN_genTracks_eta_pt')
outDict['Id'+'GH']=fileDict['Id'+'GH'].Get('NUM_MediumID_DEN_genTracks_eta_pt')
outDict['Trigger'+'BCDEF']=fileDict['Trigger'+'BCDEF'].Get('IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio')
outDict['Trigger'+'GH']=fileDict['Trigger'+'GH'].Get('IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio')
for sf in SFlist :
    outDict[sf+'POG'+'unreb'] =  outDict[sf+'BCDEF'].Clone(sf+'POG'+'unreb')
    outDict[sf+'POG'+'unreb'].Add(outDict[sf+'BCDEF'],outDict[sf+'GH'],era_ratios[0],era_ratios[1])
    
#----------------- rebin POG SF ---------------#
nEtaBins = outDict['Reco'].GetXaxis().GetNbins()
nPtBins = outDict['Reco'].GetYaxis().GetNbins()
for sf in SFlist :
    outDict[sf+'POG'] = outDict['Reco'].Clone(sf+'POG')
    for eta in range(1,nEtaBins+1) :
        for pt in range(1,nPtBins+1) :
            ptBin =  outDict[sf+'POG'+'unreb'].GetYaxis().FindBin(outDict[sf+'POG'].GetYaxis().GetBinCenter(pt))
            if sf!='Trigger' : #trigger defined for abseta
                etaBin = outDict[sf+'POG'+'unreb'].GetXaxis().FindBin(outDict[sf+'POG'].GetXaxis().GetBinCenter(eta))
            else :
                etaBin = outDict[sf+'POG'+'unreb'].GetXaxis().FindBin(abs(outDict[sf+'POG'].GetXaxis().GetBinCenter(eta)))
            val = outDict[sf+'POG'+'unreb'].GetBinContent(etaBin,ptBin)
            err = outDict[sf+'POG'+'unreb'].GetBinError(etaBin,ptBin)
            outDict[sf+'POG'].SetBinContent(eta,pt,val)
            outDict[sf+'POG'].SetBinError(eta,pt,err)
outDict['Reco'+'POG'] = outDict['Id'+'POG'].Clone('RecoPOG')
outDict['Reco'+'POG'].Multiply(outDict['Iso'+'POG'])

#----------------- ratios Whelicity/POG ---------------#

outDict['Ratio'+'Trigger'+'Plus'+'POG'+'Trigger'] = outDict['Trigger'+'Plus'].Clone('WH_triggerPlus_divided_POG_trigger')
outDict['Ratio'+'Trigger'+'Plus'+'POG'+'Trigger'].Divide(outDict['Trigger'+'POG'])

outDict['Ratio'+'Trigger'+'Minus'+'POG'+'Trigger'] = outDict['Trigger'+'Minus'].Clone('WH_triggerMinus_divided_POG_trigger')
outDict['Ratio'+'Trigger'+'Minus'+'POG'+'Trigger'].Divide(outDict['Trigger'+'POG'])

outDict['Ratio'+'Reco'+'POG'+'Reco'] = outDict['Reco'].Clone('WH_Reco_divided_POG_Reco')
outDict['Ratio'+'Reco'+'POG'+'Reco'].Divide(outDict['Reco'+'POG'])

outDict['Ratio'+'Reco'+'POG'+'Id'] = outDict['Reco'].Clone('WH_Reco_divided_POG_Id')
outDict['Ratio'+'Reco'+'POG'+'Id'].Divide(outDict['Id'+'POG'])

outDict['Ratio'+'Reco'+'POG'+'Iso'] = outDict['Reco'].Clone('WH_Reco_divided_POG_Iso')
outDict['Ratio'+'Reco'+'POG'+'Iso'].Divide(outDict['Iso'+'POG'])

#----------------- canvas with both histos ---------------#
tempDict = {}
legDict={}
for pt in range(1,nPtBins+1) :
    outDict[str(pt)+'_c_WH_triggerPlus_POG_trigger'] = ROOT.TCanvas(str(pt)+"_WH_triggerPlus_POG_trigger",str(pt)+"_WH_triggerPlus_POG_trigger",800,600)
    outDict[str(pt)+'_c_WH_triggerPlus_POG_trigger'].cd()
    tempDict[str(pt)+'_Trigger'+'Plus'] = outDict['Trigger'+'Plus'].ProjectionX(str(pt)+'_Trigger'+'Plus',pt,pt)
    tempDict[str(pt)+'_Trigger'+'POG'] = outDict['Trigger'+'POG'].ProjectionX(str(pt)+'_Trigger'+'POG',pt,pt)
    tempDict[str(pt)+'_Trigger'+'Plus'].SetLineColor(1)
    tempDict[str(pt)+'_Trigger'+'Plus'].SetMarkerColor(1)
    tempDict[str(pt)+'_Trigger'+'POG'].SetLineColor(2)
    tempDict[str(pt)+'_Trigger'+'POG'].SetMarkerColor(2)
    tempDict[str(pt)+'_Trigger'+'Plus'].DrawCopy()
    tempDict[str(pt)+'_Trigger'+'POG'].DrawCopy("SAME")
    legDict[str(pt)+'WH_triggerPlus_POG_trigger'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
    legDict[str(pt)+'WH_triggerPlus_POG_trigger'].AddEntry(tempDict[str(pt)+'_Trigger'+'Plus'],"WHelicity Trigger Plus")
    legDict[str(pt)+'WH_triggerPlus_POG_trigger'].AddEntry(tempDict[str(pt)+'_Trigger'+'POG'], "POG Trigger")
    legDict[str(pt)+'WH_triggerPlus_POG_trigger'].Draw("SAME")

    outDict[str(pt)+'_c_WH_triggerMinus_POG_trigger'] = ROOT.TCanvas(str(pt)+"_WH_triggerMinus_POG_trigger",str(pt)+"_WH_triggerMinus_POG_trigger",800,600)
    outDict[str(pt)+'_c_WH_triggerMinus_POG_trigger'].cd()
    tempDict[str(pt)+'_Trigger'+'Minus'] = outDict['Trigger'+'Minus'].ProjectionX(str(pt)+'_Trigger'+'Minus',pt,pt)
    tempDict[str(pt)+'_Trigger'+'POG'] = outDict['Trigger'+'POG'].ProjectionX(str(pt)+'_Trigger'+'POG',pt,pt)
    tempDict[str(pt)+'_Trigger'+'Minus'].SetLineColor(1)
    tempDict[str(pt)+'_Trigger'+'Minus'].SetMarkerColor(1)
    tempDict[str(pt)+'_Trigger'+'POG'].SetLineColor(2)
    tempDict[str(pt)+'_Trigger'+'POG'].SetMarkerColor(2)
    tempDict[str(pt)+'_Trigger'+'Minus'].DrawCopy()
    tempDict[str(pt)+'_Trigger'+'POG'].DrawCopy("SAME")
    legDict[str(pt)+'WH_triggerMinus_POG_trigger'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
    legDict[str(pt)+'WH_triggerMinus_POG_trigger'].AddEntry(tempDict[str(pt)+'_Trigger'+'Minus'],"WHelicity Trigger Minus")
    legDict[str(pt)+'WH_triggerMinus_POG_trigger'].AddEntry(tempDict[str(pt)+'_Trigger'+'POG'], "POG Trigger")
    legDict[str(pt)+'WH_triggerMinus_POG_trigger'].Draw("SAME")

    outDict[str(pt)+'_c_WH_Reco_divided_POG_Reco'] = ROOT.TCanvas(str(pt)+"_WH_Reco_divided_POG_Reco",str(pt)+"_WH_Reco_divided_POG_Reco",800,600)
    outDict[str(pt)+'_c_WH_Reco_divided_POG_Reco'].cd()
    tempDict[str(pt)+'_Reco'] = outDict['Reco'].ProjectionX(str(pt)+'_Reco',pt,pt)
    tempDict[str(pt)+'_Reco'+'POG']=outDict['Reco'+'POG'].ProjectionX(str(pt)+'_Reco'+'POG',pt,pt)
    tempDict[str(pt)+'_Reco'].SetLineColor(1)
    tempDict[str(pt)+'_Reco'].SetMarkerColor(1)
    tempDict[str(pt)+'_Reco'+'POG'].SetLineColor(2)
    tempDict[str(pt)+'_Reco'+'POG'].SetMarkerColor(2)
    tempDict[str(pt)+'_Reco'].DrawCopy()
    tempDict[str(pt)+'_Reco'+'POG'].DrawCopy("SAME")
    legDict[str(pt)+'WH_Reco_divided_POG_Reco'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
    legDict[str(pt)+'WH_Reco_divided_POG_Reco'].AddEntry(tempDict[str(pt)+'_Reco'],"WHelicity Reco")
    legDict[str(pt)+'WH_Reco_divided_POG_Reco'].AddEntry(tempDict[str(pt)+'_Reco'+'POG'], "POG Reco")
    legDict[str(pt)+'WH_Reco_divided_POG_Reco'].Draw("SAME")

    outDict[str(pt)+'_c_WH_Reco_divided_POG_Id'] = ROOT.TCanvas(str(pt)+"_WH_Reco_divided_POG_Id",str(pt)+"_WH_Reco_divided_POG_Id",800,600)
    outDict[str(pt)+'_c_WH_Reco_divided_POG_Id'].cd()
    tempDict[str(pt)+'_Reco'] = outDict['Reco'].ProjectionX(str(pt)+'_Reco',pt,pt)
    tempDict[str(pt)+'_Id'+'POG'] = outDict['Id'+'POG'].ProjectionX(str(pt)+'_Id'+'POG',pt,pt)
    tempDict[str(pt)+'_Reco'].SetLineColor(1)
    tempDict[str(pt)+'_Reco'].SetMarkerColor(1)
    tempDict[str(pt)+'_Id'+'POG'].SetLineColor(2)
    tempDict[str(pt)+'_Id'+'POG'].SetMarkerColor(2)
    tempDict[str(pt)+'_Reco'].DrawCopy()
    tempDict[str(pt)+'_Id'+'POG'].DrawCopy("SAME")
    legDict[str(pt)+'WH_Reco_divided_POG_Id'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
    legDict[str(pt)+'WH_Reco_divided_POG_Id'].AddEntry(tempDict[str(pt)+'_Reco'],"WHelicity Reco")
    legDict[str(pt)+'WH_Reco_divided_POG_Id'].AddEntry(tempDict[str(pt)+'_Id'+'POG'], "POG Id")
    legDict[str(pt)+'WH_Reco_divided_POG_Id'].Draw("SAME") 

    outDict[str(pt)+'_c_WH_Reco_divided_POG_Iso'] = ROOT.TCanvas(str(pt)+"_WH_Reco_divided_POG_Iso",str(pt)+"_WH_Reco_divided_POG_Iso",800,600)
    outDict[str(pt)+'_c_WH_Reco_divided_POG_Iso'].cd()
    tempDict[str(pt)+'_Reco'] = outDict['Reco'].ProjectionX(str(pt)+'_Reco',pt,pt)
    tempDict[str(pt)+'_Iso'+'POG']= outDict['Id'+'POG'].ProjectionX(str(pt)+'_Iso'+'POG',pt,pt)
    tempDict[str(pt)+'_Reco'].SetLineColor(1)
    tempDict[str(pt)+'_Reco'].SetMarkerColor(1)
    tempDict[str(pt)+'_Iso'+'POG'].SetLineColor(2)
    tempDict[str(pt)+'_Iso'+'POG'].SetMarkerColor(2)
    tempDict[str(pt)+'_Reco'].DrawCopy()
    tempDict[str(pt)+'_Iso'+'POG'].DrawCopy("SAME")
    legDict[str(pt)+'WH_Reco_divided_POG_Iso'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
    legDict[str(pt)+'WH_Reco_divided_POG_Iso'].AddEntry(tempDict[str(pt)+'_Reco'],"WHelicity Reco")
    legDict[str(pt)+'WH_Reco_divided_POG_Iso'].AddEntry(tempDict[str(pt)+'_Iso'+'POG'], "POG Iso")
    legDict[str(pt)+'WH_Reco_divided_POG_Iso'].Draw("SAME")       


        
    

#output saving
output = ROOT.TFile(OUTPUT+".root", "recreate")
output.cd()
for ind, histo in outDict.iteritems() :
    if 'Ratio' in ind : continue
    if 'GH' in ind or 'BCDEF' in ind : continue
    if 'c_' in ind: continue
    histo.SetName(ind)
    histo.Write()
    
ratio = output.mkdir('ratios')
ratio.cd()
for ind, histo in outDict.iteritems() :
    if not 'Ratio' in ind : continue
    histo.Write()

canvas = output.mkdir('canvas')
canvas.cd()
ratiolist = ['_c_WH_triggerPlus_POG_trigger','_c_WH_triggerMinus_POG_trigger','_c_WH_Reco_divided_POG_Reco','_c_WH_Reco_divided_POG_Id','_c_WH_Reco_divided_POG_Iso']
for r in ratiolist :
    for pt in range(1,nPtBins+1) :
         outDict[str(pt)+r].Write()

    # if not 'c_' in ind : continue
    
    # histo.Write()

 