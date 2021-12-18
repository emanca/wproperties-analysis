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
fileInput = ROOT.TFile.Open('ScaleFactors.root')

outDict = {}

#----------------- GET WHelicity SF ---------------#

signList = ['Minus','Plus']
direcList = ['Up','Down']
systList = ['0','1','2','Flat']
for s in signList :
    outDict[s] = fileInput.Get("WHelicity/WHelicity"+s)
    for d in direcList :
        for sy in systList :
            outDict[s+sy+d] = fileInput.Get("WHelicity/WHelicity"+s+'Syst'+sy+d)

#chech if some nom is not in between syst
nEtaBins = outDict['Plus'].GetXaxis().GetNbins()
nPtBins = outDict['Plus'].GetYaxis().GetNbins()
for s in signList :
    for sy in systList :
        for eta in range(1,nEtaBins+1) :
            for pt in range(1,nPtBins+1) :
                valUp = outDict[s+sy+'Up'].GetBinContent(eta,pt)
                valDown = outDict[s+sy+'Down'].GetBinContent(eta,pt)
                val = outDict[s].GetBinContent(eta,pt)
                if val < min(valUp,valDown) or val > max(valUp,valDown) :
                    print "PROBLEM at bin (s,eta,pt)=", s,eta,pt, "for syst", sy, ", (down, val, up)=", valDown,val,valUp
                # elif val>valUp or val<valDown:
                #     print "Wrong order, at bin (s,eta,pt)=", s,eta,pt, "for syst", sy, ", (down, val, up)=", valDown,val,valUp



# #----------------- ratios Whelicity/POG ---------------#

# outDict['Ratio'+'Trigger'+'Plus'+'POG'+'Trigger'] = outDict['Trigger'+'Plus'].Clone('WH_triggerPlus_divided_POG_trigger')
# outDict['Ratio'+'Trigger'+'Plus'+'POG'+'Trigger'].Divide(outDict['Trigger'+'POG'])

# outDict['Ratio'+'Trigger'+'Minus'+'POG'+'Trigger'] = outDict['Trigger'+'Minus'].Clone('WH_triggerMinus_divided_POG_trigger')
# outDict['Ratio'+'Trigger'+'Minus'+'POG'+'Trigger'].Divide(outDict['Trigger'+'POG'])

# outDict['Ratio'+'Reco'+'POG'+'Reco'] = outDict['Reco'].Clone('WH_Reco_divided_POG_Reco')
# outDict['Ratio'+'Reco'+'POG'+'Reco'].Divide(outDict['Reco'+'POG'])

# outDict['Ratio'+'Reco'+'POG'+'Id'] = outDict['Reco'].Clone('WH_Reco_divided_POG_Id')
# outDict['Ratio'+'Reco'+'POG'+'Id'].Divide(outDict['Id'+'POG'])

# outDict['Ratio'+'Reco'+'POG'+'Iso'] = outDict['Reco'].Clone('WH_Reco_divided_POG_Iso')
# outDict['Ratio'+'Reco'+'POG'+'Iso'].Divide(outDict['Iso'+'POG'])

# #----------------- canvas with both histos ---------------#
# tempDict = {}
# legDict={}
# for pt in range(1,nPtBins+1) :
#     outDict[str(pt)+'_c_WH_triggerPlus_POG_trigger'] = ROOT.TCanvas(str(pt)+"_WH_triggerPlus_POG_trigger",str(pt)+"_WH_triggerPlus_POG_trigger",800,600)
#     outDict[str(pt)+'_c_WH_triggerPlus_POG_trigger'].cd()
#     tempDict[str(pt)+'_Trigger'+'Plus'] = outDict['Trigger'+'Plus'].ProjectionX(str(pt)+'_Trigger'+'Plus',pt,pt)
#     tempDict[str(pt)+'_Trigger'+'POG'] = outDict['Trigger'+'POG'].ProjectionX(str(pt)+'_Trigger'+'POG',pt,pt)
#     tempDict[str(pt)+'_Trigger'+'Plus'].SetLineColor(1)
#     tempDict[str(pt)+'_Trigger'+'Plus'].SetMarkerColor(1)
#     tempDict[str(pt)+'_Trigger'+'POG'].SetLineColor(2)
#     tempDict[str(pt)+'_Trigger'+'POG'].SetMarkerColor(2)
#     tempDict[str(pt)+'_Trigger'+'Plus'].DrawCopy()
#     tempDict[str(pt)+'_Trigger'+'POG'].DrawCopy("SAME")
#     legDict[str(pt)+'WH_triggerPlus_POG_trigger'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
#     legDict[str(pt)+'WH_triggerPlus_POG_trigger'].AddEntry(tempDict[str(pt)+'_Trigger'+'Plus'],"WHelicity Trigger Plus")
#     legDict[str(pt)+'WH_triggerPlus_POG_trigger'].AddEntry(tempDict[str(pt)+'_Trigger'+'POG'], "POG Trigger")
#     legDict[str(pt)+'WH_triggerPlus_POG_trigger'].Draw("SAME")

#     outDict[str(pt)+'_c_WH_triggerMinus_POG_trigger'] = ROOT.TCanvas(str(pt)+"_WH_triggerMinus_POG_trigger",str(pt)+"_WH_triggerMinus_POG_trigger",800,600)
#     outDict[str(pt)+'_c_WH_triggerMinus_POG_trigger'].cd()
#     tempDict[str(pt)+'_Trigger'+'Minus'] = outDict['Trigger'+'Minus'].ProjectionX(str(pt)+'_Trigger'+'Minus',pt,pt)
#     tempDict[str(pt)+'_Trigger'+'POG'] = outDict['Trigger'+'POG'].ProjectionX(str(pt)+'_Trigger'+'POG',pt,pt)
#     tempDict[str(pt)+'_Trigger'+'Minus'].SetLineColor(1)
#     tempDict[str(pt)+'_Trigger'+'Minus'].SetMarkerColor(1)
#     tempDict[str(pt)+'_Trigger'+'POG'].SetLineColor(2)
#     tempDict[str(pt)+'_Trigger'+'POG'].SetMarkerColor(2)
#     tempDict[str(pt)+'_Trigger'+'Minus'].DrawCopy()
#     tempDict[str(pt)+'_Trigger'+'POG'].DrawCopy("SAME")
#     legDict[str(pt)+'WH_triggerMinus_POG_trigger'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
#     legDict[str(pt)+'WH_triggerMinus_POG_trigger'].AddEntry(tempDict[str(pt)+'_Trigger'+'Minus'],"WHelicity Trigger Minus")
#     legDict[str(pt)+'WH_triggerMinus_POG_trigger'].AddEntry(tempDict[str(pt)+'_Trigger'+'POG'], "POG Trigger")
#     legDict[str(pt)+'WH_triggerMinus_POG_trigger'].Draw("SAME")

#     outDict[str(pt)+'_c_WH_Reco_divided_POG_Reco'] = ROOT.TCanvas(str(pt)+"_WH_Reco_divided_POG_Reco",str(pt)+"_WH_Reco_divided_POG_Reco",800,600)
#     outDict[str(pt)+'_c_WH_Reco_divided_POG_Reco'].cd()
#     tempDict[str(pt)+'_Reco'] = outDict['Reco'].ProjectionX(str(pt)+'_Reco',pt,pt)
#     tempDict[str(pt)+'_Reco'+'POG']=outDict['Reco'+'POG'].ProjectionX(str(pt)+'_Reco'+'POG',pt,pt)
#     tempDict[str(pt)+'_Reco'].SetLineColor(1)
#     tempDict[str(pt)+'_Reco'].SetMarkerColor(1)
#     tempDict[str(pt)+'_Reco'+'POG'].SetLineColor(2)
#     tempDict[str(pt)+'_Reco'+'POG'].SetMarkerColor(2)
#     tempDict[str(pt)+'_Reco'].DrawCopy()
#     tempDict[str(pt)+'_Reco'+'POG'].DrawCopy("SAME")
#     legDict[str(pt)+'WH_Reco_divided_POG_Reco'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
#     legDict[str(pt)+'WH_Reco_divided_POG_Reco'].AddEntry(tempDict[str(pt)+'_Reco'],"WHelicity Reco")
#     legDict[str(pt)+'WH_Reco_divided_POG_Reco'].AddEntry(tempDict[str(pt)+'_Reco'+'POG'], "POG Reco")
#     legDict[str(pt)+'WH_Reco_divided_POG_Reco'].Draw("SAME")

#     outDict[str(pt)+'_c_WH_Reco_divided_POG_Id'] = ROOT.TCanvas(str(pt)+"_WH_Reco_divided_POG_Id",str(pt)+"_WH_Reco_divided_POG_Id",800,600)
#     outDict[str(pt)+'_c_WH_Reco_divided_POG_Id'].cd()
#     tempDict[str(pt)+'_Reco'] = outDict['Reco'].ProjectionX(str(pt)+'_Reco',pt,pt)
#     tempDict[str(pt)+'_Id'+'POG'] = outDict['Id'+'POG'].ProjectionX(str(pt)+'_Id'+'POG',pt,pt)
#     tempDict[str(pt)+'_Reco'].SetLineColor(1)
#     tempDict[str(pt)+'_Reco'].SetMarkerColor(1)
#     tempDict[str(pt)+'_Id'+'POG'].SetLineColor(2)
#     tempDict[str(pt)+'_Id'+'POG'].SetMarkerColor(2)
#     tempDict[str(pt)+'_Reco'].DrawCopy()
#     tempDict[str(pt)+'_Id'+'POG'].DrawCopy("SAME")
#     legDict[str(pt)+'WH_Reco_divided_POG_Id'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
#     legDict[str(pt)+'WH_Reco_divided_POG_Id'].AddEntry(tempDict[str(pt)+'_Reco'],"WHelicity Reco")
#     legDict[str(pt)+'WH_Reco_divided_POG_Id'].AddEntry(tempDict[str(pt)+'_Id'+'POG'], "POG Id")
#     legDict[str(pt)+'WH_Reco_divided_POG_Id'].Draw("SAME") 

#     outDict[str(pt)+'_c_WH_Reco_divided_POG_Iso'] = ROOT.TCanvas(str(pt)+"_WH_Reco_divided_POG_Iso",str(pt)+"_WH_Reco_divided_POG_Iso",800,600)
#     outDict[str(pt)+'_c_WH_Reco_divided_POG_Iso'].cd()
#     tempDict[str(pt)+'_Reco'] = outDict['Reco'].ProjectionX(str(pt)+'_Reco',pt,pt)
#     tempDict[str(pt)+'_Iso'+'POG']= outDict['Id'+'POG'].ProjectionX(str(pt)+'_Iso'+'POG',pt,pt)
#     tempDict[str(pt)+'_Reco'].SetLineColor(1)
#     tempDict[str(pt)+'_Reco'].SetMarkerColor(1)
#     tempDict[str(pt)+'_Iso'+'POG'].SetLineColor(2)
#     tempDict[str(pt)+'_Iso'+'POG'].SetMarkerColor(2)
#     tempDict[str(pt)+'_Reco'].DrawCopy()
#     tempDict[str(pt)+'_Iso'+'POG'].DrawCopy("SAME")
#     legDict[str(pt)+'WH_Reco_divided_POG_Iso'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
#     legDict[str(pt)+'WH_Reco_divided_POG_Iso'].AddEntry(tempDict[str(pt)+'_Reco'],"WHelicity Reco")
#     legDict[str(pt)+'WH_Reco_divided_POG_Iso'].AddEntry(tempDict[str(pt)+'_Iso'+'POG'], "POG Iso")
#     legDict[str(pt)+'WH_Reco_divided_POG_Iso'].Draw("SAME")       


        
    

# #output saving
# output = ROOT.TFile(OUTPUT+".root", "recreate")
# output.cd()
# for ind, histo in outDict.iteritems() :
#     if 'Ratio' in ind : continue
#     if 'GH' in ind or 'BCDEF' in ind : continue
#     if 'c_' in ind: continue
#     histo.SetName(ind)
#     histo.Write()
    
# ratio = output.mkdir('ratios')
# ratio.cd()
# for ind, histo in outDict.iteritems() :
#     if not 'Ratio' in ind : continue
#     histo.Write()

# canvas = output.mkdir('canvas')
# canvas.cd()
# ratiolist = ['_c_WH_triggerPlus_POG_trigger','_c_WH_triggerMinus_POG_trigger','_c_WH_Reco_divided_POG_Reco','_c_WH_Reco_divided_POG_Id','_c_WH_Reco_divided_POG_Iso']
# for r in ratiolist :
#     for pt in range(1,nPtBins+1) :
#          outDict[str(pt)+r].Write()

#     # if not 'c_' in ind : continue
    
#     # histo.Write()

 