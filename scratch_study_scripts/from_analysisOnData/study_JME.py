import ROOT
import math
from array import array
import copy

ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

fileDict = {}
hDict = {}
cDict = {}
legDict = {}

path = '../../../analysisOnData/test_JMEstudy/'

fileDict['nocut'] = ROOT.TFile.Open(path+'WToMu_plots_noMtCut_finebinning.root')
fileDict['cut'] = ROOT.TFile.Open(path+'WToMu_plots_finebinning.root')

#check mt cut effect (confronto histo con e senza taglio)
for k in ['nocut','cut'] :
    
    #nominal--------
    
    h2_nom_pt = fileDict[k].Get('prefit_Signal/Nominal/Mu1_pt')
    hDict[k+'pt'+'nom'] = h2_nom_pt.ProjectionX(k+'_nom_pt_plus',2,2)
    
    h2_nom_MT = fileDict[k].Get('prefit_Signal/Nominal/MT')
    hDict[k+'MT'+'nom'] = h2_nom_MT.ProjectionX(k+'_nom_MT_plus',2,2)
    
    h3_nom = fileDict[k].Get('prefit_Signal/Nominal/PtvsMT')    
    h3_nom.GetZaxis().SetRange(2,2)#plus charge
    hDict[k+'ptMT'+'nom'] = h3_nom.Project3D('yx')
    hDict[k+'ptMT'+'nom'].SetName(k+'_nom_ptMT_plus')
    
    
    #jes up--------
    
    h2_jesUp_pt = fileDict[k].Get('prefit_Signal/jme/Mu1_pt_jesTotalUp')
    hDict[k+'pt'+'jesUp'] = h2_jesUp_pt.ProjectionX(k+'_jesUp_pt_plus',2,2)
    
    h2_jesUp_pt = fileDict[k].Get('prefit_Signal/jme/MT_jesTotalUp')
    hDict[k+'MT'+'jesUp'] = h2_jesUp_pt.ProjectionX(k+'_jesUp_MT_plus',2,2)
    
    h3_jesUp = fileDict[k].Get('prefit_Signal/jme/PtvsMT_jesTotalUp')    
    h3_jesUp.GetZaxis().SetRange(2,2)#plus charge
    hDict[k+'ptMT'+'jesUp'] = h3_jesUp.Project3D('yx')
    hDict[k+'ptMT'+'jesUp'].SetName(k+'_jesUp_ptMT_plus')
    
    #unclust en up -------
    h2_EuUp_pt = fileDict[k].Get('prefit_Signal/jme/Mu1_pt_unclustEnUp')
    hDict[k+'pt'+'EuUp'] = h2_EuUp_pt.ProjectionX(k+'_EuUp_pt_plus',2,2)
    
    h2_EuUp_pt = fileDict[k].Get('prefit_Signal/jme/MT_unclustEnUp')
    hDict[k+'MT'+'EuUp'] = h2_EuUp_pt.ProjectionX(k+'_EuUp_MT_plus',2,2)
    
    h3_EuUp = fileDict[k].Get('prefit_Signal/jme/PtvsMT_unclustEnUp')    
    h3_EuUp.GetZaxis().SetRange(2,2)#plus charge
    hDict[k+'ptMT'+'EuUp'] = h3_EuUp.Project3D('yx')
    hDict[k+'ptMT'+'EuUp'].SetName(k+'_EuUp_ptMT_plus')

    for metv in ['jesUp','EuUp'] :
        for var in ['pt', 'MT'] :
            cDict[k+metv+var] = ROOT.TCanvas(k+'_'+metv+'_'+var,k+'_'+metv+'_'+var,800,600)
            cDict[k+metv+var].cd()
            hDict[k+var+'nom'].SetLineColor(2)
            hDict[k+var+'nom'].SetLineWidth(2)
            hDict[k+var+'nom'].GetXaxis().SetTitle(var)
            hDict[k+var+'nom'].Draw()
            hDict[k+var+metv].SetLineWidth(2)
            hDict[k+var+metv].Draw("same")
            hDict[k+var+metv].SetStats(0)
            hDict[k+var+'nom'].SetStats(0)
            
            legDict[k+metv+var] = ROOT.TLegend(0.8,0.8,1,1)
            legDict[k+metv+var].AddEntry(hDict[k+var+'nom'],'nominal')
            legDict[k+metv+var].AddEntry(hDict[k+var+metv],metv)
            legDict[k+metv+var].Draw("same")
            
        for var in ['pt', 'MT', 'ptMT'] :
            hDict[k+metv+var+'nom'+'ratio'] = hDict[k+var+metv].Clone(k+'_'+metv+'_'+var+'_'+'nom'+'_'+'ratio')
            hDict[k+metv+var+'nom'+'ratio'].Divide(hDict[k+var+'nom'])
            if var!='ptMT' :
                hDict[k+metv+var+'nom'+'ratio'].GetYaxis().SetTitle(metv+'/nom')
            else :
                hDict[k+metv+var+'nom'+'ratio'].GetZaxis().SetTitle(metv+'/nom')
                hDict[k+metv+var+'nom'+'ratio'].GetXaxis().SetTitle('pT')
                hDict[k+metv+var+'nom'+'ratio'].GetYaxis().SetTitle('MT')

outFile = ROOT.TFile('study_JME.root', "recreate")
outFile.cd()
for h,hval in hDict.items() :
    hval.Write()
    
for c,cval in cDict.items() :
    cval.Write()
    




