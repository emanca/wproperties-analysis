import ROOT
import math

systDict = {
            "_LHEScaleWeight" : ["_LHEScaleWeight_muR0p5_muF0p5", "_LHEScaleWeight_muR0p5_muF1p0","_LHEScaleWeight_muR1p0_muF0p5","_LHEScaleWeight_muR1p0_muF2p0","_LHEScaleWeight_muR2p0_muF1p0", "_LHEScaleWeight_muR2p0_muF2p0"],
            "_LHEPdfWeight" : ["_LHEPdfWeightHess" + str(i)  for i in range(1, 61)],
        }
varList = ['Y','Pt']

file_norew = ROOT.TFile.Open("genInfo_norew.root")
file_rew = ROOT.TFile.Open("genInfo_rewFix.root")

histo = {}

#get histos
for v in varList :
    histo['rew'+v] = file_rew.Get('angularCoefficients/'+v)
    histo['norew'+v] = file_norew.Get('angularCoefficients/'+v)
    
    for sKind, sList in systDict.items():
        for sName in sList :
            histo['rew'+v+sName] = file_rew.Get('angularCoefficients'+sKind+'/'+v+sName)
            histo['norew'+v+sName] = file_norew.Get('angularCoefficients'+sKind+'/'+v+sName)

#build bands:
for v in varList :
    histo['rew'+v+'band'] = histo['rew'+v].Clone('rew'+v+'band')
    histo['norew'+v+'band'] = histo['norew'+v].Clone('norew'+v+'band')
    for i in range(1,histo['rew'+v+'band'].GetNbinsX()+1) :
        err = 0.
        err_norew = 0.
        for sKind, sList in systDict.items():
            if 'Pdf' in sKind :
                errPdf=0.
                errPdf_norew=0.
                for sName in sList :
                    errPdf+= (histo['rew'+v+sName].GetBinContent(i)-histo['rew'+v].GetBinContent(i))**2
                    errPdf_norew+= (histo['norew'+v+sName].GetBinContent(i)-histo['norew'+v].GetBinContent(i))**2
            if 'Scale' in sKind :
                errScale = 0.
                errScale_norew = 0.
                for sName in sList :
                    errTemp = (histo['rew'+v+sName].GetBinContent(i)-histo['rew'+v].GetBinContent(i))**2
                    errTemp_norew= (histo['norew'+v+sName].GetBinContent(i)-histo['norew'+v].GetBinContent(i))**2
                    if errTemp>errScale : 
                        errScale = errTemp
                    if errTemp_norew>errScale_norew : 
                        errScale_norew = errTemp_norew
        err = math.sqrt(errScale+errPdf)
        err_norew = math.sqrt(errScale_norew+errPdf_norew)
        histo['rew'+v+'band'].SetBinError(i,err)
        histo['norew'+v+'band'].SetBinError(i,err_norew)

canDict = {}
outFile = ROOT.TFile("study_comparison_reweighting.root", "recreate")
for v in varList :
    canDict[v] = ROOT.TCanvas('comparison'+v,'comparison'+v,800,600)
    canDict[v].cd()
    canDict[v].SetGridx()
    canDict[v].SetGridy()
    canDict['l'+v] = ROOT.TLegend(0.6,0.75,0.9,0.9)
    
    histo['rew'+v].SetLineWidth(2)
    histo['rew'+v].SetLineColor(ROOT.kBlue)
    histo['rew'+v].SetMarkerColor(ROOT.kBlue)
    histo['rew'+v+'band'].SetLineColor(ROOT.kBlue-7)
    histo['rew'+v+'band'].SetFillColor(ROOT.kBlue-7)
    histo['rew'+v+'band'].SetFillStyle(0)
    
    histo['norew'+v].SetLineWidth(2)
    histo['norew'+v].SetLineColor(ROOT.kRed)
    histo['norew'+v].SetMarkerColor(ROOT.kRed)
    histo['norew'+v+'band'].SetLineColor(ROOT.kRed-7)
    histo['norew'+v+'band'].SetFillColor(ROOT.kRed-7)
    histo['norew'+v+'band'].SetFillStyle(0)
    
    histo['norew'+v].Draw()
    histo['norew'+v+'band'].DrawCopy("E2 same")
    histo['norew'+v+'band'].SetFillStyle(3002)
    histo['norew'+v+'band'].Draw("E2 same")
    histo['rew'+v].Draw("same")
    histo['rew'+v+'band'].DrawCopy("E2 same")
    histo['rew'+v+'band'].SetFillStyle(3002)
    histo['rew'+v+'band'].Draw("E2 same")
    
    canDict['l'+v].AddEntry(histo['rew'+v],'reweighted')
    canDict['l'+v].AddEntry(histo['norew'+v],'not reweighted')
    canDict['l'+v].Draw("same")
    
    outFile.cd()
    canDict[v].Write()
    
                       
                
                

    
        
    
    


