import os
import ROOT
import copy
import sys
import argparse
import math

sys.path.append('../../bkgAnalysis')
import bkg_utils
ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)

class plotter:
    
    def __init__(self, outDir, inDir = ''):

        self.indir = inDir # indir containig the various outputs
        self.outdir = outDir
        
        self.sampleDict = {
            'WToMu' : ['WToMu_plots.root', 'templates_Signal', 'W^{+}#rightarrow #mu^{+}#nu_{#mu}'],
            'SIGNAL_Fake' :['FakeFromData_plots.root','templates_fakes', 'QCD']
        }
        
        self.extSyst = copy.deepcopy(bkg_utils.bkg_systematics)
        self.extSyst['Nominal'] = ['']
        
        self.signDict = {
            'minus' : [1,  "W^{-}#rightarrow #mu^{-}#bar{#nu}_{#mu}"],
            'plus' :  [2, "W^{+}#rightarrow #mu^{+}#nu_{#mu}"]
        }
        
        self.histoDict ={} 
        
        self.varName = 'template'
    
    def getHistos(self) :
        for f,fileInfo in self.sampleDict.items() :
            inFile = ROOT.TFile.Open(self.indir+'/'+fileInfo[0])
            # inFile = ROOT.TFile.Open(self.outdir+'/hadded/'+fileInfo[0])
            for sKind, sList in self.extSyst.items():
                varname = self.varName
                if f=='WToMu' : #patch 
                    if sKind=='jmeVars' or sKind=='ptScaleVars' :
                        varname = varname+'s'
                for sName in sList :
                    inFile.cd()
                    if ROOT.gDirectory.Get(fileInfo[1]+'/'+sKind+'/'+varname+'_'+sName)==None : #this syst is missing -. take the nominal
                        if sName!='' or f!='Data': print("missing syst:", sName, " for file", f)
                        h2 = inFile.Get(fileInfo[1]+'/Nominal/template')
                    else : 
                        h2 = inFile.Get(fileInfo[1]+'/'+sKind+'/'+varname+'_'+sName)
                    for s,sInfo in self.signDict.items() :
                        h2.GetZaxis().SetRange(sInfo[0],sInfo[0])
                        self.histoDict[f+s+sName] = h2.Project3D(f+'_'+s+'_'+sName+'_yx')
                        # self.varBinWidth_modifier(self.histoDict[f+s+var+sName])
        
    def plotEtaPtInt(self, skipSyst) :
        self.getHistos()
        fname = "{dir}/EtaPtPlots.root".format(dir=self.outdir)
        outFile =  ROOT.TFile(fname, "RECREATE")
        
        for f,fileInfo in self.sampleDict.items() :
            for s,sInfo in self.signDict.items() :
                for sKind, sList in bkg_utils.bkg_systematics.items():
                    if f=='SIGNAL_Fake' and sKind=='LHEPdfWeightVars' : continue #PATCH
                    if sKind in skipSyst : continue #skipped systs
                    for sName in sList :
                        hRatio = self.histoDict[f+s+sName].Clone(f+'_'+s+'_'+sName+'_Syst_over_Nom')
                        hRatio.Divide(self.histoDict[f+s])
                        hRatio.SetTitle(f+'_'+s+'_'+sName+'_Syst/Nom')
                        
                        can = ROOT.TCanvas(hRatio.GetName(),hRatio.GetName(),800,600)
                        can.cd()
                        hRatio.GetXaxis().SetTitle('muon #eta')
                        hRatio.GetYaxis().SetTitle('muon p_{T} [GeV]')
                        hRatio.Draw("colz")
                        can.Update() 
                        
                        palette =hRatio.GetListOfFunctions().FindObject("palette")
                        palette.SetX1NDC(0.905)
                        palette.SetX2NDC(0.925)
                        can.Modified()
                        can.Update()
                        
                        can.SaveAs("{dir}/{c}.pdf".format(dir=self.outdir,c=can.GetName()))
                        can.SaveAs("{dir}/{c}.png".format(dir=self.outdir,c=can.GetName()))
                        
                        outFile.cd()
                        can.Write()
            
                #build the band/nom plots
                hBandRatio = self.histoDict[f+s+sName].Clone(f+'_'+s+'_Band_over_Nom')
                for xx in range(1,hBandRatio.GetNbinsX()+1) :
                    for yy in range(1,hBandRatio.GetNbinsY()+1) :
                        delta=0
                        for sKind, sList in bkg_utils.bkg_systematics.items():
                            if sKind in skipSyst : continue #skipped systs
                            if f=='SIGNAL_Fake' and sKind=='LHEPdfWeightVars' : continue #PATCH
                            for sName in sList :
                                if 'Down' in sName : continue
                                if 'LHE' in sName : continue 
                                if 'Up' in sName :
                                    systDown =  sName.replace("Up","Down")
                                    delta += ((self.histoDict[f+s+sName].GetBinContent(xx,yy)-self.histoDict[f+s+systDown].GetBinContent(xx,yy))/self.histoDict[f+s].GetBinContent(xx,yy))**2
                        delta = 0.5*math.sqrt(delta)
                        
                        deltaLHE=0 #LHES variations
                        for sKind, sList in bkg_utils.bkg_systematics.items():
                            if sKind in skipSyst : continue #skipped systs
                            if f=='SIGNAL_Fake' and sKind=='LHEPdfWeightVars' : continue #PATCH
                            for sName in sList :
                                if not 'LHE' in sName: continue 
                                deltaLHE+= (self.histoDict[f+s+sName].GetBinContent(xx,yy)/self.histoDict[f+s].GetBinContent(xx,yy)-1)**2
                        deltaLHE = math.sqrt(deltaLHE)
                        delta= delta+deltaLHE
                                            
                        hBandRatio.SetBinContent(xx,yy, delta)
                
                hBandRatio.SetTitle(f+'_'+s+'_Band/Nom')
                canBand = ROOT.TCanvas(hBandRatio.GetName(),hBandRatio.GetName(),800,600)
                canBand.cd()
                hBandRatio.GetXaxis().SetTitle('muon #eta')
                hBandRatio.GetYaxis().SetTitle('muon p_{T} [GeV]')
                hBandRatio.Draw("colz")
                canBand.Update() 
                
                palette =hBandRatio.GetListOfFunctions().FindObject("palette")
                palette.SetX1NDC(0.905)
                palette.SetX2NDC(0.925)
                canBand.Modified()
                canBand.Update()
                
                print("BAND NAME", canBand.GetName())
                canBand.SaveAs("{dir}/{c}.pdf".format(dir=self.outdir,c=canBand.GetName()))
                canBand.SaveAs("{dir}/{c}.png".format(dir=self.outdir,c=canBand.GetName()))
                
                outFile.cd()
                can.Write()
                
                    
        outFile.Close()
                    

print("FEATUREs NOT IMPLEMENTED:")
print("1) hadd")
print("2) comparison signal - bkg")

parser = argparse.ArgumentParser("")
parser.add_argument('-hadd','--hadd', type=int, default=True,help="hadd of the output of RDF")
parser.add_argument('-o','--output', type=str, default='TEST',help="name of the output directory")
parser.add_argument('-i','--input', type=str, default='TEST',help="name of the input direcory root file")
parser.add_argument('-s','--skipSyst', type=str, default='',nargs='*', help="list of skipped syst class as in bkgAnalysis/bkg_utils.py, separated by space")

args = parser.parse_args()
HADD = args.hadd
OUTPUT = args.output
INPUT = args.input
skippedSyst =args.skipSyst

p=plotter(outDir=OUTPUT, inDir = INPUT)
p.plotEtaPtInt(skipSyst=skippedSyst)

