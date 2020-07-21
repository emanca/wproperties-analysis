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
    
    def __init__(self, outDir, inDir = '', norm = 1):

        self.indir = inDir # indir containig the various outputs
        self.outdir = outDir
        self.norm = norm
                
        self.sampleDict = {
            "WToMu"      :  ['WToMu_plots.root',       'prefit_Signal',         ROOT.kRed+2,       "W^{+}#rightarrow #mu^{+}#nu_{#mu}"],         
            "DYJets"      : ['DYJets_plots.root',      'prefit_Signal',         ROOT.kAzure+2,     "DYJets"],              
            "WtoTau"      : ['WToTau_plots.root',      'prefit_Signal',         ROOT.kSpring+9,    "W^{#pm}#rightarrow #tau^{#pm}#nu_{#tau}"],   
            "TW"          : ['TW_plots.root',          'prefit_Signal',         ROOT.kGreen+3,     "Top"],    
            "TTbar"       : ['TTJets_plots.root',      'prefit_Signal',         ROOT.kGreen+3,     "Top"],    
            "ST"          : ['SingleTop_plots.root',   'prefit_Signal',         ROOT.kGreen+3,     "Top"],    
            "DiBoson"     : ['Diboson_plots.root',     'prefit_Signal',         ROOT.kViolet+2,    "di-boson"],     
            "SIGNAL_Fake" : ['Data_plots.root',        'prefit_fakes',          ROOT.kGray,        "QCD"],     
            "Data"        : ['Data_plots.root',        'prefit_Signal',         1,                 "Data"]
            # "WToMu"      :  ['WJets_plots.root',       'prefit_SignalWToMu',    ROOT.kRed+2,       "W^{+}#rightarrow #mu^{+}#nu_{#mu}"],         
            # "WtoTau"      : ['WJets_plots.root',       'prefit_SignalWToTau',   ROOT.kSpring+9,    "W^{#pm}#rightarrow #tau^{#pm}#nu_{#tau}"],         
            }   
        
        self.variableDict = {
            "Mu1_pt_plus"   :  ["Mu1_pt",   "p_{T} (#mu^{+})",    "Events /binWidth"],
            "Mu1_pt_minus"  :  ["Mu1_pt",   "p_{T} (#mu^{-})",    "Events /binWidth"],
            "Mu1_eta_plus"  :  ["Mu1_eta",  "#eta (#mu^{+})",     "Events /binWidth"],
            "Mu1_eta_minus" :  ["Mu1_eta",  "#eta (#mu^{-})",     "Events /binWidth"],
            "MT_plus"       :  ["MT",       "M_{T} (#mu^{+})",    "Events /binWidth"],      
            "MT_minus"      :  ["MT",       "M_{T} (#mu^{-})",    "Events /binWidth"],      
        }
        
        self.signDict = {
            'minus' : [1, ROOT.kRed+2, "W^{-}#rightarrow #mu^{-}#bar{#nu}_{#mu}"],
            'plus' :  [2, ROOT.kRed+1, "W^{+}#rightarrow #mu^{+}#nu_{#mu}"]
        }
        
        self.sampleOrder = ['DiBoson','WtoTau','ST','TTbar','TW','DYJets','SIGNAL_Fake','WToMu']
        self.histoDict ={} 
        
        self.extSyst = copy.deepcopy(bkg_utils.bkg_systematics)
        self.extSyst['Nominal'] = ['']
        
        self.LHEdict = {
            'Down' : ["LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0", "LHEScaleWeight_muR1p0_muF0p5"],
            'Up' : ["LHEScaleWeight_muR2p0_muF2p0", "LHEScaleWeight_muR2p0_muF1p0","LHEScaleWeight_muR1p0_muF2p0"]   
        }


        if not os.path.exists(self.outdir):
            os.system("mkdir -p " + self.outdir)

    def varBinWidth_modifier(self,histo):
        variable_width = False
        bin_width = histo.GetXaxis().GetBinWidth(1)
        for i in range(2, histo.GetNbinsX()+1):
            if histo.GetXaxis().GetBinWidth(i)!=bin_width:
                variable_width = True
                break
        if variable_width :
            for i in range(1, histo.GetNbinsX()+1):
                old = histo.GetBinContent(i)
                bin_width = histo.GetXaxis().GetBinWidth(i)
                histo.SetBinContent(i, old/bin_width)
        
    def CloneEmpty(self, histo,name) :
        hout = histo.Clone(name)
        for i in range(0, histo.GetNbinsX()+2) :
            hout.SetBinContent(i,0)
            hout.SetBinError(i,0)
        return hout 
            
    def getHistos(self):
        for f,fileInfo in self.sampleDict.iteritems() :
            # inFile = ROOT.TFile.Open(self.indir+'/'+fileInfo[0])
            inFile = ROOT.TFile.Open(self.outdir+'/hadded/'+fileInfo[0])
            for sKind, sList in self.extSyst.iteritems():
                for sName in sList :
                    for var, varInfo in self.variableDict.iteritems() :
                        inFile.cd()
                        if ROOT.gDirectory.Get(fileInfo[1]+'/'+sKind+'/'+varInfo[0]+'_'+sName)==None : #this syst is missing --> take the nominal
                            if sName!='' or f!='Data': print "missing syst:", sName, " for file", f
                            h2 = inFile.Get(fileInfo[1]+'/Nominal/'+varInfo[0])
                        else : 
                            h2 = inFile.Get(fileInfo[1]+'/'+sKind+'/'+varInfo[0]+'_'+sName)
                        for s,sInfo in self.signDict.iteritems() :
                            self.histoDict[f+s+var+sName] = h2.ProjectionX(h2.GetName() + s, sInfo[0],sInfo[0])
                            self.varBinWidth_modifier(self.histoDict[f+s+var+sName])
   
                            
    def plotStack(self,skipSyst=[]):

        self.getHistos()
        fname = "{dir}/stackplots.root".format(dir=self.outdir)
        outFile =  ROOT.TFile(fname, "RECREATE")
        
        for var,varInfo in self.variableDict.iteritems() :
            
            if var.endswith('minus') : s='minus'  
            else : s='plus'  
            
            legend = ROOT.TLegend(0.5, 0.6, 0.90, 0.85)            
            hStack = ROOT.THStack('stack_'+var,varInfo[1])
            hSum = self.CloneEmpty(self.histoDict[self.sampleOrder[0]+s+var],'hsum_'+var) #sum of MC for ratio evaluation
            hData = self.histoDict['Data'+s+var].Clone('hData_'+var)

            #build stacked plot
            for sample in self.sampleOrder :
                self.histoDict[sample+s+var].SetFillStyle(1001)
                self.histoDict[sample+s+var].SetLineWidth(1)
                self.histoDict[sample+s+var].SetFillColor(self.sampleDict[sample][2])
                hStack.Add(self.histoDict[sample+s+var])
                hSum.Add(self.histoDict[sample+s+var])
                if sample == 'TTbar' or sample == 'TW' : continue #all top together, only one legend
                if sample =='WToMu' : #signal has different color and legend entry
                    self.histoDict[sample+s+var].SetFillColor(self.signDict[s][1])
                    legTag = self.signDict[s][2]
                else : legTag = self.sampleDict[sample][3]
                legend.AddEntry(self.histoDict[sample+s+var],legTag,  "f")                    
            legend.AddEntry(hData, 'Data', "PE1")
            
            #build ratio plot Data/pred
            hRatio = hData.Clone('hRatio_'+var)
            hRatio.Divide(hSum) #nominal ratio
            hRatioDict = {} #syst ratio
            for sKind, sList in bkg_utils.bkg_systematics.iteritems():
                if sKind in skipSyst : continue #skipped systs
                for sName in sList :
                    hSumSyst= self.CloneEmpty(self.histoDict[self.sampleOrder[0]+s+var],'hsum_'+var+'_'+sName)
                    for sample in self.sampleOrder :
                        hSumSyst.Add(self.histoDict[sample+s+var+sName])
                    hRatioDict[sName] = hData.Clone('hRatio_'+var+'_'+sName)
                    hRatioDict[sName].Divide(hSumSyst)
                    
            hRatioBand = hRatio.Clone('hRatioBand') #systband
            for i in range(1,hRatioBand.GetNbinsX()+1) :
                delta = 0
                for syst, hsyst in hRatioDict.iteritems() :
                    if 'Down' in syst : continue
                    if syst in self.LHEdict['Down']: continue
                    if 'Up' in syst :
                        systDown =  syst.replace("Up","Down")
                    else :
                        for jj in range(len(self.LHEdict['Up'])) :
                            if syst == self.LHEdict['Up'][jj] :
                                systDown = self.LHEdict['Down'][jj] 
                    
                    delta += (hsyst.GetBinContent(i)-hRatioDict[systDown].GetBinContent(i))**2
                    # delta + = (hsyst.GetBinContent(i)-hRatio.GetBinContent(i))**2
                    if (hRatioDict[systDown].GetBinContent(i)<hRatio.GetBinContent(i) and hRatioDict[syst].GetBinContent(i)<hRatio.GetBinContent(i)) or (hRatioDict[systDown].GetBinContent(i)>hRatio.GetBinContent(i) and hRatioDict[syst].GetBinContent(i)>hRatio.GetBinContent(i)) : #nominal not in between systs
                        print var,"WARNING: systematic", syst," up/down not around nominal in bin", i, hRatioDict[systDown].GetBinContent(i), hRatio.GetBinContent(i), hRatioDict[syst].GetBinContent(i)
                    
                delta = 0.5*math.sqrt(delta)
                hRatioBand.SetBinError(i, delta)
    
            #build the canvas
            can = ROOT.TCanvas(var,var,800,700)
            can.cd()
            pad_histo = ROOT.TPad("pad_histo_"+var, "pad_histo_"+var,0,0.2,1,1)
            pad_ratio = ROOT.TPad("pad_ratio_"+var, "pad_ratio_"+var,0,0,1,0.2)
            pad_histo.SetBottomMargin(0.02)
            pad_histo.Draw()
            pad_ratio.SetTopMargin(0)
            pad_ratio.SetBottomMargin(0.25)
            pad_ratio.Draw()
            
            pad_histo.cd()
            pad_histo.SetGridx()
            pad_histo.SetGridy()
            hStack.Draw("HIST")
            hData.Draw("PE1SAME")
            legend.Draw("SAME")
            
            pad_ratio.cd()
            pad_ratio.SetGridx()
            pad_ratio.SetGridy()
            hRatioBand.Draw('E3')
            hRatio.Draw("PE1SAME")
            
            #aesthetic features
            legend.AddEntry(hRatioBand, "Tot. Syst. Unc.")
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            legend.SetNColumns(2)
            
            hStack.SetMaximum(1.5*max(hData.GetMaximum(),hStack.GetMaximum())) 
            hStack.GetYaxis().SetTitle(varInfo[2])
            hStack.GetYaxis().SetTitleOffset(1.3)
            hStack.GetXaxis().SetTitle(varInfo[1])
            hStack.GetXaxis().SetTitleOffset(3)
            hStack.GetXaxis().SetLabelOffset(3)
            
            hData.SetMarkerStyle(20)
            hData.SetMarkerColor(self.sampleDict['Data'][2])
           
            hRatio.SetMarkerStyle(1)
            hRatio.SetLineColor(1)
            hRatio.SetLineWidth(2)
            
            hRatioBand.SetLineWidth(0)
            hRatioBand.SetFillColor(ROOT.kCyan-4)
            hRatioBand.SetFillStyle(3001)
            hRatioBand.SetMarkerStyle(1)
            hRatioBand.SetTitle("")
            hRatioBand.GetYaxis().SetTitle('Data/Pred')
            hRatioBand.GetYaxis().SetTitleOffset(0.25)
            hRatioBand.GetYaxis().SetNdivisions(506)
            hRatioBand.SetTitleSize(0.15,'y')
            hRatioBand.SetLabelSize(0.12,'y')
            hRatioBand.GetXaxis().SetTitle(varInfo[1])
            hRatioBand.GetXaxis().SetTitleOffset(0.6)
            hRatioBand.SetTitleSize(0.18,'x')
            hRatioBand.SetLabelSize(0.15,'x')
            if var.startswith('MT') :
                hRatioBand.GetYaxis().SetRangeUser(0.8,1.2)
            
            # can.Update()
            can.SaveAs("{dir}/{c}.pdf".format(dir=self.outdir,c=can.GetName()))
            can.SaveAs("{dir}/{c}.png".format(dir=self.outdir,c=can.GetName()))

            outFile.cd()
            can.Write()
            
            #debug:
            hRatio.Write()
            for syst, hsyst in hRatioDict.iteritems() :
                hsyst.Write()
            
        outFile.Close()


def prepareHistos(inDir,outDir) :
    cmdList = []
    if not os.path.exists(outDir): os.system("mkdir -p " + outDir)
    if not os.path.isdir(outDir+'/hadded'): os.system('mkdir '+outDir+'/hadded')
        
    # cmdList.append('cp  '+inDir+'SingleMuonData_plots.root '+outDir+'/hadded/Data_plots.root')
    # cmdList.append('cp  '+inDir+'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_plots.root '+outDir+'/hadded/WJets_plots.root')
    cmdList.append('cp  '+inDir+'SingleMuonData_plots.root '+outDir+'/hadded/Data_plots.root')
    cmdList.append('cp  '+inDir+'WToMu_plots.root '+outDir+'/hadded/WToMu_plots.root')
    cmdList.append('cp  '+inDir+'WToTau_plots.root '+outDir+'/hadded/WToTau_plots.root')
    cmdList.append('hadd -f '+outDir+'/hadded/DYJets_plots.root '+inDir+'DYJetsToLL_M-*')
    cmdList.append('hadd -f '+outDir+'/hadded/TTJets_plots.root '+inDir+'TTJets*')
    cmdList.append('hadd -f '+outDir+'/hadded/SingleTop_plots.root  '+inDir+'ST_t-channel_* '+inDir+'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_plots.root ')
    cmdList.append('hadd -f '+outDir+'/hadded/TW_plots.root  '+inDir+'ST_tW_*')
    cmdList.append('hadd -f '+outDir+'/hadded/Diboson_plots.root '+inDir+'WW_TuneCUETP8M1_13TeV-pythia8_plots.root '+inDir+'WZ_TuneCUETP8M1_13TeV-pythia8_plots.root '+inDir+'ZZ_TuneCUETP8M1_13TeV-pythia8_plots.root')
            
    for i in cmdList :
        os.system(i)

    
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

if HADD :
    prepareHistos(inDir=INPUT,outDir=OUTPUT)
p=plotter(outDir=OUTPUT, inDir = INPUT)
p.plotStack(skipSyst=skippedSyst)