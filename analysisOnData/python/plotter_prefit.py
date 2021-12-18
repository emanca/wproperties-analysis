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
    
    def __init__(self, outDir, inDir = '',SBana=False,bkgOnlySyst=False):

        self.indir = inDir # indir containig the various outputs
        self.outdir = outDir
        self.PDFvar = 'LHEPdfWeight'
        self.SBana = SBana
        self.bkgOnlySyst = bkgOnlySyst
        
        if not self.SBana :
            dirMC = 'prefit_Signal'
            dirFakes = 'prefit_fakes'
            self.SBsuff = ''
        else :
            dirMC = 'prefit_Sideband'
            dirFakes = 'prefit_fakes_SideBand'
            self.SBsuff = '_SideBand'
            
                
        self.sampleDict = {
            "WToMu"      :  ['WToMu_plots.root',        dirMC,         ROOT.kRed+2,       "W^{+}#rightarrow #mu^{+}#nu_{#mu}"],         
            "DYJets"      : ['DYJets_plots.root',       dirMC,         ROOT.kAzure+2,     "Drell-Yan"],              
            "WtoTau"      : ['WToTau_plots.root',       dirMC,         ROOT.kSpring+9,    "W^{#pm}#rightarrow #tau^{#pm}#nu_{#tau}"],   
            "Top"       : ['Top_plots.root',            dirMC,         ROOT.kGreen+3,     "Top"],    
            "DiBoson"     : ['Diboson_plots.root',      dirMC,         ROOT.kViolet+2,    "Diboson"],     
            "SIGNAL_Fake" : ['FakeFromData_plots.root', dirFakes,      ROOT.kGray,        "QCD"],     
            "Data"        : ['Data_plots.root',         dirMC,         1,                 "Data"]
            }    
        
        # self.variableDict = {
        #     "Mu1_pt_plus"   :  ["Mu1_pt",   "p_{T}^{#mu} distribution, W^{+}",    "dN/dp_{T}^{#mu} [GeV^{-1}]", "p_{T}^{#mu} [GeV]"],
        #     "Mu1_pt_minus"  :  ["Mu1_pt",   "p_{T}^{#mu} distribution, W^{-}",    "dN/dp_{T}^{#mu} [GeV^{-1}]", "p_{T}^{#mu} [GeV]"],
        #     "Mu1_eta_plus"  :  ["Mu1_eta",  "#eta distribution, W^{+}",            "dN/d#eta", "#eta^{#mu}"],
        #     "Mu1_eta_minus" :  ["Mu1_eta",  "#eta distribution, W^{-}",            "dN/d#eta", "#eta^{#mu}"],
        #     "MT_plus"       :  ["MT",       "m_{T} distribution, W^{+}",           "dN/dm_{T} [GeV^{-1}]", "m_{T} [GeV]"],      
        #     "MT_minus"      :  ["MT",       "m_{T} distribution, W^{-}",           "dN/dm_{T} [GeV^{-1}]", "m_{T} [GeV]"],      
        # }
        
        self.variableDict = {
            "Mu1_pt_plus"   :  ["Mu1_pt",   "p_{T}^{#mu} distribution, W^{+}",    "Events / ", "p_{T}^{#mu} [GeV]", " GeV"],
            "Mu1_pt_minus"  :  ["Mu1_pt",   "p_{T}^{#mu} distribution, W^{-}",    "Events / ", "p_{T}^{#mu} [GeV]", " GeV"],
            "Mu1_eta_plus"  :  ["Mu1_eta",  "#eta distribution, W^{+}",            "Events / ", "#eta^{#mu}", ""],
            "Mu1_eta_minus" :  ["Mu1_eta",  "#eta distribution, W^{-}",            "Events / ", "#eta^{#mu}", ""],
            "MT_plus"       :  ["MT",       "m_{T} distribution, W^{+}",           "Events / ", "m_{T} [GeV]", " GeV"],      
            "MT_minus"      :  ["MT",       "m_{T} distribution, W^{-}",           "Events / ", "m_{T} [GeV]", " GeV"],      
        }
        
        self.signDict = {
            'minus' : [1, ROOT.kRed+1, "W^{-}#rightarrow #mu^{-}#bar{#nu}_{#mu}"],
            'plus' :  [2, ROOT.kRed+2, "W^{+}#rightarrow #mu^{+}#nu_{#mu}"]
        }
        
        self.sampleOrder = ['DiBoson','WtoTau','Top','DYJets','SIGNAL_Fake','WToMu']
        self.histoDict ={} 
        
        self.extSyst = copy.deepcopy(bkg_utils.bkg_systematics)
        self.extSyst['Nominal'] = ['']
        
        # self.LHEdict = {
        #     'Down' : ["LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0", "LHEScaleWeight_muR1p0_muF0p5"],
        #     'Up' : ["LHEScaleWeight_muR2p0_muF2p0", "LHEScaleWeight_muR2p0_muF1p0","LHEScaleWeight_muR1p0_muF2p0"]   
        # }
        
        self.groupedSystColors = {
            "WHSF"  : [ROOT.kGreen+1, 'Scale Factors'],
            "LHEScaleWeight" : [ROOT.kViolet-2, 'MC Scale'],
            "ptScale" : [ROOT.kYellow+2, 'p_{T} Scale'],
            "jme" : [ROOT.kAzure+10, 'MET'],
            "LHEPdfWeight" : [ROOT.kRed+1, 'PDF+#alpha_{s}'],
            "Nominal" : [1, 'Stat. Unc.'],
            "PrefireWeight" : [ROOT.kSpring+10, 'Prefire'],
            "alphaS" : [ROOT.kOrange-3, '#alpha_{s}'],
            "LHEScaleWeight_WQTlow" : [ROOT.kViolet+7, "MC Scale Wqt<5"],
            "LHEScaleWeight_WQTmid" : [ROOT.kViolet+7, "MC Scale 5<Wqt<15"],
            "LHEScaleWeight_WQThigh" : [ROOT.kViolet+7, "q_{T}^{V}"], #it contains also two previous lines (not plotted)
            "mass" : [ROOT.kBlue-4, 'm_{W}', 35],
            "lumi" : [ROOT.kOrange-7,"Lumi"],
            "topXSec" : [ROOT.kCyan-6,"#sigma_{t}"],
            "dibosonXSec" : [ROOT.kCyan-1,"#sigma_{diboson}"],
            "tauXSec" : [ROOT.kTeal,"#sigma_{W#rightarrow#tau#nu}+#sigma_{t}+#sigma_{diboson}"],
            "lepVeto" : [ROOT.kMagenta-7,"Lepton veto"],
        }
        
        self.flatDict = {
            "lumi": {
                "procs": ["WToMu", 'Top','DiBoson',"WtoTau","DYJets","SIGNAL_Fake"],
                "weight" : 1.025 
            },
            "topXSec":{
                "procs": ['Top'],
                "weight" : 1.060 
            },
            "dibosonXSec":{
                "procs": ['DiBoson'],
                "weight" : 1.160 
            },
            "tauXSec":{
                "procs": ["WtoTau"],
                "weight" : 1.040 
            },
            "lepVeto":{
                "procs": ["DYJets"],
                "weight" : 1.020 
            }
        }
        
        self.doNotPlotGroup = ["LHEScaleWeight_WQTlow", "LHEScaleWeight_WQTmid","LHEScaleWeight", "topXSec", "dibosonXSec", "alphaS" ] #already included in other groups (summed)



        if not os.path.exists(self.outdir):
            os.system("mkdir -p " + self.outdir)

    def varBinWidth_modifier(self,histo):
        variable_width = False
        bin_width = histo.GetXaxis().GetBinWidth(1)
        for i in range(2, histo.GetNbinsX()+1):
            if histo.GetXaxis().GetBinWidth(i)!=bin_width:
                variable_width = True
                break
        if variable_width or 1 : #all the histograms expressed as a d N/dX
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
        for f,fileInfo in self.sampleDict.items() :
            inFile = ROOT.TFile.Open(self.indir+'/hadded/'+fileInfo[0])
            for sKind, sList in self.extSyst.items():
                if sKind in self.flatDict and f!='SIGNAL_Fake' and not (self.bkgOnlySyst and f=='WToMu'):
                        continue
                for sName in sList :
                    for var, varInfo in self.variableDict.items() :
                        inFile.cd()
                        if '_WQT' in sKind and f!='SIGNAL_Fake' :
                            sName = sName.replace(sKind.replace('LHEScaleWeight',''),'') #replace for instance: LHEScaleWeight_muR0p5_muF0p5_WQTlow-> LHEScaleWeight_muR0p5_muF0p5
                        if ROOT.gDirectory.Get(fileInfo[1]+'/'+sKind+'/'+varInfo[0]+'_'+sName)==None or (self.bkgOnlySyst and f=='WToMu') : #this syst is missing --> take the nominal or (skip syst because bkgOnlySyst)
                            if sName!='' or f!='Data': print("missing syst:", sName, " for file", f)
                            h2 = inFile.Get(fileInfo[1]+'/Nominal/'+varInfo[0])
                        else : 
                            h2 = inFile.Get(fileInfo[1]+'/'+sKind+'/'+varInfo[0]+'_'+sName)
                            
                            if f!='SIGNAL_Fake' : #WtoTau should be the only with WQT syst, except fake (WJets skipped, not needed)
                                if sKind=='LHEScaleWeight_WQTlow' :
                                    h2.Add(inFile.Get(fileInfo[1]+'/Nom_WQTmid/'+varInfo[0]))
                                    h2.Add(inFile.Get(fileInfo[1]+'/Nom_WQThigh/'+varInfo[0]))
                                if sKind=='LHEScaleWeight_WQTmid' :
                                    h2.Add(inFile.Get(fileInfo[1]+'/Nom_WQTlow/'+varInfo[0]))
                                    h2.Add(inFile.Get(fileInfo[1]+'/Nom_WQThigh/'+varInfo[0]))
                                if sKind=='LHEScaleWeight_WQThigh' :
                                    h2.Add(inFile.Get(fileInfo[1]+'/Nom_WQTlow/'+varInfo[0]))
                                    h2.Add(inFile.Get(fileInfo[1]+'/Nom_WQTmid/'+varInfo[0]))
                            
                        if '_WQT' in sKind and f!='SIGNAL_Fake' :
                            if f=='WToMu' : 
                                h2 = inFile.Get(fileInfo[1]+'/Nominal/'+varInfo[0]) #get the nominal for the signal, wqt should not be used for it
                            sName=sName+sKind.replace('LHEScaleWeight','')
                        
                        if sKind=='LHEScaleWeight' and f!='DYJets' :
                            h2 = inFile.Get(fileInfo[1]+'/Nominal/'+varInfo[0])#use the LHEscale only for the Z (this should be automatically included in analysisOnMC/WJets step) 
                        
                        for s,sInfo in self.signDict.items() :
                            self.histoDict[f+s+var+sName] = h2.ProjectionX(h2.GetName() + s, sInfo[0],sInfo[0])
                            # self.varBinWidth_modifier(self.histoDict[f+s+var+sName])
            
            for sKind, sList in self.extSyst.items():  #flat syst building
                if sKind not in self.flatDict : continue 
                if f=='SIGNAL_Fake' : continue
                if self.bkgOnlySyst and f=='WToMu' : continue 
                for sName in sList : 
                    for var, varInfo in self.variableDict.items() :
                        h2 = inFile.Get(fileInfo[1]+'/Nominal/'+varInfo[0])    
                        procFlag = False 
                        for proc in self.flatDict[sKind]['procs'] :
                            if proc in f : #proc are part of the full filename of f
                                procFlag = True 
                        if procFlag :
                            if 'Up' in sName :
                                h2.Scale(self.flatDict[sKind]['weight'])
                            if 'Down' in sName :
                                h2.Scale(1/self.flatDict[sKind]['weight'])
                        for s,sInfo in self.signDict.items() :
                            self.histoDict[f+s+var+sName] = h2.ProjectionX(h2.GetName() + s, sInfo[0],sInfo[0])
                            # self.varBinWidth_modifier(self.histoDict[f+s+var+sName])         
                                            
    def plotStack(self,skipSyst=[]):

        self.getHistos()
        fname = "{dir}/stackPlots{suff}.root".format(dir=self.outdir,suff=self.SBsuff)
        outFile =  ROOT.TFile(fname, "RECREATE")
        
        cmslab = "#bf{CMS} #scale[0.7]{#it{Work in progress}}"
        lumilab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
        cmsLatex = ROOT.TLatex()
        
        for var,varInfo in self.variableDict.items() :
            
            if var.endswith('minus') : s='minus'  
            else : s='plus'  
            
            legend = ROOT.TLegend(0.5, 0.55, 0.90, 0.85)            
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
            
            legend.AddEntry(hData, 'Data', "PE1")
            for sample in self.sampleOrder[::-1]: #reverse order
                if sample =='WToMu' : #signal has different color and legend entry
                    self.histoDict[sample+s+var].SetFillColor(self.signDict[s][1])
                    legTag = self.signDict[s][2]
                else : legTag = self.sampleDict[sample][3]
                legend.AddEntry(self.histoDict[sample+s+var],legTag,  "f")                    
            
            
            #build ratio plot Data/pred
            hRatio = hData.Clone('hRatio_'+var)
            hRatio.Divide(hSum) #nominal ratio
            hRatioDict = {} #syst ratio
            for sKind, sList in bkg_utils.bkg_systematics.items():
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
                for syst, hsyst in hRatioDict.items() :
                    if 'Down' in syst : continue
                    # if syst in self.LHEdict['Down']: continue
                    if 'LHE' in syst : continue 
                    if 'Up' in syst :
                        systDown =  syst.replace("Up","Down")
                    # else :
                    #     for jj in range(len(self.LHEdict['Up'])) :
                    #         if syst == self.LHEdict['Up'][jj] :
                    #             systDown = self.LHEdict['Down'][jj] 
                    
                        delta += (hsyst.GetBinContent(i)-hRatioDict[systDown].GetBinContent(i))**2
                        # delta + = (hsyst.GetBinContent(i)-hRatio.GetBinContent(i))**2
                        if (hRatioDict[systDown].GetBinContent(i)<hRatio.GetBinContent(i) and hRatioDict[syst].GetBinContent(i)<hRatio.GetBinContent(i)) or (hRatioDict[systDown].GetBinContent(i)>hRatio.GetBinContent(i) and hRatioDict[syst].GetBinContent(i)>hRatio.GetBinContent(i)) : #nominal not in between systs
                            print(var,"WARNING: systematic", syst," up/down not around nominal in bin", i, hRatioDict[systDown].GetBinContent(i), hRatio.GetBinContent(i), hRatioDict[syst].GetBinContent(i))
                            
                delta = 0.25*delta #1/2 for up/down, 1/2 for sqrt
                
                deltaPDF=0 #LHE PDF variations (wrt nominal)
                for syst, hsyst in hRatioDict.items() : 
                    if not 'LHEPdf' in syst: continue 
                    Nrepl = 1.
                    # if 'LHEPdfWeight' in syst :
                    #     Nrepl = float(len(bkg_utils.bkg_systematics[self.PDFvar]))
                    deltaPDF+= (1/Nrepl)*(hsyst.GetBinContent(i)-hRatio.GetBinContent(i))**2
                
                deltaScale=0 #LHE Scale variations (envelope)
                for syst, hsyst in hRatioDict.items() : 
                    if not 'LHEScale' in syst: continue 
                    # if 'WQT' in syst: continue
                    deltaScale_temp= (hsyst.GetBinContent(i)-hRatio.GetBinContent(i))**2
                    if deltaScale_temp>deltaScale : 
                        deltaScale = deltaScale_temp
                
                # deltaWQT = 0 #if i want to sum wqt linearly and not **2
                # for Bin in ['low', 'mid', 'high'] :
                #     deltaWQT_Bin = 0
                #     for syst, hsyst in hRatioDict.items() : 
                #         if 'WQT' in syst: continue
                #         if Bin not in syst : continue 
                #         deltaWQT_temp= abs(hsyst.GetBinContent(i)-hRatio.GetBinContent(i))
                #         if deltaWQT_temp>deltaWQT : 
                #             deltaWQT_Bin = deltaWQT_temp
                #     deltaWQT+=deltaWQT_bin
                        
                deltaSum = math.sqrt(delta+deltaPDF+deltaScale)
                hRatioBand.SetBinError(i, deltaSum)
    
            #build the canvas
            can = ROOT.TCanvas('prefit_'+var,'prefit_'+var,800,700)
            can.cd()
            pad_histo = ROOT.TPad("pad_histo_"+var, "pad_histo_"+var,0,0.2,1,1)
            pad_ratio = ROOT.TPad("pad_ratio_"+var, "pad_ratio_"+var,0,0,1,0.2)
            pad_histo.SetBottomMargin(0.02)
            pad_histo.Draw()
            pad_ratio.SetTopMargin(0)
            pad_ratio.SetBottomMargin(0.32)
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
            if 'MT' in var :
                hRatioBand.Draw('E2')
            else :
                hRatioBand.Draw('E3')
            hRatio.Draw("PE1SAME")
            
            #aesthetic features
            legend.AddEntry(hRatioBand, "Tot. Syst. Unc.")
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            legend.SetNColumns(2)
            
            hStack.SetMaximum(1.7*max(hData.GetMaximum(),hStack.GetMaximum())) 
            # hStack.GetYaxis().SetTitle(varInfo[2])
            if 'MT' in var :
                hStack.GetYaxis().SetTitle(varInfo[2]+str(int(hStack.GetXaxis().GetBinWidth(1)))+varInfo[4])
            else :
                hStack.GetYaxis().SetTitle(varInfo[2]+str(round(hStack.GetXaxis().GetBinWidth(1),1))+varInfo[4])
            hStack.GetYaxis().SetTitleOffset(1.3)
            hStack.GetXaxis().SetTitle(varInfo[3])
            hStack.GetXaxis().SetTitleOffset(3)
            hStack.GetXaxis().SetLabelOffset(3)
            hStack.SetTitle("")
            
            hData.SetMarkerStyle(20)
            hData.SetMarkerColor(self.sampleDict['Data'][2])
           
            hRatio.SetMarkerStyle(1)#20
            hRatio.SetLineColor(1)
            hRatio.SetLineWidth(2)
            
            hRatioBand.SetLineWidth(0)
            hRatioBand.SetFillColor(ROOT.kOrange)#(ROOT.kCyan-4)
            hRatioBand.SetFillStyle(3001)
            hRatioBand.SetMarkerStyle(1)
            hRatioBand.SetTitle("")
            hRatioBand.GetYaxis().SetTitle('Data/Pred')
            hRatioBand.GetYaxis().SetTitleOffset(0.25)
            hRatioBand.GetYaxis().SetNdivisions(506)
            hRatioBand.SetTitleSize(0.15,'y')
            hRatioBand.SetLabelSize(0.12,'y')
            hRatioBand.GetXaxis().SetTitle(varInfo[3])
            hRatioBand.GetXaxis().SetTitleOffset(0.8)
            hRatioBand.SetTitleSize(0.18,'x')
            hRatioBand.SetLabelSize(0.15,'x')
            if var.startswith('MT') :
                hRatioBand.GetYaxis().SetRangeUser(0.8,1.23)
            elif var.startswith('Mu1_pt'):
                hRatioBand.GetYaxis().SetRangeUser(0.88,1.18)
            else :
                hRatioBand.GetYaxis().SetRangeUser(0.9,1.13)
            
            pad_histo.cd()
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.SetTextSize(cmsLatex.GetTextSize())
            cmsLatex.DrawLatex(1-pad_histo.GetRightMargin(),1-0.8*pad_histo.GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11)
            cmsLatex.DrawLatex(0.07+pad_histo.GetLeftMargin(),1-0.8*pad_histo.GetTopMargin(),cmslab)
            
            # can.Update()
            can.SaveAs("{dir}/{c}{suff}.pdf".format(dir=self.outdir,c=can.GetName(),suff=self.SBsuff))
            can.SaveAs("{dir}/{c}{suff}.png".format(dir=self.outdir,c=can.GetName(),suff=self.SBsuff))

            outFile.cd()
            can.Write()
            
            #debug:
            # hRatio.Write()
            # for syst, hsyst in hRatioDict.iteritems() :
            #     hsyst.Write()
            
        outFile.Close()

    def plotSyst(self,skipSyst=[]) :
        fname = "{dir}/systPlots{suff}.root".format(dir=self.outdir,suff=self.SBsuff)
        outFile =  ROOT.TFile(fname, "RECREATE")
        
        cmslab = "#bf{CMS} #scale[0.7]{#it{Work in progress}}"
        lumilab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
        cmsLatex = ROOT.TLatex()
        
        for var,varInfo in self.variableDict.items() :
            if var.endswith('minus') : s='minus'  
            else : s='plus' 
            
            hdict = {}
            
            #completely broken syst 
            legend = ROOT.TLegend(0.15, 0.6, 0.85, 0.85) 
            hData = self.histoDict['Data'+s+var].Clone('hData_'+var) 
            for sKind, sList in self.extSyst.items(): #summing MCs
                for sName in sList :    
                    hdict[sName] = self.CloneEmpty(self.histoDict[self.sampleOrder[0]+s+var+sName],'hsum_'+var+'_'+sName)
                    for sample in self.sampleOrder :
                        hdict[sName].Add(self.histoDict[sample+s+var+sName])
            
            hSumNOM =hdict[''].Clone("sumNom"+var)          
            for sKind, sList in self.extSyst.items(): #ratios
                for sName in sList :    
                    hdict[sName].Divide(hSumNOM)
                    
            for x in range(1,hdict[''].GetNbinsX()+1) : #stat error using data
                if hdict[''].GetBinContent(x)==0 :
                    hdict[''].SetBinError(x,0)
                else :
                    hdict[''].SetBinError(x,hData.GetBinError(x)/ hSumNOM.GetBinContent(x))
            
            #LHEPDF variation band
            hdict['LHEPdfUp'] = hdict[''].Clone('hsum_'+var+'LHEPdfUp')
            hdict['LHEPdfDown'] = hdict[''].Clone('hsum_'+var+'LHEPdfDown')
            for x in range(1,hdict['LHEPdfUp'].GetNbinsX()+1) :
                stdVar = 0
                for sName in bkg_utils.bkg_systematics[self.PDFvar] :
                    stdVar+= (hdict[sName].GetBinContent(x)-hdict[''].GetBinContent(x))**2
                # Nrepl = float(len(bkg_utils.bkg_systematics[self.PDFvar])) 
                Nrepl=1.
                stdVar = math.sqrt(stdVar/Nrepl)
                hdict['LHEPdfUp'].SetBinContent(x, hdict['LHEPdfUp'].GetBinContent(x)+stdVar)
                hdict['LHEPdfDown'].SetBinContent(x, hdict['LHEPdfDown'].GetBinContent(x)-stdVar)
            
            #build the canvas
            can = ROOT.TCanvas('prefit_'+var+'_systBreakdown','prefit_'+var+'_systBreakdown',800,600)
            can.cd()
            can.SetGridx()
            can.SetGridy()
            
            hdict[''].SetLineWidth(0)
            hdict[''].SetLineColor(1)
            hdict[''].SetFillColor(1)
            hdict[''].SetFillStyle(3001)
            hdict[''].Draw()
            hdict[''].Draw("same E2")# same 0P5
            hdict[''].SetTitle(varInfo[1]+', systematic breakdown')
            hdict[''].GetYaxis().SetTitle('Syst/Nom')
            hdict[''].GetXaxis().SetTitle(varInfo[3])
            if 'MT' in var :
                hdict[''].GetYaxis().SetRangeUser(0.8,1.3)
            else :
                 hdict[''].GetYaxis().SetRangeUser(0.95,1.1)
            legend.AddEntry(hdict[''], 'Stat. from Data')
            # legend.SetFillStyle(0)
            # legend.SetBorderSize(0)
            legend.SetNColumns(3)
            
            colorList = [600,616,416,632,432,800,900,880,840,820,920,800,850,950,840,920,640,690]
            colorNumber = 1
            colorCounter = 0
            modSyst = copy.deepcopy(bkg_utils.bkg_systematics)
            modSyst['LHEPdfUpDown'] = ['LHEPdfUp','LHEPdfDown']
            for sKind, sList in modSyst.items():
                if sKind==self.PDFvar : continue #skip 100 lines on the plot!
                colorNumber = colorList[colorCounter]
                colorCounter = colorCounter+1
                for sName in sList : 
                    colorNumber = colorNumber-2
                    if colorNumber < colorList[colorCounter-1]-10 :
                        colorNumber = colorList[colorCounter]+2
                        
                    if sKind in skipSyst : continue #skipped systs
                    if self.PDFvar in skipSyst : 
                        if sKind=='LHEPdfUpDown' : continue 
    
                    hdict[sName].SetLineWidth(3)
                    hdict[sName].SetLineColor(1)  
                    hdict[sName].SetLineColor(colorNumber)
                    hdict[sName].Draw("hist SAME")
                    
                    hdict[sName].SetFillStyle(0)
                    hdict[sName].SetFillColor(0)
    
                    legend.AddEntry(hdict[sName], sName.replace('LHEScaleWeight','Scale'))
            legend.Draw("SAME")
        
            can.SaveAs("{dir}/{c}{suff}.pdf".format(dir=self.outdir,c=can.GetName(),suff=self.SBsuff))
            can.SaveAs("{dir}/{c}{suff}.png".format(dir=self.outdir,c=can.GetName(),suff=self.SBsuff))
            
            outFile.cd()
            can.Write()  
            
            
            
            #grouped syst
            legend2 = ROOT.TLegend(0.15, 0.7, 0.7, 0.9) 
            for sKind, sList in self.extSyst.items(): #ratios
                hdict[sKind] = self.CloneEmpty(self.histoDict[self.sampleOrder[0]+s+var],'hsum_'+var+'_'+sKind)
                for i in range(1,hdict[sKind].GetNbinsX()+1) :
                    delta = 0
                    for sName in sList :
                        if 'Down' in sName : continue
                        # if sName in self.LHEdict['Down']: continue
                        if 'LHE' in sName : continue
                        if sName=='' : continue
                        if 'Up' in sName :
                            systDown =  sName.replace("Up","Down")
                        # else :
                        #     for jj in range(len(self.LHEdict['Up'])) :
                        #         if sName == self.LHEdict['Up'][jj] :
                        #             systDown = self.LHEdict['Down'][jj] 
                        
                            delta += (hdict[sName].GetBinContent(i)-hdict[systDown].GetBinContent(i))**2    
                    delta = 0.5*math.sqrt(delta)
                    if sKind=='Nominal' :
                        delta = hdict[''].GetBinError(i)
                    if 'LHEPdf' in sKind :  
                    # if sKind=='LHEScaleWeightVars' :  
                        delta=0
                        for sName in sList :                            
                            delta+= (hdict[sName].GetBinContent(i)-hdict[''].GetBinContent(i))**2
                        # if sKind==self.PDFvar :
                        #     Nrepl = float(len(sList)) 
                        # else : 
                        Nrepl=1.#hessian approach
                        delta = math.sqrt(delta/Nrepl)
                    if 'LHEScale' in sKind :
                        delta = 0
                        for sName in sList :
                            delta_temp = (hdict[sName].GetBinContent(i)-hdict[''].GetBinContent(i))**2
                            if delta_temp> delta :
                                delta=delta_temp
                        delta = math.sqrt(delta)
                    hdict[sKind].SetBinContent(i, delta) 
                hdict[sKind].SetFillStyle(0)
                hdict[sKind].SetFillColor(0)
            sKind = 'LHEScaleWeight_WQThigh'  
            for ipt in range(1,hdict[sKind].GetNbinsX()+1) :
                deltaSumWQT = hdict[sKind].GetBinContent(ipt)**2
                deltaSumWQT += hdict[sKind.replace('high','mid')].GetBinContent(ipt)**2
                deltaSumWQT += hdict[sKind.replace('high','low')].GetBinContent(ipt)**2
                deltaSumWQT += hdict['LHEScaleWeight'].GetBinContent(ipt)**2
                deltaSumWQT = math.sqrt(deltaSumWQT)
                hdict[sKind].SetBinContent(ipt , deltaSumWQT)
            
            sKind = 'LHEPdfWeight'  
            for ipt in range(1,hdict[sKind].GetNbinsX()+1) :
                deltaSumPDF = hdict[sKind].GetBinContent(ipt)**2
                deltaSumPDF += hdict["alphaS"].GetBinContent(ipt)**2
                deltaSumPDF = math.sqrt(deltaSumPDF)
                hdict[sKind].SetBinContent(ipt , deltaSumPDF)
            
            sKind = 'tauXSec'  
            for ipt in range(1,hdict[sKind].GetNbinsX()+1) :
                deltaSumXsec = hdict[sKind].GetBinContent(ipt)**2
                deltaSumXsec += hdict["dibosonXSec"].GetBinContent(ipt)**2
                deltaSumXsec += hdict["topXSec"].GetBinContent(ipt)**2
                deltaSumXsec = math.sqrt(deltaSumXsec)
                hdict[sKind].SetBinContent(ipt , deltaSumXsec)
  
            #build the canvas
            can2 = ROOT.TCanvas('prefit_'+var+'_systComparison','prefit_'+var+'_systComparison',800,600)
            can2.cd()
            can2.SetGridx()
            can2.SetGridy()
            can2.SetLogy()
            
            hdict['Nominal'].SetLineWidth(3)
            hdict['Nominal'].SetLineColor(1)
            # hdict['Nominal'].Draw()
            # hdict['Nominal'].Draw()# same 0P5
            # if not self.bkgOnlySyst :
            #     hdict['Nominal'].SetTitle(varInfo[1]+', systematics breakdown')
            # else :
            #     hdict['Nominal'].SetTitle(varInfo[1]+', systematics breakdown (QCD syst only)')
            hdict['Nominal'].SetTitle('')
            hdict['Nominal'].GetYaxis().SetTitle('|Var-Nom| / Nom (Event yield)')
            if self.bkgOnlySyst :
                hdict['Nominal'].GetYaxis().SetTitle('|Var-Nom| / Nom (Event yield, QCD var. only)')
            hdict['Nominal'].GetXaxis().SetTitle(varInfo[3])
            # if 'MT' in var :
            hdict['Nominal'].GetYaxis().SetRangeUser(0.0001,1)
            hdict['Nominal'].SetFillColor(1)
            hdict['Nominal'].SetFillStyle(3002)
            # else :
            #      hdict['Nominal'].GetYaxis().SetRangeUser(0.95,1.1)
            
            legend2.AddEntry(hdict['Nominal'], 'Stat. from Data')
            legend2.SetFillStyle(0)
            legend2.SetBorderSize(0)
            legend2.SetNColumns(2)
            
            hdict['sum'] = self.CloneEmpty(hdict['Nominal'],'hsum_'+var+'_'+'sum')                            
            hdict['sum'].SetFillColor(ROOT.kOrange)
            hdict['sum'].SetLineColor(ROOT.kOrange)
            hdict['sum'].SetFillStyle(3003)
            for i in range(1,hdict['sum'].GetNbinsX()+1) :
                sumOfDelta=0
                for sKind, sList in bkg_utils.bkg_systematics.items():
                    if sKind in skipSyst : continue #skipped systs
                    # if 'WQTlow' in sKind or 'WQTmid' in sKind: continue #high is the sum
                    if sKind in self.doNotPlotGroup: continue
                    sumOfDelta+= hdict[sKind].GetBinContent(i)**2
                sumOfDelta = math.sqrt(sumOfDelta)
                hdict['sum'].SetBinContent(i,sumOfDelta) 
            legend2.AddEntry(hdict['sum'], 'Squared Sum of Syst.')
            hdict['sum'].Draw('hist')
            hdict['Nominal'].Draw('hist SAME')
            
            # colorNumber=2
            for sKind, sList in bkg_utils.bkg_systematics.items():
                    # if colorNumber==5 : colorNumber+=1

                    if sKind in skipSyst : continue #skipped systs
                    # if 'WQTlow' in sKind or 'WQTmid' in sKind: continue #high is the sum
                    if sKind in self.doNotPlotGroup: continue
                    
                    hdict[sKind].SetLineWidth(3)
                    hdict[sKind].SetLineColor(self.groupedSystColors[sKind][0])
                    hdict[sKind].Draw("hist SAME")
                    # colorNumber+=1
                    # legend2.AddEntry(hdict[sKind], sKind.replace('Vars',''))
                    legend2.AddEntry(hdict[sKind],self.groupedSystColors[sKind][1])
                    
            legend2.Draw("SAME")
            
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.SetTextSize(cmsLatex.GetTextSize())
            cmsLatex.DrawLatex(1-can2.GetRightMargin(),1-0.8*can2.GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11)
            cmsLatex.DrawLatex(can2.GetLeftMargin(),1-0.8*can2.GetTopMargin(),cmslab)
            if s=='plus' :     
                proclab = "W^{+}#rightarrow #mu^{+}#nu_{#mu}"
            else :
                proclab = "W^{-}#rightarrow #mu^{-}#bar#nu_{#mu}"
            cmsLatex.SetTextAlign(31)
            cmsLatex.DrawLatex(1-can2.GetRightMargin()-0.02,1-can2.GetTopMargin()-0.05,proclab)
            
            can2.SaveAs("{dir}/{c}{suff}.pdf".format(dir=self.outdir,c=can2.GetName(),suff=self.SBsuff))
            can2.SaveAs("{dir}/{c}{suff}.png".format(dir=self.outdir,c=can2.GetName(),suff=self.SBsuff))
            
            outFile.cd()
            can2.Write()  
                    
        outFile.Close()

def prepareHistos(inDir,outDir) :
    cmdList = []
    if not os.path.exists(outDir): os.system("mkdir -p " + outDir)
    if not os.path.isdir(outDir+'/hadded'): os.system('mkdir '+ outDir+'/hadded')
    cmdList.append('cp  '+inDir+'/SingleMuonData_plots.root '+outDir+'/hadded/Data_plots.root')
    cmdList.append('cp  '+inDir+'/FakeFromData_plots.root '+outDir+'/hadded/FakeFromData_plots.root')
    cmdList.append('cp  '+inDir+'/WToMu_plots.root '+outDir+'/hadded/WToMu_plots.root')
    cmdList.append('cp  '+inDir+'/WToTau_plots.root '+outDir+'/hadded/WToTau_plots.root')
    cmdList.append('hadd -f '+outDir+'/hadded/DYJets_plots.root '+inDir+'/DYJetsToLL_M-*')
    cmdList.append('hadd -f '+outDir+'/hadded/Top_plots.root '+inDir+'/TTJets* ' +inDir+'/ST* ')
    cmdList.append('hadd -f '+outDir+'/hadded/Diboson_plots.root '+inDir+'/WW_TuneCUETP8M1_13TeV-pythia8_plots.root '+inDir+'/WZ_TuneCUETP8M1_13TeV-pythia8_plots.root '+inDir+'/ZZ_TuneCUETP8M1_13TeV-pythia8_plots.root')

    for i in cmdList :
        os.system(i)

#########################################################################
#
#  usage: python plotter_prefit.py --hadd 1 --output OUT --input output
#
#  - the --input should be the output of runAnalysisOn*.py 0 0
#  - the --hadd 1 is needed only the first time the input is processed
#    and it prepare the input for the stacked plots
#  - skipSyst allow to skip some syst category in the band plotting
#  - systComp do the comparsion between various syst
#########################################################################

    
parser = argparse.ArgumentParser("")
parser.add_argument('-a','--hadd', type=int, default=True,help="hadd of the output of RDF")
parser.add_argument('-o','--output', type=str, default='TEST',help="name of the output directory")
parser.add_argument('-i','--input', type=str, default='TEST',help="name of the input direcory with HADDED root files/if used with HADD=1, will be set to output dir")
parser.add_argument('-s','--skipSyst', type=str, default='',nargs='*', help="list of skipped syst class as in bkgAnalysis/bkg_utils.py, separated by space")
parser.add_argument('-c','--systComp', type=int, default=True,help="systematic uncertainity comparison plots")
parser.add_argument('-sb', '--SBana',type=int, default=False, help="run also on the sideband (clousure test)")
parser.add_argument('-bos', '--bkgOnlySyst',type=int, default=False, help="run with syst only for Bkg")


args = parser.parse_args()
HADD = args.hadd
OUTPUT = args.output
INPUT = args.input
skippedSyst =args.skipSyst
SYSTCOMP = args.systComp
SBana = args.SBana
BKGONLYSYST = args.bkgOnlySyst

if HADD :
    prepareHistos(inDir=INPUT,outDir=OUTPUT)
    INPUT=OUTPUT

p=plotter(outDir=OUTPUT, inDir = INPUT,bkgOnlySyst=BKGONLYSYST)
p.plotStack(skipSyst=skippedSyst)
if SYSTCOMP :
    p.plotSyst(skipSyst=skippedSyst)

if SBana :
    pSB=plotter(outDir=OUTPUT, inDir = INPUT,SBana=SBana,bkgOnlySyst=BKGONLYSYST)
    pSB.plotStack(skipSyst=skippedSyst)
    if SYSTCOMP :
        pSB.plotSyst(skipSyst=skippedSyst)
