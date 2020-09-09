import ROOT
import pickle
from termcolor import colored
import math
from HiggsAnalysis.CombinedLimit.DatacardParser import *
from collections import OrderedDict
import copy

class fitUtils:
    def __init__(self, fsig, fmap, fbkg = {}, channel ="WPlus"):
        
        self.templates2D = {}
        self.processes = []
        self.signals = []

        #combine utils
        self.channel = channel
        self.shapeMap = {}
        self.helGroups = OrderedDict()
        self.sumGroups = OrderedDict()
        self.helMetaGroups = OrderedDict()
        
        self.templSystematics = {
            "Nominal" : [""],
            "mass" : ["mass"],
            #"WHSFVars"  : ["WHSFSyst0", "WHSFSyst1","WHSFSyst2","WHSFSystFlat"],
            #"LHEScaleWeightVars" : ["LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0","LHEScaleWeight_muR1p0_muF0p5","LHEScaleWeight_muR1p0_muF2p0","LHEScaleWeight_muR2p0_muF1p0", "LHEScaleWeight_muR2p0_muF2p0"],
            #"ptScaleVars" : [ "corrected"], 
            #"jmeVars" : ["jesTotal", "unclustEn"],
            #"LHEPdfWeightVars" : ["LHEPdfWeightHess{}".format(i+1) for i in range(60)]
        }
        
        #all the files that are needed
        self.fmap = ROOT.TFile.Open(fmap) #file containing the angular coefficient values and inclusive pt-y map
        self.fbkg = fbkg
        self.fsig = ROOT.TFile.Open(fsig)
        
        #get the inclusive pt-y map to unfold
        self.imap = self.fmap.Get("accMaps/mapTot")
        self.xsec = self.fmap.Get("accMaps/sumw")

    def getTemplates(self):
        
        for key in self.fsig.Get("Nominal").GetListOfKeys():
            if not 'mass' in key.GetName():
                self.templates2D[key.GetName()] = {}
                self.templates2D[key.GetName()]['Nominal']=[]
                self.templates2D[key.GetName()]['mass']=[]
                self.processes.append(key.GetName())
                if not "helXsecs7" in key.GetName() and not "helXsecs8" in key.GetName() and not "helXsecs9" in key.GetName():
                    self.signals.append(key.GetName())
                self.templates2D[key.GetName()]['Nominal'].append(self.fsig.Get("Nominal").Get(key.GetName()))

        for proc in self.processes:
            self.templates2D[proc]['mass'].append(self.fsig.Get("Nominal").Get(proc+'_massUp'))
            self.templates2D[proc]['mass'].append(self.fsig.Get("Nominal").Get(proc+'_massDown'))

        for proc in self.processes:
            for syst,variations in self.templSystematics.iteritems():
                if 'Nominal' in syst or 'mass' in syst: continue
                self.templates2D[proc][syst] = []
                #dive into file and add the relevant histograms
                if self.fsig.GetDirectory(syst):
                    for var in variations:
                        temp = self.fsig.Get(syst).Get(proc+'_'+var+'Up')
                        self.templates2D[proc][syst].append(temp)
                        temp = self.fsig.Get(syst).Get(proc+'_'+var+'Down')
                        self.templates2D[proc][syst].append(temp)


        bkg_list = ["DY","Diboson","Top","Fake","Tau","LowAcc","data_obs"]
        for proc in bkg_list:
            if not self.fbkg[proc]=="":
                aux = ROOT.TFile.Open(self.fbkg[proc])
                self.templates2D[proc] = {}
                for syst,variations in self.templSystematics.iteritems():
                    self.templates2D[proc][syst] = []
                    if 'mass' in syst: continue
                    for var in variations:
                        #dive into file and add the relevant histograms
                        if 'Nominal' in syst:
                            temp = aux.Get(syst+"/"+"templates")
                        else:
                            if aux.Get(syst).GetListOfKeys().Contains(syst+"/"+"templates_"+var+'Up'):
                                temp = aux.Get(syst+"/"+"templates_"+var+'Up')
                                self.templates2D[proc][syst].append(copy.deepcopy(temp))
                                temp = aux.Get(syst+"/"+"templates_"+var+'Down')
                                self.templates2D[proc][syst].append(copy.deepcopy(temp))
                                
        self.processes.extend(bkg_list)
    def shapeFile(self):

        shapeOut = ROOT.TFile(self.channel+'.root', 'recreate')
                
        for proc in self.processes:
            for syst in self.templSystematics:
                for temp in self.templates2D[proc][syst]:
                               
                    nbins = temp.GetNbinsX()*temp.GetNbinsY()
                    temp.Sumw2()
                    new = temp.GetName()
                    old = new + '_roll'
                    temp.SetName(old)
                    unrolledtemp = ROOT.TH1F(new, '', nbins, 0., nbins+1)
        
                    for ibin in range(1, temp.GetNbinsX()+1):
                        for jbin in range(1, temp.GetNbinsY()+1):
                
                            bin1D = temp.GetBin(ibin,jbin)
                            unrolledtemp.SetBinContent(bin1D, temp.GetBinContent(ibin,jbin))
                            unrolledtemp.SetBinError(bin1D, temp.GetBinError(ibin,jbin))
                
                    shapeOut.cd()
                    unrolledtemp.Write()
                    
                
        shapeOutxsec = ROOT.TFile(self.channel+'_xsec.root', 'recreate')

        self.xsec.Scale(61526.7*1000.) #xsec in fb
        self.xsec.Write()
        
        #just a bunch of useful factors
        self.factors = {}
        self.factors["A0"]= 2.
        self.factors["A1"]=2.*math.sqrt(2)
        self.factors["A2"]=4.
        self.factors["A3"]=4.*math.sqrt(2)
        self.factors["A4"]=2.
        self.factors["A5"]=2.
        self.factors["A6"]=2.*math.sqrt(2)
        self.factors["A7"]=4.*math.sqrt(2)

        self.helXsecs = OrderedDict()
        self.helXsecs["L"] = "A0"
        self.helXsecs["I"] = "A1" 
        self.helXsecs["T"] = "A2" 
        self.helXsecs["A"] = "A3" 
        self.helXsecs["P"] = "A4" 
        self.helXsecs["7"] = "A5" 
        self.helXsecs["8"] = "A6" 
        self.helXsecs["9"] = "A7" 
        self.helXsecs["UL"] = "AUL" 

        for proc in self.processes:
            if proc in self.signals: #give the correct xsec to unfold
                #helXsecsUL_y_6_qt_1_massUp  
                
                iY = int(proc.split('_')[2])
                iQt = int(proc.split('_')[4])
                coeff = proc.split('_')[0].replace('helXsecs','')

                tmp = ROOT.TH1D(proc,proc, 1, 0., 1.)
                cont = self.xsec.GetBinContent(iY,iQt)
                tmp.SetBinContent(1, cont)
                nsum = (3./16./math.pi)
                
                if not "UL" in proc: #rescale for the releative xsec
                    hAC = self.fmap.Get("angularCoefficients/harmonics{}".format(self.helXsecs[coeff]))
                    nsum = nsum*hAC.GetBinContent(iY,iQt)/self.factors[self.helXsecs[coeff]]
                tmp.Scale(nsum)
                shapeOutxsec.cd()
                tmp.Write()
            else:
                tmp = ROOT.TH1D(proc, proc, 1, 0.,1.)
                if proc == "data_obs": tmp.SetBinContent(1, 1.)
                else: tmp.SetBinContent(1, 0.)
                shapeOutxsec.cd()
                tmp.Write()
    
    def fillHelGroup(self):

        for i in range(1, self.imap.GetNbinsX()+1):
            for j in range(1, self.imap.GetNbinsY()+1):

                s = 'y_{i}_pt_{j}'.format(i=i,j=j)
                self.helGroups[s] = []
                
                for hel in self.helXsecs:
                    if 'helXsecs_'+hel+'_'+s in self.signals:

                        self.helGroups[s].append('helXsecs_'+hel+'_'+s)
                                
                if self.helGroups[s] == []:
                    del self.helGroups[s]
    def fillHelMetaGroup(self):

        for i in range(1, self.imap.GetNbinsX()+1):
            s = 'y_{i}'.format(i=i)
            self.helMetaGroups[s] = []
            for key in self.sumGroups:
                if s in key:
                    self.helMetaGroups[s].append(key)
            
            if self.helMetaGroups[s] == []:
                    del self.helMetaGroups[s]
        
        for j in range(1, self.imap.GetNbinsY()+1):
            s = 'pt_{j}'.format(j=j)
            self.helMetaGroups[s] = []
            for key in self.sumGroups:
                if 'pt' in key and key.split('_')[3]==str(j):
                    self.helMetaGroups[s].append(key)
        
            if self.helMetaGroups[s] == []:
                    del self.helMetaGroups[s]
    def fillSumGroup(self):

        for i in range(1, self.imap.GetNbinsX()+1):
            s = 'y_{i}'.format(i=i)
            for hel in self.helXsecs:
                for signal in self.signals:
                    if 'helXsecs_'+hel+'_'+s in signal:
                        self.sumGroups['helXsecs_'+hel+'_'+s] = []
                        for j in range(1, self.imap.GetNbinsY()+1):
                            if 'helXsecs_'+hel+'_'+'y_{i}_pt_{j}'.format(i=i,j=j) in self.signals:
                                self.sumGroups['helXsecs_'+hel+'_'+s].append('helXsecs_'+hel+'_'+s+'_pt_{j}'.format(j=j))
        
        for j in range(1, self.imap.GetNbinsY()+1):
            s = 'pt_{j}'.format(j=j)
            for hel in self.helXsecs:
                for signal in self.signals:
                    if signal.split('_')[1]==hel and signal.split('_')[5]==str(j):
                        self.sumGroups['helXsecs_'+hel+'_'+s] = []
                        for i in range(1, self.imap.GetNbinsX()+1):
                            if 'helXsecs_'+hel+'_'+'y_{i}_pt_{j}'.format(i=i,j=j) in self.signals:
                            #print i, signal, 'helXsecs_'+hel+'_'+'y_{i}_pt_{j}'.format(i=i,j=j)
                            #print 'append', 'helXsecs_'+hel+'_y_{i}_'.format(i=i)+s, 'to', 'helXsecs_'+hel+'_'+s
                                self.sumGroups['helXsecs_'+hel+'_'+s].append('helXsecs_'+hel+'_y_{i}_'.format(i=i)+s)
        #print self.sumGroups
    def makeDatacard(self):

        self.DC = Datacard()

        ############## Setup the datacard (must be filled in) ###########################

        self.DC.bins =   [self.channel, self.channel+'_xsec'] # <type 'list'>
        self.DC.obs =    {} # <type 'dict'>
        self.DC.processes =  self.processes # <type 'list'>
        self.DC.signals =    self.signals # <type 'list'>
        self.DC.isSignal =   {} # <type 'dict'>
        for proc in self.processes:
            if proc in self.signals:
                self.DC.isSignal[proc] = True
            else:
                self.DC.isSignal[proc] = False
        self.DC.keyline = [] # <type 'list'> # not used by combine-tf
        self.DC.exp =    {} # <type 'dict'>
        self.DC.exp[self.channel] = {}
        self.DC.exp[self.channel+'_xsec'] = {}
        for proc in self.processes:
            self.DC.exp[self.channel][proc] = -1.00
            self.DC.exp[self.channel+'_xsec'][proc] = -1.00
        self.DC.systs =  [] # <type 'list'>
        
        ## list of [{bin : {process : [input file, path to shape, path to shape for uncertainty]}}]
        aux = {}
        aux[self.channel] = {}
        aux[self.channel+'_xsec'] = {}
        for proc in self.processes:
            if "lowAcc" in proc:
                aux[self.channel][proc] = 1.0
                aux[self.channel+'_xsec'][proc] = 0.0
            else:
                aux[self.channel][proc] = 0.0
                aux[self.channel+'_xsec'][proc] = 0.0
        
        for i in range(60):
            self.DC.systs.append(('pdf{}'.format(i+1), False, 'shape', [], aux))
        
        aux2 = {}
        aux2[self.channel] = {}
        aux2[self.channel+'_xsec'] = {}
        for proc in self.processes:
            aux2[self.channel][proc] = 1.0
            aux2[self.channel+'_xsec'][proc] = 0.0
                    
        self.DC.systs.append(('mass', False, 'shape', [], aux2))
        
        self.DC.shapeMap = 	{self.channel: {'*': [self.channel+'.root', '$PROCESS', '']},\
        self.channel+'_xsec': {'*': [self.channel+'_xsec.root', '$PROCESS', '']}} # <type 'dict'>
        self.DC.hasShapes =  True # <type 'bool'>
        self.DC.flatParamNuisances =  {} # <type 'dict'>
        self.DC.rateParams =  {} # <type 'dict'>
        self.DC.extArgs =    {} # <type 'dict'>
        self.DC.rateParamsOrder  =  set([]) # <type 'set'>
        self.DC.frozenNuisances  =  set([]) # <type 'set'>
        self.DC.systematicsShapeMap =  {} # <type 'dict'>
        self.DC.nuisanceEditLines    =  [] # <type 'list'>
        self.DC.groups   =  {'pdfs': set(['pdf{}'.format(i+1) for i in range(60)])} # <type 'dict'>
        self.DC.discretes    =  [] # <type 'list'>
        self.DC.helGroups = self.helGroups
        self.DC.sumGroups = self.sumGroups
        self.DC.helMetaGroups = self.helMetaGroups

        filehandler = open('{}.pkl'.format(self.channel), 'w')
        pickle.dump(self.DC, filehandler)
