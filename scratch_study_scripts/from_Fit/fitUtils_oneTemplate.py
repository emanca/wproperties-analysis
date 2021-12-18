import ROOT
import pickle
from termcolor import colored
import math
from HiggsAnalysis.CombinedLimit.DatacardParser import *
from collections import OrderedDict
import copy
from systToapply import systematicsDict
import numpy as np

class fitUtils:
    def __init__(self, fmap, channel ="WPlus", doSyst=False):
        
        self.nBinsQt = 1#9
        self.nBinsY = 1#6
        
        self.doSyst = doSyst
        self.processes = []
        self.signals = []

        #combine utils
        self.channel = channel
        self.shapeMap = {}
        self.helGroups = OrderedDict()
        self.sumGroups = OrderedDict()
        self.helMetaGroups = OrderedDict()
        self.poly1DRegGroups = OrderedDict()
        self.poly2DRegGroups = OrderedDict()
        
        self.templSystematics = systematicsDict

        #all the files that are needed
        self.fmap = ROOT.TFile.Open(fmap) #file containing the angular coefficient values and inclusive pt-y map
        #get the inclusive pt-y map to unfold
        self.imap = self.fmap.Get("accMaps/mapTot")
        self.xsec = self.fmap.Get("accMaps/sumw")
        #just a bunch of useful factors
        self.factors = {}
        self.factors["A0"] = 2.
        self.factors["A1"] = 2.*math.sqrt(2)
        self.factors["A2"] = 4.
        self.factors["A3"] = 4.*math.sqrt(2)
        self.factors["A4"] = 2.
        self.factors["A5"] = 2.
        self.factors["A6"] = 2.*math.sqrt(2)
        self.factors["A7"] = 4.*math.sqrt(2)

        self.helXsecs = OrderedDict()
        self.helXsecs["UL"] = "AUL"
    def fillProcessList(self):
        for hel in self.helXsecs:
            for i in range(1, self.nBinsY+1): #binsY
                for j in range(1, self.nBinsQt+1): #binsPt
                    proc = 'helXsecs' + hel + '_y_{}'.format(i)+'_qt_{}'.format(j)
                    self.processes.append(proc)
                    if not "helXsecs7" in proc and not "helXsecs8" in proc and not "helXsecs9" in proc:
                        self.signals.append(proc)
        # bkg_list = ["DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"]  # bkg_list = ["Fake","LowAcc"] # bkg_list = ["LowAcc"]
        bkg_list = []
        self.processes.extend(bkg_list)
    def shapeFile(self):
        
        if self.doSyst :
            xsecVarDict = {
            '_mass': ['_massUp','_massDown'],
            # '_LHEPdfWeight' : set(['_LHEPdfWeightHess{}Up'.format(i+1) for i in range(60)]+['_LHEPdfWeightHess{}Down'.format(i+1) for i in range(60)]),
            # '_alphaS' : ['_alphaSUp', '_alphaSDown'],
            '' : ['']
            }
        else :
            xsecVarDict = {'' : ['']}

        shapeOutxsec = ROOT.TFile(self.channel+'_xsec.root', 'recreate')

        # self.xsec.Scale(61526.7*1000.*35.9) #xsec in fb x integrated luminosity
        self.xsec = self.fmap.Get("accMaps/mapTot") #xsec in fb x integrated luminosity
        self.xsec.Write()
        
        NevTot ={}
        for sKind,sList in xsecVarDict.items() :
            if sKind=='' : continue 
            for sName in sList :
                NevTot[sKind+sName] = self.fmap.Get("angularCoefficients"+sKind+"/mapTot"+sName)
                # NevTot[sKind+sName].Scale(61526.7*1000.*35.9) #xsec in fb x integrated luminosity
                NevTot[sKind+sName].Write()
                # NevTot[sKind+sName] = self.xsec#temporary to test
        print("PROCESSES:", self.processes)
        print("PROCESSES:", self.signals)
        for proc in self.processes:
            if proc in self.signals: #give the correct xsec to unfold
                for sKind,sList in xsecVarDict.items() :
                    for sName in sList :
                        # for ud in ['Up','Down'] :
                        #     if sName =='nom' and ud=='Down' : continue
                        #     if sName =='nom' : ud=''
                            
                        iY = int(proc.split('_')[2])
                        iQt = int(proc.split('_')[4])
                        coeff = proc.split('_')[0].replace('helXsecs','')
                        # print('coeff=',coeff, "----", proc, "---",sName)

                        # tmp = ROOT.TH1D(proc,proc, 1, 0., 1.)
                        # if 'Up' in sName : sName_h = sName.replace('Up', 'Down') #josh debug, flip up/down in masked only
                        # elif 'Down' in sName :  sName_h = sName.replace('Down', 'Up') #josh debug
                        # else : sName_h = sName #josh debug
                        # tmp = ROOT.TH1D(proc+sName_h,proc+sName_h, 1, 0., 1.) #josh debug
                        tmp = ROOT.TH1D(proc+sName,proc+sName, 1, 0., 1.)
                        # if sName =='' :
                        cont = self.xsec.GetBinContent(iY,iQt)
                        # hsww = self.fmap.Get("accMaps/sumw")
                        # print("DEBUGGG: "+sName+", mapTot",cont, hsww.GetBinContent(iY,iQt), hsww.GetBinContent(iY,iQt)*61526.7*1000.*35.9 )
                        # else :
                        #     cont = NevTot[sKind+sName].GetBinContent(iY,iQt)
                        #     print("DEBUGGG: "+sName+", mapTot",cont)

                        tmp.SetBinContent(1, cont)
                        nsum = (3./16./math.pi)
                        
                        if not "UL" in proc: #rescale for the releative xsec
                            print('wtf!!!')
                        else : #josh debug
                            if 'Up' in sName :
                                nsum = nsum*(1.1)
                            if 'Down' in sName :
                                nsum = nsum*(1./1.1)  
                        print("debug xsec value--->",nsum)
                        tmp.Scale(nsum)
                        shapeOutxsec.cd()
                        tmp.Write()
            else:
                for sKind,sList in xsecVarDict.items() :
                    for sName in sList :
                        tmp = ROOT.TH1D(proc+sName, proc+sName, 1, 0.,1.)
                        if proc == "data_obs": tmp.SetBinContent(1, 1.)
                        else: tmp.SetBinContent(1, 0.)
                        shapeOutxsec.cd()
                        tmp.Write()
    def fillHelGroup(self):

        for i in range(1, self.imap.GetNbinsX()+1):
            for j in range(1, self.imap.GetNbinsY()+1):

                s = 'y_{i}_qt_{j}'.format(i=i,j=j)
                self.helGroups[s] = []
                
                for hel in self.helXsecs:
                    if 'helXsecs'+hel+'_'+s in self.signals:

                        self.helGroups[s].append('helXsecs'+hel+'_'+s)
                                
                if self.helGroups[s] == []:
                    del self.helGroups[s]
    def fillHelMetaGroup(self):

        for i in range(1, self.nBinsY+1):
            s = 'y_{i}'.format(i=i)
            self.helMetaGroups[s] = []
            for key in self.sumGroups:
                if s in key:
                    self.helMetaGroups[s].append(key)
            
            if self.helMetaGroups[s] == []:
                    del self.helMetaGroups[s]
        
        for j in range(1, self.nBinsQt+1):
            s = 'qt_{j}'.format(j=j)
            self.helMetaGroups[s] = []
            for key in self.sumGroups:
                if 'qt' in key and key.split('_')[2]==str(j):
                    self.helMetaGroups[s].append(key)
        
            if self.helMetaGroups[s] == []:
                    del self.helMetaGroups[s]
        #print self.helMetaGroups
    def fillSumGroup(self):

        for i in range(1, self.nBinsY+1):
            s = 'y_{i}'.format(i=i)
            for hel in self.helXsecs:
                for signal in self.signals:
                    if 'helXsecs'+hel+'_'+s in signal:
                        self.sumGroups['helXsecs'+hel+'_'+s] = []
                        for j in range(1, self.nBinsQt+1):
                            if 'helXsecs'+hel+'_'+'y_{i}_qt_{j}'.format(i=i,j=j) in self.signals:
                                self.sumGroups['helXsecs'+hel+'_'+s].append('helXsecs'+hel+'_'+s+'_qt_{j}'.format(j=j))
        
        for j in range(1, self.nBinsQt+1):
            s = 'qt_{j}'.format(j=j)
            for hel in self.helXsecs:
                for signal in self.signals:
                    if signal.split('_')[0] == 'helXsecs'+hel and signal.split('_')[4] == str(j):
                        self.sumGroups['helXsecs'+hel+'_'+s] = []
                        for i in range(1, self.nBinsY+1):
                            if 'helXsecs'+hel+'_'+'y_{i}_qt_{j}'.format(i=i,j=j) in self.signals:
                            #print i, signal, 'helXsecs'+hel+'_'+'y_{i}_pt_{j}'.format(i=i,j=j)
                            #print 'append', 'helXsecs'+hel+'_y_{i}_'.format(i=i)+s, 'to', 'helXsecs'+hel+'_'+s
                                self.sumGroups['helXsecs'+hel+'_'+s].append('helXsecs'+hel+'_y_{i}_'.format(i=i)+s)
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
        if self.doSyst:
            for syst in self.templSystematics: #loop over systematics
                if 'Nominal' in syst: continue
                for var in self.templSystematics[syst]["vars"]:
                    aux = {} #each sys will have a separate aux dict
                    aux[self.channel] = {}
                    aux[self.channel+'_xsec'] = {}
                    for proc in self.processes: 
                        if proc in self.templSystematics[syst]["procs"]:
                            aux[self.channel][proc] = self.templSystematics[syst]["weight"]
                            aux[self.channel+'_xsec'][proc] = 0.0
                        else:
                            if "Signal" in self.templSystematics[syst]["procs"] and "hel" in proc:
                                aux[self.channel][proc] = self.templSystematics[syst]["weight"]
                                if syst in ["alphaS", "LHEPdfWeight", "mass"] : #theo nuis. applied to signal
                                    print("xsec variation")
                                    aux[self.channel+'_xsec'][proc] = self.templSystematics[syst]["weight"]
                                else :
                                    aux[self.channel+'_xsec'][proc] = 0.0
                            else:
                                print("chemminchia e: ",syst, proc, self.templSystematics[syst]["procs"])
                                aux[self.channel][proc] = 0.0
                                aux[self.channel+'_xsec'][proc] = 0.0

                    self.DC.systs.append((var, False, self.templSystematics[syst]["type"], [], aux))
        
        self.DC.groups = {'mass': ['mass'],
                        #  'pdfs': set(['LHEPdfWeightHess{}'.format(i+1) for i in range(60)]+['alphaS']),
                        #  'WHSFStat': set(["WHSFSyst0Eta{}".format(i) for i in range(1, 49)]+["WHSFSyst1Eta{}".format(i) for i in range(1, 49)]+["WHSFSyst2Eta{}".format(i) for i in range(1, 49)]),
                        #  'WHSFSyst': ['WHSFSystFlat'],
                        #  'jme': set(['jesTotal', 'unclustEn']),
                        #  'PrefireWeight':['PrefireWeight'],
                        #   'CMSlumi' :['CMSlumi'],
                        #   "DYxsec" :["DYxsec"],
                        #   "Topxsec" :["Topxsec"],
                        #   "Dibosonxsec" :["Dibosonxsec"],
                        #   "Tauxsec" :["Tauxsec"],
                        #   "LeptonVeto" : ["LeptonVeto"],
                        #   "FakeNorm" :['FakeNorm'],
                        #   "LHEScaleWeight" : ["LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0","LHEScaleWeight_muR1p0_muF0p5","LHEScaleWeight_muR1p0_muF2p0","LHEScaleWeight_muR2p0_muF1p0","LHEScaleWeight_muR2p0_muF2p0"],
                        #   "WQT" : ["LHEScaleWeight_muR0p5_muF0p5_WQTlow", "LHEScaleWeight_muR0p5_muF1p0_WQTlow","LHEScaleWeight_muR1p0_muF0p5_WQTlow","LHEScaleWeight_muR1p0_muF2p0_WQTlow","LHEScaleWeight_muR2p0_muF1p0_WQTlow","LHEScaleWeight_muR2p0_muF2p0_WQTlow", 
                        #         "LHEScaleWeight_muR0p5_muF0p5_WQTmid", "LHEScaleWeight_muR0p5_muF1p0_WQTmid","LHEScaleWeight_muR1p0_muF0p5_WQTmid","LHEScaleWeight_muR1p0_muF2p0_WQTmid","LHEScaleWeight_muR2p0_muF1p0_WQTmid","LHEScaleWeight_muR2p0_muF2p0_WQTmid", 
                        #         "LHEScaleWeight_muR0p5_muF0p5_WQThigh", "LHEScaleWeight_muR0p5_muF1p0_WQThigh","LHEScaleWeight_muR1p0_muF0p5_WQThigh","LHEScaleWeight_muR1p0_muF2p0_WQThigh","LHEScaleWeight_muR2p0_muF1p0_WQThigh","LHEScaleWeight_muR2p0_muF2p0_WQThigh"],
                        #   'ptScale': set(["Eta{}zptsyst".format(j) for j in range(1, 5)] + ["Eta{}Ewksyst".format(j) for j in range(1, 5)] + ["Eta{}deltaMsyst".format(j) for j in range(1, 5)]+["Eta{}stateig{}".format(j, i) for i in range(0, 99) for j in range(1, 5)]),
                        #  'pdfs': ['LHEPdfWeightHess3'],
                         }  # <type 'dict'>
        
        self.DC.shapeMap = 	{self.channel: {'*': [self.channel+'.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']},\
        self.channel+'_xsec': {'*': [self.channel+'_xsec.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']}} # <type 'dict'>
        self.DC.hasShapes =  True # <type 'bool'>
        self.DC.flatParamNuisances =  {} # <type 'dict'>
        self.DC.rateParams =  {} # <type 'dict'>
        self.DC.extArgs =    {} # <type 'dict'>
        self.DC.rateParamsOrder  =  set([]) # <type 'set'>
        self.DC.frozenNuisances  =  set([]) # <type 'set'>
        self.DC.systematicsShapeMap =  {} # <type 'dict'>
        self.DC.nuisanceEditLines    =  [] # <type 'list'>
        self.DC.discretes    =  [] # <type 'list'>
        self.DC.helGroups = self.helGroups
        self.DC.sumGroups = self.sumGroups
        self.DC.helMetaGroups = self.helMetaGroups
        self.DC.noiGroups = {'mass':['mass']}

        coeff = [0,2]
        for j in coeff:
            testnames = []
            for i in range(1, self.nBinsY+1):
                testnames.append("y_{}_helmeta_A{}".format(i,j))

            etas = [0.2, 0.6, 1.0, 1.4, 1.8, 2.2]


            #self.poly1DRegGroups["poly1dyA{}".format(j)] = {"names": testnames, "bincenters": etas, "firstorder": 0, "lastorder": 2}

            testnames = []
            for i in range(1, self.nBinsQt+1):
                testnames.append("qt_{}_helmeta_A{}".format(i, j))

            pts = [2., 6., 10., 14., 18., 22., 26., 30.]


            #self.poly1DRegGroups["poly1dqtA{}".format(j)] = {"names": testnames, "bincenters": pts, "firstorder": 1, "lastorder": 3}

        coeff = [1, 3, 4]
        for j in coeff:
            testnames = []
            for i in range(1, self.nBinsY+1):
                testnames.append("y_{}_helmeta_A{}".format(i, j))

            etas = [0.2, 0.6, 1.0, 1.4, 1.8, 2.2]


            #self.poly1DRegGroups["poly1dyA{}".format(j)] = {"names": testnames, "bincenters": etas, "firstorder": 1, "lastorder": 2}
        
        coeff = [1, 3]
        for j in coeff:
            testnames = []
            for i in range(1, self.nBinsQt+1):
                testnames.append("qt_{}_helmeta_A{}".format(i, j))

            pts = [2., 6., 10., 14., 18., 22., 26., 30.]


            #self.poly1DRegGroups["poly1dqtA{}".format(j)] = {"names": testnames, "bincenters": pts, "firstorder": 1, "lastorder": 3}

        testnames = []
        for i in range(1, self.nBinsQt+1):
            testnames.append("qt_{}_helmeta_A4".format(i))

        pts = [2., 6., 10., 14., 18., 22., 26., 30.]

        #self.poly1DRegGroups["poly1dqtA4"] = {"names": testnames, "bincenters": pts, "firstorder": 0, "lastorder": 3}
        
        self.DC.poly1DRegGroups = self.poly1DRegGroups

        #etas = np.array([0.2/2.4, 0.6/2.4, 1.0/2.4, 1.4/2.4, 1.8/2.4, 2.2/2.4])
        #pts = np.array([2./32., 6./32., 10./32., 14./32., 18./32., 22./32., 26./32., 30./32.])
        etas = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2])
        # etas = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2])
        pts = np.array([1., 3., 5., 7., 9., 11., 14., 19., 27.])

        #etas = etas/2.4
        #pts = pts/32.

        testnames = []
        bincenters = []
        for i in range(self.nBinsY):
            for j in range(self.nBinsQt):
                testnames.append("y_%i_qt_%i_A4" % (i+1, j+1))
                bincenters.append([etas[i], pts[j]])

        self.poly2DRegGroups["poly2dA4"] = {"names": testnames, "bincenters": bincenters, "firstorder": (1, 0), "lastorder": (3, 4), "fullorder": (5, 7)}

        coeff = [1, 3]
        for c in coeff:
            testnames = []
            for i in range(self.nBinsY):
                for j in range(self.nBinsQt):
                    testnames.append("y_%i_qt_%i_A%i" % (i+1, j+1,c))

        self.poly2DRegGroups["poly2dA%i"%c] = {"names": testnames, "bincenters": bincenters, "firstorder": (1, 1), "lastorder": (3, 4), "fullorder": (5, 7)}
        
        coeff = [0, 2]
        for c in coeff:
            testnames = []
            for i in range(self.nBinsY):
                for j in range(self.nBinsQt):
                    testnames.append("y_%i_qt_%i_A%i" % (i+1, j+1, c))

        self.poly2DRegGroups["poly2dA%i" % c] = {"names": testnames, "bincenters": bincenters, "firstorder": (0, 1), "lastorder": (3, 4), "fullorder": (5, 7)}
        
        self.DC.poly2DRegGroups = self.poly2DRegGroups
        filehandler = open('{}.pkl'.format(self.channel), 'w')
        pickle.dump(self.DC, filehandler)
