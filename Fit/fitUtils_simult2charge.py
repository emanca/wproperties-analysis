import ROOT
import pickle
from termcolor import colored
import math
from HiggsAnalysis.CombinedLimit.DatacardParser import *
from collections import OrderedDict
import copy
from systToapply_simult2charge import systematicsDict
import numpy as np

class fitUtils:
    def __init__(self, fmap, channel =["Wplus","Wminus"], doSyst=False, systDict=''):
        
        self.nBinsQt = 11#13#9#8#29
        self.nBinsY = 6#7
        
        self.doSyst = doSyst
        self.processes = []
        self.signals = []

        #combine utils
        self.channels = channel
        self.shapeMap = {}
        self.helGroups = OrderedDict()
        self.sumGroups = OrderedDict()
        self.helMetaGroups = OrderedDict()
        self.poly1DRegGroups = OrderedDict()
        self.poly2DRegGroups = OrderedDict()
        
        if systDict == '' :
            self.templSystematics = systematicsDict
        else : 
            self.templSystematics = systDict 
            
        #all the files that are needed
        self.fmap = {}
        for c in self.channels :
            self.fmap[c] = ROOT.TFile.Open(fmap[c]) #file containing the angular coefficient values and inclusive pt-y map
        #get the inclusive pt-y map to unfold
        self.imap = self.fmap[self.channels[0]].Get("accMaps/mapTot") #not charge dependent
        self.xsec = {}
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
        self.helXsecs["L"] = "A0"
        self.helXsecs["I"] = "A1"
        self.helXsecs["T"] = "A2"
        self.helXsecs["A"] = "A3"
        self.helXsecs["P"] = "A4"
        self.helXsecs["7"] = "A5"
        self.helXsecs["8"] = "A6"
        self.helXsecs["9"] = "A7"
        self.helXsecs["UL"] = "AUL"
    def fillProcessList(self):
        for c in self.channels :
            for hel in self.helXsecs:
                for i in range(1, self.nBinsY+1): #binsY
                    for j in range(1, self.nBinsQt+1): #binsPt
                        proc = c+'helXsecs' + hel + '_y_{}'.format(i)+'_qt_{}'.format(j)
                        self.processes.append(proc)
                        if not "helXsecs7" in proc and not "helXsecs8" in proc and not "helXsecs9" in proc:
                            self.signals.append(proc)
            bkg_list = [c+"DYJets",c+"DiBoson",c+"Top",c+"Fake",c+"WtoTau",c+"LowAcc"] 
            # bkg_list = ["DiBoson","Top","Fake","WtoTau"] 
            self.processes.extend(bkg_list)
    def shapeFile(self):
        self.shapeOutxsec = {}
        
        if self.doSyst :
            xsecVarDict = {
            '_mass': ['_massUp','_massDown'],
            '_LHEPdfWeight' : set(['_LHEPdfWeightHess{}Up'.format(i+1) for i in range(60)]+['_LHEPdfWeightHess{}Down'.format(i+1) for i in range(60)]),
            '_alphaS' : ['_alphaSUp', '_alphaSDown'],
            '' : ['']
            }
        else :
            xsecVarDict = {'' : ['']}
        
        for c in self.channels :
            self.shapeOutxsec[c] = ROOT.TFile(c+'_xsec.root', 'recreate')

            # self.xsec.Scale(61526.7*1000.*35.9) #xsec in fb x integrated luminosity
            self.xsec[c] = self.fmap[c].Get("accMaps/mapTot") #xsec in fb x integrated luminosity
            self.xsec[c].Write()
        
            NevTot ={}
            for sKind,sList in xsecVarDict.items() :
                if sKind=='' : continue 
                for sName in sList :
                    NevTot[sKind+sName] = self.fmap[c].Get("angularCoefficients"+sKind+"/mapTot"+sName)
                    # NevTot[sKind+sName].Scale(61526.7*1000.*35.9) #xsec in fb x integrated luminosity
                    NevTot[sKind+sName].Write()
                    # NevTot[sKind+sName] = self.xsec#temporary to test
            
            for proc in self.processes:
                if c not in proc : continue
                if proc in self.signals: #give the correct xsec to unfold
                    for sKind,sList in xsecVarDict.items() :
                        for sName in sList :
                            # for ud in ['Up','Down'] :
                            #     if sName =='nom' and ud=='Down' : continue
                            #     if sName =='nom' : ud=''
                            # print(proc)
                            iY = int(proc.split('_')[2])
                            iQt = int(proc.split('_')[4])
                            coeff = proc.split('_')[0].replace('helXsecs','').replace(c,'')
                            # print('coeff=',coeff, "----", proc, "---",sName)

                            # tmp = ROOT.TH1D(proc,proc, 1, 0., 1.)
                            tmp = ROOT.TH1D(proc+sName,proc+sName, 1, 0., 1.)
                            if sName =='' :
                                cont = self.xsec[c].GetBinContent(iY,iQt)
                            else :
                                cont = NevTot[sKind+sName].GetBinContent(iY,iQt)
                            tmp.SetBinContent(1, cont)
                            nsum = (3./16./math.pi)
                            
                            if not "UL" in proc: #rescale for the releative xsec
                                if sKind!='': 
                                    if 'alpha' not in sName and 'mass' not in sName:
                                        sNameMod = sName.replace('Up','').replace('Down','')
                                    else :
                                        sNameMod = sName
                                    hAC = self.fmap[c].Get("angularCoefficients{}/harmonics{}{}{}".format(sKind,self.helXsecs[coeff],sNameMod,sName))
                                    nsum = nsum*hAC.GetBinContent(iY,iQt)/self.factors[self.helXsecs[coeff]]
                                else : 
                                    hAC = self.fmap[c].Get("angularCoefficients/harmonics{}_nom_nom".format(self.helXsecs[coeff]))
                                    nsum = nsum*hAC.GetBinContent(iY,iQt)/self.factors[self.helXsecs[coeff]]
                    
                            tmp.Scale(nsum)
                            self.shapeOutxsec[c].cd()
                            tmp.Write()
                else:
                    for sKind,sList in xsecVarDict.items() :
                        for sName in sList :
                            tmp = ROOT.TH1D(proc+sName, proc+sName, 1, 0.,1.)
                            if proc == "data_obs": tmp.SetBinContent(1, 1.)
                            else: tmp.SetBinContent(1, 0.)
                            self.shapeOutxsec[c].cd()
                            tmp.Write()
            self.shapeOutxsec[c].Close()
            
    def fillHelGroup(self):
        for c in self.channels :
            for i in range(1, self.imap.GetNbinsX()+1):
                for j in range(1, self.imap.GetNbinsY()+1):

                    s = 'y_{i}_qt_{j}'.format(i=i,j=j)
                    self.helGroups[c+s] = []
                    
                    for hel in self.helXsecs:
                        if c+'helXsecs'+hel+'_'+s in self.signals:

                            self.helGroups[c+s].append(c+'helXsecs'+hel+'_'+s)
                                    
                    if self.helGroups[c+s] == []:
                        del self.helGroups[c+s]
    def fillHelMetaGroup(self):
        for c in self.channels :
            for i in range(1, self.nBinsY+1):
                s = 'y_{i}'.format(i=i)
                self.helMetaGroups[c+s] = []
                for key in self.sumGroups:
                    if s in key and c in key:
                        self.helMetaGroups[c+s].append(key)
                
                if self.helMetaGroups[c+s] == []:
                        del self.helMetaGroups[c+s]
            
            for j in range(1, self.nBinsQt+1):
                s = 'qt_{j}'.format(j=j)
                self.helMetaGroups[c+s] = []
                for key in self.sumGroups:
                    if 'qt' in key and key.split('_')[2]==str(j) and c in key:
                        self.helMetaGroups[c+s].append(key)
            
                if self.helMetaGroups[c+s] == []:
                        del self.helMetaGroups[c+s]
        #print self.helMetaGroups
    def fillSumGroup(self):
        for c in self.channels :
            for i in range(1, self.nBinsY+1):
                s = 'y_{i}'.format(i=i)
                for hel in self.helXsecs:
                    for signal in self.signals:
                        if c+'helXsecs'+hel+'_'+s in signal:
                            self.sumGroups[c+'helXsecs'+hel+'_'+s] = []
                            for j in range(1, self.nBinsQt+1):
                                if c+'helXsecs'+hel+'_'+'y_{i}_qt_{j}'.format(i=i,j=j) in self.signals:
                                    self.sumGroups[c+'helXsecs'+hel+'_'+s].append(c+'helXsecs'+hel+'_'+s+'_qt_{j}'.format(j=j))
            
            for j in range(1, self.nBinsQt+1):
                s = 'qt_{j}'.format(j=j)
                for hel in self.helXsecs:
                    for signal in self.signals:
                        if signal.split('_')[0] == c+'helXsecs'+hel and signal.split('_')[4] == str(j):
                            self.sumGroups[c+'helXsecs'+hel+'_'+s] = []
                            for i in range(1, self.nBinsY+1):
                                if c+'helXsecs'+hel+'_'+'y_{i}_qt_{j}'.format(i=i,j=j) in self.signals:
                                #print i, signal, 'helXsecs'+hel+'_'+'y_{i}_pt_{j}'.format(i=i,j=j)
                                #print 'append', 'helXsecs'+hel+'_y_{i}_'.format(i=i)+s, 'to', 'helXsecs'+hel+'_'+s
                                    self.sumGroups[c+'helXsecs'+hel+'_'+s].append(c+'helXsecs'+hel+'_y_{i}_'.format(i=i)+s)
            #print self.sumGroups
    def makeDatacard(self):

        self.DC = Datacard()

        ############## Setup the datacard (must be filled in) ###########################

        self.DC.bins =   [self.channels[0], self.channels[0]+'_xsec', self.channels[1], self.channels[1]+'_xsec'] # <type 'list'>
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
        for c in self.channels : 
            self.DC.exp[c] = {}
            self.DC.exp[c+'_xsec'] = {}
            for proc in self.processes:
                if c not in proc : continue 
                self.DC.exp[c][proc] = -1.00
                self.DC.exp[c+'_xsec'][proc] = -1.00
        self.DC.systs =  [] # <type 'list'>
        ## list of [{bin : {process : [input file, path to shape, path to shape for uncertainty]}}]
        if self.doSyst:
            for syst in self.templSystematics: #loop over systematics
                if 'Nominal' in syst: continue
                for var in self.templSystematics[syst]["vars"]:
                    aux = {} #each sys will have a separate aux dict
                    for c in self.channels : 
                        aux[c] = {}
                        aux[c+'_xsec'] = {}
                        for proc in self.processes: 
                            if c not in proc : continue
                            if proc in self.templSystematics[syst]["procs"] or proc.replace(c,"") in self.templSystematics[syst]["procs"]: #bkg two or one charge 
                                aux[c][proc] = self.templSystematics[syst]["weight"]
                                aux[c+'_xsec'][proc] = 0.0
                            else:
                                if (("Signal" in self.templSystematics[syst]["procs"] or c+'Signal' in self.templSystematics[syst]["procs"]) and "hel" in proc): #sig two or one charge
                                    aux[c][proc] = self.templSystematics[syst]["weight"]
                                    if syst in ["alphaS", "LHEPdfWeight", "mass"] : #theo nuis. applied to signal
                                        aux[c+'_xsec'][proc] = self.templSystematics[syst]["weight"]
                                    else :
                                        aux[c+'_xsec'][proc] = 0.0                                    
                                else:
                                    aux[c][proc] = 0.0
                                    aux[c+'_xsec'][proc] = 0.0

                    self.DC.systs.append((var, False, self.templSystematics[syst]["type"], [], aux))

        self.DC.groups = {'mass': set(['mass']),
                         'pdfs': set(['LHEPdfWeightHess{}'.format(i+1) for i in range(60)]+['alphaS']),
                         'WHSFStat': set(["WplusWHSFSyst0Eta{}".format(i) for i in range(1, 49)]+["WplusWHSFSyst1Eta{}".format(i) for i in range(1, 49)]+["WplusWHSFSyst2Eta{}".format(i) for i in range(1, 49)]+
                                    ["WminusWHSFSyst0Eta{}".format(i) for i in range(1, 49)]+["WminusWHSFSyst1Eta{}".format(i) for i in range(1, 49)]+["WminusWHSFSyst2Eta{}".format(i) for i in range(1, 49)]),
                         'WHSFSyst': set(['WplusWHSFSystFlat','WminusWHSFSystFlat']),
                         'jme': set(['WplusjesTotal', 'WplusunclustEn','WminusjesTotal', 'WminusunclustEn']),
                         'PrefireWeight':['WplusPrefireWeight','WminusPrefireWeight'],
                          'CMSlumi' :set(['CMSlumi','lumi']),
                          "ewkXsec" : set(["Topxsec","Dibosonxsec","Tauxsec"]),
                          "LeptonVeto" : set(["WplusLeptonVeto","WminusLeptonVeto"]),
                          "WQT" : set(["LHEScaleWeight_muR0p5_muF0p5_WQTlow", "LHEScaleWeight_muR0p5_muF1p0_WQTlow","LHEScaleWeight_muR1p0_muF0p5_WQTlow","LHEScaleWeight_muR1p0_muF2p0_WQTlow","LHEScaleWeight_muR2p0_muF1p0_WQTlow","LHEScaleWeight_muR2p0_muF2p0_WQTlow", 
                                "LHEScaleWeight_muR0p5_muF0p5_WQTmid", "LHEScaleWeight_muR0p5_muF1p0_WQTmid","LHEScaleWeight_muR1p0_muF0p5_WQTmid","LHEScaleWeight_muR1p0_muF2p0_WQTmid","LHEScaleWeight_muR2p0_muF1p0_WQTmid","LHEScaleWeight_muR2p0_muF2p0_WQTmid", 
                                "LHEScaleWeight_muR0p5_muF0p5_WQThigh", "LHEScaleWeight_muR0p5_muF1p0_WQThigh","LHEScaleWeight_muR1p0_muF0p5_WQThigh","LHEScaleWeight_muR1p0_muF2p0_WQThigh","LHEScaleWeight_muR2p0_muF1p0_WQThigh","LHEScaleWeight_muR2p0_muF2p0_WQThigh",
                                "LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0","LHEScaleWeight_muR1p0_muF0p5","LHEScaleWeight_muR1p0_muF2p0","LHEScaleWeight_muR2p0_muF1p0","LHEScaleWeight_muR2p0_muF2p0"]),
                          'ptScale':set(["Wpluscorrected","Wminuscorrected"]), 
                          'QCDnorm' : set(['WplusQCDnorm','WminusQCDnorm'])
                         }  # <type 'dict'>
        self.DC.groups['other'] = set()
        self.DC.groups['allSyst'] = set()
        self.DC.groups['otherNoMass'] = set()
        for k,vals in self.DC.groups.items() :
            if k=='otherNoMass' or k=='other' or k=='allSyst' : continue 
            for val in vals  :
                    self.DC.groups['allSyst'].add(val)
            if k=='pdfs' or k=='WHSFStat' or k=='CMSlumi' or k=='WQT' : 
                continue
            else : 
                for val in vals  :
                    self.DC.groups['other'].add(val)
                if k!='mass' :
                    for val in vals  :
                        self.DC.groups['otherNoMass'].add(val)
        
        self.DC.shapeMap = 	{
            self.channels[0]: {'*': [self.channels[0]+'.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']},\
            self.channels[0]+'_xsec': {'*': [self.channels[0]+'_xsec.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']},\
            self.channels[1]: {'*': [self.channels[1]+'.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']},\
            self.channels[1]+'_xsec': {'*': [self.channels[1]+'_xsec.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']}
            } # <type 'dict'>
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

            pts = [1., 3., 5., 7., 9., 11., 14., 19., 27.]


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

            pts = [1., 3., 5., 7., 9., 11., 14., 19., 27.]


            #self.poly1DRegGroups["poly1dqtA{}".format(j)] = {"names": testnames, "bincenters": pts, "firstorder": 1, "lastorder": 3}

        testnames = []
        for i in range(1, self.nBinsQt+1):
            testnames.append("qt_{}_helmeta_A4".format(i))

        pts = [1., 3., 5., 7., 9., 11., 14., 19., 27.]

        #self.poly1DRegGroups["poly1dqtA4"] = {"names": testnames, "bincenters": pts, "firstorder": 0, "lastorder": 3}

        ################################################################################################        
        # self.DC.poly1DRegGroups = self.poly1DRegGroups

        # #etas = np.array([0.2/2.4, 0.6/2.4, 1.0/2.4, 1.4/2.4, 1.8/2.4, 2.2/2.4])
        # #pts = np.array([2./32., 6./32., 10./32., 14./32., 18./32., 22./32., 26./32., 30./32.])
        # etas = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2])
        # # etas = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2])
        # pts = np.array([1., 3., 5., 7., 9., 11., 14., 19., 27.])

        # #etas = etas/2.4
        # #pts = pts/32.

        # testnames = []
        # bincenters = []
        # for i in range(self.nBinsY):
        #     for j in range(self.nBinsQt):
        #         testnames.append("y_%i_qt_%i_A4" % (i+1, j+1))
        #         bincenters.append([etas[i], pts[j]])

        # self.poly2DRegGroups["poly2dA4"] = {"names": testnames, "bincenters": bincenters, "firstorder": (1, 0), "lastorder": (3, 4), "fullorder": (5, 7)}

        # coeff = [1, 3]
        # for c in coeff:
        #     testnames = []
        #     for i in range(self.nBinsY):
        #         for j in range(self.nBinsQt):
        #             testnames.append("y_%i_qt_%i_A%i" % (i+1, j+1,c))

        # self.poly2DRegGroups["poly2dA%i"%c] = {"names": testnames, "bincenters": bincenters, "firstorder": (1, 1), "lastorder": (3, 4), "fullorder": (5, 7)}
        
        # coeff = [0, 2]
        # for c in coeff:
        #     testnames = []
        #     for i in range(self.nBinsY):
        #         for j in range(self.nBinsQt):
        #             testnames.append("y_%i_qt_%i_A%i" % (i+1, j+1, c))

        # self.poly2DRegGroups["poly2dA%i" % c] = {"names": testnames, "bincenters": bincenters, "firstorder": (0, 1), "lastorder": (3, 4), "fullorder": (5, 7)}
        
        # self.DC.poly2DRegGroups = self.poly2DRegGroups
        # filehandler = open('{}.pkl'.format(self.channel), 'w')
        # pickle.dump(self.DC, filehandler)
        ################################################################################################        

        
        
        
        
        self.DC.poly1DRegGroups = self.poly1DRegGroups
        etas = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2])
        # etas = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.6])
        # pts = np.array([1., 3., 5., 7., 9., 11., 14., 19., 27.])
        # pts = np.array([1., 3., 5., 7., 9., 11., 13., 15., 18., 22., 28., 38., 60.])#Extended
        pts = np.array([1., 3., 5., 7., 9., 11., 14., 18., 23., 30., 46.])#Extended-reduced
        # pts =np.array([0.5, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.5, 13.5, 14.5, 15.5,19.,27.])   # len=32-->29 fitted bin

        
        # etas = etas/2.4
        # pts = pts/32.
        
                
        testnames = []
        bincenters = []
        for i in range(self.nBinsY):
            for j in range(self.nBinsQt):
                testnames.append("y_%i_qt_%i_A0" % (i+1, j+1))
                bincenters.append([etas[i], pts[j]])

        self.poly2DRegGroups["poly2dA0"] = {"names": testnames, "bincenters": bincenters, "firstorder": (0, 1), "lastorder": (2, 3), "fullorder": (len(etas)-1, len(pts)-1)}#full=max number of pars=nbins-1
        # self.poly2DRegGroups["poly2dA0"] = {"names": testnames, "bincenters": bincenters, "firstorder": (0, 1), "lastorder": (2, 3), "fullorder": (5, 8)}#full=max number of pars=nbins-1


        
        testnames = []
        bincenters = []
        for i in range(self.nBinsY):
            for j in range(self.nBinsQt):
                testnames.append("y_%i_qt_%i_A1" % (i+1, j+1))
                bincenters.append([etas[i], pts[j]])

        self.poly2DRegGroups["poly2dA1"] = {"names": testnames, "bincenters": bincenters, "firstorder": (1, 1), "lastorder": (2, 5), "fullorder": (len(etas)-1, len(pts)-1)} 
        # self.poly2DRegGroups["poly2dA1"] = {"names": testnames, "bincenters": bincenters, "firstorder": (1, 1), "lastorder": (3, 3), "fullorder": (5, 8)} 
        
        testnames = []
        bincenters = []
        for i in range(self.nBinsY):
            for j in range(self.nBinsQt):
                testnames.append("y_%i_qt_%i_A2" % (i+1, j+1))
                bincenters.append([etas[i], pts[j]])

        self.poly2DRegGroups["poly2dA2"] = {"names": testnames, "bincenters": bincenters, "firstorder": (0, 1), "lastorder": (1, 4), "fullorder": (len(etas)-1, len(pts)-1)} 
        # self.poly2DRegGroups["poly2dA2"] = {"names": testnames, "bincenters": bincenters, "firstorder": (0, 1), "lastorder": (2, 4), "fullorder": (5, 8)} 
        
        testnames = []
        bincenters = []
        for i in range(self.nBinsY):
            for j in range(self.nBinsQt):
                testnames.append("y_%i_qt_%i_A3" % (i+1, j+1))
                bincenters.append([etas[i], pts[j]])

        self.poly2DRegGroups["poly2dA3"] = {"names": testnames, "bincenters": bincenters, "firstorder": (0, 1), "lastorder": (2, 4), "fullorder": (len(etas)-1, len(pts)-1)} 
        # self.poly2DRegGroups["poly2dA3"] = {"names": testnames, "bincenters": bincenters, "firstorder": (0, 1), "lastorder": (2, 3), "fullorder": (5, 8)} 

        testnames = []
        bincenters = []
        for i in range(self.nBinsY):
            for j in range(self.nBinsQt):
                testnames.append("y_%i_qt_%i_A4" % (i+1, j+1))
                bincenters.append([etas[i], pts[j]])

        self.poly2DRegGroups["poly2dA4"] = {"names": testnames, "bincenters": bincenters, "firstorder": (1, 0), "lastorder": (3, 5), "fullorder": (len(etas)-1, len(pts)-1)} 
        # self.poly2DRegGroups["poly2dA4"] = {"names": testnames, "bincenters": bincenters, "firstorder": (1, 0), "lastorder": (3, 4), "fullorder": (5, 8)} 

        self.DC.poly2DRegGroups = self.poly2DRegGroups
        filehandler = open('Wall.pkl', 'w')
        pickle.dump(self.DC, filehandler)
