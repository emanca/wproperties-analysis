import ROOT
import pickle
from termcolor import colored
import math
from HiggsAnalysis.CombinedLimit.DatacardParser import *
#import CombineHarvester.CombineTools.ch as ch
from collections import OrderedDict
import copy

class fitUtils:
    def __init__(self, ftemplates, fmap, fbkg=[]):
        
        self.ftemplates = ROOT.TFile.Open(ftemplates) #file containing the N-dimensional templates
        self.fmap = ROOT.TFile.Open(fmap) #file containing the angular coefficient values and inclusive pt-y map
        self.fbkg = fbkg
        self.templates2D = {}
        
        if not fbkg == []:
            
            self.templates2D["DY"] = []
            self.templates2D["Diboson"] = []
            self.templates2D["Top"] = []
            self.templates2D["Fake"] = []
            self.templates2D["Tau"] = []
            
            for f in fbkg:
                aux = ROOT.TFile.Open(f)
                if "Fake" in f:
                    faux = aux.Get('templates_fakes/Nominal')
                else:
                    faux = aux.Get('templates_Signal/Nominal')
                templ = ROOT.TH3D
                templ = faux.Get('template')
                if "DY" in f:
                    templ.GetZaxis().SetRange(2,2)
                    taux = templ.Project3D("yxe")
                    taux.SetName("DY")
                    self.templates2D["DY"].append(copy.deepcopy(taux))
                elif "Diboson" in f:
                    templ.GetZaxis().SetRange(2,2)
                    taux = templ.Project3D("yxe")
                    taux.SetName("Diboson")
                    self.templates2D["Diboson"].append(copy.deepcopy(taux))
                elif "TT" in f:
                    templ.GetZaxis().SetRange(2,2)
                    taux = templ.Project3D("yxe")
                    taux.SetName("Top")
                    self.templates2D["Top"].append(copy.deepcopy(taux))
                elif "Fake" in f:
                    templ.GetZaxis().SetRange(2,2)
                    taux = templ.Project3D("yxe")
                    taux.SetName("Fake")
                    self.templates2D["Fake"].append(copy.deepcopy(taux))
                elif "Tau" in f:
                    templ.GetZaxis().SetRange(2,2)
                    taux = templ.Project3D("yxe")
                    taux.SetName("Tau")
                    self.templates2D["Tau"].append(copy.deepcopy(taux))

        #get a list of the 3D templates
        self.templates3D = []
        for key in self.ftemplates.Get('templates').GetListOfKeys():
            if not 'TH3' in key.GetClassName(): continue
            #print key.GetName()
            self.templates3D.append(self.ftemplates.Get('templates').Get(key.GetName()))
        
        #get the inclusive pt-y map
        self.imap = ROOT.TH2D
        self.imap = self.fmap.Get("accMaps/mapTot")
        #use this histogram to check that histograms are correctly normalised
        self.closure = copy.deepcopy(ROOT.TH2D("clos", "clos", self.imap.GetXaxis().GetNbins(),self.imap.GetXaxis().GetXbins().GetArray(),self.imap.GetYaxis().GetNbins(),self.imap.GetYaxis().GetXbins().GetArray()))
        self.xsec = ROOT.TH2D
        self.xsec = self.fmap.Get("accMaps/sumw")

        self.templates1D = {}
        #organise the projected templates into helXsecs
        self.clist = ["L", "I", "T", "A", "P", "7", "8", "9", "UL"]
        for c in self.clist:
            self.templates2D[c] = []
            self.templates1D[c] = []
        #add low acceptance templates
        self.templates2D["lowAcc"] = []
        self.templates2D["lowAcc"].append(copy.deepcopy(self.ftemplates.Get('dataObs').Get('lowAcc')))
        #add data
        self.templates2D["data_obs"] = []
        self.templates2D["data_obs"].append(copy.deepcopy(self.ftemplates.Get('dataObs').Get('data_obs')))
        
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
        
        #combine utils
        self.shapeFile='testbkg'
        self.processes = []
        self.signals = []
        self.shapeMap = {}
        self.helGroups = OrderedDict()
        self.sumGroups = OrderedDict()
        self.helMetaGroups = OrderedDict()
        
    def project3Dto2D(self):
        # returns a list of th2 ordered by rapidity bin
        for th3 in self.templates3D:
            th3.Sumw2()
        
            for ibin in range(1, th3.GetNbinsZ()+1):

                lowEdgeY = th3.GetZaxis().GetBinLowEdge(ibin)
                upEdgeY = th3.GetZaxis().GetBinUpEdge(ibin)

                th3.GetZaxis().SetRange(ibin, ibin)
                proj = th3.Project3D("y_{ibin}_yxe".format(ibin=ibin))
                
                name = proj.GetName()
                name = name.replace('_yx', '')

                coeff = name.split('_')[3]
                jbin = int(name.split('_')[1])
                ibin = int(name.split('_')[5])
                
                new = 'helXsecs_'+coeff+'_y_{}'.format(ibin)+'_pt_{}'.format(jbin)
                #print colored(new,'blue')
                proj.SetName(new)

                if proj.Integral()==0.0: #debug line
                    #print colored(proj.GetName(),'red')
                    continue 

                self.templates2D[coeff].append(proj)
    
    def symmetrisePDF(self):
        for i in range(60):
            th2var = ROOT.TH2D()
            th2var = self.ftemplates.Get('dataObs_LHEPdfWeight').Get('lowAcc_LHEPdfWeightHess{}'.format(i+1))
            th2central = self.ftemplates.Get('dataObs').Get('lowAcc')
            th2c = th2central.Clone()
            th2varD = th2var.Clone()
            th2var.Divide(th2c)
            th2c.Divide(th2varD)
            
            th2Up = ROOT.TH2D("up","up",th2central.GetXaxis().GetNbins(),th2central.GetXaxis().GetXbins().GetArray(),th2central.GetYaxis().GetNbins(),th2central.GetYaxis().GetXbins().GetArray())
            th2Down = ROOT.TH2D("down","down",th2central.GetXaxis().GetNbins(),th2central.GetXaxis().GetXbins().GetArray(),th2central.GetYaxis().GetNbins(),th2central.GetYaxis().GetXbins().GetArray())
            
            for j in range(1,th2central.GetNbinsX()+1):
                for k in range(1,th2central.GetNbinsY()+1):
                
                    th2Up.SetBinContent(j,k,th2central.GetBinContent(j,k)*th2var.GetBinContent(j,k))
                    th2Down.SetBinContent(j,k,th2central.GetBinContent(j,k)*th2c.GetBinContent(j,k))
                
            th2Up.SetName('lowAcc_pdf{}Up'.format(i+1))
            self.templates2D["lowAcc"].append(copy.deepcopy(th2Up))
            
            th2Down.SetName('lowAcc_pdf{}Down'.format(i+1))
            self.templates2D["lowAcc"].append(copy.deepcopy(th2Down))
            
    def unrollTemplates(self):
        # returns a th1 out of a th2
        shapeOut = ROOT.TFile(self.shapeFile+'.root', 'recreate')
  
        for kind,templList in self.templates2D.iteritems():
            for templ in templList:
                nbins = (templ.GetNbinsX()+2)*(templ.GetNbinsY()+2)
                templ.Sumw2()
                new = templ.GetName()
                old = new + '_roll'
                templ.SetName(old)
                
                unrolledtempl = ROOT.TH1F(new, '', nbins-1, 0., nbins-1)
        
                for ibin in range(1, templ.GetNbinsX()+1):
                    for ybin in range(1, templ.GetNbinsY()+1):
                
                        bin1D = templ.GetBin(ibin,ybin)
                        unrolledtempl.SetBinContent(bin1D, templ.GetBinContent(ibin,ybin))
                        unrolledtempl.SetBinError(bin1D, templ.GetBinError(ibin,ybin))
                
                #if not "data_obs" in kind and not "DY" in kind and not "Fake" in kind and not "Diboson" in kind and not "Top" in kind and not "Tau" in kind:
                if not "data_obs" in kind and not "lowAcc" in kind:
                    self.processes.append(unrolledtempl.GetName())
                if kind in self.clist:
                    self.normTempl(unrolledtempl)
                    self.templates1D[kind].append(copy.deepcopy(unrolledtempl))
                    if not "helXsecs_7" in unrolledtempl.GetName() and not "helXsecs_8" in unrolledtempl.GetName() and not "helXsecs_9" in unrolledtempl.GetName() and not unrolledtempl.Integral() == 0.0:
                    #if not "helXsecs_7" in unrolledtempl.GetName() and not "helXsecs_8" in unrolledtempl.GetName() and not "helXsecs_9" in unrolledtempl.GetName() and unrolledtempl.Integral() > 100000:
                        self.signals.append(unrolledtempl.GetName())
                
                shapeOut.cd()
                unrolledtempl.Write()
        
    def normTempl(self, t):

        name = t.GetName()
        #print name, "normalising"
        coeff = name.split('_')[1]
        jbin = int(name.split('_')[5])
        ibin = int(name.split('_')[3])

        nsum = (3./16./math.pi)*self.imap.GetBinContent(ibin,jbin)

        if not "UL" in coeff:
            hAC = ROOT.TH2D
            hAC = self.fmap.Get("angularCoefficients/harmonics{}".format(self.helXsecs[coeff]))
            nsum = nsum*hAC.GetBinContent(ibin,jbin)/self.factors[self.helXsecs[coeff]]
        t.Scale(nsum)

    def xsecMap(self):

        shapeOut = ROOT.TFile(self.shapeFile+'_xsec.root', 'recreate')
        shapeOut.cd()

        out = ROOT.TFile('xsec.root', 'recreate')
        out.cd()

        self.xsec.Scale(61526.7*1000.) #xsec in fb
        self.xsec.Write()

        for kind,templList in self.templates2D.iteritems():
            if kind in self.clist:
                for t in templList:
                        
                    name = t.GetName()
                    if "_roll" in name:
                        name = name.replace("_roll","")

                    jbin = int(name.split('_')[5])
                    ibin = int(name.split('_')[3])
                    tmp = ROOT.TH1D(name,name, 1, 0., 1.)
                    cont = self.xsec.GetBinContent(ibin,jbin)
                    tmp.SetBinContent(1, cont)

                    nsum = (3./16./math.pi)
                    if not "UL" in kind:
                        hAC = ROOT.TH2D
                        hAC = self.fmap.Get("angularCoefficients/harmonics{}".format(self.helXsecs[kind]))

                        nsum = nsum*hAC.GetBinContent(ibin,jbin)/self.factors[self.helXsecs[kind]]

                    tmp.Scale(nsum)
                    shapeOut.cd()
                    tmp.Write()
            else:
                tmp = ROOT.TH1D(kind, kind, 1, 0.,1.)
                if kind == "data_obs": tmp.SetBinContent(1, 1.)
                else: tmp.SetBinContent(1, 0.)
                tmp.Write()
                    
        #build closure map
        for i in range(1, self.imap.GetNbinsX()+1):
            for j in range(1, self.imap.GetNbinsY()+1):

                hlist = []

                for kind,templList in self.templates1D.iteritems():
                    for t in templList:
                        name = t.GetName()
                        s = 'y_{i}_pt_{j}'.format(i=i,j=j)
                            
                        if name.replace('helXsecs_'+kind+'_', '') == s:
                            hlist.append(t)
                yields = 0
                for h in hlist:
                    #print colored(h.GetName(), "red"), i,j
                    #print colored(h.Integral(), 'yellow')
                    yields=yields+h.Integral(0, h.GetNbinsX()+2,)
                #print i,j,yields
                self.closure.SetBinContent(i,j,yields)


        self.closure.Divide(self.imap)
        out = ROOT.TFile('clos.root', 'recreate')
        out.cd()
        self.closure.Write()

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

        self.DC.bins =   ['bin1', 'Wplus'] # <type 'list'>
        self.DC.obs =    {} # <type 'dict'>
        self.processes.append("lowAcc")
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
        self.DC.exp['bin1'] = {}
        self.DC.exp['Wplus'] = {}
        for proc in self.processes:
            self.DC.exp['bin1'][proc] = -1.00
            self.DC.exp['Wplus'][proc] = -1.00
        self.DC.systs =  [] # <type 'list'>
        
        ## list of [{bin : {process : [input file, path to shape, path to shape for uncertainty]}}]
        aux = {}
        aux['bin1'] = {}
        aux['Wplus'] = {}
        for proc in self.processes:
            if "lowAcc" in proc:
                aux['bin1'][proc] = 1.0
                aux['Wplus'][proc] = 0.0
            else:
                aux['bin1'][proc] = 0.0
                aux['Wplus'][proc] = 0.0
        
        for i in range(60):
            self.DC.systs.append(('pdf{}'.format(i+1), False, 'shape', [], aux))
        
        aux2 = {}
        aux2['bin1'] = {}
        aux2['Wplus'] = {}
        for proc in self.processes:
            if "lowAcc" in proc:
                aux2['bin1'][proc] = 1.0
                aux2['Wplus'][proc] = 1.0
            else:
                aux2['bin1'][proc] = 0.0
                aux2['Wplus'][proc] = 0.0
        
        self.DC.systs.append(('mass', False, 'shape', [], aux2))
        
        self.DC.shapeMap = 	{'bin1': {'*': [self.shapeFile+'.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']},\
        'Wplus': {'*': [self.shapeFile+'_xsec.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']}} # <type 'dict'>
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


        filehandler = open('{}.pkl'.format(self.shapeFile), 'w')
        pickle.dump(self.DC, filehandler)
