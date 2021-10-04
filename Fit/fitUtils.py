import ROOT
import pickle
from termcolor import colored
import math
from HiggsAnalysis.CombinedLimit.DatacardParser import *
from collections import OrderedDict
import copy
from systToapply import systematicsDict
import numpy as np
import sys
sys.path.append('../Common/data')
import h5py

class fitUtils:
    def __init__(self, channels =["Wplus_preVFP","Wplus_postVFP"], doSyst=False):
        
        self.doSyst = doSyst
        self.processes = []
        self.signals = []

        #combine utils
        self.channels = channels
        self.shapeMap = {}
        self.helGroups = OrderedDict()
        self.sumGroups = OrderedDict()
        self.helMetaGroups = OrderedDict()
        
        # self.templSystematics = systematicsDict
        self.templSystematics = {}

        # load files
        self.ftempl = {}
        self.data = {}
        self.templ = {}
        self.templw2 = {}
        self.gen = {}
        self.lowacc = {}
        self.lowaccw2 = {}
        self.Wtau = {}
        self.Wtauw2 = {}
        self.DY = {}
        self.DYw2 = {}
        self.Top = {}
        self.Topw2 = {}
        self.Diboson = {}
        self.Dibosonw2 = {}
        self.fakeslow = {}
        self.fakesloww2 = {}
        self.fakeshigh = {}
        self.fakeshighw2 = {}
        for chan in self.channels:
            self.ftempl[chan] = h5py.File('../Common/shapes{}.hdf5'.format(chan), mode='r+')
            self.data[chan] = self.ftempl[chan]['data_obs'][:]
            self.templ[chan] = self.ftempl[chan]['template'][:]
            self.templw2[chan] = self.ftempl[chan]['template_sumw2'][:]
            self.gen[chan] = self.ftempl[chan]['helicity'][:]
            self.lowacc[chan] = self.ftempl[chan]['lowacc'][:]
            self.lowaccw2[chan] = self.ftempl[chan]['lowacc_sumw2'][:]
            self.Wtau[chan] = self.ftempl[chan]['Wtau'][:]
            self.Wtauw2[chan] = self.ftempl[chan]['Wtau_sumw2'][:]
            self.DY[chan] = self.ftempl[chan]['DY'][:]
            self.DYw2[chan] = self.ftempl[chan]['DY_sumw2'][:]
            self.Top[chan] = self.ftempl[chan]['Top'][:]
            self.Topw2[chan] = self.ftempl[chan]['Top_sumw2'][:]
            self.Diboson[chan] = self.ftempl[chan]['Diboson'][:]
            self.Dibosonw2[chan] = self.ftempl[chan]['Diboson_sumw2'][:]
            self.fakeslow[chan] = self.ftempl[chan]['fakesLowMt'][:]
            self.fakesloww2[chan] = self.ftempl[chan]['fakesLowMt_sumw2'][:]
            self.fakeshigh[chan] = self.ftempl[chan]['fakesHighMt'][:]
            self.fakeshighw2[chan] = self.ftempl[chan]['fakesHighMt_sumw2'][:]

        self.helXsecs = ['L', 'I', 'T', 'A', 'P','7','8', '9','UL']

        # reduce bins to acceptance
        from binning import ptBins, etaBins, mTBins, etaBins, isoBins, chargeBins, yBins, qtBins, cosThetaBins, phiBins
        threshold_y = np.digitize(2.4,yBins)-1
        threshold_qt = np.digitize(60.,qtBins)-1
        self.yBins = np.array(yBins[:threshold_y+1])
        self.qtBins = np.array(qtBins[:threshold_qt+1])
        self.yBinsC = 0.5*(self.yBins[1:]+self.yBins[:-1])
        self.qtBinsC = 0.5*(self.qtBins[1:]+self.qtBins[:-1])
        print(self.yBins,self.qtBins)
        print(self.yBinsC,self.qtBinsC)
        
    def fillProcessList(self):
        for hel in self.helXsecs:
            for i in range(len(self.yBinsC)):
                for j in range(len(self.qtBinsC)):
                    proc = 'helXsecs' + hel + '_y_{}'.format(i)+'_qt_{}'.format(j)
                    self.processes.append(proc)
                    #if "helXsecsUL" in proc:
                    if not "helXsecs7" in proc and not "helXsecs8" in proc and not "helXsecs9" in proc:
                        self.signals.append(proc)
        # bkg_list = ["DY","Diboson","Top","fakesLowMt","fakesHighMt", "Wtau","LowAcc"]
        # bkg_list = ["DY","Diboson","Top","Wtau","LowAcc"]

        bkg_list = []
        self.processes.extend(bkg_list)
    
    def shapeFile(self):
        dtype = 'float64'
        for chan in self.channels:
            with h5py.File('{}.hdf5'.format(chan), mode="w") as f:
                for proc in self.processes:
                    if "helXsecs" in proc:
                        coeff = self.helXsecs.index(proc.split('_')[0].replace('helXsecs',''))
                        iY = int(proc.split('_')[2])
                        iQt = int(proc.split('_')[4])
                        # print(chan,proc,self.templ[chan].shape)
                        dset_templ = f.create_dataset(proc, self.templ[chan][...,-1,0,iY,iQt,coeff].ravel().shape, dtype=dtype)
                        dset_templ[...] = self.templ[chan][...,-1,0,iY,iQt,coeff].ravel()
                        dset_templw2 = f.create_dataset(proc+'_sumw2', self.templw2[chan][...,-1,0,iY,iQt,coeff].ravel().shape, dtype=dtype)
                        dset_templw2[...] = self.templw2[chan][...,-1,0,iY,iQt,coeff].ravel()

                # dset_bkg = f.create_dataset("DY", self.DY.ravel().shape, dtype=dtype)
                # dset_bkg[...] = self.DY.ravel()
                # dset_bkgw2 = f.create_dataset("DY_sumw2", self.DYw2.ravel().shape, dtype=dtype)
                # dset_bkgw2[...] = self.DYw2.ravel()

                # dset_bkg = f.create_dataset("Diboson", self.Diboson.ravel().shape, dtype=dtype)
                # dset_bkg[...] = self.Diboson.ravel()
                # dset_bkgw2 = f.create_dataset("Diboson_sumw2", self.Dibosonw2.ravel().shape, dtype=dtype)
                # dset_bkgw2[...] = self.Dibosonw2.ravel()

                # dset_bkg = f.create_dataset("Top", self.Top.ravel().shape, dtype=dtype)
                # dset_bkg[...] = self.Top.ravel()
                # dset_bkgw2 = f.create_dataset("Top_sumw2", self.Topw2.ravel().shape, dtype=dtype)
                # dset_bkgw2[...] = self.Topw2.ravel()

                # dset_bkg = f.create_dataset("Wtau", self.Wtau.ravel().shape, dtype=dtype)
                # dset_bkg[...] = self.Wtau.ravel()
                # dset_bkgw2 = f.create_dataset("Wtau_sumw2", self.Wtauw2.ravel().shape, dtype=dtype)
                # dset_bkgw2[...] = self.Wtauw2.ravel()

                # dset_bkg = f.create_dataset("fakesLowMt", self.fakeslow.ravel().shape, dtype=dtype)
                # dset_bkg[...] = self.fakeslow.ravel()
                # dset_bkgw2 = f.create_dataset("fakesLowMt_sumw2", self.fakesloww2.ravel().shape, dtype=dtype)
                # dset_bkgw2[...] = self.fakesloww2.ravel()

                # dset_bkg = f.create_dataset("fakesHighMt", self.fakeshigh.ravel().shape, dtype=dtype)
                # dset_bkg[...] = self.fakeshigh.ravel()
                # dset_bkgw2 = f.create_dataset("fakesHighMt_sumw2", self.fakeshighw2.ravel().shape, dtype=dtype)
                # dset_bkgw2[...] = self.fakeshighw2.ravel()

                # dset_bkg = f.create_dataset("LowAcc", self.lowacc[...,-1,0].ravel().shape, dtype=dtype)
                # dset_bkg[...] = self.lowacc[...,-1,0].ravel()
                # dset_bkgw2 = f.create_dataset("LowAcc_sumw2", self.lowacc[...,-1,0].ravel().shape, dtype=dtype)
                # dset_bkgw2[...] = self.lowaccw2[...,-1,0].ravel()

                dset_data = f.create_dataset('data_obs', self.data[chan][...,-1,0].ravel().shape, dtype=dtype)
                dset_data[...] = self.data[chan][...,-1,0].ravel()

    
                # # copy fakes nuisances in shape file
                # for type in ["Up","Down"]:
                #     for i in range(48*30):
                #         histo = self.ftempl['fakesLowMt_fakeNormLowMtBin{}{}'.format(i,type)][:]
                #         dset = f.create_dataset(name='fakesLowMt_fakeNormLowMtBin{}{}'.format(i,type), shape=histo.ravel().shape, dtype='float64')
                #         dset[...] = histo.ravel()

                #         histo = self.ftempl['fakesHighMt_fakeNormHighMtBin{}{}'.format(i,type)][:]
                #         dset = f.create_dataset(name='fakesHighMt_fakeNormHighMtBin{}{}'.format(i,type), shape=histo.ravel().shape, dtype='float64')
                #         dset[...] = histo.ravel()

                #         histo = self.ftempl['fakesLowMt_fakeShapeBin{}{}'.format(i,type)][:]
                #         dset = f.create_dataset(name='fakesLowMt_fakeShapeBin{}{}'.format(i,type), shape=histo.ravel().shape, dtype='float64')
                #         dset[...] = histo.ravel()

                #         histo = self.ftempl['fakesHighMt_fakeShapeBin{}{}'.format(i,type)][:]
                #         dset = f.create_dataset(name='fakesHighMt_fakeShapeBin{}{}'.format(i,type), shape=histo.ravel().shape, dtype='float64')
                #         dset[...] = histo.ravel()
                #         pass

    def maskedChannels(self):
        dtype = 'float64'
        for chan in self.channels:
            with h5py.File('{}_xsec.hdf5'.format(chan), mode="w") as f:
                for proc in self.processes:
                    if "helXsecs" in proc: #give the correct xsec to unfold
                        coeff = self.helXsecs.index(proc.split('_')[0].replace('helXsecs',''))
                        iY = int(proc.split('_')[2])
                        iQt = int(proc.split('_')[4])
                        dset_masked = f.create_dataset(proc, [1], dtype=dtype)
                        dset_masked[...] = self.gen[chan][iY,iQt,coeff]
                    else:
                        dset_masked = f.create_dataset(proc, [1], dtype=dtype)
                        dset_masked[...] = 0.
                dset_masked = f.create_dataset("data_obs", [1], dtype=dtype)
                dset_masked[...] = 1.

    def setPreconditionVec(self):
        f=h5py.File('fitresults_asimov_onlysig.hdf5', 'r+')
        hessian = f['hess'][:]
        eig, U = np.linalg.eigh(hessian)
        M1 = np.matmul(np.diag(1./np.sqrt(eig)),U.T)
        # print(M1,np.linalg.inv(np.linalg.inv(M1)))
        self.preconditioner = M1
        self.invpreconditioner = np.linalg.inv(self.preconditioner)
        # print(np.matmul(self.preconditioner,np.linalg.inv(self.preconditioner)))
        # test = np.matmul(self.preconditioner,np.linalg.inv(self.preconditioner)) - np.identity(M1.shape[0])
        # print(np.max(np.abs(test)))
        # self.preconditioner = np.identity(len(self.signals))
        # self.invpreconditioner = np.identity(len(self.signals))

    def fillHelGroup(self):

        for i in range(len(self.yBinsC)):
            for j in range(len(self.qtBinsC)):

                s = 'y_{i}_qt_{j}'.format(i=i,j=j)
                self.helGroups[s] = []
                
                for hel in self.helXsecs:
                    if 'helXsecs'+hel+'_'+s in self.signals:

                        self.helGroups[s].append('helXsecs'+hel+'_'+s)
                                
                if self.helGroups[s] == []:
                    del self.helGroups[s]
    def fillHelMetaGroup(self):

        for i in range(len(self.yBinsC)):
            s = 'y_{i}'.format(i=i)
            self.helMetaGroups[s] = []
            for key in self.sumGroups:
                if s in key:
                    self.helMetaGroups[s].append(key)
            
            if self.helMetaGroups[s] == []:
                    del self.helMetaGroups[s]
        
        for j in range(len(self.qtBinsC)):
            s = 'qt_{j}'.format(j=j)
            self.helMetaGroups[s] = []
            for key in self.sumGroups:
                if 'qt' in key and key.split('_')[2]==str(j):
                    self.helMetaGroups[s].append(key)
        
            if self.helMetaGroups[s] == []:
                    del self.helMetaGroups[s]
        # print self.helMetaGroups
    def fillSumGroup(self):

        for i in range(len(self.yBinsC)):
            s = 'y_{i}'.format(i=i)
            for hel in self.helXsecs:
                for signal in self.signals:
                    if 'helXsecs'+hel+'_'+s in signal:
                        self.sumGroups['helXsecs'+hel+'_'+s] = []
                        for j in range(len(self.qtBinsC)):
                            if 'helXsecs'+hel+'_'+'y_{i}_qt_{j}'.format(i=i,j=j) in self.signals:
                                self.sumGroups['helXsecs'+hel+'_'+s].append('helXsecs'+hel+'_'+s+'_qt_{j}'.format(j=j))
        
        for j in range(len(self.qtBinsC)):
            s = 'qt_{j}'.format(j=j)
            for hel in self.helXsecs:
                for signal in self.signals:
                    if signal.split('_')[0] == 'helXsecs'+hel and signal.split('_')[4] == str(j):
                        self.sumGroups['helXsecs'+hel+'_'+s] = []
                        for i in range(len(self.yBinsC)):
                            if 'helXsecs'+hel+'_'+'y_{i}_qt_{j}'.format(i=i,j=j) in self.signals:
                            #print i, signal, 'helXsecs'+hel+'_'+'y_{i}_pt_{j}'.format(i=i,j=j)
                            #print 'append', 'helXsecs'+hel+'_y_{i}_'.format(i=i)+s, 'to', 'helXsecs'+hel+'_'+s
                                self.sumGroups['helXsecs'+hel+'_'+s].append('helXsecs'+hel+'_y_{i}_'.format(i=i)+s)
    def makeDatacard(self):

        self.DC = Datacard()

        ############## Setup the datacard (must be filled in) ###########################

        self.DC.bins =   self.channels
        self.DC.bins.extend([c+'_xsec' for c in self.channels]) # <type 'list'>
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
        for chan in self.channels:
            self.DC.exp[chan] = {}
            self.DC.exp[chan+'_xsec'] = {}
            for proc in self.processes:
                self.DC.exp[chan][proc] = -1.00
                self.DC.exp[chan+'_xsec'][proc] = -1.00
        self.DC.systs =  [] # <type 'list'>
        # aux = {} #each sys will have a separate aux dict
        # aux[self.channel] = {}
        # aux[self.channel+'_xsec'] = {}
        # for i in range(48*30):
        #     for proc in self.processes:
        #         aux[self.channel][proc] = 0.
        #         aux[self.channel+'_xsec'][proc] = 0.
        #     aux[self.channel]['fakesLowMt'] = 1.
        #     aux[self.channel]['fakesHighMt'] = 1.
        #     self.DC.systs.append(('fakeShapeBin{}'.format(i), False, 'shapeNoConstraint', [], aux))
        # aux = {} #each sys will have a separate aux dict
        # aux[self.channel] = {}
        # aux[self.channel+'_xsec'] = {}
        # for proc in self.processes:
        #         aux[self.channel][proc] = 0.
        #         aux[self.channel+'_xsec'][proc] = 0.
        # aux[self.channel]['fakesLowMt'] = 1.5
        # aux[self.channel]['fakesHighMt'] = 0
        # self.DC.systs.append(('fakesNormLowMt', False, 'lnNNoConstraint', [], aux))
        # aux = {} #each sys will have a separate aux dict
        # aux[self.channel] = {}
        # aux[self.channel+'_xsec'] = {}
        # for proc in self.processes:
        #         aux[self.channel][proc] = 0.
        #         aux[self.channel+'_xsec'][proc] = 0.
        # aux[self.channel]['fakesLowMt'] = 0.
        # aux[self.channel]['fakesHighMt'] = 1.5
        # self.DC.systs.append(('fakesNormHighMt', False, 'lnNNoConstraint', [], aux))
        # aux = {} #each sys will have a separate aux dict
        # aux[self.channel] = {}
        # aux[self.channel+'_xsec'] = {}
        # for i in range(48*30):
        #     for proc in self.processes:
        #         aux[self.channel][proc] = 0.
        #         aux[self.channel+'_xsec'][proc] = 0.
        #     aux[self.channel]['fakesLowMt'] = 1.
        #     aux[self.channel]['fakesHighMt'] = 0.
        #     self.DC.systs.append(('fakeNormLowMtBin{}'.format(i), False, 'shapeNoConstraint', [], aux))
        # aux = {} #each sys will have a separate aux dict
        # aux[self.channel] = {}
        # aux[self.channel+'_xsec'] = {}
        # for i in range(48*30):
        #     for proc in self.processes:
        #         aux[self.channel][proc] = 0.
        #         aux[self.channel+'_xsec'][proc] = 0.
        #     aux[self.channel]['fakesLowMt'] = 0.
        #     aux[self.channel]['fakesHighMt'] = 1.
        #     self.DC.systs.append(('fakeNormHighMtBin{}'.format(i), False, 'shapeNoConstraint', [], aux))
        # list of [{bin : {process : [input file, path to shape, path to shape for uncertainty]}}]
        # if self.doSyst:
        #     for syst in self.templSystematics: #loop over systematics
        #         for var in self.templSystematics[syst]["vars"]:
        #             aux = {} #each sys will have a separate aux dict
        #             aux[self.channel] = {}
        #             aux[self.channel+'_xsec'] = {}
        #             for proc in self.processes: 
        #                 if proc in self.templSystematics[syst]["procs"]:
        #                     aux[self.channel][proc] = self.templSystematics[syst]["weight"]
        #                     aux[self.channel+'_xsec'][proc] = 0.0
        #                 else:
        #                     if "Signal" in self.templSystematics[syst]["procs"] and "hel" in proc:
        #                         aux[self.channel][proc] = self.templSystematics[syst]["weight"]
        #                         if syst in ["alphaS", "LHEPdfWeight", "mass"] : #theo nuis. applied to signal
        #                             aux[self.channel+'_xsec'][proc] = self.templSystematics[syst]["weight"]
        #                         else :
        #                             aux[self.channel+'_xsec'][proc] = 0.0
        #                     else:
        #                         aux[self.channel][proc] = 0.0
        #                         aux[self.channel+'_xsec'][proc] = 0.0

        #             self.DC.systs.append((var, False, self.templSystematics[syst]["type"], [], aux))
        self.DC.groups = {}
        # self.DC.groups = {'mass': ['mass'],
        #                  'pdfs': set(['LHEPdfWeightHess{}'.format(i) for i in range(103)]),
        #                   #'jme': set(['jesTotal', 'unclustEn']),
        #                   # 'PrefireWeight':['PrefireWeight'],
        #                  }  # <type 'dict'>
        for chan in self.channels:
            self.DC.shapeMap = 	{chan: {'*': [chan+'.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']},\
            chan+'_xsec': {'*': [chan+'_xsec.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']}} # <type 'dict'>
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
        # self.DC.noiGroups = {'mass':['mass']}
        # self.DC.noiGroups = {}

        self.DC.preconditioner  = self.preconditioner 
        self.DC.invpreconditioner  = self.invpreconditioner 
        
        filehandler = open('{}.pkl'.format("WPlus"), 'w')
        pickle.dump(self.DC, filehandler)
