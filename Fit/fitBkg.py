import ROOT
import pickle
from termcolor import colored
import math
from HiggsAnalysis.CombinedLimit.DatacardParser import *
from collections import OrderedDict
import copy
import numpy as np

class fitUtils:
    def __init__(self, channel ="bkgWplus"):
        
        self.processes = ['ewk','fakesLowMt','fakesHighMt']
        self.signals = ['ewk']

        #combine utils
        self.channel = channel

    def makeDatacard(self):

        self.DC = Datacard()

        ############## Setup the datacard (must be filled in) ###########################

        self.DC.bins =   [self.channel] # <type 'list'>
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
        for proc in self.processes:
            self.DC.exp[self.channel][proc] = -1.00
        self.DC.systs =  [] # <type 'list'>
        # list of [{bin : {process : [input file, path to shape, path to shape for uncertainty]}}]
        aux = {} #each sys will have a separate aux dict
        aux[self.channel] = {}
        for i in range(48*30):
            aux[self.channel]['fakesLowMt'] = 1.
            aux[self.channel]['fakesHighMt'] = 1.
            aux[self.channel]['ewk'] = 0.
            self.DC.systs.append(('fakeShapeBin{}'.format(i), False, 'shapeNoConstraint', [], aux))
        aux = {} #each sys will have a separate aux dict
        aux[self.channel] = {}
        for i in range(48*30):
            aux[self.channel]['fakesLowMt'] = 1.
            aux[self.channel]['fakesHighMt'] = 0.
            aux[self.channel]['ewk'] = 0.
            self.DC.systs.append(('fakeNormLowMtBin{}'.format(i), False, 'shapeNoConstraint', [], aux))
        aux = {} #each sys will have a separate aux dict
        aux[self.channel] = {}
        for i in range(48*30):
            aux[self.channel]['fakesLowMt'] = 0.
            aux[self.channel]['fakesHighMt'] = 1.
            aux[self.channel]['ewk'] = 0.
            self.DC.systs.append(('fakeNormHighMtBin{}'.format(i), False, 'shapeNoConstraint', [], aux))
        # aux = {} #each sys will have a separate aux dict
        # aux[self.channel] = {}
        # aux[self.channel]['fakesLowMt'] = 1.5
        # aux[self.channel]['fakesHighMt'] = 0
        # aux[self.channel]['ewk'] = 0.
        # self.DC.systs.append(('fakesNormLowMt', False, 'lnNNoConstraint', [], aux))
        # aux = {} #each sys will have a separate aux dict
        # aux[self.channel] = {}
        # aux[self.channel]['fakesLowMt'] = 0.
        # aux[self.channel]['fakesHighMt'] = 1.5
        # aux[self.channel]['ewk'] = 0.
        # self.DC.systs.append(('fakesNormHighMt', False, 'lnNNoConstraint', [], aux))
        # self.DC.groups = {'fakes': set(['fakeShapeBin{}'.format(i) for i in range(48*30)])}
        self.DC.shapeMap = 	{self.channel: {'*': [self.channel+'.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']}}
        self.DC.hasShapes =  True # <type 'bool'>
        self.DC.flatParamNuisances =  {} # <type 'dict'>
        self.DC.rateParams =  {} # <type 'dict'>
        self.DC.extArgs =    {} # <type 'dict'>
        self.DC.rateParamsOrder  =  set([]) # <type 'set'>
        self.DC.frozenNuisances  =  set([]) # <type 'set'>
        self.DC.systematicsShapeMap =  {} # <type 'dict'>
        self.DC.nuisanceEditLines    =  [] # <type 'list'>
        self.DC.discretes    =  [] # <type 'list'>

        filehandler = open('{}.pkl'.format(self.channel), 'w')
        pickle.dump(self.DC, filehandler)

fit = fitUtils()
fit.makeDatacard()

