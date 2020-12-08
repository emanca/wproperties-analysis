import sys
import ROOT
import os
import copy
import json
import argparse

import bkg_utils

class bkg_prepareHistos:
    def __init__(self, outDir, inputDir,extrap=True) :
        self.outDir = outDir
        self.inputDir = inputDir
        self.extrap = extrap
        
        self.systDict = copy.deepcopy(bkg_utils.bkg_systematics)
        self.sNameNom = ''
        self.sKindNom = 'Nominal'
        self.systDict.update({self.sKindNom:[self.sNameNom]})
    
    
    def fileListBuilder(self) :
        self.fileList = []
        with open('../analysisOnData/data/samples_2016.json') as f:
            samples = json.load(f)
        for sample in samples :
            if samples[sample]['datatype']=='DATA':  continue
            if sample.startswith('WJetsToLNu') : continue 
            self.fileList.append(sample)
        self.fileList.append('WToMu')#not in the json list
        self.fileList.append('WToTau')#not in the json list
        print(self.fileList)
        
        
    def prepare(self) :
        varName = 'templates'
        suffString = '_plots'
        regList = ["Signal", "Signal_aiso","Sideband","Sideband_aiso"]
        
        if self.extrap :
            for lcut, lbin in bkg_utils.looseCutDict.items() :
                regList.append('Sideband'+lcut)
                regList.append('Sideband_aiso'+lcut)
                
        for f in self.fileList : 
            inFile =  ROOT.TFile.Open(self.inputDir+f+suffString+".root")
            outFile =  ROOT.TFile(self.outDir+f+".root", "recreate")
            for r in regList :
                regExtrapFlag = self.isExtrapReg(r)
                for sKind, sList in self.systDict.items():  
                    if sKind!='Nominal' and regExtrapFlag : #extrap region only for nominal
                        continue
                    if not outFile.GetDirectory('templates_'+r+'/'+sKind):                           
                        outFile.mkdir('templates_'+r+'/'+sKind)
                    for sName in sList : 
                        inFile.cd()
                        if sKind!=self.sKindNom :
                            systName = '_'+sName
                        else :
                            systName = sName
                        #print "Trying to read:", 'templates_'+r+'/'+sKind+'/'+varName+systName, ' in ', f
                        if ROOT.gDirectory.Get('templates_'+r+'/'+sKind+'/'+varName+systName)==None : #this syst is not present
                            print("no syst in:", f, r, sKind, systName)
                            h = inFile.Get('templates_'+r+'/'+self.sKindNom+'/'+varName+self.sNameNom)
                            h.SetName(varName+systName)
                        else :
                            h = inFile.Get('templates_'+r+'/'+sKind+'/'+varName+systName)
                        outFile.cd('templates_'+r+'/'+sKind)
                        h.Write()
                        
    def h_add(self) :
        cmdList = []
        if not os.path.isdir(self.outDir+'/hadded'): os.system('mkdir '+self.outDir+'/hadded')
        
        cmdList.append('hadd ./'+self.outDir+'hadded/WToMuNu.root ./'+self.outDir+'*.root')
        cmdList.append('cp  ./'+self.inputDir+'SingleMuonData_plots.root ./'+self.outDir+'hadded/Data.root')
                
        for i in cmdList :
            os.system(i)

    def isExtrapReg(self, reg) :
        FlagOut = False
        for lcut, lbin in bkg_utils.looseCutDict.items() :
                if reg == 'Sideband_aiso'+lcut:
                    FlagOut = True
                if reg == 'Sideband'+lcut:
                    FlagOut = True
        return  FlagOut 


parser = argparse.ArgumentParser("")
parser.add_argument('-inputDir', '--inputDir',type=str, default='./data/', help="input dir name")
parser.add_argument('-outputDir', '--outputDir',type=str, default='./bkg/', help="output dir name")

args = parser.parse_args()
inputDir = args.inputDir
outputDir = args.outputDir

preparator = bkg_prepareHistos(outDir=outputDir+'/', inputDir=inputDir+'/',extrap=True)
preparator.fileListBuilder()
preparator.prepare()
preparator.h_add()

