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
        
        self.flatDict = {
            "lumi": {
                "procs": ["WToMu", "ST_","TTJets_","WW_","WZ_","ZZ_","WToTau","DYJets"],
                "weight" : 1.025 
            },
            "topXSec":{
                "procs": ["ST_","TTJets_"],
                "weight" : 1.060 
            },
            "dibosonXSec":{
                "procs": ["WW_","WZ_","ZZ_"],
                "weight" : 1.160 
            },
            "tauXSec":{
                "procs": ["WToTau"],
                "weight" : 1.040 
            },
            "lepVeto":{
                "procs": ["DYJets"],
                "weight" : 1.020 
            }
        }
        if not os.path.isdir(self.outDir+'/hadded'): os.system('mkdir '+self.outDir+'/hadded')
    
    
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
                    if sKind in self.flatDict :
                        continue
                    if sKind!='Nominal' and regExtrapFlag : #extrap region only for nominal
                        continue
                    # if 'Nom_WQT' in sKind: continue # used in LHEScaleWeight_WQt case, not useful as standalone. Commented because has been removed from the bkg_syst_dict
                    if not outFile.GetDirectory('templates_'+r+'/'+sKind):                           
                        outFile.mkdir('templates_'+r+'/'+sKind)
                    for sName in sList : 
                        if '_WQT' in sKind :
                            #  sName = sName+'_'+sKind.replace('LHEScaleWeight_','')
                             sName = sName.replace(sKind.replace('LHEScaleWeight',''),'') #replace for instance: LHEScaleWeight_muR0p5_muF0p5_WQTlow-> LHEScaleWeight_muR0p5_muF0p5
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
                            
                            #Wpt shape systematics block
                            if sKind=='LHEScaleWeight_WQTlow' :
                                h.Add(inFile.Get('templates_'+r+'/Nom_WQTmid/'+varName+self.sNameNom))
                                h.Add(inFile.Get('templates_'+r+'/Nom_WQThigh/'+varName+self.sNameNom))
                            if sKind=='LHEScaleWeight_WQTmid' :
                                h.Add(inFile.Get('templates_'+r+'/Nom_WQTlow/'+varName+self.sNameNom))
                                h.Add(inFile.Get('templates_'+r+'/Nom_WQThigh/'+varName+self.sNameNom))
                            if sKind=='LHEScaleWeight_WQThigh' :
                                h.Add(inFile.Get('templates_'+r+'/Nom_WQTlow/'+varName+self.sNameNom))
                                h.Add(inFile.Get('templates_'+r+'/Nom_WQTmid/'+varName+self.sNameNom))
                        if '_WQT' in sKind : 
                            h.SetName(h.GetName()+sKind.replace('LHEScaleWeight',''))   
                        outFile.cd('templates_'+r+'/'+sKind)
                        # h.RebinY(4)#rebinned pt to have 2 GeV binning
                        h.Write()
                        
                for sKind, sList in self.systDict.items():  #flat syst building
                    if sKind not in self.flatDict : continue 
                    if regExtrapFlag : continue
                    if not outFile.GetDirectory('templates_'+r+'/'+sKind):                           
                        outFile.mkdir('templates_'+r+'/'+sKind)
                    for sName in sList : 
                        h = inFile.Get('templates_'+r+'/'+self.sKindNom+'/'+varName+self.sNameNom)
                        h = h.Clone(varName+'_'+sName)
                        procFlag = False 
                        for proc in self.flatDict[sKind]['procs'] :
                            if proc in f : #proc are part of the full filename of f
                                procFlag = True 
                        if procFlag :
                            if 'Up' in sName :
                                h.Scale(self.flatDict[sKind]['weight'])
                            if 'Down' in sName :
                                h.Scale(1/self.flatDict[sKind]['weight'])
                        outFile.cd('templates_'+r+'/'+sKind)
                        # h.RebinY(4)#rebinned pt to have 2 GeV binning (not decomment this line)
                        h.Write()
        
        #data rebin
        inFile =  ROOT.TFile.Open(self.inputDir+'SingleMuonData_plots.root')
        outFile =  ROOT.TFile(self.outDir+"/hadded/Data.root", "recreate")
        print("Data (rebin)")
        for r in regList :
            sKind = self.sKindNom
            sName = self.sNameNom
            regExtrapFlag = self.isExtrapReg(r)
            if not outFile.GetDirectory('templates_'+r+'/'+sKind): 
                outFile.mkdir('templates_'+r+'/'+sKind)
            inFile.cd()
            h = inFile.Get('templates_'+r+'/'+sKind+'/'+varName)  
            # h.RebinY(4)#rebinned pt to have 2 GeV binning   
            outFile.cd('templates_'+r+'/'+sKind)   
            h.Write()         
                
                        
    def h_add(self) :
        cmdList = []
        
        cmdList.append('hadd -f ./'+self.outDir+'hadded/WToMuNu.root ./'+self.outDir+'*.root')
        # cmdList.append('cp  ./'+self.inputDir+'SingleMuonData_plots.root ./'+self.outDir+'hadded/Data.root')
                
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

