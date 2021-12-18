import os
import sys
import ROOT
import json
import argparse
import copy
from RDFtree import RDFtree
from multiprocessing import Process

sys.path.append('python/')
sys.path.append('data/')
from systematics import systematics
from selections import selections, selectionVars, selections_bkg

from getLumiWeight import getLumiWeight

ROOT.gSystem.Load('bin/libAnalysisOnData.so')
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;");

def RDFprocessMC(fvec, outputDir, sample, xsec, fileSF, ncores, systType, pretendJob=True,SBana=False):
    ROOT.ROOT.EnableImplicitMT(ncores)
    print("running with {} cores".format(ncores))
    p = RDFtree(outputDir = outputDir, inputFile = fvec, outputFile="{}_plots.root".format(sample), pretend=pretendJob)
    p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.baseDefinitions(True, False),ROOT.weightDefinitions(fileSF),getLumiWeight(xsec=xsec, inputFile=fvec)])
    #selections bkg also includes Signal
    QCDreg = {
        "iso":         "Vtype==0 && HLT_SingleMu24 && MET_filters==1 && nVetoElectrons==0 && 1",
        "aiso":    "Vtype==1 && HLT_SingleMu24 && MET_filters==1 && nVetoElectrons==0 && 1",
    }
    for region,cut in QCDreg.items():
        print("running in region {}".format(region))
                    
        if 'aiso' in region:
            weight = 'float(puWeight*PrefireWeight*lumiweight)'
        else:
            weight = 'float(puWeight*PrefireWeight*lumiweight*WHSF)'
        
        nom = ROOT.vector('string')()
        nom.push_back("")
        #last argument refers to histo category - 0 = Nominal, 1 = Pt scale , 2 = MET scale
        print("branching nominal")

        p.branch(nodeToStart = 'defs', nodeToEnd = 'QCDhisto_{}/Nominal'.format(region), modules = [ROOT.QCDhistos(cut, weight, nom,"Nom",0)])    
    
    p.getOutput()
    p.saveGraph()


def main():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-p', '--pretend',type=int, default=False, help="run over a small number of event")
    parser.add_argument('-i', '--inputDir',type=str, default='/scratchnvme/wmass/NanoAOD2016-V2/', help="input dir name")    
    parser.add_argument('-o', '--outputDir',type=str, default='./output/', help="output dir name")
    parser.add_argument('-c', '--ncores',type=int, default=64, help="number of cores used")
    parser.add_argument('-sb', '--SBana',type=int, default=False, help="run also on the sideband (clousure test)")
    args = parser.parse_args()
    pretendJob = args.pretend
    inDir = args.inputDir
    outputDir = args.outputDir
    ncores = args.ncores

    SBana = args.SBana
    if pretendJob:
        print("Running a test job over a few events")
    else:
        print("Running on full dataset")
    samples={}
    with open('data/samples_2016_QCD.json') as f:
        samples = json.load(f)
    for sample in samples:
        if not samples[sample]['datatype']=='MC': continue
        #if not samples[sample]['multiprocessing']: continue
        #WJets to be run by a separate config
        if 'WJets' in sample : continue
        print(sample)
        direc = samples[sample]['dir']
        xsec = samples[sample]['xsec']
        cores = 0 if samples[sample]['multiprocessing'] else ncores 		
        fvec=ROOT.vector('string')()
        for dirname,fname in direc.items():
            ##check if file exists or not
            inputFile = '{}/{}/tree.root'.format(inDir, dirname)
            isFile = os.path.isfile(inputFile)  
            if not isFile:
                print(inputFile, " does not exist")
                continue
            fvec.push_back(inputFile)

        if fvec.empty():
            print("No files found for directory:", samples[sample], " SKIPPING processing")
            continue
        print(fvec) 
        fileSF = ROOT.TFile.Open("data/ScaleFactors_OnTheFly.root")
        systType = samples[sample]['systematics']
        RDFprocessMC(fvec, outputDir, sample, xsec, fileSF, cores, systType, pretendJob,SBana)

if __name__ == "__main__":
    main()
