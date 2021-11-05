import os
import sys
import ROOT
import argparse
import copy
import time
from datetime import datetime

FWKBASE=os.getenv('FWK_BASE')
print('Adding path: {}/RDFprocessor/framework'.format(FWKBASE))
sys.path.append('{}/RDFprocessor/framework'.format(FWKBASE))
from RDFtree import RDFtree

sys.path.append('{}/Common/data'.format(FWKBASE))
from samples_2016_ulV2 import samplespreVFP, samplespostVFP
from genSumW import sumwDictpreVFP, sumwDictpostVFP
#from genSumWClipped import sumwClippedDictpreVFP, sumwClippedDictpostVFP


ROOT.gSystem.Load('{}/nanotools/bin/libNanoTools.so'.format(FWKBASE))
sys.path.append('../nanotools')
from nanoSequence import nanoSequence

ROOT.gSystem.Load('{}/nanotools/bin/libNanoTools.so'.format(FWKBASE))
sys.path.append('{}/nanotools'.format(FWKBASE))
from nanoSequence import nanoSequence

sys.path.append('{}/templateMaker/'.format(FWKBASE))
from dySequence import dySelectionSequence

ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")
eras = ["preVFP","postVFP"]

def RDFprocess(fvec, outputDir, sample, xsec, systType, sumw, era, pretendJob):
    print("processing ", sample)
    p = RDFtree(outputDir = outputDir, inputFile = fvec, outputFile="{}.root".format(sample), pretend=pretendJob)
    postnano, endNode=nanoSequence(p, systType, sample, xsec, sumw, era)
    print("Post nano node name: ", endNode)
    #return postnano
    resultNode=dySelectionSequence(postnano, xsec, systType, sumw, endNode, era)
    return resultNode


def main():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-p', '--pretend',type=bool, default=False, help="run over a small number of event")
    parser.add_argument('-r', '--report',type=bool, default=False, help="Prints the cut flow report for all named filters")
    parser.add_argument('-o', '--outputDir',type=str, default='outputDY', help="output dir name")
    parser.add_argument('-i', '--inputDir',type=str, default='/scratchnvme/wmass/NANOJEC/', help="input dir name")    

    args = parser.parse_args()
    pretendJob = args.pretend
    inDir = args.inputDir
    RDFtrees = {}
    for era in eras:
        outputDir = args.outputDir + '_' + era
        ##Add era to input dir
        fullinDir=inDir+era
        if pretendJob:
            print("Running a test job over a few events")
        else:
            print("Running on full dataset")
        ROOT.ROOT.EnableImplicitMT(48)
        RDFtrees[era] = {}
        samples = samplespreVFP
        sumwClippedDict=sumwDictpreVFP
        if era == 'postVFP': 
            samples = samplespostVFP
            sumwClippedDict=sumwDictpostVFP

        for sample in samples:
            #print('analysing sample: %s'%sample)
            if 'WPlus' in sample or 'WMinus' in sample: continue
            checkS= 'DYJetsToMuMu' in sample or 'data' in sample
            #print('CheckSample={}'.format(checkS))
            #if not checkS : continue
            direc = samples[sample]['dir']
            xsec = samples[sample]['xsec']
            fvec=ROOT.vector('string')()
            for d in direc:
                targetDir='{}/{}/'.format(fullinDir, d)
                for f in os.listdir(targetDir):#check the directory
                    if not f.endswith('.root'): continue
                    inputFile=targetDir+f
                    #print(f)
                    fvec.push_back(inputFile)
            if fvec.empty():
                print("No files found for directory:", samples[sample], " SKIPPING processing")
                continue
            print(fvec)         
            systType = samples[sample]['nsyst']
            sumw=1.
            if not 'data' in sample:
                sumw=sumwClippedDict[sample]        
            RDFtrees[era][sample] = RDFprocess(fvec, outputDir, sample, xsec, systType, sumw, era, pretendJob)
    #sys.exit(0)
    #now trigger all the event loops at the same time:
    objList = []
    cutFlowreportDict = {}
    for era in eras:
        print(era, "merging objects")
        samples = samplespreVFP
        sumwClippedDict=sumwDictpreVFP
        if era == 'postVFP': 
            # samples = wsignalNLO_postVFP
            samples = samplespostVFP
            sumwClippedDict= sumwDictpostVFP
        for sample in samples:
            for sample in samples:
                if 'WPlus' in sample or 'WMinus' in sample: continue
                checkS= 'DYJetsToMuMu' in sample or 'data' in sample
                #if not checkS : continue
                print(sample)
                RDFtreeDict = RDFtrees[era][sample].getObjects()
                if args.report: cutFlowreportDict[sample] = RDFtrees[era][sample].getCutFlowReport()
                for node in RDFtreeDict:
                    objList.extend(RDFtreeDict[node])
    print("end merging objects")
    #magic happens here
    start = time.time()
    ROOT.RDF.RunGraphs(objList)
    #now write the histograms:
    for era in eras:
        for sample in samples:    
            if 'WPlus' in sample or 'WMinus' in sample: continue
            checkS= 'DYJetsToMuMu' in sample or 'data' in sample
            #if not checkS : continue
            print(sample)
            #RDFtrees[sample].getOutput()
            RDFtrees[era][sample].gethdf5Output()
            if args.report: cutFlowreportDict[sample].Print()
            #RDFtrees[sample].saveGraph()

    print('all samples processed in {} s'.format(time.time()-start))

#At some point this should be what we want
#p.Histogram(columns = ["dimuonMass", "Mu1_eta", "Mu1_pt", "Mu2_eta", "Mu2_pt", "dimuonPt", "dimuonY", "MET_T1_pt"], types = ['float']*8,node='defs',histoname=ROOT.string('data_muons'),bins = [zmassBins, etaBins, ptBins, etaBins, qtBins, etaBins, metBins], variations = [])
if __name__ == "__main__":
    main()
