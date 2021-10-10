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
from genSumW import sumwDictpreVFP, sumwDictpostVFP
# from samples_2016_ulCentral import samplespreVFP, samplespostVFP
from samples_2016_ulV2 import samplespreVFP, samplespostVFP, wsignalNLO_preVFP, wsignalNLO_postVFP
#from samples_2016_ul import samplespreVFP, samplespostVFP

ROOT.gSystem.Load('{}/nanotools/bin/libNanoTools.so'.format(FWKBASE))
sys.path.append('{}/nanotools'.format(FWKBASE))
from nanoSequence import nanoSequence

sys.path.append('{}/templateMaker/'.format(FWKBASE))
from wSequence import wSelectionSequence, wSelectionHelWeightsSequence, wSelectionDifferentialSequence

ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")

eras = ["preVFP","postVFP"]

def RDFprocess(fvec, outputDir, sample, xsec, systType, sumw, era, pretendJob, helWeights=False):
    print("processing ", sample)
    p = RDFtree(outputDir = outputDir, inputFile = fvec, outputFile="{}.root".format(sample if not helWeights else sample+'_helweights'), pretend=pretendJob)
    postnano, endNode=nanoSequence(p, systType, sample, xsec, sumw, era)
    print("Post nano node name: ", endNode)
    #return postnano
    if not helWeights: 
        resultNode = wSelectionSequence(postnano, systType, endNode, era)
        # if sample =="WJetsToLNu_0J" or sample=='WJetsToLNu_1J' or sample=='WJetsToLNu_2J':
        if sample =="WPlusJetsToMuNu" or sample=='WMinusJetsToMuNu':
            resultNode = wSelectionDifferentialSequence(resultNode,era,sample)
    else: resultNode = wSelectionHelWeightsSequence(postnano, endNode,era)
    
    return resultNode


def main():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-p', '--pretend',type=bool, default=False, help="run over a small number of event")
    parser.add_argument('-r', '--report',type=bool, default=False, help="Prints the cut flow report for all named filters")
    parser.add_argument('-o', '--outputDir',type=str, default='outputW', help="output dir name")
    parser.add_argument('-i', '--inputDir',type=str, default='/scratchnvme/wmass/NANOJEC/', help="input dir name")    
    parser.add_argument('-helWeights', '--helWeights',type=bool, default=False, help="derive helicity weights for reweighting")    

    RDFtrees = {}
    args = parser.parse_args()
    for era in eras:
        pretendJob = args.pretend
        inDir = args.inputDir
        outputDir = args.outputDir+"_"+era
        helWeights = args.helWeights
        ##Add era to input dir
        inDir+=era
        if pretendJob:
            print("Running a test job over a few events")
        else:
            print("Running on full dataset")
        ROOT.ROOT.EnableImplicitMT(48)
        RDFtrees[era] = {}

        samples = samplespreVFP
        # samples = wsignalNLO_preVFP
        sumwClippedDict=sumwDictpreVFP
        if era == 'postVFP': 
            # samples = wsignalNLO_postVFP
            samples = samplespostVFP
            sumwClippedDict= sumwDictpostVFP

        for sample in samples:
            if helWeights:
                if not 'WPlusJetsToMuNu' in sample or 'WMinusJetsToMuNu' in sample: continue
                # if not 'WJetsToLNu_0J' in sample and not 'WJetsToLNu_1J' in sample and not 'WJetsToLNu_2J' in sample: continue
            systType = samples[sample]['nsyst']
            # if systType==0: continue #only process samples with pdfs
            if  'WPlusJetsToMuNu' not in sample and 'data' not in sample: continue
            # if not 'WJetsToLNu_0J' in sample and not 'WJetsToLNu_1J' in sample and not 'WJetsToLNu_2J' in sample: continue
            direc = samples[sample]['dir']
            xsec = samples[sample]['xsec']
            fvec=ROOT.vector('string')()
            for d in direc:
                targetDir='{}/{}/'.format(inDir, d)
                # mergedDir=targetDir+'merged/'
                # if os.path.isdir(mergedDir):
                #     targetDir=mergedDir
                #for now  let's run on unmerged file.
                for f in os.listdir(targetDir):#check the directory
                    if not f.endswith('.root'): continue
                    inputFile=targetDir+f
                    #print(f)
                    fvec.push_back(inputFile)
            if fvec.empty():
                print("No files found for directory:", samples[sample], " SKIPPING processing")
                continue
            sumw=1.
            if not 'data' in sample:
                sumw=sumwClippedDict[sample]
            RDFtrees[era][sample] = RDFprocess(fvec, outputDir, sample, xsec, systType, sumw, era, pretendJob, helWeights)

    #now trigger all the event loops at the same time:
    objList = []
    cutFlowreportDict = {}
    for era in eras:
        print(era, "merging objects")
        sumwClippedDict=sumwDictpreVFP
        
        samples = samplespreVFP
        # samples = wsignalNLO_preVFP
        sumwClippedDict=sumwDictpreVFP
        if era == 'postVFP': 
            # samples = wsignalNLO_postVFP
            samples = samplespostVFP
            sumwClippedDict= sumwDictpostVFP
        for sample in samples:
            if helWeights:
                if not 'WPlusJetsToMuNu' in sample or 'WMinusJetsToMuNu' in sample: continue
                # if not 'WJetsToLNu_0J' in sample and not 'WJetsToLNu_1J' in sample and not 'WJetsToLNu_2J' in sample: continue
            systType = samples[sample]['nsyst']
            # if systType==0: continue #only process samples with pdfs
            if 'WPlusJetsToMuNu' not in sample and 'data' not in sample: continue
            # if not 'WJetsToLNu_0J' in sample and not 'WJetsToLNu_1J' in sample and not 'WJetsToLNu_2J' in sample: continue
            RDFtreeDict = RDFtrees[era][sample].getObjects()
            if args.report: cutFlowreportDict[sample] = RDFtrees[era][sample].getCutFlowReport('defs')
            for node in RDFtreeDict:
                objList.extend(RDFtreeDict[node])
    print("end merging objects")
    #magic happens here
    start = time.time()
    print("triggering event loop...")
    ROOT.RDF.RunGraphs(objList)
    #now write the histograms:
    for era in eras:
        for sample in samples:
            if helWeights:
                if not 'WPlusJetsToMuNu' in sample or 'WMinusJetsToMuNu' in sample: continue
                # if not 'WJetsToLNu_0J' in sample and not 'WJetsToLNu_1J' in sample and not 'WJetsToLNu_2J' in sample: continue
            systType = samples[sample]['nsyst']
            # if systType==0: continue #only process samples with pdfs
            if 'WPlusJetsToMuNu' not in sample and 'data' not in sample: continue
            # if not 'WJetsToLNu_0J' in sample and not 'WJetsToLNu_1J' in sample and not 'WJetsToLNu_2J' in sample: continue
            RDFtrees[era][sample].gethdf5Output()
            if args.report: cutFlowreportDict[sample].Print()
            # RDFtrees[era][sample].saveGraph()

    print('all samples processed in {} s'.format(time.time()-start))

if __name__ == "__main__":
    main()

