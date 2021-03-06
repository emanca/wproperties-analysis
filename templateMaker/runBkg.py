import os
import sys
import ROOT
import argparse
import copy
import time
from datetime import datetime
sys.path.append('../RDFprocessor/framework')
sys.path.append('../Common/data')
from RDFtree import RDFtree
from samples_2016 import samples
sys.path.append('python/')
from getLumiWeight import getLumiWeight
from binning import ptBins, etaBins, mTBins, isoBins, chargeBins
from externals import filePt, fileY, fileSF

ROOT.gSystem.Load('bin/libAnalysisOnData.so')
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")

def RDFprocess(fvec, outputDir, sample, xsec, systType, pretendJob):
    print("processing ", sample)
    p = RDFtree(outputDir = outputDir, inputFile = fvec, outputFile="{}.root".format(sample), pretend=pretendJob)

    if systType == 0: #this is data
        p.branch(nodeToStart='input', nodeToEnd='defs', modules=[ROOT.recoDefinitions(False, False)])
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="HLT_SingleMu24", filtername="{:20s}".format("Pass HLT"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(Vtype==0 || Vtype==1)", filtername="{:20s}".format("Vtype selection"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="HLT_SingleMu24", filtername="{:20s}".format("Pass HLT"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="MET_filters==1", filtername="{:20s}".format("Pass MET filter"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="nVetoElectrons==0", filtername="{:20s}".format("Electron veto"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Idx_mu1>-1", filtername="{:20s}".format("Atleast 1 mu"))

        p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso"], types = ['float']*5,node='defs',histoname=ROOT.string('data_obs'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        return p
    elif systType < 2: #this is MC with no PDF variations
        p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.recoDefinitions(True, False),ROOT.recoWeightDefinitions(fileSF),getLumiWeight(xsec=xsec, inputFile=fvec, genEvsbranch = "genEventSumw_")])
    else:
        if 'DY' in sample: #reweight full Z kinematics
            p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.recoDefinitions(True, False),ROOT.recoWeightDefinitions(fileSF),getLumiWeight(xsec=xsec, inputFile=fvec, genEvsbranch = "genEventSumw_"),ROOT.Replica2Hessian(),ROOT.reweightFromZ(filePt,fileY,False,False)])
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="HLT_SingleMu24", filtername="{:20s}".format("Pass HLT"))
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(Vtype==0 || Vtype==1)", filtername="{:20s}".format("Vtype selection"))
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="HLT_SingleMu24", filtername="{:20s}".format("Pass HLT"))
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="MET_filters==1", filtername="{:20s}".format("Pass MET filter"))
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="nVetoElectrons==0", filtername="{:20s}".format("Electron veto"))
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Idx_mu1>-1", filtername="{:20s}".format("Aleast 1 mu"))
            p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MTVars","Mu1_relIso", "lumiweight", "weightPt", "weightY", "PrefireWeight", "puWeight", "WHSFVars"], types = ['float']*11,node='defs',histoname=ROOT.string('ewk'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        else:
            p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.recoDefinitions(True, False),ROOT.recoWeightDefinitions(fileSF),getLumiWeight(xsec=xsec, inputFile=fvec, genEvsbranch = "genEventSumw_"),ROOT.Replica2Hessian(),ROOT.reweightFromZ(filePt,fileY,True,False)])
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="HLT_SingleMu24", filtername="{:20s}".format("Pass HLT"))
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(Vtype==0 || Vtype==1)", filtername="{:20s}".format("Vtype selection"))
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="HLT_SingleMu24", filtername="{:20s}".format("Pass HLT"))
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="MET_filters==1", filtername="{:20s}".format("Pass MET filter"))
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="nVetoElectrons==0", filtername="{:20s}".format("Electron veto"))
            p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Idx_mu1>-1", filtername="{:20s}".format("Aleast 1 mu"))
            p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MTVars","Mu1_relIso", "lumiweight", "weightPt", "PrefireWeight", "puWeight", "WHSFVars"], types = ['float']*10,node='defs',histoname=ROOT.string('ewk'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
    return p

def main():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-p', '--pretend',type=bool, default=False, help="run over a small number of event")
    parser.add_argument('-r', '--report',type=bool, default=False, help="Prints the cut flow report for all named filters")
    parser.add_argument('-o', '--outputDir',type=str, default='output', help="output dir name")
    parser.add_argument('-i', '--inputDir',type=str, default='/scratchnvme/wmass/NanoAOD2016-V2/', help="input dir name")    

    args = parser.parse_args()
    pretendJob = args.pretend
    now = datetime.now()
    dt_string = now.strftime("_%d_%m_%Y_%H_%M_%S")
    outputDir = args.outputDir + dt_string
    inDir = args.inputDir
    if pretendJob:
        print("Running a test job over a few events")
    else:
        print("Running on full dataset")
    ROOT.ROOT.EnableImplicitMT(128)
    RDFtrees = {}
    
    for sample in samples:
        #print('analysing sample: %s'%sample)

        direc = samples[sample]['dir']
        xsec = samples[sample]['xsec']
        fvec=ROOT.vector('string')()
        for d in direc:
            ##check if file exists or not
            inputFile = '{}/{}/tree.root'.format(inDir, d)
            isFile = os.path.isfile(inputFile)
            if not isFile:
                print(inputFile, " does not exist")
                continue
            fvec.push_back(inputFile)
            #f = ROOT.TFile(inputFile)
            #t = f.Get('Events')
            #it = t.GetClusterIterator(0)
            #counter = 0
            #n_entries = t.GetEntries()
            #while it.Next() != n_entries:
            #    counter += 1
            #print('number of clusters:',counter, sample)
        if fvec.empty():
            print("No files found for directory:", samples[sample], " SKIPPING processing")
            continue
        #print(fvec) 
        systType = samples[sample]['nsyst']
        RDFtrees[sample] = RDFprocess(fvec, outputDir, sample, xsec, systType, pretendJob)

    #now trigger all the event loops at the same time:
    objList = []
    for sample in samples:
        RDFtreeDict = RDFtrees[sample].getObjects()
        for node in RDFtreeDict:
            objList.extend(RDFtreeDict[node])
    #magic happens here
    start = time.time()
    ROOT.RDF.RunGraphs(objList)
    #now write the histograms:
    
    for sample in samples:
        print(sample)
        #RDFtrees[sample].getOutput()
        RDFtrees[sample].gethdf5Output()
        if args.report: RDFtrees[sample].getCutFlowReport()

    print('all samples processed in {} s'.format(time.time()-start))
if __name__ == "__main__":
    main()
