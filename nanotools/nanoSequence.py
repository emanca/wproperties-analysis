import os
import sys
import ROOT
import argparse
import copy
import time
from datetime import datetime

FWKBASE=os.getenv('FWK_BASE')
sys.path.append('{}/RDFprocessor/framework'.format(FWKBASE))
from RDFtree import RDFtree

sys.path.append('{}/Common/data/'.format(FWKBASE))
from externals import pufile_mc_UL2016, pufile_data_UL2016_allData, pufile_data_UL2016_preVFP, pufile_data_UL2016_postVFP
from externals import datajson,filemuPrefire
from dataluminosity import lumi_preVFP,lumi_postVFP,lumi_total2016
from genSumW import sumwClipSamples

JMEBASE="/scratchnvme/wmass/CMSJMECalculators"
sys.path.append("{}/include/".format(JMEBASE))
sys.path.append("{}/interface/".format(JMEBASE))
sys.path.append("{}/lib/".format(JMEBASE))
ROOT.gSystem.Load('{}/lib/libCMSJMECalculators.so'.format(JMEBASE))
ROOT.gSystem.Load('{}/lib/libCMSJMECalculatorsDict.so'.format(JMEBASE))

ROOT.gSystem.Load('{}/nanotools/bin/libNanoTools.so'.format(FWKBASE))
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")

#pass a RDFtree object
def nanoSequence(rdftree, systType, sample, xsec, sumw, era):
    endNode="input" #this is to protect if no postnano sequence is run
    if systType == 0: #this is data
        endNode='postnano'
        rdftree.branch(nodeToStart='input', nodeToEnd='postnano', modules=[ROOT.isGoodLumi(datajson), ROOT.trigObjMatchProducer()])
        rdftree.EventFilter(nodeToStart='postnano', nodeToEnd='postnano', evfilter="isGoodLumi==true", filtername="{:20s}".format("good Lumi"))
    else: #this is mc
        luminosityN = lumi_total2016
        if era == 'preVFP' :
            luminosityN = lumi_preVFP
        else: luminosityN = lumi_postVFP
        
        clip=False
        if sample in sumwClipSamples: #add other relevant sam
            clip=True

        endNode='postnano'
        datapu = pufile_data_UL2016_preVFP
        mcprofName="Pileup_nTrueInt_Wplus_preVFP"
        eraCode=1#1 for preVFP, 2 for postVFP
            
        if era == 'postVFP' : 
            datapu=pufile_data_UL2016_postVFP
            mcprofName="Pileup_nTrueInt_Wplus_preVFP"
            eraCode=2
        if not datapu :
            print("DATA pu hist doesn't exist")
            sys.exit(2)
        if not pufile_mc_UL2016 :
            print("MC pu hist doesn't exist")
            sys.exit(2)   
        print('Using puWeight producer era code=',eraCode)
        #TFile *puMC, TFile *puData, TString hmcName, TString hdataName, bool dosyst, booql fixlargeW = true, bool normtoArea = true
        #rdftree.branch(nodeToStart='input', nodeToEnd='postnano', modules=[ROOT.puWeightProducer(mcprofName, "pileup", False, True, True), ROOT.trigObjMatchProducer()])
        #FOR ROOT HISTOS
        print(sumwClipSamples)
        print(sample,'clipping:',clip)
        rdftree.branch(nodeToStart='input', nodeToEnd='postnano', modules=[ROOT.lumiWeight(xsec=xsec, sumw=sumw, targetLumi = luminosityN, clip=clip),ROOT.puWeightProducer(eraCode), ROOT.trigObjMatchProducer(),ROOT.muonPrefireWeightProducer(filemuPrefire, eraCode)])
        if systType == 2:#for signal MC
            rdftree.branch(nodeToStart='postnano', nodeToEnd='postnano', modules=[ROOT.genLeptonSelector(), ROOT.CSvariableProducer(), ROOT.genVProducer()])
            if sample in ['WPlusJetsToMuNu', 'WMinusJetsToMuNu']:
                print('Will run jme producer for sample:' + sample)
                rdftree.branch(nodeToStart='postnano', nodeToEnd='postnano', modules=[ROOT.jmeProducer(era=eraCode)])
    return rdftree,endNode
