import os
import sys
import ROOT
import argparse
sys.path.append('../RDFprocessor/framework')

from RDFtree import RDFtree
sys.path.append('python/')
sys.path.append('utils/')

from systematics import systematics
from getLumiWeight import getLumiWeight

ROOT.gSystem.Load('bin/libSignalAnalysis.so')

filePt = ROOT.TFile.Open("../analysisOnData/data/histoUnfoldingSystPt_nsel2_dy3_rebin1_default.root")
fileY = ROOT.TFile.Open("../analysisOnData/data/histoUnfoldingSystRap_nsel2_dy3_rebin1_default.root")

parser = argparse.ArgumentParser('')
parser.add_argument('-c', '--ncores',type=int, default=128, help="number of cores used")

args = parser.parse_args()
ncores = args.ncores

ROOT.ROOT.EnableImplicitMT(ncores)
print("running with {} cores".format(ncores))

inputFile = '/scratchnvme/wmass/WJetsNoCUT_v2/tree_*_*.root'

p = RDFtree(outputDir = 'GenInfo', inputFile = inputFile, outputFile="study_runAcceptanceMap.root")
p.branch(nodeToStart='input', nodeToEnd='basicSelection', modules=[getLumiWeight(xsec=61526.7, inputFile=inputFile), ROOT.reweightFromZ(filePt, fileY),ROOT.baseDefinitions()])

Wcharge = {"Wplus":"GenPart_pdgId[GenPart_preFSRMuonIdx]<0","Wminus":"GenPart_pdgId[GenPart_preFSRMuonIdx]>0"}
for charge,filter in Wcharge.items():
    p.branch(nodeToStart='basicSelection', nodeToEnd='angularCoefficients_{}'.format(charge), modules=[ROOT.accMap(filter)])
p.getOutput()

