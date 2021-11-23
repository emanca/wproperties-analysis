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
parser.add_argument('-a', '--runAC', default=False, action='store_true', help='Use to run the Angular Coefficients with all the variations')
parser.add_argument('-t', '--runTemplates', default=False, action='store_true', help='Use to run the templates')
parser.add_argument('-c', '--ncores',type=int, default=128, help="number of cores used")

args = parser.parse_args()
runAC = args.runAC
runTemplates = args.runTemplates
ncores = args.ncores

ROOT.ROOT.EnableImplicitMT(ncores)
print("running with {} cores".format(ncores))

inputFile = '/scratchnvme/wmass/WJetsNoCUT_v2/tree_*_*.root'

p = RDFtree(outputDir = 'GenInfo', inputFile = inputFile, outputFile="genInfo.root")
p.branch(nodeToStart='input', nodeToEnd='basicSelection', modules=[getLumiWeight(xsec=61526.7, inputFile=inputFile), ROOT.reweightFromZ(filePt, fileY),ROOT.baseDefinitions(), ROOT.defineHarmonics(), ROOT.Replica2Hessian()])

# fileACplus = ROOT.TFile.Open("./GenInfo_3D/genInput_v7_syst_Wplus.root")
# fileACminus = ROOT.TFile.Open("./GenInfo_3D/genInput_v7_syst_Wminus.root")
# p = RDFtree(outputDir = 'GenInfo_3D', inputFile = inputFile, outputFile="genInfo_v8.root")
# p.branch(nodeToStart='input', nodeToEnd='basicSelection', modules=[getLumiWeight(xsec=61526.7, inputFile=inputFile), ROOT.reweightFromZ(filePt, fileY),ROOT.baseDefinitions(), ROOT.getACValues(fileACplus,fileACminus), ROOT.defineHarmonics(), ROOT.Replica2Hessian()])

Wcharge = {"Wplus":"GenPart_pdgId[GenPart_preFSRMuonIdx]<0","Wminus":"GenPart_pdgId[GenPart_preFSRMuonIdx]>0"}
# Wcharge = {"Wplus":"GenPart_pdgId[GenPart_preFSRMuonIdx]<0 && abs(genVtype)==14","Wminus":"GenPart_pdgId[GenPart_preFSRMuonIdx]>0 && abs(genVtype)==14"}
if runAC:
    for charge,filter in Wcharge.items():
        p.branch(nodeToStart='basicSelection', nodeToEnd='angularCoefficients_{}'.format(charge), modules=[ROOT.accMap(filter), ROOT.AngCoeff(filter)])
        # p.branch(nodeToStart='basicSelection', nodeToEnd='angularCoefficients_{}'.format(charge), modules=[ROOT.accMap(filter),ROOT.getWeights(), ROOT.AngCoeff(filter)])

        #weight variations
        for s,variations in systematics.items():
            print("branching weight variations", s)
            vars_vec = ROOT.vector('string')()
            for var in variations[0]:
                vars_vec.push_back(var)
        
            p.branch(nodeToStart='basicSelection', nodeToEnd='angularCoefficients_{}_{}'.format(charge, s), modules=[ROOT.accMap(filter),ROOT.getMassWeights(), ROOT.AngCoeff(filter, vars_vec, variations[1])])

    p.getOutput()

if runTemplates:
    
    mass = ROOT.vector('string')()
    mass.push_back("_massDown")
    mass.push_back("")
    mass.push_back("_massUp")

    fileAC = ROOT.TFile.Open("genInput.root")
    p.branch(nodeToStart = 'basicSelection', nodeToEnd = 'harmonicsWeights',modules = [ROOT.getACValues(fileAC)])
    p.branch(nodeToStart = 'harmonicsWeights', nodeToEnd = 'accMap', modules =[ROOT.getAccMap(fileAC),ROOT.getMassWeights()])
    p.branch(nodeToStart = 'accMap', nodeToEnd = 'templates', modules =[ROOT.getWeights(), ROOT.templateBuilder()])
    p.branch(nodeToStart = 'accMap', nodeToEnd = 'dataObs', modules =[ROOT.dataObs(mass,"massWeights")])

    #weight variations
    for s,variations in systematics.items():
        print("branching weight variations", s)
        vars_vec = ROOT.vector('string')()
        for var in variations[0]:
            vars_vec.push_back(var)
        p.branch(nodeToStart = 'accMap', nodeToEnd = 'dataObs_{}'.format(s), modules =[ROOT.dataObs(vars_vec,variations[1])])

    p.getOutput()
    # p.saveGraph()
