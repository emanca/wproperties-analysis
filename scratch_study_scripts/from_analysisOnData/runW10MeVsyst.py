import os
import sys
import ROOT
import json
import argparse
import subprocess
sys.path.append('../RDFprocessor/framework')
from RDFtree import RDFtree
sys.path.append('python/')
sys.path.append('data/')
from systematics import systematics
from selections import selections, selectionVars, selections_bkg

from getLumiWeight import getLumiWeight

ROOT.gSystem.Load('bin/libAnalysisOnData.so')
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")

syst10MeV = {
    "LHEPdfWeight" : (["_LHEPdfWeightHess{}".format(i+1) for i in range(60)], "LHEPdfWeightHess"),
    "alphaS"       : (["_alphaSUp", "_alphaSDown"], "alphaSVars"),
    "mass"         : (["_massUp","_massDown"], "massWeights"),
    # "Qt"          : (["_QtUp","_QtDown"], "Qtweights"),
    "Qt"          : (["_QtTheoryUp","_QtTheoryDown","_Qt8Up","_Qt8Down","_Qt3Up","_Qt3Down","_Qt1Up","_Qt1Down","_QtCMSUp","_QtCMSDown"   ], "Qtweights")

}


def RDFprocessW10MeVSyst(fvec, outputDir, sample, xsec, fileSF, fileScale, ncores, pretendJob,SBana=False):
    ROOT.ROOT.EnableImplicitMT(ncores)
    print("running with {} cores for sample:{}".format(ncores, sample)) 
    wdecayselections = { 
        'WToMu'  : ' && (abs(genVtype) == 14)',
        'WToTau' : ' && (abs(genVtype) == 16)'
    }
    p = RDFtree(outputDir = outputDir, inputFile = fvec, outputFile="{}_plots.root".format(sample), pretend=pretendJob)
    filePt = ROOT.TFile.Open("data/histoUnfoldingSystPt_nsel2_dy3_rebin1_default.root")
    fileY = ROOT.TFile.Open("data/histoUnfoldingSystRap_nsel2_dy3_rebin1_default.root")
    p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.reweightFromZ(filePt,fileY),ROOT.baseDefinitions(True, True),ROOT.rochesterVariations(fileScale), ROOT.weightDefinitions(fileSF),getLumiWeight(xsec=xsec, inputFile=fvec, genEvsbranch = "genEventSumw"),ROOT.Replica2Hessian(),ROOT.getMassWeights()])
    for region,cut in selections_bkg.items():
        if region!='Signal' : continue 
        weight = 'float(puWeight*PrefireWeight*lumiweight*WHSF*weightPt)'
        print(weight, "NOMINAL WEIGHT")
        nom = ROOT.vector('string')()
        nom.push_back("")
        wtomu_cut = cut + wdecayselections['WToMu']
        print("Region=", region)
        print("WToMu cut=", wtomu_cut)
        #last argument refers to histo category - 0 = Nominal, 1 = Pt scale , 2 = MET scale
        print("branching nominal for region:", region) 
        #Nominal templates
        p.branch(nodeToStart = 'defs', nodeToEnd = '{}/prefit_{}/Nominal'.format('WToMu', region), modules = [ROOT.muonHistos(wtomu_cut, weight, nom,"Nom",0)])     
        #weight variations
        for s,variations in syst10MeV.items():
            
            print("branching weight variations", s)            
            #  if 'Wqt' in s : # Wqt syst using scale
            #     wtomu_modcut = wtomu_cut + variations[2]
            # else :
            #     wtomu_modcut = wtomu_cut
            var_weight = weight.replace(s, "1.") #do something only if s is in the weight 
            vars_vec = ROOT.vector('string')()
            for var in variations[0]:
                vars_vec.push_back(var)
            print(weight,"\t",var_weight, "MODIFIED WEIGHT")
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = '{}/prefit_{}/{}'.format('WToMu', region,s), modules = [ROOT.muonHistos(wtomu_cut, var_weight,vars_vec,variations[1], 0)])
    p.getOutput()
    p.saveGraph()
    #if not pretendJob:
    #split the output file into decay modes
    os.chdir(outputDir)
    os.system("ls -ltrh")
    #os.system("root -b -q ../python/splitWJets.C")
    os.chdir("../")




def main():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-p', '--pretend',type=int, default=False, help="run over a small number of event")
    parser.add_argument('-i', '--inputDir',type=str, default="/scratchnvme/wmass/NanoAOD2016-V2/", help="input dir name")
    parser.add_argument('-o', '--outputDir',type=str, default='./output/', help="output dir name")
    parser.add_argument('-c', '--ncores',type=int, default=128, help="number of cores used")
    parser.add_argument('-sb', '--SBana', type=int, default=False, help="run also on the sideband (clousure test)")

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
    with open('data/samples_2016.json') as f:
        samples = json.load(f)
    sample='WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'
    print(sample)
    direc = samples[sample]['dir']
    xsec = samples[sample]['xsec']
    fvec=ROOT.vector('string')()
    for dirname,fname in direc.items():
        inputFile = '/scratchnvme/wmass/WJetsNoCUT_v2/tree_*_*.root'
        isFile = os.path.isfile(inputFile)  
        if not isFile:
            print(inputFile, " does not exist")
            #continue
        fvec.push_back(inputFile)
        break
    if fvec.empty():
        print("No files found for directory:", samples[sample], " SKIPPING processing")
        sys.exit(1)
    print(fvec) 
    
    fileSF = ROOT.TFile.Open("data/ScaleFactors_OnTheFly.root")
    
    fileScale = ROOT.TFile.Open("data/muscales_extended.root")
    print("REMINDER: you have to modify: massWeights-->10 MeV")
    RDFprocessW10MeVSyst(fvec, outputDir, sample, xsec, fileSF, fileScale, ncores, pretendJob,SBana)


if __name__ == "__main__":
    main()
