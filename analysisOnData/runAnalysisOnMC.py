import os
import sys
import ROOT
import json
import argparse
import copy
sys.path.append('../RDFprocessor/framework')
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
    filePt = ROOT.TFile.Open("data/histoUnfoldingSystPt_nsel2_dy3_rebin1_default.root")
    fileY = ROOT.TFile.Open("data/histoUnfoldingSystRap_nsel2_dy3_rebin1_default.root")
    #sample specific systematics
    systematicsFinal=copy.deepcopy(systematics)
    if systType != 2:
        p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.baseDefinitions(True, False),ROOT.weightDefinitions(fileSF),getLumiWeight(xsec=xsec, inputFile=fvec)])
    else: #this run for DY sample only (only DY has PDF and need reweighting)
        p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.reweightFromZ(filePt,fileY,False),ROOT.baseDefinitions(True, False),ROOT.weightDefinitions(fileSF),getLumiWeight(xsec=xsec, inputFile=fvec),ROOT.Replica2Hessian()])
    #selections bkg also includes Signal
    for region,cut in selections_bkg.items():
        print("running in region {}".format(region))
                    
        if 'aiso' in region:
            weight = 'float(puWeight*PrefireWeight*lumiweight)'
        else:
            weight = 'float(puWeight*PrefireWeight*lumiweight*WHSF)'
        
        if systType == 2: #this run for DY sample only (only DY has PDF and need reweighting)
            # weight = weight + '*float(weightY)'
            weight = weight + '*float(weightPt*weightY)'
        
        print(weight, "NOMINAL WEIGHT")

        nom = ROOT.vector('string')()
        nom.push_back("")
        #last argument refers to histo category - 0 = Nominal, 1 = Pt scale , 2 = MET scale
        print("branching nominal")

        if region == "Signal" or (region=='Sideband' and SBana):
            p.branch(nodeToStart = 'defs', nodeToEnd = 'prefit_{}/Nominal'.format(region), modules = [ROOT.muonHistos(cut, weight, nom,"Nom",0)])  
        p.branch(nodeToStart = 'defs', nodeToEnd = 'templates_{}/Nominal'.format(region), modules = [ROOT.templates(cut, weight, nom,"Nom",0)])    
  
       #weight variations
        for s,variations in systematicsFinal.items():
            if 'WQT' in s : continue 
            if "LHEScaleWeight" in s and systType != 2 :  continue #==2 --> DY (and W)
            if "LHEPdfWeight"   in s and systType != 2 :  continue
            if "alphaS"         in s and systType != 2 :  continue
            if 'WHSF' in s and 'aiso' in region: continue
            if 'mass' in s : continue
            
            # if s not in weight :  #qui entra aiso-WHSF e non deve!
            #     var_weight = weight
            # else:
                # if 'WHSF' in s and 'aiso' in region: continue
                # var_weight = weight.replace(s, "1.")
            var_weight = weight.replace(s, "1.")

            vars_vec = ROOT.vector('string')()
            for var in variations[0]:
                vars_vec.push_back(var)
            print("branching weight variations", s)
            print(weight,var_weight, "MODIFIED WEIGHT")
            if region == "Signal" or (region=='Sideband' and SBana): 
                p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = 'prefit_{}/{}'.format(region,s), modules = [ROOT.muonHistos(cut,var_weight,vars_vec,variations[1], 0)])
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = 'templates_{}/{}'.format(region,s), modules = [ROOT.templates(cut,var_weight,vars_vec,variations[1], 0)])

        #column variations#weight will be nominal, cut will vary
        for vartype, vardict in selectionVars.items():
            cut_vec = ROOT.vector('string')()
            var_vec = ROOT.vector('string')()
            for selvar, hcat in vardict.items() :
                newcut = cut.replace('MT', 'MT'+selvar)
                if 'corrected' in selvar:
                    newcut = newcut.replace('Mu1_pt', 'Mu1_pt'+selvar)
                cut_vec.push_back(newcut)
                var_vec.push_back(selvar)
            if region == "Signal" or (region=='Sideband' and SBana):
                p.branch(nodeToStart = 'defs', nodeToEnd = 'prefit_{}/{}'.format(region,vartype), modules = [ROOT.muonHistos(cut_vec, weight, nom,"Nom",hcat,var_vec)])  
            p.branch(nodeToStart = 'defs', nodeToEnd = 'templates_{}/{}'.format(region,vartype), modules = [ROOT.templates(cut_vec, weight, nom,"Nom",hcat,var_vec)])  

    p.getOutput()
    # p.saveGraph()


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
    with open('data/samples_2016.json') as f:
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
