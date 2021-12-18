import os
import copy
import math
import argparse
import sys
import ROOT
import multiprocessing
import numpy as np
import root_numpy as rootnp

import bkg_utils
import bkg_fakerateAnalyzer

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(True)


#############################################################################################################################################
#
#  usage: python bkg_config.py --mainAna 1 --CFAna 1  --syst 1 --compAna 1 --inputDir NAME/hadded/ --outputDir OUTNAME/ --ncores 64
#
#  the full path of this config is the following
#  1) run the nominal analysis (mainAna)
#  2) run the analysis for each systematic variation (from bkg_utils dictionary)
#  3) run again the nominal analysis but wiht the correlatedFit to fit in pt the fakerate_analysis
#  4) run each systematic variation with the correlatedFit enabled
#  5) run the comparison plots between the various syst. This step is not needed to produce the bkg output, but is a validation step only
#  6) build the output file using the output from step 1-4
#
#  Internal parameters:
#  - STATANA: if true run a proper estimation of the statistical uncertainity, if the correlatedFit is used
#  - CORRFITFINAL: if true make the comparison of the step 5) using the correlatedFit results
#  - TEMPLATE: if true plots also the QCD templates and not only the fakerate and propmptrate
#  - NOM: nominal name in the input
#  - EXTRAP: 3 additional runs of the nominal analysis to estimate the extrapolation uncertainity
#  - NCORES: number of cores used, if >1 activate multiprocessing
#
#  Code structure: 
#   - bkg_config.py use bkg_utils.py to load some common dictionaries and lists
#   - bkg_fakerateAnalyzer.py is the class which is used to perform each step of the analysis. 
#     Each time the main analysis is run is a different class object 
#     (nom + one for each syst, for main and correlatedFit analyses + one final for comparison)
#  - several class-object calls are realized with an helper function: fakerate_analysis, which configure the object and run some methods
#
#############################################################################################################################################

parser = argparse.ArgumentParser("")
parser.add_argument('-systAna', '--syst',type=int, default=False, help="enable systematics analysis")
parser.add_argument('-mainAna', '--mainAna',type=int, default=False, help="main bkg analysis")
parser.add_argument('-correlatedFitAna', '--CFAna',type=int, default=False, help="second run of main ana, with correlated Fit")
parser.add_argument('-comparisonAna', '--compAna',type=int, default=False, help="comparison plots of systematic analysis")
parser.add_argument('-inputDir', '--inputDir',type=str, default='./data/', help="input dir name")
parser.add_argument('-outputDir', '--outputDir',type=str, default='./bkg_V2/', help="output dir name")
parser.add_argument('-ncores', '--ncores',type=int, default=64, help="number of cores used")
parser.add_argument('-s', '--SBana',type=int, default=False, help="run also on the sideband (clousure test)")

args = parser.parse_args()
mainAna = args.mainAna
comparisonAna = args.compAna
correlatedFitAna = args.CFAna
systAna = args.syst
inputDir = args.inputDir
outputDir = args.outputDir
ncores = args.ncores
SBana = args.SBana


#internal parameters:
CORRFITFINAL= False #correlated fit in the final plots
STATANA = CORRFITFINAL  #or False or CORRFITFINAL
TEMPLATE = True
NOM = ['Nominal','']
EXTRAP = False #extrapolation syst
EXTRAPCORR = False #extrapolation correction. not applied by default 
NCORES=ncores
if NCORES>1 :
    MULTICORE=True
else :
    MULTICORE=False

def getPtEtaBinning(fileName, histoPathName='templates_Signal/Nominal/templates') :
    openFile = ROOT.TFile.Open(fileName) 
    refHisto = openFile.Get(histoPathName)
    hArr,binArr = rootnp.hist2array(hist=refHisto,return_edges=True)
    # print("etaBinning", binArr[0])
    # print("ptBinning", binArr[1])
    
    etaBins = binArr[0].tolist()
    ptBins = binArr[1].tolist()
    etaBins = [round(el,2) for el in etaBins]
    ptBins = [round(el,2) for el in ptBins]
    
    return etaBins,ptBins
    
etaBins,ptBins = getPtEtaBinning(fileName=inputDir+'/WToMuNu.root')
   
def fakerate_analysis(systKind, systName,correlatedFit,statAna, template, ptBinning=ptBins, etaBinning=etaBins, outdir=outputDir, indir=inputDir, extrapSuff='',regionSuff='') :
    outdirBkg = outdir+'/bkg_'+systName+extrapSuff        
    if not os.path.isdir(outdirBkg): os.system('mkdir '+outdirBkg)
    
    fake = bkg_fakerateAnalyzer.bkg_analyzer(systKind=systKind,systName=systName,correlatedFit=correlatedFit,statAna=statAna, ptBinning=ptBinning, etaBinning=etaBinning, outdir=outdirBkg, inputDir=indir,extrapCorr=EXTRAPCORR, nameSuff=regionSuff)
    fake.main_analysis(correlatedFit=correlatedFit,template=template,output4Plots=template,produce_ext_output=True,extrapSuff=extrapSuff)
    if template :
        fake.bkg_plots()

def runSingleAna(par) :
    fakerate_analysis_call = fakerate_analysis(systKind=par[0], systName=par[1], correlatedFit=par[2], statAna=par[3], template = par[4], regionSuff=par[5])


print("----> Background analyzer:")
if not os.path.isdir(outputDir): os.system('mkdir '+outputDir)

regionList = ['']
if SBana :
    regionList.append('SideBand')

for reg in regionList : 
    if mainAna :
        print("--> Not correlated analysis path:")
        fakerate_analysis(systKind=NOM[0], systName=NOM[1], correlatedFit=False, statAna=False, template = TEMPLATE, regionSuff=reg)
        if systAna :
            processesMain = []
            for sKind, sList in bkg_utils.bkg_systematics.items():
                # if sKind=='LHEScaleWeight' : continue 
                for sName in sList :
                    print("-> systematic:", sKind,sName)
                    if MULTICORE :
                        processesMain.append((sKind,sName,False,False,TEMPLATE,reg))
                    else :
                        fakerate_analysis(systKind=sKind, systName=sName, correlatedFit=False, statAna=False, template = TEMPLATE, regionSuff=reg)
                if MULTICORE :
                    poolMain = multiprocessing.Pool(NCORES)
                    poolMain.map(runSingleAna,processesMain)

    if correlatedFitAna :
        print("--> Correlated analysis path:")
        fakerate_analysis(systKind=NOM[0], systName=NOM[1], correlatedFit=True, statAna=False, template = TEMPLATE, regionSuff=reg)
        if STATANA :
            print("-> statistical uncertainity analysis:")
            fakerate_analysis(systKind=NOM[0], systName=NOM[1], correlatedFit=True, statAna=STATANA, template = TEMPLATE, regionSuff=reg)
        if systAna :
            processesCF = []
            for sKind, sList in bkg_utils.bkg_systematics.items():
                # if sKind=='LHEScaleWeight' : continue 
                for sName in sList : 
                    print("-> systematatic:", sKind,sName)
                    if MULTICORE :
                        processesCF.append((sKind,sName,True,False,TEMPLATE,reg))
                    else :
                        fakerate_analysis(systKind=sKind, systName=sName, correlatedFit=True, statAna=False, template = TEMPLATE, regionSuff=reg)
            if MULTICORE :
                poolCF = multiprocessing.Pool(NCORES)
                poolCF.map(runSingleAna,processesCF)

    if EXTRAP and reg!='SideBand' :
        print("--> Extrapolation syst evaluation:")
        CFstring = ''
        if CORRFITFINAL :
            CFstring+='_CF'
        if STATANA :
            CFstring+='statAna'
        if (mainAna or correlatedFitAna):
            for lcut, lbin in bkg_utils.looseCutDict.items() :    
                print("-> Mt bin:", lbin)
                fakerate_analysis(systKind=NOM[0], systName=NOM[1], correlatedFit=False, statAna=False, template = TEMPLATE, extrapSuff=lcut)
        fakeExtrap = bkg_fakerateAnalyzer.bkg_analyzer(systKind=NOM[0],systName=NOM[1],correlatedFit=CORRFITFINAL,statAna=False, ptBinning=ptBins, etaBinning=etaBins, outdir=outputDir, inputDir=inputDir)
        if not os.path.isdir(outputDir+'/extrapolation_syst'): os.system('mkdir '+outputDir+'/extrapolation_syst')
        fakeExtrap.extrapolationSyst(extrapDict=bkg_utils.looseCutDict,linearFit=True, CFstring=CFstring)
        fakeExtrap.extrapolationSyst(extrapDict=bkg_utils.looseCutDict,linearFit=False, CFstring=CFstring)
        if EXTRAPCORR :
            print("extrapCorr")
            fakerate_analysis(systKind=NOM[0], systName=NOM[1], correlatedFit=CORRFITFINAL, statAna=STATANA, template = TEMPLATE)
                    
    if comparisonAna :
        print("--> Syst comparison plots:")
        fakeFinal = bkg_fakerateAnalyzer.bkg_analyzer(systKind=NOM[0],systName=NOM[1],correlatedFit=CORRFITFINAL,statAna=False, ptBinning=ptBins, etaBinning=etaBins, outdir=outputDir+'/bkg_'+NOM[1], inputDir=inputDir, nameSuff=reg)
        fakeFinal.syst_comparison(systDict=bkg_utils.bkg_systematics, SymBands=True, outDir=outputDir, noratio=False, statAna=STATANA)
        fakeFinal.buildOutput(outputDir=outputDir,statAna=STATANA)
