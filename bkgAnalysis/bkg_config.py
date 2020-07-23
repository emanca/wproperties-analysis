import os
import copy
import math
import argparse
import sys
import ROOT
import multiprocessing

import bkg_utils
import bkg_fakerateAnalyzer

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(True)


#############################################################################################################################################
#
#  usage: python bkg_config.py --mainAna 1 --CFAna 1  --syst 1 --compAna 1 --inputDir NAME/hadded/ --outputDir OUTNAME/
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

args = parser.parse_args()
mainAna = args.mainAna
comparisonAna = args.compAna
correlatedFitAna = args.CFAna
systAna = args.syst
inputDir = args.inputDir
outputDir = args.outputDir

#internal parameters:
STATANA = True
CORRFITFINAL= True #correlated fit in the final plots
TEMPLATE = True
NOM = ['Nominal','']
EXTRAP = False #extrapolation syst
NCORES=1
if NCORES>1 :
    MULTICORE=True
else :
    MULTICORE=False


def fakerate_analysis(systKind, systName,correlatedFit,statAna, template, ptBinning=bkg_utils.ptBinning, etaBinning=bkg_utils.etaBinning, outdir=outputDir, indir=inputDir, extrapSuff='') :
    outdirBkg = outdir+'/bkg_'+systName+extrapSuff        
    if not os.path.isdir(outdirBkg): os.system('mkdir '+outdirBkg)
    
    fake = bkg_fakerateAnalyzer.bkg_analyzer(systKind=systKind,systName=systName,correlatedFit=correlatedFit,statAna=statAna, ptBinning=ptBinning, etaBinning=etaBinning, outdir=outdirBkg, inputDir=indir)
    fake.main_analysis(correlatedFit=correlatedFit,template=template,output4Plots=template,produce_ext_output=True,extrapSuff=extrapSuff)
    if template :
        fake.bkg_plots()

def runSingleAna(par) :
    fakerate_analysis_call = fakerate_analysis(systKind=par[0], systName=par[1], correlatedFit=par[2], statAna=par[3], template = par[4])


print "----> Background analyzer:"
if not os.path.isdir(outputDir): os.system('mkdir '+outputDir)

if mainAna :
    print "--> Not correlated analysis path:"
    fakerate_analysis(systKind=NOM[0], systName=NOM[1], correlatedFit=False, statAna=False, template = TEMPLATE)
    if systAna :
        processesMain = []
        for sKind, sList in bkg_utils.bkg_systematics.iteritems():
            for sName in sList :
                print "-> systematic:", sKind,sName
                if MULTICORE :
                    processesMain.append((sKind,sName,False,False,TEMPLATE))
                else :
                    fakerate_analysis(systKind=sKind, systName=sName, correlatedFit=False, statAna=False, template = TEMPLATE)
            if MULTICORE :
                poolMain = multiprocessing.Pool(NCORES)
                poolMain.map(runSingleAna,processesMain)
    if EXTRAP:
        for lcut, lbin in bkg_utils.looseCutDict.iteritems() :    
            fakerate_analysis(systKind=NOM[0], systName=NOM[1], correlatedFit=False, statAna=False, template = TEMPLATE, extrapSuff=lcut)
        fakeExtrap = bkg_fakerateAnalyzer.bkg_analyzer(systKind=NOM[0],systName=NOM[1],correlatedFit=CORRFITFINAL,statAna=False, ptBinning=bkg_utils.ptBinning, etaBinning=bkg_utils.etaBinning, outdir=outputDir+'/bkg_', inputDir=inputDir)
        if not os.path.isdir(outputDir+'/extrapolation_plots'): os.system('mkdir '+outputDir+'/extrapolation_plots')
        fakeExtrap.extrapolationSyst(extrapDict=bkg_utils.looseCutDict,linearFit=True)
        fakeExtrap.extrapolationSyst(extrapDict=bkg_utils.looseCutDict,linearFit=False)
    
if correlatedFitAna :
    print "--> Correlated analysis path:"
    fakerate_analysis(systKind=NOM[0], systName=NOM[1], correlatedFit=True, statAna=False, template = TEMPLATE)
    if STATANA :
        print "-> statistical uncertainity analysis:"
        fakerate_analysis(systKind=NOM[0], systName=NOM[1], correlatedFit=True, statAna=STATANA, template = TEMPLATE)
    if systAna :
        processesCF = []
        for sKind, sList in bkg_utils.bkg_systematics.iteritems():
            for sName in sList : 
                print "-> systematatic:", sKind,sName
                if MULTICORE :
                    processesCF.append((sKind,sName,True,False,TEMPLATE))
                else :
                    fakerate_analysis(systKind=sKind, systName=sName, correlatedFit=True, statAna=False, template = TEMPLATE)
        if MULTICORE :
            poolCF = multiprocessing.Pool(NCORES)
            poolCF.map(runSingleAna,processesCF)
            
if comparisonAna :
    print "--> Syst comparison plots..."
    fakeFinal = bkg_fakerateAnalyzer.bkg_analyzer(systKind=NOM[0],systName=NOM[1],correlatedFit=CORRFITFINAL,statAna=False, ptBinning=bkg_utils.ptBinning, etaBinning=bkg_utils.etaBinning, outdir=outputDir+'/bkg_'+NOM[1], inputDir=inputDir)
    fakeFinal.syst_comparison(systDict=bkg_utils.bkg_systematics, SymBands=True, outDir=outputDir, noratio=False, statAna=STATANA)
    fakeFinal.buildOutput(outputDir=outputDir,statAna=STATANA)
