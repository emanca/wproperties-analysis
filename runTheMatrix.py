import sys
import ROOT
import os
import copy
import argparse
import time
# from contextlib import contextmanager
################################################################################################################################################
#
#  usage: python runTheMatrix.py --inputDir INPUTDIR --outputDir OUTPUT --bkgOutput BKGOUT --ncores 64 --bkgFile MYBKG --bkgPrep 1 --bkgAna 1 --prefit 1 --plotter 1
#
#  each parameters is described below in the argparse (INPUTDIR can be omitted, OUTPUT can be the same of BKGOUT)
#
#  NB: if --bkgFile myBKG --> special string to use the file that runTheMatrix produce at step2. 
#
#################################################################################################################################################
# @contextmanager
# def cd(newdir):
#     prevdir = os.getcwd()
#     os.chdir(os.path.expanduser(newdir))
#     try: 
#         yield

parser = argparse.ArgumentParser("")
parser.add_argument('-i',   '--inputDir',           type=str, default='/scratchnvme/wmass/NanoAOD2016-V2/',    help="input dir name with the trees (data+bkg)")
parser.add_argument('-o',   '--outputDir',          type=str, default='output/',    help="output dir name of step1 and step3 (inside analysisOnData/)")
parser.add_argument('-bo',  '--bkgOutput',          type=str, default='bkg/',       help="output dir name for bkgAna (inside bkgAnalysis/)")
parser.add_argument('-f',   '--bkgFile',            type=str, default='MYBKG',      help="bkg parameters file path/name.root, or the special 'MYBKG' for the one produced in the same loop")
parser.add_argument('-c',   '--ncores',             type=int, default=128,          help="number of cores used")
parser.add_argument('-sb',  '--SBana',              type=int, default=0,            help="run also on the sideband (clousure test)")
parser.add_argument('-r',   '--regFit',             type=int, default=0,            help="apply regularization to the fit")
parser.add_argument('-t',   '--toy',                type=int, default=-1,           help="number of toy in the fit, -1=asimov")
parser.add_argument('-step0', '--genAna',           type=int, default=False,        help="run the analysisOnGen, the result are saved in analysisOngen/")
parser.add_argument('-step1', '--bkgPrep',          type=int, default=False,        help="run the bkg input preparation")
parser.add_argument('-step2', '--bkgAna',           type=int, default=False,        help="run the bkg analysis")
parser.add_argument('-step3', '--prefit',           type=int, default=False,        help="run the prefit hitograms building")
parser.add_argument('-step4', '--plotter',          type=int, default=False,        help="run the prefit plotter, the result are saved in outputDir/plot/")
parser.add_argument('-step5', '--fitPrep',          type=int, default=False,        help="run the fit input preparation (the shapefile) saved in Fit/outputDir/")
parser.add_argument('-step6', '--fit',              type=int, default=False,        help="run the fit, result saved in Fit/outputDir/, must be runned as standalone step (used cmsenv)")
parser.add_argument('-step7', '--postfit',          type=int, default=False,        help="run the postfit plotting")
parser.add_argument('-step8', '--aposteriori',      type=int, default=False,        help="run the aposteriori fit, must be runned as a standalone step (used cmsenv)")
parser.add_argument('-step10', '--regularize',      type=int, default=False,        help="regularization study (not needed in the main analysis path")

args = parser.parse_args()
inputDir = args.inputDir
outputDir = args.outputDir
bkgOutput = args.bkgOutput
bkgFile = args.bkgFile
ncores = str(args.ncores)
SBana = str(args.SBana)
regFit = str(args.regFit)
toy = str(args.toy)
step0 = args.genAna
step1 = args.bkgPrep
step2 = args.bkgAna
step3 = args.prefit
step4 = args.plotter
step5 = args.fitPrep
step6 = args.fit
step7 = args.postfit
step8 = args.aposteriori
step10 = args.regularize

if step6 and (step0 or step1 or step3) :
    print("step6 (Fit) cannot be runned in the same configuration of RDF related steps 0,1,3 (RDF)")
    assert(0)
if not step6 :
    os.system( 'source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc9-opt/setup.sh')

tic=time.time()
runTimes=[]
if bkgFile == 'MYBKG' :
    # bkgFile = '../bkgAnalysis/'+bkgOutput+'/bkg_parameters_CFstatAna.root'
    bkgFile = '../bkgAnalysis/'+bkgOutput+'/bkg_parameters.root'
    #if not step2 :
    #    raw_input('You are using self-produced bkg file but you are not running step2, are you sure? [press Enter to continue]')

if step0 :
    s0start=time.time()
    print("step5: gen analysis...")
    os.chdir('./analysisOnGen')
    os.system('python runAnalysisOnGen.py --runAC --ncores 128')
    os.system('python scripts/prepareAngularCoeff.py -o genInput -i GenInfo/genInfo.root')
    os.chdir('../')
    s0end=time.time()
    runTimes.append(s0end - s0start)
else :   runTimes.append(0.)

if step1 :
    s1start=time.time()
    print("step1: bkg input preparation... ")
    os.chdir('./analysisOnData')
    if not os.path.isdir(outputDir): os.system('mkdir '+ outputDir)
    os.system('python runAnalysisOnMC.py        -i='+inputDir+' -o='+outputDir+' -c='+ncores+' -sb='+SBana)
    os.system('python runAnalysisOnWJetsMC.py   -i='+inputDir+' -o='+outputDir+' -c='+ncores+' -sb='+SBana+' -b=1')
    os.system('python runAnalysisOnData.py      -i='+inputDir+' -o='+outputDir+' -c='+ncores+' -sb='+SBana+' -b=1')
    os.chdir('../')
    s1end=time.time()
    runTimes.append(s1end - s1start)
else :   runTimes.append(0.)

if step2 :
    s2start=time.time()
    print("step2: bkg analysis...")
    if not os.path.isdir('bkgAnalysis/'+bkgOutput): os.system('mkdir bkgAnalysis/'+bkgOutput)
    if not os.path.isdir('bkgAnalysis/'+bkgOutput+'/bkgInput/'): os.system('mkdir bkgAnalysis/'+bkgOutput+'/bkgInput/')
    os.chdir('bkgAnalysis')
    #inputDir for bkgAnalysis is the outputDir of step1
    os.system('python bkg_prepareHistos.py --inputDir ../analysisOnData/'+outputDir+' --outputDir '+bkgOutput+'/bkgInput/')
    os.system('python bkg_config.py --mainAna 1 --CFAna 0 --inputDir '+bkgOutput+'/bkgInput/hadded/ --outputDir '+bkgOutput+' --syst 1 --compAna 1 --ncores '+ncores+ ' --SBana '+SBana)
    os.chdir('../')
    s2end=time.time()
    runTimes.append(s2end - s2start)
else :   runTimes.append(0.)

if step3 :
    s3start=time.time()
    print("step3: Fake from Data...")
    os.chdir('analysisOnData')
    if not os.path.isdir(outputDir): os.system('mkdir '+ outputDir)
    #ncores is optimized and set in the config itself, so no need to pass here
    os.system('python runAnalysisOnData.py      -i='+inputDir+' -o='+outputDir+' -c='+ncores+' -sb='+SBana+' -b=0 -f'+bkgFile)
    # os.system('python runAnalysisOnWJetsMC.py   -i='+inputDir+' -o='+outputDir+' -c='+ncores+' -sb='+SBana+' -b=0')
    os.chdir('../')
    s3end=time.time()
    runTimes.append(s3end - s3start)
else :   runTimes.append(0.)

if step4 :
    s4start=time.time()
    print("step4: prefit plotter...")
    os.chdir('analysisOnData/python')
    if not os.path.isdir('../'+outputDir): os.system('mkdir ../'+outputDir)
    os.system('python plotter_prefit.py --hadd 1 --output ../'+outputDir+'/plot/ --input ../'+outputDir+' --systComp 1'+' -sb='+SBana)

    sys.path.append('../../bkgAnalysis')
    import bkg_utils
    for sKind,sList in bkg_utils.bkg_systematics.items() : 
        skipList = ""
        for sKindInt,sListInt in bkg_utils.bkg_systematics.items() :
            if sKindInt==sKind : continue
            else : skipList+= ' '+str(sKindInt)
        print("Skipped systematics:", skipList) 
        os.system('python plotter_prefit.py --hadd 0 --output ../'+outputDir+'/plot_only_'+str(sKind)+' --input ../'+outputDir+'/plot/  --systComp 1 --skipSyst '+skipList+' -sb='+SBana)
    os.system('python plotter_prefit.py --hadd 0 --output ../'+outputDir+'/plot_bkgOnlySyst/ --input ../'+outputDir+'/plot/  --systComp 1'+' -sb='+SBana+' --bkgOnlySyst 1')
    os.chdir('../../')
    s4end=time.time()
    runTimes.append(s4end - s4start)
else :   runTimes.append(0.)

if step5 :
    s5start=time.time()
    print("step5: fit shape file preparation...")
    os.chdir('analysisOnData')
    os.system('python runAnalysisOnWJetsMC.py   -i='+inputDir+' -o='+outputDir+' -c='+ncores+' -sb='+SBana+' -b=0')
    os.chdir('../')
    os.chdir('analysisOnData/python/')
    if not os.path.isdir('../'+outputDir+'/plot/hadded/'): 
        print('missing input directory ../'+outputDir+'/plot/hadded/')
        assert(0)
    os.system('cp ../'+outputDir+'/WToMu_AC_plots.root ../'+outputDir+'/plot/hadded/')
    os.system('python plotter_SignalACtemplates.py --input ../' + outputDir+'/plot/hadded/') 
    os.system('python plotter_template2D.py --input ../' + outputDir+'/plot/hadded/') 
    if not os.path.isdir('../../Fit/'+outputDir): os.system('mkdir ../../Fit/'+outputDir)
    os.system('mv Wminus_reco.root Wplus_reco.root accMap.root ../../Fit/'+outputDir)
    os.chdir('../../')
    s5end=time.time()
    runTimes.append(s5end - s5start)
else :   runTimes.append(0.)

if step6 :
    s6start=time.time()
    print("step6: run the fit...")
    if not os.path.isdir('Fit/'+outputDir): os.system('mkdir Fit/'+outputDir)
    if os.path.isdir('Fit/CMSSW_10_6_18/CMSSW_10_6_18/src/HiggsAnalysis/CombinedLimit/src/'): 
        os.chdir('Fit/CMSSW_10_6_18/CMSSW_10_6_18/src/HiggsAnalysis/CombinedLimit/src/')
        os.system('cmsenv')
        os.chdir('../../../../../'+outputDir)
    else :
        print('missing Fit/CMSSW_10_6_18/CMSSW_10_6_18/src/HiggsAnalysis/CombinedLimit/src/')
        assert(0)
    os.system('python ../runFit.py --impact 1 --postfit 1 --regularize '+regFit+' --tau=1e3 --toy '+toy+ ' --cores '+ncores)
    os.chdir('../../')
    s6end=time.time()
    runTimes.append(s6end - s6start)
else :   runTimes.append(0.)
    
if step7 :
    s7start=time.time()
    print("step7: fit result plotting...")
    os.chdir('Fit/'+outputDir)
    os.system("python ../plotter_fitResult.py --output fitPlots_Wplus --fitFile fit_Wplus_reco.root --suff reco --input ../../analysisOnGen/genInput_Wplus.root --impact 1 --postfit 1")
    if regFit :
        os.system('python ../rebuild_regulFunc.py -o fitRegulFunc_Wplus -f fit_Wplus_reco.root -i fitPlots_Wplus_angCoeff.root')
        os.system('python ../plotter_fitResult.py --output fitPlots_Wplus_rebuild --fitFile fit_Wplus_reco.root --suff reco --input ../../analysisOnGen/genInput_Wplus.root --impact 1 --aposteriori fitRegulFunc_Wplus.root')
    os.chdir('../../')
    s7end=time.time()
    runTimes.append(s7end - s7start)
    print("note: to have resonable impact for the mass additional procedure must be done. see comments in runTheMatrix")
    # - the mass nuisance type must be changed to "shapeNoConstraint" and the fit must be repeated --> [FIT1]
    # - to obtain the statistical uncertainity on the mass another run of the fit is needed, without all the syst but "mass" --> [FIT2]
    # - then the plotter_fitResult.py must be runned on [FIT2] --> [OUT2]
    # - then the output of this plotter must be used as "--mass [OUT2]" parameter to plot [FIT1]
else :   runTimes.append(0.)

if step8 :
    s8start=time.time()
    print("step8: aposteriori fit...")
    # if not os.path.isdir('Fit/aposterioriFit/'+outputDir): os.system('mkdir Fit/aposterioriFit/'+outputDir)
    os.chdir('Fit/aposterioriFit/')
    if os.path.isdir('CMSSW_11_2_0_pre8/src/'): 
        os.chdir('CMSSW_11_2_0_pre8/src/')
        os.system('cmsenv')
        os.chdir('../../')
    os.system('python aposterioriFit.py --fitInput ../'+outputDir+'/fitPlots_Wplus_angCoeff.root --output ../'+outputDir+'/aposterioriFit_Wplus')
    os.chdir('../'+outputDir)
    os.system('python ../plotter_fitResult.py --output fitPlots_Wplus_aposteriori --fitFile fit_Wplus_reco.root --suff reco --input ../../analysisOnGen/genInput_Wplus.root --impact 1 --postfit 1 --aposteriori aposterioriFit_Wplus.root')
    os.chdir('../../')
    s8end=time.time()
    runTimes.append(s8end - s8start)
else :   runTimes.append(0.)

if step10 :
    s10start=time.time()
    print("step10: regularization study...")
    # if not os.path.isdir('Fit/aposterioriFit/'+outputDir): os.system('mkdir Fit/aposterioriFit/'+outputDir)
    os.chdir('regularization/')
    os.system('python regularizationFit.py --validation_only 0 -o '+outputDir+' --map01 0 -i ../analysisOnGen/GenInfo_backup/genInput_fineBinned_4regularization_fix_Wplus.root ./analysisOnGen/GenInfo_backup/genInput_fineBinned_4regularization_fix_Wminus.root')
    os.system('python regFit_jobLauncher.py -o '+outputDir+' -i ../analysisOnGen/GenInfo_backup/genInput_fineBinned_4regularization_fix_Wplus.root ./analysisOnGen/GenInfo_backup/genInput_fineBinned_4regularization_fix_Wminus.root')
    os.system('python regularizationFit.py --validation_only 0 -o '+outputDir+' --map01 0 --syst_val 1 -i ../analysisOnGen/GenInfo_backup/genInput_fineBinned_4regularization_fix_Wplus.root ./analysisOnGen/GenInfo_backup/genInput_fineBinned_4regularization_fix_Wminus.root')
    os.chdir('../')
    s10end=time.time()
    runTimes.append(s10end - s10start)
else :   runTimes.append(0.)
    
        
toc=time.time()
for i in range(10) :
    print("Step"+str(i)+" completed in:", runTimes[i], " seconds")
print("Total runtime:", toc - tic, " seconds")
