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
parser.add_argument('-i', '--inputDir',type=str, default='/scratchnvme/wmass/NanoAOD2016-V2/', help="input dir name with the trees")
parser.add_argument('-o', '--outputDir',type=str, default='output/', help="output dir name of step1 and step3 (inside analysisOnData/)")
parser.add_argument('-t', '--bkgOutput',type=str, default='bkg/', help="output dir name for bkgAna (inside bkgAnalysis/)")
parser.add_argument('-f', '--bkgFile',type=str, default='/scratch/bertacch/wmass/wproperties-analysis/bkgAnalysis/TEST_runTheMatrix/bkg_parameters_CFstatAna.root', help="bkg parameters file path/name.root, or the special 'MYBKG' for the one produced in the same loop")
parser.add_argument('-c', '--ncores',type=int, default=128, help="number of cores used")
parser.add_argument('-q', '--bkgPrep',type=int, default=False, help="run the bkg input preparation")
parser.add_argument('-w', '--bkgAna',type=int, default=False, help="run the bkg analysis")
parser.add_argument('-e', '--prefit',type=int, default=False, help="run the prefit hitograms building")
parser.add_argument('-r', '--plotter',type=int, default=False, help="run the prefit plotter, the result are saved in outputDir/plot/")
parser.add_argument('-sb', '--SBana',type=int, default=0, help="run also on the sideband (clousure test)")


args = parser.parse_args()
inputDir = args.inputDir
outputDir = args.outputDir
bkgOutput = args.bkgOutput
bkgFile = args.bkgFile
ncores = str(args.ncores)
SBana = str(args.SBana)
step1 = args.bkgPrep
step2 = args.bkgAna
step3 = args.prefit
step4 = args.plotter

tic=time.time()
runTimes=[]
if bkgFile == 'MYBKG' :
    bkgFile = '../bkgAnalysis/'+bkgOutput+'/bkg_parameters_CFstatAna.root'
    #if not step2 :
    #    raw_input('You are using self-produced bkg file but you are not running step2, are you sure? [press Enter to continue]')

if step1 :
    s1start=time.time()
    print("step1: bkg input preparation... ")
    os.chdir('./analysisOnData')
    if not os.path.isdir(outputDir): os.system('mkdir '+ outputDir)
    os.system('python runAnalysis.py -p=0 -b=1 -i='+inputDir+' -c='+ncores+' -o='+outputDir + ' -sb='+SBana)
    # os.system('python runAnalysisOnMC.py      -i='+inputDir+' -c='+ncores+' -o='+outputDir + ' -sb='+SBana)
    # os.system('python runAnalysisOnWJetsMC.py -i='+inputDir+' -c='+ncores+' -o='+outputDir + ' -sb='+SBana +' -b=1')
    # os.system('python runAnalysisOnData.py    -i='+inputDir+' -c='+ncores+' -o='+outputDir + ' -sb='+SBana +' -b=1')
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
    os.system('python bkg_config.py --mainAna 1 --CFAna 1 --inputDir '+bkgOutput+'/bkgInput/hadded/ --outputDir '+bkgOutput+' --syst 1 --compAna 1 --ncores '+ncores+ ' --SBana '+SBana)
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
    os.system('python runAnalysis.py -b=0 -i='+inputDir+ ' -o=' +outputDir+ ' -f='+bkgFile + ' -sb='+SBana)
    # os.system('python runAnalysisOnData.py    -i='+inputDir+' -c='+ncores+' -o='+outputDir + ' -sb='+SBana +' -b=0' + ' -f='+bkgFile)
    os.chdir('../')
    s3end=time.time()
    runTimes.append(s3end - s3start)
else :   runTimes.append(0.)

if step4 :
    s4start=time.time()
    print("step4: plotter...")
    os.chdir('analysisOnData/python')
    if not os.path.isdir('../'+outputDir): os.system('mkdir ../'+outputDir)
    os.system('python plotter_prefit.py --hadd 1 --output ../'+outputDir+'/plot/ --input ../'+outputDir+' --systComp 1'+' -sb='+SBana)

    if not os.path.isdir('../'+outputDir+'/plot/hadded/template2D/'): os.system('mkdir -p ../'+outputDir+'/plot/hadded/template2D/')
    os.system('python plotter_template2D.py -o=../'+outputDir+'/plot/hadded/template2D/ -i=../'+outputDir+'/plot/hadded')

    sys.path.append('../../bkgAnalysis')
    import bkg_utils
    for sKind,sList in bkg_utils.bkg_systematics.items() : 
        skipList = ""
        for sKindInt,sListInt in bkg_utils.bkg_systematics.items() :
            if sKindInt==sKind : continue
            else : skipList+= ' '+str(sKindInt)
        print("Skipped systematics:", skipList) 
        os.system('python plotter_prefit.py --hadd 0 --output ../'+outputDir+'/plot_only_'+str(sKind)+' --input ../'+outputDir+'/plot/  --systComp 1 --skipSyst '+skipList+' -sb='+SBana)
    s4end=time.time()
    runTimes.append(s4end - s4start)
else :   runTimes.append(0.)
    
toc=time.time()
print("Step1 completed in:", runTimes[0], " seconds")
print("Step2 completed in:", runTimes[1], " seconds")
print("Step3 completed in:", runTimes[2], " seconds")
print("Step4 completed in:", runTimes[3], " seconds")
print("Total runtime:", toc - tic, " seconds")
