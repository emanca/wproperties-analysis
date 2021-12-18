import os
import sys
import ROOT
import json
import argparse
import copy
from RDFtree import RDFtree
from multiprocessing import Process, cpu_count
from runAnalysisOnMC import RDFprocessMC
from runAnalysisOnWJetsMC import RDFprocessWJetsMC
from runAnalysisOnData import RDFprocessData, RDFprocessfakefromData

sys.path.append('python/')
sys.path.append('data/')
from systematics import systematics
from selections import selections, selectionVars, selections_bkg
from getLumiWeight import getLumiWeight

ROOT.gSystem.Load('bin/libAnalysisOnData.so')
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;");

parser = argparse.ArgumentParser("")
parser.add_argument('-p', '--pretend',type=int, default=0, help="run over a small number of event")
parser.add_argument('-c', '--ncores',type=int, default=64, help="number of cores used")
parser.add_argument('-o', '--outputDir',type=str, default='./output/', help="output dir name")
parser.add_argument('-i', '--inputDir',type=str, default='/scratchnvme/wmass/NanoAOD2016-V2/', help="input dir name")
parser.add_argument('-b', '--runBKG',type=int, default=False, help="prepare the input of the bkg analysis, if =false run the prefit Plots")
parser.add_argument('-f', '--bkgFile',type=str, default='/scratch/bertacch/wmass/wproperties-analysis/bkgAnalysis/TEST_runTheMatrix/bkg_parameters_C\
FstatAna.root', help="bkg parameters file path/name.root")
parser.add_argument('-s', '--SBana',type=int, default=False, help="run also on the sideband (clousure test)")
args = parser.parse_args()
pretendJob = args.pretend
ncores = args.ncores
outputDir = args.outputDir
inDir = args.inputDir
runBKG = args.runBKG
bkgFile = args.bkgFile
SBana = args.SBana

#in the new machine, optimal performance with 128 cores is seen.
# ncmax = cpu_count()/2 if cpu_count() > 64 else cpu_count()

samples={}
with open('data/samples_2016.json') as f:
  samples = json.load(f)

multiprocs=[]
if runBKG : #produces templates for all regions and prefit for signal
  for sample in samples:
      if not samples[sample]['datatype']=='MC': continue
      print("Sample:", sample)
      #if not samples[sample]['multiprocessing']: continue
      #WJets to be run by a separate config
      #if 'WJets' in sample : continue
    #   print(sample)
      direc = samples[sample]['dir']
      xsec = samples[sample]['xsec']
      cores = 0 if samples[sample]['multiprocessing'] else ncores
      fvec=ROOT.vector('string')()
      for dirname,fname in direc.items():
          ##check if file exists or not
          #inputFile = '/scratchssd/sroychow/NanoAOD2016-V2/{}/tree.root'.format(dirname)
          inputFile = '{}/{}/tree.root'.format(inDir, dirname)
          isFile = os.path.isfile(inputFile)  
          if not isFile:
              print(inputFile, " does not exist")
              continue
          fvec.push_back(inputFile)
          print(inputFile)
      if fvec.empty():
         print("No files found for directory:", samples[sample], " SKIPPING processing")
         continue
      print(fvec) 
      fileSF = ROOT.TFile.Open("data/ScaleFactors_OnTheFly.root")
      fileScale = ROOT.TFile.Open("data/muscales_extended.root")
      systType = samples[sample]['systematics']

      if 'WJets' in sample : #nc = ncpus/2
          wjfvec=ROOT.vector('string')()
          wjfvec.push_back('/scratchnvme/wmass/WJetsNoCUT_v2/tree_*_*.root')
          print("WJETS is run with:", wjfvec)
          p = Process(target=RDFprocessWJetsMC, args=(wjfvec, outputDir, sample, xsec, fileSF, fileScale, ncores, pretendJob, runBKG,SBana))
          p.start()
          multiprocs.append(p)
      else:
          p = Process(target=RDFprocessMC, args=(fvec, outputDir, sample, xsec, fileSF, cores, systType, pretendJob,SBana))	
          p.start()
          multiprocs.append(p)


###For Data
print("Sample:DATA")
fvec=ROOT.vector('string')()
for sample in samples:
    if not samples[sample]['datatype']=='DATA': continue
    direc = samples[sample]['dir']
    for dirname,fname in direc.items():
     #inputFile = '/scratchssd/sroychow/NanoAOD2016-V2/{}/tree.root'.format(dirname)   
        inputFile = '{}/{}/tree.root'.format(inDir, dirname)
        isFile = os.path.isfile(inputFile)
        if not isFile:
            print(inputFile, " does not exist")
            continue
        fvec.push_back(inputFile)
if fvec.empty():
    print("No files found for json provided\n")
    sys.exit(1)
print(fvec)
if runBKG : #produces templates for all regions and prefit for signal
    p = Process(target=RDFprocessData, args=(fvec, outputDir, ncores, pretendJob, SBana))
    p.start()
    multiprocs.append(p)
else : #produces Fake contribution to prefit plots computed from data
    p = Process(target=RDFprocessfakefromData, args=(fvec, outputDir, bkgFile, ncores, pretendJob, SBana))
    p.start()
    multiprocs.append(p)

for p in multiprocs:  
    #p.start()
    p.join()
