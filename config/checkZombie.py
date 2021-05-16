import os
import sys
import ROOT
import argparse
import copy
import time
from datetime import datetime

FWKBASE=os.getenv('FWK_BASE')

sys.path.append('{}/Common/data'.format(FWKBASE))
from samples_2016_ulV2 import samplespreVFP, samplespostVFP


ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")

def main():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-i', '--inputDir',type=str, default='/scratchnvme/wmass/NANOMAY2021/', help="input dir name")    
    parser.add_argument('-e', '--era',type=str, default='preVFP', help="either (preVFP|postVFP)")    

    args = parser.parse_args()
    now = datetime.now()
    inDir = args.inputDir
    era=args.era
    ##Add era to input dir
    inDir+=era
    
    samples = samplespreVFP
    if era == 'postVFP': 
        samples = samplespostVFP

    for sample in samples:
        print('analysing sample: %s'%sample)
        # if not ("WPlusJetsToMuNu" in sample or "WPlusJetsToTauNu" in sample): continue
        direc = samples[sample]['dir']
        for d in direc:
            targetDir='{}/{}/merged/'.format(inDir, d)
            for f in os.listdir(targetDir):#check the directory
                if not f.endswith('.root'): continue
                inputFile=targetDir+f
                tf=ROOT.TFile.Open(inputFile)
                if not tf or tf.IsZombie():
                    print("{} is zombie".format(inputFile))
                    continue
                evt=ROOT.TTree()
                evt=tf.Get("Events")
                if not evt or evt.GetEntries() == 0:
                    print("Events tree is empty in {} ".format(inputFile))
                runtree=ROOT.TTree()
                runtree=tf.Get("Runs")
                if not runtree  or runtree.GetEntries() == 0:
                    print("Run tree is empty in {} ".format(inputFile))


if __name__ == "__main__":
    main()
