import ROOT
import json
import sys
import os

sys.path.append('data')

ROOT.gSystem.CompileMacro("getClippedSumW.cpp", "kO")
#from samples_2016_ulCentral import samplespreVFP, samplespostVFP
from samples_2016_ulV2 import samplespreVFP, samplespostVFP, wsignalNLO_preVFP, wsignalNLO_postVFP
eras = ["preVFP", "postVFP"]

def genwtsum(sample, fvec) :
    print("processing ", sample)

ROOT.ROOT.EnableImplicitMT(64)

for era in eras:
    inDir='/scratchnvme/wmass/NANOJEC/'
    inDir+=era
    samples = samplespreVFP if era=="preVFP" else samplespostVFP
    # samples = wsignalNLO_preVFP if era=="preVFP" else wsignalNLO_postVFP
    sumwProc={}
    sumwRunProc = {}
    objList = []
    for sample in samples:
        if 'data' in sample: continue
        if not "WPlusJetsToMuNu" in sample and not "WMinusJetsToMuNu" in sample: continue
        computeYes= 'DY' in sample or 'WPlus' in sample or 'WMinus' in sample
        fvec=ROOT.vector('string')()
        direc = samples[sample]['dir']
        for d in direc:
            targetDir='{}/{}/'.format(inDir, d)
            # mergedDir=targetDir+'merged/'
            # if os.path.isdir(mergedDir):
            #     targetDir=mergedDir
            for f in os.listdir(targetDir):
                if not f.endswith('.root'): continue
                inputFile=targetDir+f
                fvec.push_back(inputFile)
        if fvec.empty():
            print("No files found for directory:", samples[sample], " SKIPPING processing")
            continue
        # print(fvec)
        if computeYes: sumwProc[sample] = ROOT.getClippedSumW(fvec)
        else:
            print('using Runs')
            RDF = ROOT.ROOT.RDataFrame
            runs = RDF('Runs', fvec)
            sumwProc[sample] = runs.Sum("genEventSumw")
        RDF = ROOT.ROOT.RDataFrame
        runs = RDF('Runs', fvec)
        runs = runs.Define("effEntries", "genEventSumw_*genEventSumw_/genEventSumw2_")
        effectiveEntries = runs.Sum("effEntries").GetValue()
        lumieq_w = effectiveEntries/ (samples[sample]["xsec"]*1000) 
        print(sample,samples[sample]["xsec"],lumieq_w)
        objList.append(sumwProc[sample])

    ROOT.RDF.RunGraphs(objList)

    sumwClipped = {} 

    for sample in samples:
        if sample not in sumwProc.keys() : continue
        sumwClipped[sample]=sumwProc[sample].GetValue()

    with open('genSumWOutput{}.py'.format(era), 'w') as fp:
        json.dump(sumwClipped, fp, indent=4)
