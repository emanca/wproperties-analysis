import os
import ROOT
import copy
import sys
import argparse
import math

parser = argparse.ArgumentParser("")
parser.add_argument('-i','--input', type=str, default='../analysisOnData/test/',help="name of the input direcory")

args = parser.parse_args()
inputDir = args.input
HADD = False

sampleList = ['Data','WToMu', 'WToTau', 'DY', 'DiBoson','Top','QCD']
outDict = {}
signList = ['Plus','Minus']

sampleDict = {
            "WToMu"       : ['WToMu_plots.root',       'prefit_Signal',         "WToMu",        ['WToMu_plots.root']],         
            "DY"          : ['DYJets_plots.root',      'prefit_Signal',          "DY",          ['DYJetsToLL_M-*']],
            "WToTau"      : ['WToTau_plots.root',      'prefit_Signal',          "WToTau",     ['WToTau_plots.root']],   
            "TW"          : ['TW_plots.root',          'prefit_Signal',           "Top",       ['ST_tW_*']],
            "TTbar"       : ['TTJets_plots.root',      'prefit_Signal',           "Top",       ['TTJets*']],
            "ST"          : ['SingleTop_plots.root',   'prefit_Signal',           "Top",       ['ST_t-*', 'ST_s-*']],
            "DiBoson"     : ['Diboson_plots.root',     'prefit_Signal',           "di-boson",  ['WW_TuneCUETP8M1_13TeV-pythia8_plots.root','WZ_TuneCUETP8M1_13TeV-pythia8_plots.root', 'ZZ_TuneCUETP8M1_13TeV-pythia8_plots.root']],    
            "QCD"         : ['FakeFromData_plots.root', 'prefit_fakes',           "QCD",        ['FakeFromData_plots.root']],     
            "Data"        : ['Data_plots.root',        'prefit_Signal',            "Data",      ['SingleMuonData_plots.root']]
}

if HADD :
    cmdList = []
    for sample, info in sampleDict.items() :
        haddStr = ' '
        for inn in info[3] :
            haddStr = haddStr+' '+inputDir+'/'+inn
        cmdList.append('hadd -f '+info[0]+ haddStr )
        print('hadd -f '+info[0]+ haddStr)
    # cmdList.append('cp  '+inDir+'/SingleMuonData_plots.root Data_plots.root')
    # cmdList.append('cp  '+inDir+'/FakeFromData_plots.root FakeFromData_plots.root')
    # cmdList.append('cp  '+inDir+'/WToMu_plots.root WToMu_plots.root')
    # cmdList.append('cp  '+inDir+'/WToTau_plots.root WToTau_plots.root')
    # cmdList.append('hadd -f DYJets_plots.root '+inDir+'/DYJetsToLL_M-*')
    # cmdList.append('hadd -f Top_plots.root '+inDir+'/TTJets* ' +inDir+'/ST* '+inDir+'/ST_tW_* ')
    # cmdList.append('hadd -f Diboson_plots.root '+inDir+'/WW_TuneCUETP8M1_13TeV-pythia8_plots.root '+inDir+'/WZ_TuneCUETP8M1_13TeV-pythia8_plots.root '+inDir+'/ZZ_TuneCUETP8M1_13TeV-pythia8_plots.root')
    
    for i in cmdList :
        os.system(i)





for sample, info in sampleDict.items() :
    # inFile = ROOT.TFile.Open(inputDir+'/'+info[0])
    inFile = ROOT.TFile.Open(info[0])
    h2 = inFile.Get(info[1]+'/Nominal/Mu1_eta')
    h2.GetYaxis().SetRangeUser(0,2)
    outDict[sample+'Plus'] = h2.Integral()
    h2.GetYaxis().SetRangeUser(-2,0)
    outDict[sample+'Minus'] = h2.Integral()

topPlus = 0
topMinus = 0
for sample, info in sampleDict.items() :
    if 'Top' in info[2] :
        topPlus+=outDict[sample+'Plus']
        topMinus+=outDict[sample+'Minus'] 
outDict['Top'+'Plus'] = topPlus
outDict['Top'+'Minus'] = topMinus

print('-------------------------------------------')
print('-------------------------------------------')
for sample in sampleList :
    for s in signList :
        print(sample, s, ":", outDict[sample+s], '    (', outDict[sample+s]/outDict['WToMu'+s]*100,'percent of signal MC)')
print('-------------------------------------------')
print('-------------------------------------------')