import ROOT
from termcolor import colored
import math 
from fitUtils import fitUtils
import os
import argparse
import systToapply
import copy

parser = argparse.ArgumentParser("")
parser.add_argument('-i', '--impact',      type=int, default=False, help="make impact plots (doImpact)")
parser.add_argument('-pf', '--postfit',      type=int, default=False, help="save postfit plots (saveHists)")
parser.add_argument('-r', '--regularize',   type=int, default=False, help="apply regularization (doRegularization))")
parser.add_argument('-s', '--tau',          type=str, default='1e4', help="set strenght of regularization (regularizationTau)")
parser.add_argument('-t', '--toy',          type=str, default='-1', help="number of toy, -1=asimov")
parser.add_argument('-c', '--cores',          type=str, default='-1', help="number of cores, -1=all")
# parser.add_argument('-pdf', '--PDFprefit',          type=int, default=False, help="run likelihood scan to obtain prefit PDF impacts from mw shift")
args = parser.parse_args()
impact = args.impact 
postfit = args.postfit
regularize = args.regularize
tau = args.tau
toy = args.toy
cores = args.cores
# PDF = args.PDFprefit


CTFmodifier = ''
if impact :     CTFmodifier+= ' --doImpacts '
if postfit :    CTFmodifier+= ' --saveHists --computeHistErrors'
if regularize : CTFmodifier+= ' --doRegularization '
if regularize : CTFmodifier+= ' --regularizationTau='+tau

excludeList = [
    # 'LHEPdfWeightHess44',
    # 'LHEPdfWeightHess45',
    # 'LHEPdfWeightHess46',
    # 'LHEPdfWeightHess47',
    # 'LHEPdfWeightHess40',
    # 'LHEPdfWeightHess41',
    # 'LHEPdfWeightHess42',
    # 'LHEPdfWeightHess43',
    # 'LHEPdfWeightHess26',
    # 'LHEPdfWeightHess27',
    # 'LHEPdfWeightHess24',
    # 'LHEPdfWeightHess50',
    # 'LHEPdfWeightHess54',
    # 'LHEPdfWeightHess60',
    # 'LHEPdfWeightHess53',
    # 'LHEPdfWeightHess32',
    # 'LHEPdfWeightHess48',
    # 'LHEPdfWeightHess49',
    # 'LHEPdfWeightHess52',
    # 'LHEPdfWeightHess25',
    
    # 'LHEPdfWeightHess17',
    # 'LHEPdfWeightHess18',
    # 'LHEPdfWeightHess19',
    # 'LHEPdfWeightHess30',
    # 'LHEPdfWeightHess31',
    # 'LHEPdfWeightHess36',
    # 'LHEPdfWeightHess37',
    # 'LHEPdfWeightHess34',
    # 'LHEPdfWeightHess35',
    # 'LHEPdfWeightHess38',
    # 'LHEPdfWeightHess51',
    # 'LHEPdfWeightHess39',
    # 'LHEPdfWeightHess58',
    # 'LHEPdfWeightHess55',
    # 'LHEPdfWeightHess59',
    # 'LHEPdfWeightHess1',
    # 'LHEPdfWeightHess2',
    # 'alphaS',
    # 'LHEPdfWeightHess4',
    # 'LHEPdfWeightHess5',
    # 'LHEPdfWeightHess6',
    # 'LHEPdfWeightHess7',
    # 'LHEPdfWeightHess8',
    # 'LHEPdfWeightHess9',
    # 'LHEPdfWeightHess33',
    # 'LHEPdfWeightHess56',
    # 'LHEPdfWeightHess29',
    # 'LHEPdfWeightHess28',
    # 'LHEPdfWeightHess57',
    # 'LHEPdfWeightHess21',
    # 'LHEPdfWeightHess20',
    # 'LHEPdfWeightHess23',
    # 'LHEPdfWeightHess25',
    # 'LHEPdfWeightHess24',
    # 'LHEPdfWeightHess27',
    # 'LHEPdfWeightHess26',
    # 'LHEPdfWeightHess43',
    # 'LHEPdfWeightHess42',
    # 'LHEPdfWeightHess41',
    # 'LHEPdfWeightHess40',
    # 'LHEPdfWeightHess47',
    # 'LHEPdfWeightHess46',
    # 'LHEPdfWeightHess45',
    # 'LHEPdfWeightHess44',
    # 'LHEPdfWeightHess52',
    # 'LHEPdfWeightHess49',
    # 'LHEPdfWeightHess48',
    # 'LHEPdfWeightHess32',
    # 'LHEPdfWeightHess53',
    # 'LHEPdfWeightHess60',
    # 'LHEPdfWeightHess54',
    # 'LHEPdfWeightHess50',
    
    
    # 'LHEPdfWeightHess59', #skipped, doesn't work! 
    
    'LHEPdfWeightHess22', #skipped, doesn't work! 
    # 'LHEPdfWeightHess3', #skipped, doesn't work!  
    # 'LHEPdfWeightHess16' #skipped, doesn't work! 
]


charges = ["Wplus","Wminus"]
# charges = ["Wminus"]
pdfList = set(['LHEPdfWeightHess{}'.format(i+1) for i in range(60)]+['alphaS'])
# pdfList = set(['LHEPdfWeightHess{}'.format(i+1) for i in range(2)])
print(pdfList)
for charge in charges:
    for pdf in pdfList :
        if pdf not in excludeList : continue
        print("processing pdf:", pdf)
 
        
        systKind = 'LHEPdfWeight'
        if 'alpha' in pdf :
            systKind='alphaS'
        systDict = copy.deepcopy(systToapply.systematicsDict)
        del systDict['LHEPdfWeight']
        del systDict['alphaS']
        systDict[systKind] = {
                    "vars":[pdf],
                    "procs": ["Signal", "DYJets", "Fake", "WtoTau", "LowAcc"],
                    "type": "shape",
                    "weight" : 1.
                    }
        
        # systDict = {
        #     "Nominal": {},
        #     "mass" : {
        #             "vars":["mass"],
        #             "procs": ["Signal", "LowAcc", "Fake"],
        #             "type": "shapeNoConstraint",
        #             "weight" : 1.
        #             },
        #     systKind : {
        #             "vars":[pdf],
        #             "procs": ["Signal", "DYJets", "Fake", "WtoTau", "LowAcc"],
        #             "type": "shape",
        #             "weight" : 1.
        #             },
        # }
        
        fmap = '../../analysisOnGen/genInput_{}.root'.format(charge)
        f = fitUtils(fmap, channel=charge+"_reco", doSyst=True, systDict =systDict )
        f.fillProcessList()
        f.shapeFile()
        f.fillHelGroup()
        f.fillSumGroup()
        f.fillHelMetaGroup()
        f.makeDatacard()    
        
        text2hd5f = 'text2hdf5.py --allowNegativeExpectation --doSystematics 1 --maskedChan={}_xsec {}.pkl'.format(f.channel,f.channel)
        print('executing', text2hd5f) 
        os.system(text2hd5f) 
        
        scanConf = '--scanRangeUsePrefit 1 --scanPoint 2 --scanRange 1 --scan {}'.format(pdf)
        combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} {} --POIMode none'.format(toy, f.channel, f.channel, CTFmodifier,cores,scanConf )
        print('executing', combinetf)
        os.system(combinetf)
        
        os.system('mkdir PDFscan_extra/'+pdf+'_'+charge)
        os.system('mv '+f.channel+'.pkl PDFscan_extra/'+pdf+'_'+charge)
        os.system('mv '+f.channel+'.pkl.hdf5 PDFscan_extra/'+pdf+'_'+charge)
        os.system('mv '+f.channel+'_xsec.root PDFscan_extra/'+pdf+'_'+charge)
        os.system('mv fit_'+f.channel+'.root PDFscan_extra/'+pdf+'_'+charge)
    
    assert(0) #no minus

