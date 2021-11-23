import ROOT
from termcolor import colored
import math 
from fitUtils_simult2charge import fitUtils
import os
import argparse

parser = argparse.ArgumentParser("")
parser.add_argument('-i', '--impact',      type=int, default=False, help="make impact plots (doImpact)")
parser.add_argument('-pf', '--postfit',      type=int, default=False, help="save postfit plots (saveHists)")
parser.add_argument('-r', '--regularize',   type=int, default=False, help="apply regularization (doRegularization))")
parser.add_argument('-s', '--tau',          type=str, default='1e4', help="set strenght of regularization (regularizationTau)")
parser.add_argument('-t', '--toy',          type=str, default='-1', help="number of toy, -1=asimov")
parser.add_argument('-c', '--cores',          type=str, default='-1', help="number of cores, -1=all")
args = parser.parse_args()
impact = args.impact 
postfit = args.postfit
regularize = args.regularize
tau = args.tau
toy = args.toy
cores = args.cores


CTFmodifier = ''
if impact :     CTFmodifier+= ' --doImpacts '
if postfit :    CTFmodifier+= ' --saveHists --computeHistErrors'
if regularize : CTFmodifier+= ' --doRegularization '
if regularize : CTFmodifier+= ' --regularizationTau='+tau

charges = ["Wplus","Wminus"]
nameOut = "Wall"
fmap = {
    "Wplus" : '../../analysisOnGen/genInput_Wplus.root',
    "Wminus" : '../../analysisOnGen/genInput_Wminus.root',
}

channelList = ["Wplus","Wminus"]

f = fitUtils(fmap, channel=channelList, doSyst=True)
f.fillProcessList()
f.shapeFile()
f.fillHelGroup()
f.fillSumGroup()
f.fillHelMetaGroup()
f.makeDatacard()    

text2hd5f = 'text2hdf5.py --allowNegativeExpectation --doSystematics 1 --maskedChan={}_xsec --maskedChan={}_xsec {}.pkl --sparse'.format(charges[0],charges[1], nameOut)
print('executing', text2hd5f) 
os.system(text2hd5f) 


combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {}'.format(toy, nameOut, nameOut, CTFmodifier,cores)

print('executing', combinetf)
os.system(combinetf)
