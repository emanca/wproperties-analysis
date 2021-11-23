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




charges = ["Wplus","Wminus"]
Nrand = 1000 # number of template-toys produced
# charges = ["Wminus"]
for charge in charges:
    for n in range(1,Nrand+1) :
        print("processing toy template randomized:", n)
        # os.system('mkdir templateRand/rand'+n+'_'+charge)
        os.chdir('templateRand/rand'+str(n)+'_'+charge.replace('W',''))

        fmap = '../../../../analysisOnGen/genInput_{}.root'.format(charge)
        f = fitUtils(fmap, channel=charge+"_reco", doSyst=True)
        f.fillProcessList()
        f.shapeFile()
        f.fillHelGroup()
        f.fillSumGroup()
        f.fillHelMetaGroup()
        f.makeDatacard()    
        
        text2hd5f = 'text2hdf5.py --allowNegativeExpectation --doSystematics 1 --maskedChan={}_xsec {}.pkl'.format(f.channel,f.channel)
        print('executing', text2hd5f) 
        os.system(text2hd5f) 
        
        combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --binByBinStat --correlateXsecStat'.format(toy, f.channel, f.channel, CTFmodifier,cores)
        # combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {}'.format(toy, f.channel, f.channel, CTFmodifier,cores)
        print('executing', combinetf)
        os.system(combinetf)
        os.chdir('../../')
    
    assert(0) #no minus

