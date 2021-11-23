import ROOT
from termcolor import colored
import math 
from fitUtils import fitUtils
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
# if postfit :    CTFmodifier+= ' --saveHists'
if regularize : CTFmodifier+= ' --doRegularization '
if regularize : CTFmodifier+= ' --regularizationTau='+tau

charges = ["Wplus","Wminus"]
# charges = ["Wminus"]
for charge in charges:
    fmap = '../../analysisOnGen/genInput_{}.root'.format(charge)
    # fmap = '../../analysisOnGen/genInput_{}_2GeV_qt80_y28.root'.format(charge)
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
    
    #OLD
    # combinetf = 'combinetf.py --allowNegativePOI --binByBinStat --correlateXsecStat --doRegularization --regularizationTau=1e1 -t-1 {}.pkl.hdf5 -o fit_{}.root'.format(
    # combinetf = 'combinetf.py --allowNegativePOI --binByBinStat -t -1 {}.pkl.hdf5 -o fit_{}.root --doImpacts --saveHists'.format(
    # combinetf = 'combinetf.py --allowNegativePOI --binByBinStat --doRegularization --regularizationTau=1e4 --doImpacts --saveHists -t -1 {}.pkl.hdf5 -o fit_{}.root'.format(
        # f.channel, f.channel)
    
    #good one (with and without BBB)
    # combinetf = 'combinetf.py --allowNegativePOI --binByBinStat -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --correlateXsecStat'.format(toy, f.channel, f.channel, CTFmodifier,cores) # --fitverbose 9
    # combinetf = 'combinetf.py --allowNegativePOI --binByBinStat -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {}'.format(toy, f.channel, f.channel, CTFmodifier,cores) # --fitverbose 9
    # combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --seed 7031993'.format(toy, f.channel, f.channel, CTFmodifier,cores) #working toy
    combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {}'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    # combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --POIMode none'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    


    #other config
    # combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --scanRangeUsePrefit 1 --scanPoint 3 --scanRange 1 --scan helXsecsUL_y_1_qt_1_mu'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    # combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --scanRangeUsePrefit 1 --scanPoint 3 --scanRange 1 --scan mass --fitverbose 9'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    # combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --scanPoint 3 --scanRange 1 --scan mass --fitverbose 9'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    # combinetf = 'combinetf.py --allowNegativePOI --binByBinStat -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --correlateXsecStat --yieldProtectionCutoff 1'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    # combinetf = 'combinetf.py --binByBinStat -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --correlateXsecStat'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    # combinetf = 'combinetf.py --binByBinStat -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --correlateXsecStat'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    # combinetf = 'combinetf.py --allowNegativePOI -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {}'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    # combinetf = 'combinetf.py -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {}'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    # combinetf = 'combinetf.py --allowNegativePOI --binByBinStat -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {}'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    # combinetf = 'combinetf.py --allowNegativePOI --binByBinStat -t {} {}.pkl.hdf5 -o fit_{}.root {} --nThreads {} --POIMode none'.format(toy, f.channel, f.channel, CTFmodifier,cores)
    print('executing', combinetf)
    os.system(combinetf)
    # assert(0)
