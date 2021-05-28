import ROOT
from termcolor import colored
import math 
from fitUtils import fitUtils
import os

charges = ["Wplus"]
for charge in charges:
    f = fitUtils(doSyst=True)
    f.fillProcessList()
    f.shapeFile()
    f.maskedChannels()
    f.fillHelGroup()
    f.setPreconditionVec()
    f.fillSumGroup()
    f.fillHelMetaGroup()
    f.makeDatacard()    
    
    text2hd5f = 'text2hdf5_npinput.py --allowNegativeExpectation --maskedChan=Wplus_preVFP_xsec --maskedChan=Wplus_postVFP_xsec {}.pkl --out {}.pkl.root'.format("WPlus","WPlus")
    print('executing', text2hd5f) 
    os.system(text2hd5f)
    # #--doRegularization --regularizationTau=1e4
    # #--allowNegativePOI
    # #--yieldProtectionCutoff 1000.
    combinetf = 'combinetf.py --fitverbose 9 -t0 --saveHists --computeHistErrors --seed 260292 --yieldProtectionCutoff 1000. --allowNegativePOI --doh5Output {}.pkl.hdf5 -o fit_{}.root'.format("WPlus", "WPlus")
    # combinetf = 'combinetf.py --fitverbose 9 -t-1 --seed 260292 --allowNegativePOI --doh5Output {}.pkl.hdf5 -o fit_{}.root'.format("WPlus", "WPlus")
    print('executing', combinetf)
    os.system(combinetf)
    assert(0)
