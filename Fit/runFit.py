import ROOT
from termcolor import colored
import math 
from fitUtils import fitUtils
import os
import argparse

parser = argparse.ArgumentParser("")
parser.add_argument('-t', '--t',type=int, default=0, help='run fit on data')

args = parser.parse_args()

toy = args.t
data=False
if toy==0: data=True
type = "asimov"
if toy>0:
    type = "toy"

charges = ["WPlus"]
for charge in charges:
    f = fitUtils(doSyst=True,channels =["{}_preVFP".format(charge),"{}_postVFP".format(charge)])
    f.fillProcessList()
    # f.shapeFile()
    # f.maskedChannels()
    f.fillHelGroup()
    f.setPreconditionVec()
    f.fillSumGroup()
    f.fillHelMetaGroup()
    f.makeDatacard()    
    text2hd5f = 'text2hdf5_npinput.py --allowNegativeExpectation --sparse --maskedChan={}_preVFP_xsec --maskedChan={}_postVFP_xsec {}.pkl --out {}.pkl.root'.format(charge,charge,charge,charge)
    print('executing', text2hd5f) 
    os.system(text2hd5f)
    #--yieldProtectionCutoff 100. --scan helXsecsA_y_5_qt_3_pmaskedexp --binByBinStat
    combinetf = 'combinetf.py --fitverbose 9 -t{} --seed 260292 --yieldProtectionCutoff 100. --allowNegativePOI --binByBinStat --doh5Output {}.pkl_sparse.hdf5 -o FitRes/fit_{}_{}_blockAi.root'.format(toy,charge, charge, "data" if data else type)
    print('executing', combinetf)
    os.system(combinetf)

# charges = ["WPlus"]
# for charge in charges:
#     f = fitUtils(doSyst=True,channels =["{}_postVFP".format(charge)])
#     f.fillProcessList()
#     f.shapeFile()
#     f.maskedChannels()
#     f.fillHelGroup()
#     f.setPreconditionVec()
#     f.fillSumGroup()
#     f.fillHelMetaGroup()
#     f.makeDatacard()    
#     text2hd5f = 'text2hdf5_npinput.py --allowNegativeExpectation --sparse --maskedChan={}_postVFP_xsec {}.pkl --out {}.pkl.root'.format(charge,charge,charge)
#     print('executing', text2hd5f) 
#     os.system(text2hd5f)
#     #--yieldProtectionCutoff 100. --scan helXsecsA_y_5_qt_3_pmaskedexp --binByBinStat
#     combinetf = 'combinetf.py --fitverbose 9 -t{} --binByBinStat --seed 260292  --yieldProtectionCutoff 100. --allowNegativePOI --doImpacts --doh5Output {}.pkl_sparse.hdf5 -o FitRes/fit_{}_{}_postVFP.root'.format(toy,charge, charge, "data" if data else type)
#     print('executing', combinetf)
#     os.system(combinetf)
#     assert(0)

# charges = ["WPlus"]
# for charge in charges:
#     f = fitUtils(doSyst=True,channels =["{}_preVFP".format(charge)])
#     f.fillProcessList()
#     # f.shapeFile()
#     # f.maskedChannels()
#     f.fillHelGroup()
#     f.setPreconditionVec()
#     f.fillSumGroup()
#     f.fillHelMetaGroup()
#     f.makeDatacard()    
#     text2hd5f = 'text2hdf5_npinput.py --allowNegativeExpectation --sparse --maskedChan={}_preVFP_xsec {}.pkl --out {}.pkl.root'.format(charge,charge,charge)
#     print('executing', text2hd5f) 
#     os.system(text2hd5f)
#     #--yieldProtectionCutoff 100. --scan helXsecsA_y_5_qt_3_pmaskedexp --binByBinStat
#     combinetf = 'combinetf.py --fitverbose 9 -t{} --binByBinStat --seed 260292  --yieldProtectionCutoff 100. --allowNegativePOI --doImpacts --doh5Output {}.pkl_sparse.hdf5 -o FitRes/fit_{}_{}_preVFP.root'.format(toy,charge, charge, "data" if data else type)
#     print('executing', combinetf)
#     os.system(combinetf)
#     assert(0)