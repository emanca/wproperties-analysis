import ROOT
from termcolor import colored
import math 
from fitUtilsZ import fitUtilsZ
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



f = fitUtilsZ(doSyst=True,channels =['Wlike_minus','Wlike_plus'])
f.fillProcessList()
f.shapeFile()
f.maskedChannels()
f.fillHelGroup()
f.setPreconditionVec()
f.fillSumGroup()
f.fillHelMetaGroup()
f.makeDatacard()    
text2hd5f = 'text2hdf5_npinput.py --allowNegativeExpectation --sparse --maskedChan=Wlike_minus_xsec --maskedChan=Wlike_plus_xsec Wlike.pkl --out Wlike.pkl.root'
print('executing', text2hd5f) 
os.system(text2hd5f)
#--yieldProtectionCutoff 100. --scan helXsecsA_y_5_qt_3_pmaskedexp --binByBinStat
combinetf = 'combinetf.py --fitverbose 9 -t{} --seed 260292 --yieldProtectionCutoff 100. --allowNegativePOI --doh5Output --binByBinStat --saveHists --doImpacts Wlike.pkl_sparse.hdf5 -o FitRes/fit_Wlike_{}.root'.format(toy, "data" if data else type)
print('executing', combinetf)
os.system(combinetf)
