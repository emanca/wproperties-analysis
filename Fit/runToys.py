import ROOT
from termcolor import colored
import math 
from fitUtils import fitUtils
import os
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count

parser = argparse.ArgumentParser("")
parser.add_argument('-t', '--t',type=int, default=1, help='number of toys to run')

args = parser.parse_args()

toys = args.t
charge = "WPlus"
seeds = np.random.randint(100000, size=(toys))
print(seeds)
def runSingleToy(toy):
    # combinetf = 'combinetf.py --nThreads=3 --bootstrapData -t25 --seed {}  --yieldProtectionCutoff 100. --allowNegativePOI  {}.pkl_sparse.hdf5 -o FitRes/fit_{}_toy{}.root'.format(seeds[toy], charge,charge, toy)
    combinetf = 'combinetf.py --nThreads=3 -t40 --bootstrapData --seed {}  --yieldProtectionCutoff 100. --allowNegativePOI  {}.pkl_sparse.hdf5 -o FitRes/fit_{}_toy{}.root'.format(seeds[toy], charge,charge, toy)
    print('executing', combinetf)
    os.system(combinetf)

f = fitUtils(doSyst=True,channels =["{}_preVFP".format(charge),"{}_postVFP".format(charge)])
f.fillProcessList()
f.shapeFile()
f.maskedChannels()
f.fillHelGroup()
f.setPreconditionVec()
f.fillSumGroup()
f.fillHelMetaGroup()
f.makeDatacard()    
text2hd5f = 'text2hdf5_npinput.py --allowNegativeExpectation --sparse --maskedChan={}_preVFP_xsec --maskedChan={}_postVFP_xsec {}.pkl --out {}.pkl.root'.format(charge,charge,charge, charge)
print('executing', text2hd5f) 
os.system(text2hd5f)

toy = [i for i in range(toys)]
pool=Pool()
pool.map(runSingleToy, toy)
