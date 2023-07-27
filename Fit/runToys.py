import ROOT
from termcolor import colored
import math 
from fitUtils import fitUtils
import os
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count

def generate_seeds(seed, num_seeds):
    seed_seq = np.random.SeedSequence(seed)
    seeds = seed_seq.generate_state(num_seeds)
    return seeds

def runSingleToy(toy):
    # combinetf = 'combinetf.py --nThreads=3 --bootstrapData -t25 --seed {}  --yieldProtectionCutoff 100. --allowNegativePOI  {}.pkl_sparse.hdf5 -o FitRes/fit_{}_toy{}.root'.format(seeds[toy], charge,charge, toy)
    combinetf = 'combinetf.py --nThreads=3 -t-1 --seed {}  --yieldProtectionCutoff 100. --allowNegativePOI  ../templateMaker/bstrp_inputs/Wplus_bstrp_{}.hdf5 --saveHists --doJacobian -o FitRes/fit_Wplus_btsrp{}_asimov.root'.format(seeds[toy],toy,toy)
    print('executing', combinetf)
    os.system(combinetf)

parser = argparse.ArgumentParser("")
parser.add_argument('-t', '--t',type=int, default=1, help='number of toys to run')

args = parser.parse_args()
toys = args.t

seed = 260292
n_bstrp = 100

# seeds = generate_seeds(seed, n_bstrp)
seeds = np.random.randint(100000, size=(toys))
print(seeds)


toy = [i for i in range(n_bstrp)]
pool=Pool()
pool.map(runSingleToy, toy)
