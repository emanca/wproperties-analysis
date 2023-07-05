from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt
from wremnants.correctionsTensor_helper import makeCorrectionsTensor
import hist
import narf

ROOT.gInterpreter.Declare('#include "interface/helHelper.h"')

class getHelWeights(module):
   
    def __init__(self, angFile, type):
        
        self.angFile=angFile
        self.type = type
        pass
      

    def run(self,d):

        self.d=d

        with h5py.File(self.angFile, "r") as f:
            results = narf.ioutils.pickle_load_h5py(f["angCoeffWZ"])
            hharmonics = results[f"hist_coeffs_{self.type}"]
        
        #convert to boost histogram with pyroot bindings
        
        hharmonicsConv = narf.hist_to_pyroot_boost(hharmonics)
        print(type(hharmonicsConv).__cpp_name__)
        # self.coeffhelper = ROOT.helHelper(hharmonicsConv)
        self.coeffhelper = ROOT.helHelper[type(hharmonicsConv)](ROOT.std.move(hharmonicsConv))
        
        self.d=self.d.Define("helVec",self.coeffhelper,["Vrap_preFSR_abs","Vpt_preFSR"])

        return self.d