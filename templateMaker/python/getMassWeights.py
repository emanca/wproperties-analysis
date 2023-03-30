from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt

ROOT.gInterpreter.Declare('#include "interface/helSystHelper.h"')

class getMassWeights(module):
   
    def __init__(self, isW):
        self.isW = isW
        pass
      
    def run(self,d):
        nweights = 23 if not self.isW else 21 #if it's not W it's Z
        # from -100 to 100 MeV with 10 MeV increment
        self.d = d
        self.d = self.d.Define("massWeight_tensor", f"wrem::vec_to_tensor_t<double, {nweights}>(MEParamWeight)")\
                       .Define("massWeight_tensor_wnom", "auto res = massWeight_tensor; res = nominal_weight*res; return res;")

        self.helper = ROOT.helSystHelper[nweights]()
        self.d = self.d.Define("massWeight_tensor_hel",self.helper,["helWeightTensor", "massWeight_tensor"])

        return self.d
