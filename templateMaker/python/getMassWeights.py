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
        self.d = self.d.Define("massDown","MEParamWeight[5]")\
                        .Define("massUp","MEParamWeight[15]")\
                        .Define("massVec","RVec<double>{massDown,massUp}")\
                        .Define("massWeight_tensor", f"Eigen::TensorFixedSize<double, Eigen::Sizes<1,2>>(wrem::vec_to_tensor_t<double, {2}>(massVec).reshape(std::array<Eigen::Index, 2>{{1, 2}}))")\
                        .Define("massWeight_tensor_wnom", "auto res = massWeight_tensor; res = nominal_weight*res; return res;")
        #.Define("massWeight_tensor", f"wrem::vec_to_tensor_t<double, {nweights}>(MEParamWeight)")\
        self.helper = ROOT.helMassHelper[2]()
        self.d = self.d.Define("massWeight_tensor_hel",self.helper,["helWeightTensor", "massWeight_tensor"])
        # print(self.d.GetColumnType("massWeight_tensor_hel"))
        # print(self.d.GetColumnType("massWeight_tensor"))
        # print(self.d.GetColumnType("massVec"))
        # print(self.d.GetColumnType("MEParamWeight_size"))

        return self.d
