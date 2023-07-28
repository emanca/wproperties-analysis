from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt

ROOT.gInterpreter.Declare('#include "interface/helSystHelper.h"')

class bootstrap(module):
   
    def __init__(self,dataset):
        self.dataset=dataset
        pass
      
    def run(self,d):

        self.d = d

        if "Wplus" in self.dataset.name:
            dataset_number = 260292
        elif "Wminus" in self.dataset.name:
            dataset_number = 161096
        self.d = self.d.Define("dataset_number",f"{dataset_number}").Define("rndPoisson", "bootstrapping(run, luminosityBlock, event,dataset_number)")\
                    .Define("rndPoisson_tensor", f"Eigen::TensorFixedSize<double, Eigen::Sizes<400>>(wrem::vec_to_tensor_t<double, 400>(rndPoisson))")
        self.helper = ROOT.helSystPrefHelper[400]()
        self.d = self.d.Define("rndPoisson_tensor_hel",self.helper,["helWeightTensor", "rndPoisson_tensor"])
        self.helper = ROOT.helMassPoissHelper[2,400]()
        self.d = self.d.Define("rndPoisson_tensor_hel_mass",self.helper,["massWeight_tensor_hel", "rndPoisson_tensor"])

        return self.d
