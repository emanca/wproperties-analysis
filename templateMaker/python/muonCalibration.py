from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt
from wremnants import muon_calibration

class muonCalibration(module):
   
    def __init__(self,dataset,args,data_calibration_helper,mc_calibration_helper,data_jpsi_crctn_helper,mc_jpsi_crctn_helper):
        self.dataset = dataset
        self.args = args
        self.cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
        self.jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper
      
    def run(self,d):

        self.d = d
        self.d = muon_calibration.define_corrected_muons(self.d, self.cvh_helper, self.jpsi_helper, args, dataset, smearing_helper, bias_helper)
        

        return self.d
