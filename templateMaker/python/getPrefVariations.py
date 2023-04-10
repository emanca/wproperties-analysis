from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt

ROOT.gInterpreter.Declare('#include "interface/helSystHelper.h"')

class getPrefVariations(module):
   
    def __init__(self,helper_stat, helper_syst):
        
        self.helper_stat = helper_stat
        self.helper_syst = helper_syst
      
    def run(self,d):

        self.d = d
        
        self.d = self.d.Define("muonL1PrefireStat_tensor", self.helper_stat, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId", "nominal_weight"])
    
        self.d = self.d.Define("muonL1PrefireSyst_tensor", self.helper_syst, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId", "nominal_weight"])

        self.d = self.d.Define("ecalL1Prefire_tensor", f"wrem::twoPointScaling(nominal_weight/L1PreFiringWeight_ECAL_Nom, L1PreFiringWeight_ECAL_Dn, L1PreFiringWeight_ECAL_Up)")

        print(self.d.GetColumnType("muonL1PrefireStat_tensor"))
        print(self.d.GetColumnType("muonL1PrefireSyst_tensor"))
        print(self.d.GetColumnType("ecalL1Prefire_tensor"))


        return self.d
