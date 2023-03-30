from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt
from wremnants import muon_selections,muon_calibration

class muonSelection(module):
   
    def __init__(self,dataset,args,data_calibration_helper,mc_calibration_helper,data_jpsi_crctn_helper,mc_jpsi_crctn_helper):
        self.dataset = dataset
        self.args = args
        self.cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
        self.jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper
      
    def run(self,d):

        self.d = d.Filter("HLT_IsoMu24 ||  HLT_IsoTkMu24")
        self.d = muon_selections.veto_electrons(self.d)
        self.d = muon_selections.apply_met_filters(self.d)

        self.d = muon_calibration.define_corrected_muons(self.d, self.cvh_helper, self.jpsi_helper, self.args, self.dataset)

        self.d = muon_selections.select_veto_muons(self.d, nMuons=2)
        self.d = muon_selections.select_good_muons(self.d, nMuons=2, use_trackerMuons=False, use_isolation=True)

        self.d = muon_selections.define_trigger_muons(self.d)

        self.d = muon_selections.select_z_candidate(self.d, 25., 55.)

        self.d = muon_selections.select_standalone_muons(self.d, self.dataset, self.args.trackerMuons, "trigMuons")
        self.d = muon_selections.select_standalone_muons(self.d, self.dataset, self.args.trackerMuons, "nonTrigMuons")

        self.d = muon_selections.apply_triggermatching_muon(self.d, self.dataset, "trigMuons_eta0", "trigMuons_phi0")
        

        return self.d
