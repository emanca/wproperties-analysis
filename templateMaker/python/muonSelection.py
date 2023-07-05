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

        if self.args.Wlike:
            nMuons = 2
            flag = "trigMuons"
        else:
            nMuons = 1
            flag = "goodMuons"

        self.d = d.Filter("HLT_IsoMu24 ||  HLT_IsoTkMu24")
        self.d = muon_selections.veto_electrons(self.d)
        self.d = muon_selections.apply_met_filters(self.d)

        self.d = muon_calibration.define_corrected_muons(self.d, self.cvh_helper, self.jpsi_helper, self.args, self.dataset)

        self.d = muon_selections.select_veto_muons(self.d, nMuons=nMuons)
        self.d = muon_selections.select_good_muons(self.d, nMuons=nMuons, use_trackerMuons=self.args.trackerMuons, use_isolation=True)

        if self.args.Wlike:
            self.d = muon_selections.define_trigger_muons(self.d)
        else:
            self.d = muon_calibration.define_corrected_reco_muon_kinematics(self.d)
        
        if self.args.Wlike:
            self.d = muon_selections.select_z_candidate(self.d, 25., 55.)

        self.d = muon_selections.select_standalone_muons(self.d, self.dataset, self.args.trackerMuons, f"{flag}")
        
        if self.args.Wlike:
            self.d = muon_selections.select_standalone_muons(self.d, self.dataset, self.args.trackerMuons, "nonTrigMuons")

        self.d = muon_selections.apply_triggermatching_muon(self.d, self.dataset, f"{flag}_eta0", f"{flag}_phi0")
        
        ### recoil stuff
        self.d = self.d.Alias("MET_corr_rec_pt", "MET_pt")
        self.d = self.d.Alias("MET_corr_rec_phi", "MET_phi")

        self.d = self.d.Define("transverseMass", "wrem::mt_2(goodMuons_pt0, goodMuons_phi0, MET_corr_rec_pt, MET_corr_rec_phi)")

        self.d = self.d.Define("deltaPhiMuonMet", "std::abs(wrem::deltaPhi(goodMuons_phi0,MET_corr_rec_phi))")
        self.d = self.d.Define("passMT", "transverseMass >= 40.0")
        self.d = self.d.Define("goodMuons_pfRelIso04_all0", "Muon_pfRelIso04_all[goodMuons][0]")
        self.d = self.d.Define("passIso", "goodMuons_pfRelIso04_all0 < 0.15")

        return self.d
