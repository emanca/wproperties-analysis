from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt

ROOT.gInterpreter.Declare('#include "interface/helSystHelper.h"')

class getSFVariations(module):
   
    def __init__(self,isWlike,helper_stat, helper_syst,helperPref_stat, helperPref_syst):
        
        self.isWlike = isWlike
        self.helper_stat = helper_stat
        self.helper_syst = helper_syst
        self.helperPref_stat = helperPref_stat
        self.helperPref_syst = helperPref_syst
      
    def run(self,d):

        self.d = d
        
        if self.isWlike:
            # muon_columns_stat = ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_charge0", "nonTrigMuons_pt0", "nonTrigMuons_eta0", "nonTrigMuons_charge0"]
            # muon_columns_syst = ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_SApt0", "trigMuons_SAeta0", "trigMuons_charge0",
            # "nonTrigMuons_pt0", "nonTrigMuons_eta0", "nonTrigMuons_SApt0", "nonTrigMuons_SAeta0", "nonTrigMuons_charge0"]
            muon_columns_stat = ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_charge0"]
            muon_columns_syst = ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_SApt0", "trigMuons_SAeta0", "trigMuons_charge0","passIso"]
        else:
            muon_columns_stat = ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_charge0"]
            muon_columns_syst = ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso"]

        # SF stat variations
        for key,helper in self.helper_stat.items():
            if "iso" in key and not self.isWlike:
                self.d = self.d.Define(f"effStatTnP_{key}_tensor", helper, [*muon_columns_stat, "passIso", "nominal_weight"])        
            else:
                self.d = self.d.Define(f"effStatTnP_{key}_tensor", helper, [*muon_columns_stat, "nominal_weight"])
            
            print(f"effStatTnP_{key}_tensor",self.d.GetColumnType(f"effStatTnP_{key}_tensor"))
    
        # SF syst variations
        self.d = self.d.Define("effSystTnP_weight", self.helper_syst, [*muon_columns_syst, "nominal_weight"])
        print("effSystTnP_weight",self.d.GetColumnType("effSystTnP_weight"))

        # redefine without weight for helicity multiplication
        # SF stat variations
        for key,helper in self.helper_stat.items():
            if "iso" in key and not self.isWlike:
                self.d = self.d.Define(f"effStatTnP_{key}_tensor_unweighted", helper, [*muon_columns_stat, "passIso", "unity"])        
            else:
                self.d = self.d.Define(f"effStatTnP_{key}_tensor_unweighted", helper, [*muon_columns_stat, "unity"])
        
        # SF syst variations
        self.d = self.d.Define("effSystTnP_weight_unweighted", self.helper_syst, [*muon_columns_syst, "unity"])
        
        # prefiring variations
        self.d = self.d.Define("muonL1PrefireStat_tensor", self.helperPref_stat, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId", "nominal_weight"])
        self.d = self.d.Define("muonL1PrefireStat_tensor_unweighted", self.helperPref_stat, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId", "unity"])
    
        self.d = self.d.Define("muonL1PrefireSyst_tensor", self.helperPref_syst, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId", "nominal_weight"])
        self.d = self.d.Define("muonL1PrefireSyst_tensor_unweighted", self.helperPref_syst, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId", "unity"])
    
        self.d = self.d.Define("ecalL1Prefire_tensor", f"wrem::twoPointScaling(nominal_weight/L1PreFiringWeight_ECAL_Nom, L1PreFiringWeight_ECAL_Dn, L1PreFiringWeight_ECAL_Up)")
        self.d = self.d.Define("ecalL1Prefire_tensor_unweighted", f"wrem::twoPointScaling(unity/L1PreFiringWeight_ECAL_Nom, L1PreFiringWeight_ECAL_Dn, L1PreFiringWeight_ECAL_Up)")

        print("muonL1PrefireStat_tensor",self.d.GetColumnType("muonL1PrefireStat_tensor"))
        print("muonL1PrefireSyst_tensor",self.d.GetColumnType("muonL1PrefireSyst_tensor"))
        print("ecalL1Prefire_tensor",self.d.GetColumnType("ecalL1Prefire_tensor"))

        '''
        effStatTnP_sf_reco_tensor Eigen::TensorFixedSize<double,Eigen::Sizes<1,4,2>,0,long>
        effStatTnP_sf_tracking_tensor Eigen::TensorFixedSize<double,Eigen::Sizes<1,3,2>,0,long>
        effStatTnP_sf_idip_tensor Eigen::TensorFixedSize<double,Eigen::Sizes<1,4,2>,0,long>
        effStatTnP_sf_trigger_tensor Eigen::TensorFixedSize<double,Eigen::Sizes<1,4,2>,0,long>
        effStatTnP_sf_iso_tensor Eigen::TensorFixedSize<double,Eigen::Sizes<1,4,1>,0,long>
        muonL1PrefireStat_tensor Eigen::TensorFixedSize<double,Eigen::Sizes<12,2>,0,long>
        muonL1PrefireSyst_tensor Eigen::TensorFixedSize<double,Eigen::Sizes<2>,0,long>
        ecalL1Prefire_tensor Eigen::TensorFixedSize<double,Eigen::Sizes<2>,0,long>
        '''
        nPtEig = [4,3,4,4,4]
        nCharges = [2,2,2,2,1]
        for ikey,key in enumerate(self.helper_stat):
            # multiply by helicity
            self.helper = ROOT.helSystSFHelper[nPtEig[ikey],nCharges[ikey]]()
            self.d = self.d.Define(f"effStatTnP_{key}_tensor_hel",self.helper,["helWeightTensor", f"effStatTnP_{key}_tensor_unweighted"])
            print(f"effStatTnP_{key}_tensor_hel",self.d.GetColumnType(f"effStatTnP_{key}_tensor_hel"))
        self.helper = ROOT.helSystSFSystVarHelper[5,49]()
        self.d = self.d.Define("effSystTnP_tensor_hel",self.helper,["helWeightTensor", "effSystTnP_weight_unweighted"])
        self.helper = ROOT.helSystSFSystVarHelper[12,2]()
        self.d = self.d.Define("muonL1PrefireStat_tensor_hel",self.helper,["helWeightTensor", "muonL1PrefireStat_tensor_unweighted"])
        self.helper = ROOT.helSystPrefHelper[2]()
        self.d = self.d.Define("muonL1PrefireSyst_tensor_hel",self.helper,["helWeightTensor", "muonL1PrefireSyst_tensor_unweighted"])
        self.d = self.d.Define("ecalL1Prefire_tensor_hel",self.helper,["helWeightTensor", "ecalL1Prefire_tensor_unweighted"])
        
        return self.d
