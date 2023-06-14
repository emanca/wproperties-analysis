from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt
from wremnants import muon_calibration

ROOT.gInterpreter.Declare('#include "interface/helSystHelper.h"')

class muonCalibration(module):
   
    def __init__(self,dataset,jpsi_crctn_data_unc_helper):
        
        self.jpsi_crctn_data_unc_helper=jpsi_crctn_data_unc_helper
        self.dataset=dataset
    def run(self,d):

        self.d = d
        jpsi_unc_helper = self.jpsi_crctn_data_unc_helper

        self.d = muon_calibration.define_cvh_reco_muon_kinematics(self.d)
        #FIXME: make the smearing weights work without filtering on taus
        require_prompt = "tau" not in self.dataset.name # for muon GEN-matching
        
        reco_sel = "vetoMuonsPre"
        self.d = muon_calibration.define_genFiltered_recoMuonSel(self.d, reco_sel, require_prompt)
        reco_sel_GF = muon_calibration.getColName_genFiltered_recoMuonSel(reco_sel, require_prompt)
        self.d = muon_calibration.define_covMatFiltered_recoMuonSel(self.d, reco_sel_GF)
        self.d = muon_calibration.define_matched_gen_muons_covMat(self.d, reco_sel_GF)
        self.d = muon_calibration.define_matched_gen_muons_kinematics(self.d, reco_sel_GF)
        self.d = muon_calibration.calculate_matched_gen_muon_kinematics(self.d, reco_sel_GF)
        self.d = muon_calibration.define_matched_genSmeared_muon_kinematics(self.d, reco_sel_GF)

        reco_sel = "trigMuons"
        self.d = muon_calibration.define_matched_gen_muons_kinematics(self.d, reco_sel)
        self.d = muon_calibration.calculate_matched_gen_muon_kinematics(self.d, reco_sel)
        self.d = muon_calibration.define_matched_gen_muons_covMat(self.d, reco_sel)
        self.d = muon_calibration.define_matched_genSmeared_muon_kinematics(self.d, reco_sel)

        for var in ['Pt', 'Eta', 'Charge', 'Qop']:
            self.d = self.d.Define(f"{reco_sel}_{var.lower()}0_gen", f"{reco_sel}_gen{var.capitalize()}[0]")
            self.d = self.d.Define(f"{reco_sel_GF}_{var.lower()}0_gen", f"{reco_sel_GF}_gen{var.capitalize()}[0]")

            self.d = self.d.Define(f"{reco_sel}_{var.lower()}0_gen_smeared", f"{reco_sel}_genSmeared{var.capitalize()}[0]")
            self.d = self.d.Define(f"{reco_sel_GF}_{var.lower()}0_gen_smeared", f"{reco_sel_GF}_genSmeared{var.capitalize()}[0]")
        self.d = self.d.Define(f"{reco_sel_GF}_covMat0", f"{reco_sel_GF}_covMat[0]")

        self.d = self.d.DefinePerSample("bool_false", "false")
        self.d = self.d.Define("muonScaleSyst_responseWeights_tensor_gensmear", jpsi_unc_helper,
                            [
                                f"{reco_sel_GF}_genQop",
                                f"{reco_sel_GF}_genPhi",
                                f"{reco_sel_GF}_genEta",
                                f"{reco_sel_GF}_genSmearedQop",
                                f"{reco_sel_GF}_genSmearedPhi",
                                f"{reco_sel_GF}_genSmearedEta",
                                f"{reco_sel_GF}_genSmearedCharge",
                                f"{reco_sel_GF}_genSmearedPt",
                                f"{reco_sel_GF}_covMat",
                                "unity",
                                "bool_false"
                            ]
                        )
        
        print("jpsi tensor",self.d.GetColumnType("muonScaleSyst_responseWeights_tensor_gensmear"))
        # Eigen::TensorFixedSize<double,Eigen::Sizes<144,2>,0,long>
        self.helper = ROOT.helSystMuCalHelper[144]()
        self.d = self.d.Define("jpsiWeight_tensor_hel",self.helper,["helWeightTensor", "muonScaleSyst_responseWeights_tensor_gensmear"])
        print("jpsi tensor hel",self.d.GetColumnType("jpsiWeight_tensor_hel"))
        return self.d
