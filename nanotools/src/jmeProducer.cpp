#include "jmeProducer.hpp"
RNode jmeProducer::run(RNode d) {
  auto getcorrMET = [this](const RVec<float> &Jet_pt, const RVec<float> &Jet_eta, const RVec<float> &Jet_phi, const RVec<float> &Jet_mass, const RVec<\
float> &Jet_rawFactor, const RVec<float> &Jet_area, const RVec<float> &Jet_muonSubtrFactor, const RVec<float> &Jet_neEmEF, const RVec<float> &Jet_chEm\
EF, const float fixedGridRhoFastjetAll, const RVec<float> &GenJet_pt, const RVec<float> &GenJet_eta, const RVec<float> &GenJet_phi, const RVec<float> \
&GenJet_mass, const float &RawMET_phi, const float &RawMET_pt, const float &MET_MetUnclustEnUpDeltaX, const float &MET_MetUnclustEnUpDeltaY, const RVe\
c<float> &CorrT1METJet_rawPt, const RVec<float> &CorrT1METJet_eta, const RVec<float> &CorrT1METJet_phi, const RVec<float> &CorrT1METJet_area, const RV\
ec<float> &CorrT1METJet_muonSubtrFactor) {
    const ROOT::VecOps::RVec<int> empty1;//empty jet id..only needed for  2017
    const ROOT::VecOps::RVec<float> empty2;//empty jet id..only needed for  2017
    const ROOT::VecOps::RVec<float> empty3;//empty jet id..only needed for  2017
    auto res = myJetVarCalc->produce(Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_rawFactor, Jet_area,
				     Jet_muonSubtrFactor, Jet_neEmEF, Jet_chEmEF, empty1,
				     fixedGridRhoFastjetAll, 1,
				     GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass,
				     RawMET_phi, RawMET_pt,
				     MET_MetUnclustEnUpDeltaX, MET_MetUnclustEnUpDeltaY,
				     CorrT1METJet_rawPt, CorrT1METJet_eta, CorrT1METJet_phi,
				     CorrT1METJet_area, CorrT1METJet_muonSubtrFactor, empty2, empty3);
    return res;
  };
  
  auto getcorrMETdata = [this](const RVec<float> &Jet_pt, const RVec<float> &Jet_eta, const RVec<float> &Jet_phi, const RVec<float> &Jet_mass, const RVec<float> &Jet_rawFactor, const RVec<float> &Jet_area, const RVec<float> &Jet_muonSubtrFactor, const RVec<float> &Jet_neEmEF, const RVec<float> &Jet_chEmEF, const float fixedGridRhoFastjetAll, const float &RawMET_phi, const float &RawMET_pt, const float &MET_MetUnclustEnUpDeltaX, const float &MET_MetUnclustEnUpDeltaY, const RVec<float> &CorrT1METJet_rawPt, const RVec<float> &CorrT1METJet_eta, const RVec<float> &CorrT1METJet_phi, const RVec<float> &CorrT1METJet_area, const RVec<float> &CorrT1METJet_muonSubtrFactor) {
    const ROOT::VecOps::RVec<int> empty1;//empty jet id..only needed for  2017
    const ROOT::VecOps::RVec<float> empty2;//empty jet id..only needed for  2017
    const ROOT::VecOps::RVec<float> empty3;//empty jet id..only needed for  2017
    const ROOT::VecOps::RVec<float> empty4;//for genJet_pt
    const ROOT::VecOps::RVec<float> empty5;//for genJet_eta
    const ROOT::VecOps::RVec<float> empty6;//for genJet_phi
    const ROOT::VecOps::RVec<float> empty7;//for genJet_mass
    auto res = myJetVarCalc->produce(Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_rawFactor, Jet_area,
				     Jet_muonSubtrFactor, Jet_neEmEF, Jet_chEmEF, empty1,
				     fixedGridRhoFastjetAll, 0,
				     empty4, empty5, empty6, empty7,
				     RawMET_phi, RawMET_pt,
				     MET_MetUnclustEnUpDeltaX, MET_MetUnclustEnUpDeltaY,
				     CorrT1METJet_rawPt, CorrT1METJet_eta, CorrT1METJet_phi,
				     CorrT1METJet_area, CorrT1METJet_muonSubtrFactor, empty2, empty3);
    return res;
  };
  
  auto df = d;
  if(isData_) {
    df = df.Define("metVars", getcorrMETdata,
		   {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_rawFactor",
		    "Jet_area", "Jet_muonSubtrFactor", "Jet_neEmEF", "Jet_chEmEF",
		    "fixedGridRhoFastjetAll",
		    "RawMET_phi", "RawMET_pt",
		    "MET_MetUnclustEnUpDeltaX", "MET_MetUnclustEnUpDeltaY",
		    "CorrT1METJet_rawPt", "CorrT1METJet_eta", "CorrT1METJet_phi",
		    "CorrT1METJet_area", "CorrT1METJet_muonSubtrFactor"});
    df = df.Define("MET_T1comp_pt", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.pt(0)); }, {"metVars"})
      .Define("MET_T1comp_phi", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.phi(0)); }, {"metVars"});
  } else {
    df = df.Define("metVars", getcorrMET,
		   {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_rawFactor",
		    "Jet_area", "Jet_muonSubtrFactor", "Jet_neEmEF", "Jet_chEmEF",
		    "fixedGridRhoFastjetAll", "GenJet_pt", "GenJet_eta", "GenJet_phi",
		    "GenJet_mass", "RawMET_phi", "RawMET_pt",
		    "MET_MetUnclustEnUpDeltaX", "MET_MetUnclustEnUpDeltaY",
		    "CorrT1METJet_rawPt", "CorrT1METJet_eta", "CorrT1METJet_phi",
		    "CorrT1METJet_area", "CorrT1METJet_muonSubtrFactor"});
    //Now define the relevant MET branches
    df = df.Define("MET_T1comp_pt", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.pt(0)); }, {"metVars"})
      .Define("MET_T1comp_pt_jerUp", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.pt(1)); }, {"metVars"})
      .Define("MET_T1comp_pt_jerDown", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.pt(2)); }, {"metVars"})
      .Define("MET_T1comp_pt_jesTotalUp", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.pt(3)); }, {"metVars"})
      .Define("MET_T1comp_pt_jesTotalDown", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.pt(4)); }, {"metVars"})
      .Define("MET_T1comp_pt_unclustEnUp", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.pt(5)); }, {"metVars"})
      .Define("MET_T1comp_pt_unclustEnDown", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.pt(6)); }, {"metVars"})
      .Define("MET_T1comp_phi", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.phi(0)); }, {"metVars"})
      .Define("MET_T1comp_phi_jerUp", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.phi(1)); }, {"metVars"})
      .Define("MET_T1comp_phi_jerDown", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.phi(2)); }, {"metVars"})
      .Define("MET_T1comp_phi_jesTotalUp", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.phi(3)); }, {"metVars"})
      .Define("MET_T1comp_phi_jesTotalDown", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.phi(4)); }, {"metVars"})
      .Define("MET_T1comp_phi_unclustEnUp", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.phi(5)); }, {"metVars"})
      .Define("MET_T1comp_phi_unclustEnDown", [](rdfhelpers::ModifiedMET vmet) { return float(vmet.phi(6)); }, {"metVars"});
  }
  return df;
}

