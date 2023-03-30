#include "interface/zSelection.hpp"
#include "interface/functions.hpp"
#include "TLorentzVector.h"

RNode zSelection::run(RNode d) {
  auto getdimuonP4 = [](float mu1_pt, float mu1_eta, float mu1_phi, float mu2_pt, float mu2_eta, float mu2_phi) {
    TLorentzVector mu1P;
    mu1P.SetPtEtaPhiM(mu1_pt, mu1_eta, mu1_phi, 0.1057);
    TLorentzVector mu2P;
    mu2P.SetPtEtaPhiM(mu2_pt, mu2_eta, mu2_phi, 0.1057);
    //float mass = (mu1P + mu2P).M();
    return (mu1P + mu2P);
  };

  auto hasTriggerMatch = [](const float mueta, const float muphi, const RVec<float> &TrigObj_eta, const RVec<float> &TrigObj_phi) {
   bool muhasTrigm = false;
   for (unsigned int jtrig = 0; jtrig < TrigObj_eta.size(); ++jtrig) {
     if (deltaR(mueta, muphi, TrigObj_eta[jtrig], TrigObj_phi[jtrig]) < 0.3) {
       muhasTrigm = true;
       break;
     }//if
   }
   
   return muhasTrigm;
  };

  //.Define("Mu2_hasTriggerMatch", "Muon_hasTriggerMatch[prefIdx[1]]")

  auto d1 = d.Define("vetoMuonsPre", "Muon_looseId && abs(Muon_dxybs) < 0.05 && Muon_charge != -99")
                .Define("vetoMuons", "vetoMuonsPre && Muon_pt > 10. && abs(Muon_eta) < 2.4")
                .Define("goodMuons", "vetoMuons && Muon_mediumId && Muon_isGlobal && Muon_pfRelIso04_all < 0.15 && Muon_highPurity")
                .Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz)< 0.2");

  auto d2 = d1.Define("Mu1_eta", "Muon_eta[goodMuons && Muon_charge>0][0]")
                .Define("Mu1_phi", "Muon_phi[goodMuons && Muon_charge>0][0]")
                .Define("Mu1_charge", "Muon_charge[goodMuons && Muon_charge>0][0]")
                .Define("Mu1_relIso", "Muon_pfRelIso04_all[goodMuons && Muon_charge>0][0]")
                .Define("Mu1_pt", "Muon_pt[goodMuons && Muon_charge>0][0]")
                .Define("Mu1_hasTriggerMatch", hasTriggerMatch, {"Mu1_eta", "Mu1_phi", "goodTrigObjs_eta", "goodTrigObjs_phi"})
                .Define("Mu2_eta", "Muon_eta[goodMuons && Muon_charge<0][0]")
                .Define("Mu2_phi", "Muon_phi[goodMuons && Muon_charge<0][0]")
                .Define("Mu2_charge", "Muon_charge[goodMuons && Muon_charge<0][0]")
                .Define("Mu2_relIso", "Muon_pfRelIso04_all[goodMuons && Muon_charge<0][0]")
                .Define("Mu2_pt", "Muon_pt[goodMuons && Muon_charge<0][0]")
                .Define("dimuonP4", getdimuonP4, {"Mu1_pt", "Mu1_eta", "Mu1_phi", "Mu2_pt", "Mu2_eta", "Mu2_phi"})
                .Define("dimuonMass", [this](TLorentzVector p)
                        { return float(p.M()); },
                        {"dimuonP4"})
                .Define("dimuonPt", [this](TLorentzVector p)
                        { return float(p.Pt()); },
                        {"dimuonP4"})
                .Define("dimuonY", [this](TLorentzVector p)
                        { return float(p.Rapidity()); },
                        {"dimuonP4"})
                .Define("nPV", [this](int p)
                        { return float(1. * p); },
                        {"PV_npvsGood"})
                .Define("uno", []()
                        { float uno = 1.; return uno; });
  // auto d1 = d.Define("wlikeMuons", "Muon_mediumId && abs(Muon_dxybs) < 0.2  && Muon_pt > 30. && abs(Muon_eta) < 2.4 && Muon_pfRelIso04_all<0.15")
  //               .Define("metMuons", "Muon_mediumId && abs(Muon_dxybs) < 0.2  && Muon_pt > 10. && Muon_pfRelIso04_all<0.5")

  //                   auto d2 = d1.Define("Mu1_eta", "Muon_eta[wlikeMuons && Muon_charge>0][0]")
  //                                 .Define("Mu1_phi", "Muon_phi[wlikeMuons && Muon_charge>0][0]")
  //                                 .Define("Mu1_charge", "Muon_charge[wlikeMuons && Muon_charge>0][0]")
  //                                 .Define("Mu1_relIso", "Muon_pfRelIso04_all[wlikeMuons && Muon_charge>0][0]")
  //                                 .Define("Mu1_pt", "Muon_pt[wlikeMuons && Muon_charge>0][0]")
  //                                 .Define("Mu1_hasTriggerMatch", hasTriggerMatch, {"Mu1_eta", "Mu1_phi", "goodTrigObjs_eta", "goodTrigObjs_phi"})
  //                                 .Define("Mu2_eta", "Muon_eta[metMuons && Muon_charge<0][0]")
  //                                 .Define("Mu2_phi", "Muon_phi[metMuons && Muon_charge<0][0]")
  //                                 .Define("Mu2_charge", "Muon_charge[metMuons && Muon_charge<0][0]")
  //                                 .Define("Mu2_relIso", "Muon_pfRelIso04_all[metMuons && Muon_charge<0][0]")
  //                                 .Define("Mu2_pt", "Muon_pt[metMuons && Muon_charge<0][0]")
  //                                 .Define("dimuonP4", getdimuonP4, {"Mu1_pt", "Mu1_eta", "Mu1_phi", "Mu2_pt", "Mu2_eta", "Mu2_phi"})
  //                                 .Define("dimuonMass", [this](TLorentzVector p)
  //                                         { return float(p.M()); },
  //                                         {"dimuonP4"})
  //                                 .Define("dimuonPt", [this](TLorentzVector p)
  //                                         { return float(p.Pt()); },
  //                                         {"dimuonP4"})
  //                                 .Define("dimuonY", [this](TLorentzVector p)
  //                                         { return float(p.Rapidity()); },
  //                                         {"dimuonP4"})
  //                                 .Define("nPV", [this](int p)
  //                                         { return float(1. * p); },
  //                                         {"PV_npvsGood"})
  //                                 .Define("uno", []()
  //                                         { float uno = 1.; return uno; });
  return d2;

}
