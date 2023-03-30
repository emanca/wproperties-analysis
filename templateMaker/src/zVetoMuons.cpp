#include "interface/zVetoMuons.hpp"
#include "interface/functions.hpp"
#include "TLorentzVector.h"

RNode zVetoMuons::run(RNode d) {

  auto d1 = d.Define("vetoMuonsPre", "Muon_looseId && abs(Muon_dxybs) < 0.05 && Muon_charge != -99")
                .Define("vetoMuons", "vetoMuonsPre && Muon_correctedPt > 10. && abs(Muon_correctedEta) < 2.4")
                .Define("goodMuons", "vetoMuons && Muon_mediumId && Muon_isGlobal && Muon_pfRelIso04_all < 0.15 && Muon_highPurity")
                .Define("vetoElectrons", "Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz)< 0.2");

  return d1;
}
