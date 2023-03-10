#include "interface/zVetoMuons.hpp"
#include "interface/functions.hpp"
#include "TLorentzVector.h"

RNode zVetoMuons::run(RNode d) {

  auto d1 = d.Define("vetoMuonsPre", "Muon_looseId && abs(Muon_dxybs) < 0.05 && Muon_charge != -99")
                .Define("vetoMuons", "vetoMuonsPre && Muon_pt > 10. && abs(Muon_eta) < 2.4")
                .Define("goodMuons", "vetoMuons && Muon_mediumId && Muon_isGlobal && Muon_pfRelIso04_all < 0.15");

                return d1;
}
