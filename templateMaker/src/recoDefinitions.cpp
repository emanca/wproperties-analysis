#include "interface/recoDefinitions.hpp"
#include "interface/functions.hpp"

RNode recoDefinitions::run(RNode d)
{
    //define all nominal quantities // true for data and MC


    auto d1 = d.Define("vetoMuons","Muon_pt > 10 && Muon_looseId && abs(Muon_eta) < 2.4 && abs(Muon_dxybs) < 0.05")
                    .Define("vetoElectrons"," Electron_pt > 10 && Electron_cutBased > 0 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz)< 0.2")
                    .Define("goodMuons","vetoMuons && Muon_pt > 25 && Muon_mediumId")
                    .Define("Mu1_eta", "Muon_eta[goodMuons][0]")
                    .Define("Mu1_phi", "Muon_phi[goodMuons][0]")
                    .Define("Mu1_charge", "float(Muon_charge[goodMuons][0])")
                    .Define("Mu1_relIso", "Muon_pfRelIso04_all[goodMuons][0]")
                    .Define("Mu1_dz", "Muon_dz[goodMuons][0]")
                    .Define("Mu1_pt", "Muon_pt[goodMuons][0]")
                    .Define("Mu1_sip3d", "Muon_sip3d[goodMuons][0]")
                    .Define("Mu1_dxy", "Muon_dxy[goodMuons][0]")
                    .Define("Mu1_hasTriggerMatch", hasTriggerMatch, {"Mu1_eta", "Mu1_phi", "goodTrigObjs_eta", "goodTrigObjs_phi"});

    return d1;
}
