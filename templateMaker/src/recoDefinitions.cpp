#include "recoDefinitions.hpp"
#include "functions.hpp"

RNode recoDefinitions::run(RNode d)
{
    //define all nominal quantities // true for data and MC

    auto d1 = d.Define("Mu1_eta", "Muon_eta[0]")
                  .Define("Mu1_phi", "Muon_phi[0]")
                  .Define("Mu1_charge", "float(Muon_charge[0])")
                  .Define("Mu1_relIso", "Muon_pfRelIso04_all[0]")
                  .Define("Mu1_dz", "Muon_dz[0]")
                  .Define("Mu1_pt", "Muon_pt[0]")
                  .Define("Mu1_sip3d", "Muon_sip3d[0]")
                  .Define("Mu1_dxy", "Muon_dxy[0]")
                  .Define("Mu1_hasTriggerMatch", hasTriggerMatch, {"Mu1_eta", "Mu1_phi", "goodTrigObjs_eta", "goodTrigObjs_phi"});

    return d1;
}
