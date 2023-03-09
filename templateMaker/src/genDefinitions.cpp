#include "interface/genDefinitions.hpp"

RNode genDefinitions::run(RNode d)
{
    auto reduceVec = [](ROOT::VecOps::RVec<float> LHE) {
        ROOT::VecOps::RVec<float> red;
        red.emplace_back(LHE[0]);
        red.emplace_back(LHE[1]);
        red.emplace_back(LHE[3]);
        red.emplace_back(LHE[5]);
        red.emplace_back(LHE[7]);
        red.emplace_back(LHE[8]);
        return red;
    };

    //define all nominal quantities // true for data and MC
    auto d1 = d.Define("GenCharge", "float(GenPart_pdgId[GenPart_preFSRMuonIdx]/TMath::Abs(GenPart_pdgId[GenPart_preFSRMuonIdx]))")
                  .Define("Mupt_preFSR", "GenPart_pt[GenPart_preFSRMuonIdx]")
                  .Define("Mueta_preFSR", "GenPart_eta[GenPart_preFSRMuonIdx]")
                  .Define("Wrap_preFSR_abs", "TMath::Abs(Wrap_preFSR)");
                  //.Define("LHEScaleWeightred", reduceVec, {"LHEScaleWeight"});
    return d1;
}