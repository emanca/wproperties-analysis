#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "interface/baseDefinitions.hpp"

RNode baseDefinitions::run(RNode d)
{
  auto defineNomweight = []() {
    ROOT::VecOps::RVec<float> One;
    One.emplace_back(1.);
    return One;
  };

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
  auto d1 = d.Define("Mupt_preFSR", "GenPart_pt[GenPart_preFSRMuonIdx]")
                .Define("Mueta_preFSR", "GenPart_eta[GenPart_preFSRMuonIdx]")
                .Define("Wrap_preFSR_abs", "TMath::Abs(Wrap_preFSR)")
                .Define("Nom", defineNomweight)
                .Define("LHEScaleWeightred", reduceVec, {"LHEScaleWeight"})
                // .Define("weight", "lumiweight*weightPt*weightY");
                .Define("weight", "lumiweight*weightPt");

  //at this point return the node in case of data
  return d1;
}