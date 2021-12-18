#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TMath.h"
#include "interface/PDFWeightsHelper.hpp"
#include "interface/TH1weightsHelper.hpp"
#include "interface/Replica2Hessian.hpp"


RNode Replica2Hessian::run(RNode d)
{

  auto newPDFweights = [this](ROOT::VecOps::RVec<float> replicas, float lhenom, unsigned long long ev) {
    std::vector<float> raw_weights;

    for (unsigned int i = 0; i < nPdfWeights_; i++)
    {
      raw_weights.push_back(replicas[i] * lhenom);
    }
    std::vector<float> raw_Hessweights;
    raw_Hessweights.resize(nPdfEigWeights_);

    pdfweightshelper_.DoMC2Hessian(lhenom, raw_weights.data(), raw_Hessweights.data());

    for (unsigned int i = 0; i < nPdfEigWeights_; i++)
    {
      raw_Hessweights[i] = raw_Hessweights[i] / lhenom;
    }
    return ROOT::VecOps::RVec<float>(raw_Hessweights);
  };
  
  auto getalphaSvarVec = [this](float wup, float wdown) {
    ROOT::VecOps::RVec<float> alpS;
    alpS.emplace_back(wup);
    alpS.emplace_back(wdown);
    return alpS;
  };

  auto d1 = d.Define("nLHEPdfWeightHess", [this]() { return nPdfEigWeights_ - 1; })
             .Define("LHEPdfWeightHess", newPDFweights, {"LHEPdfWeight", "LHEWeight_originalXWGTUP", "event"})
             .Define("alphaSUp", [this](ROOT::VecOps::RVec<float> pdfs, unsigned int npdfs) {return pdfs[npdfs-1];}, {"LHEPdfWeight", "nLHEPdfWeight"})
             .Define("alphaSDown", [this](ROOT::VecOps::RVec<float> pdfs, unsigned int npdfs) {return pdfs[npdfs-2];}, {"LHEPdfWeight", "nLHEPdfWeight"})
             .Define("alphaSVars", getalphaSvarVec, {"alphaSUp", "alphaSDown"});

  return d1;
}
/*
std::vector<ROOT::RDF::RResultPtr<TH1D>> Replica2Hessian::getTH1()
{
  return _h1List;
}
std::vector<ROOT::RDF::RResultPtr<TH2D>> Replica2Hessian::getTH2()
{
  return _h2List;
}
std::vector<ROOT::RDF::RResultPtr<TH3D>> Replica2Hessian::getTH3()
{
  return _h3List;
}

std::vector<ROOT::RDF::RResultPtr<std::vector<TH1D>>> Replica2Hessian::getGroupTH1()
{
  return _h1Group;
}
std::vector<ROOT::RDF::RResultPtr<std::vector<TH2D>>> Replica2Hessian::getGroupTH2()
{
  return _h2Group;
}
std::vector<ROOT::RDF::RResultPtr<std::vector<TH3D>>> Replica2Hessian::getGroupTH3()
{
  return _h3Group;
}

void Replica2Hessian::reset()
{

  _h1List.clear();
  _h2List.clear();
  _h3List.clear();

  _h1Group.clear();
  _h2Group.clear();
  _h3Group.clear();
}
*/
