#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "interface/defineHarmonics.hpp"

RNode defineHarmonics::run(RNode d)
{

  // angular coefficients as defined in https://arxiv.org/pdf/1609.02536.pdf

  auto getHarmonicsVec = [](float costheta, float phiRaw) {
    float phi = phiRaw-TMath::Pi();
    float P0 = 1. / 2. * (1. - 3. * costheta * costheta);
    float P1 = 2. * costheta * sqrt(1. - costheta * costheta) * cos(phi);
    float P2 = 1. / 2. * (1. - costheta * costheta) * cos(2. * phi);
    float P3 = sqrt(1. - costheta * costheta) * cos(phi);
    float P4 = costheta;
    float P5 = (1. - costheta * costheta) * sin(2. * phi);
    float P6 = 2. * costheta * sqrt(1. - costheta * costheta) * sin(phi);
    float P7 = sqrt(1. - costheta * costheta) * sin(phi);
    float PUL = 1 + costheta * costheta;

    ROOT::VecOps::RVec<float> harms;

    harms.push_back(P0);
    harms.push_back(P1);
    harms.push_back(P2);
    harms.push_back(P3);
    harms.push_back(P4);
    harms.push_back(P5);
    harms.push_back(P6);
    harms.push_back(P7);
    harms.push_back(PUL);

    return harms;
  };

  auto multByWeight = [](float a, const ROOT::VecOps::RVec<float> &w) { return a * w; };
  auto multSqByWeight = [](float a, const ROOT::VecOps::RVec<float> &w) -> ROOT::VecOps::RVec<float> { return a * w * w; };

  auto d1 = d.Define("harmonicsVec", getHarmonicsVec, {"CStheta_preFSR", "CSphi_preFSR"})
                .Define("harmonicsVecWeighted", multByWeight, {"lumiweight", "harmonicsVec"})
                .Define("harmonicsVecSqWeighted", multSqByWeight, {"lumiweight", "harmonicsVec"});

  return d1;
}