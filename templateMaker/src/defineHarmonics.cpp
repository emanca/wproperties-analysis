#include "defineHarmonics.hpp"
#include "functions.hpp"

RNode defineHarmonics::run(RNode d)
{

  // angular coefficients as defined in https://arxiv.org/pdf/1609.02536.pdf

  auto getHarmonicsVec = [](float costheta, float phi)
  {
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

  auto vecMultiplication = [](const ROOT::VecOps::RVec<float> &v1, const ROOT::VecOps::RVec<float> &v2) -> ROOT::VecOps::RVec<float>
  {
    ROOT::VecOps::RVec<float> products;

    products.reserve(v1.size() * v2.size());
    for (auto e1 : v1)
      for (auto e2 : v2)
        products.push_back(e1 * e2);

    return products;
  };

  auto d1 = d.Define("harmonicsVec", getHarmonicsVec, {"CStheta_preFSR", "CSphi_preFSR"})
                // .Define("harmonicsVec_BWmassWeights", vecMultiplication, {"harmonicsVec", "BWmassWeights"})
                .Define("harmonicsVec_LHEPdfWeight", vecMultiplication, {"harmonicsVec", "LHEPdfWeight"})
                .Define("harmonicsVec_LHEScaleWeight", vecMultiplication, {"harmonicsVec", "LHEScaleWeight"})
                .Define("DeltaPhi", deltaPhi, {"CSphi_preFSR","Vphi_preFSR"});
  return d1;
}