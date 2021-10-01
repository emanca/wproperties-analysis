#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "interface/getMassWeights.hpp"

ROOT::VecOps::RVec<float> vecMultiplication(const ROOT::VecOps::RVec<float> &v1, const ROOT::VecOps::RVec<float> &v2)
{
  ROOT::VecOps::RVec<float> products;

  products.reserve(v1.size() * v2.size());
  for (auto e1 : v1)
    for (auto e2 : v2)
      products.push_back(e1 * e2);

  return products;
}

RNode getMassWeights::run(RNode d)
{

  auto BreitWigner = [](float Q2, float M2, float G2) -> float
  {
    return 1. / ((Q2 - M2) * (Q2 - M2) + M2 * G2);
  };

  auto getBWVec = [&](float Q)
  {
    ROOT::VecOps::RVec<float> v;

    float w1 = BreitWigner(Q * Q, 80.369002 * 80.369002, 2.0476 * 2.0476) /
               // float w1 = BreitWigner(Q * Q, 80.409002 * 80.409002, 2.0476 * 2.0476) /      //10MeVsyst
               BreitWigner(Q * Q, 80.419002 * 80.419002, 2.0476 * 2.0476);
    // v.emplace_back(1.);
    float w2 = BreitWigner(Q * Q, 80.469002 * 80.469002, 2.0476 * 2.0476) /
               // float w2 = BreitWigner(Q * Q, 80.429002 * 80.429002, 2.0476 * 2.0476) /     //10MeVsyst
               BreitWigner(Q * Q, 80.419002 * 80.419002, 2.0476 * 2.0476);
    v.emplace_back(w2);
    v.emplace_back(w1);

    return v;
  };

  auto d1 = d.Define("BWmassWeights", getBWVec, {"Vmass_preFSR"}).Define("helBWweights", vecMultiplication, {"helWeights", "BWmassWeights"});
  return d1;
}