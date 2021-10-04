#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"

using RNode = ROOT::RDF::RNode;

ROOT::RDF::RResultPtr<double> getClippedSumW(std::vector<std::string> fileNames)
{
  auto df = RNode(ROOT::RDataFrame("Events", fileNames));

  auto clipGenWeight = [](float Gen_weight) {
    double sign = Gen_weight / std::abs(Gen_weight);
    //return sign;                                                                                                                                                                                   
    return sign;
  };
  auto d1 = df.Define("Generator_weight_clipped", clipGenWeight, {"genWeight"}).Sum<double>("Generator_weight_clipped");
  return d1;
}
