#include "interface/getACValues.hpp"

RNode getACValues::run(RNode d)
{

  auto getACValues = [this](float y, float pt, float charge) mutable {
        
    ROOT::VecOps::RVec<float> AngCoeff;

    int bin = _hA0_plus->FindBin(y, pt);

    if(charge>0){
      AngCoeff.push_back(_hA0_plus->GetBinContent(bin)/2.);
      AngCoeff.push_back(_hA1_plus->GetBinContent(bin)/(2. * TMath::Sqrt(2)));
      AngCoeff.push_back(_hA2_plus->GetBinContent(bin)/4.);
      AngCoeff.push_back(_hA3_plus->GetBinContent(bin)/(4. * TMath::Sqrt(2)));
      AngCoeff.push_back(_hA4_plus->GetBinContent(bin)/2.);
      AngCoeff.push_back(_hA5_plus->GetBinContent(bin)/2.);
      AngCoeff.push_back(_hA6_plus->GetBinContent(bin)/(2. * TMath::Sqrt(2)));
      AngCoeff.push_back(_hA7_plus->GetBinContent(bin)/(4. * TMath::Sqrt(2)));
    }
    else{
      AngCoeff.push_back(_hA0_minus->GetBinContent(bin)/2.);
      AngCoeff.push_back(_hA1_minus->GetBinContent(bin)/(2. * TMath::Sqrt(2)));
      AngCoeff.push_back(_hA2_minus->GetBinContent(bin)/4.);
      AngCoeff.push_back(_hA3_minus->GetBinContent(bin)/(4. * TMath::Sqrt(2)));
      AngCoeff.push_back(_hA4_minus->GetBinContent(bin)/2.);
      AngCoeff.push_back(_hA5_minus->GetBinContent(bin)/2.);
      AngCoeff.push_back(_hA6_minus->GetBinContent(bin)/(2. * TMath::Sqrt(2)));
      AngCoeff.push_back(_hA7_minus->GetBinContent(bin)/(4. * TMath::Sqrt(2)));
    }
    return AngCoeff;
  };

  auto getMapValue = [this](float y, float pt, float charge) mutable {
    int bin = _htotMap_plus->FindBin(y, pt);
    float totval=-99.;
    if(charge>0)
      totval = _htotMap_plus->GetBinContent(bin);
    else
      totval = _htotMap_minus->GetBinContent(bin);
    return totval;
  };

  auto d1 = d.Define("AngCoeffVec", getACValues, {"Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge"})
                .Define("totMap", getMapValue, {"Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge"});

  return d1;
}
