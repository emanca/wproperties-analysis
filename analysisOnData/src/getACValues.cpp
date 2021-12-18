#include "interface/getACValues.hpp"

RNode getACValues::run(RNode d)
{

  auto getACValues = [this](float y, float pt, float charge) mutable {
        
    ROOT::VecOps::RVec<float> AngCoeff;

    int bin = _hA0_plus->FindBin(y, pt);
    if(charge>0){
      AngCoeff.push_back(_hA0_plus->GetBinContent(bin));
      AngCoeff.push_back(_hA1_plus->GetBinContent(bin));
      AngCoeff.push_back(_hA2_plus->GetBinContent(bin));
      AngCoeff.push_back(_hA3_plus->GetBinContent(bin));
      AngCoeff.push_back(_hA4_plus->GetBinContent(bin));
      AngCoeff.push_back(_hA5_plus->GetBinContent(bin));
      AngCoeff.push_back(_hA6_plus->GetBinContent(bin));
      AngCoeff.push_back(_hA7_plus->GetBinContent(bin));
      AngCoeff.push_back(_hAUL_plus->GetBinContent(bin));
    }
    else{
      AngCoeff.push_back(_hA0_minus->GetBinContent(bin));
      AngCoeff.push_back(_hA1_minus->GetBinContent(bin));
      AngCoeff.push_back(_hA2_minus->GetBinContent(bin));
      AngCoeff.push_back(_hA3_minus->GetBinContent(bin));
      AngCoeff.push_back(_hA4_minus->GetBinContent(bin));
      AngCoeff.push_back(_hA5_minus->GetBinContent(bin));
      AngCoeff.push_back(_hA6_minus->GetBinContent(bin));
      AngCoeff.push_back(_hA7_minus->GetBinContent(bin));
      AngCoeff.push_back(_hAUL_minus->GetBinContent(bin));
    }
    return AngCoeff;
  };

  auto getMapValue = [this](float y, float pt, float charge) mutable {
    int bin = _htotMap_plus->FindBin(y, pt);
    float totval=-99.;
    if(charge>0) totval = _htotMap_plus->GetBinContent(bin);
    else totval = _htotMap_minus->GetBinContent(bin);
    return totval;
  };
  
  auto getACValuesVars = [this](float y, float pt, float charge) mutable {
        
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> AngCoeff;
    int bin = _hA0_plus_vars[0]->FindBin(y, pt);
    for(long unsigned int v=0; v!=_syst_name.size(); v++) {
      ROOT::VecOps::RVec<float> AngCoeff_var;
      if(charge>0){
        AngCoeff_var.push_back(_hA0_plus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA1_plus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA2_plus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA3_plus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA4_plus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA5_plus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA6_plus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA7_plus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hAUL_plus_vars[v]->GetBinContent(bin));
      }
      else{
        AngCoeff_var.push_back(_hA0_minus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA1_minus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA2_minus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA3_minus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA4_minus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA5_minus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA6_minus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hA7_minus_vars[v]->GetBinContent(bin));
        AngCoeff_var.push_back(_hAUL_minus_vars[v]->GetBinContent(bin));
      }
      AngCoeff.emplace_back(AngCoeff_var);
    };
    
    return AngCoeff;
  };

  auto getMapValueVars = [this](float y, float pt, float charge) mutable {
    int bin = _htotMap_plus_vars[0]->FindBin(y, pt);
    ROOT::VecOps::RVec<float> totval;
    for(long unsigned int v=0; v!=_syst_name.size(); v++) {
      if(charge>0) totval.push_back(_htotMap_plus_vars[v]->GetBinContent(bin));
      else totval.push_back(_htotMap_minus_vars[v]->GetBinContent(bin));
    }
    return totval;
  };

  

  if(_syst_kind==""){
    auto d1 = d.Define("AngCoeffVec", getACValues, {"Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge"})
                  .Define("totMap", getMapValue, {"Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge"});
    return d1;
  }
  else{
    auto d1 = d.Define("AngCoeffVec", getACValuesVars, {"Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge"})
                  .Define("totMap", getMapValueVars, {"Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge"});
    return d1;
  }
  
}
