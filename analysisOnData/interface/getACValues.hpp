#ifndef GETACVALUES_H
#define GETACVALUES_H

#include "module.hpp"

class getACValues : public Module
{

private:
  TFile *_AChistos_plus;
  TH2D *_hA0_plus;
  TH2D *_hA1_plus;
  TH2D *_hA2_plus;
  TH2D *_hA3_plus;
  TH2D *_hA4_plus;
  TH2D *_hA5_plus;
  TH2D *_hA6_plus;
  TH2D *_hA7_plus;
  TH2D *_hAUL_plus;
  TH2D *_htotMap_plus;

  TFile *_AChistos_minus;
  TH2D *_hA0_minus;
  TH2D *_hA1_minus;
  TH2D *_hA2_minus;
  TH2D *_hA3_minus;
  TH2D *_hA4_minus;
  TH2D *_hA5_minus;
  TH2D *_hA6_minus;
  TH2D *_hA7_minus;
  TH2D *_hAUL_minus;
  TH2D *_htotMap_minus;
  
  std::string _syst_kind = "";
  ROOT::VecOps::RVec<std::string> _syst_name;
  ROOT::VecOps::RVec<TH2D *>_hA0_plus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA1_plus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA2_plus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA3_plus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA4_plus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA5_plus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA6_plus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA7_plus_vars;
  ROOT::VecOps::RVec<TH2D *>_hAUL_plus_vars;
  ROOT::VecOps::RVec<TH2D *>_htotMap_plus_vars;

  ROOT::VecOps::RVec<TH2D *>_hA0_minus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA1_minus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA2_minus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA3_minus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA4_minus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA5_minus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA6_minus_vars;
  ROOT::VecOps::RVec<TH2D *>_hA7_minus_vars;
  ROOT::VecOps::RVec<TH2D *>_hAUL_minus_vars;
  ROOT::VecOps::RVec<TH2D *>_htotMap_minus_vars;

public:
  getACValues(TFile *AChistos_plus, TFile *AChistos_minus)
  {
    _AChistos_plus = AChistos_plus;
    _hA0_plus = (TH2D *)_AChistos_plus->Get("angularCoefficients/harmonicsA0_nom_nom");
    _hA1_plus = (TH2D *)_AChistos_plus->Get("angularCoefficients/harmonicsA1_nom_nom");
    _hA2_plus = (TH2D *)_AChistos_plus->Get("angularCoefficients/harmonicsA2_nom_nom");
    _hA3_plus = (TH2D *)_AChistos_plus->Get("angularCoefficients/harmonicsA3_nom_nom");
    _hA4_plus = (TH2D *)_AChistos_plus->Get("angularCoefficients/harmonicsA4_nom_nom");
    _hA5_plus = (TH2D *)_AChistos_plus->Get("angularCoefficients/harmonicsA5_nom_nom");
    _hA6_plus = (TH2D *)_AChistos_plus->Get("angularCoefficients/harmonicsA6_nom_nom");
    _hA7_plus = (TH2D *)_AChistos_plus->Get("angularCoefficients/harmonicsA7_nom_nom");
    _hAUL_plus = (TH2D *)_AChistos_plus->Get("angularCoefficients/harmonicsAUL_nom_nom");
    _htotMap_plus = (TH2D *)_AChistos_plus->Get("accMaps/mapTot");

    _AChistos_minus = AChistos_minus;
    _hA0_minus = (TH2D *)_AChistos_minus->Get("angularCoefficients/harmonicsA0_nom_nom");
    _hA1_minus = (TH2D *)_AChistos_minus->Get("angularCoefficients/harmonicsA1_nom_nom");
    _hA2_minus = (TH2D *)_AChistos_minus->Get("angularCoefficients/harmonicsA2_nom_nom");
    _hA3_minus = (TH2D *)_AChistos_minus->Get("angularCoefficients/harmonicsA3_nom_nom");
    _hA4_minus = (TH2D *)_AChistos_minus->Get("angularCoefficients/harmonicsA4_nom_nom");
    _hA5_minus = (TH2D *)_AChistos_minus->Get("angularCoefficients/harmonicsA5_nom_nom");
    _hA6_minus = (TH2D *)_AChistos_minus->Get("angularCoefficients/harmonicsA6_nom_nom");
    _hA7_minus = (TH2D *)_AChistos_minus->Get("angularCoefficients/harmonicsA7_nom_nom");
    _hAUL_minus = (TH2D *)_AChistos_minus->Get("angularCoefficients/harmonicsAUL_nom_nom");
    _htotMap_minus = (TH2D *)_AChistos_minus->Get("accMaps/mapTot");
  };
  
  getACValues(TFile *AChistos_plus, TFile *AChistos_minus, std::string syst_kind, std::vector<std::string> syst_name)
  {
    _syst_kind = syst_kind;
    _syst_name = syst_name;
    _AChistos_plus = AChistos_plus;
    _AChistos_minus = AChistos_minus;
    
    for(auto const& sName : _syst_name) {
      
      
      _hA0_plus_vars.emplace_back((TH2D *)_AChistos_plus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA0"+sName+""+sName)));
      _hA1_plus_vars.emplace_back((TH2D *)_AChistos_plus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA1"+sName+""+sName)));
      _hA2_plus_vars.emplace_back((TH2D *)_AChistos_plus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA2"+sName+""+sName)));
      _hA3_plus_vars.emplace_back((TH2D *)_AChistos_plus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA3"+sName+""+sName)));
      _hA4_plus_vars.emplace_back((TH2D *)_AChistos_plus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA4"+sName+""+sName)));
      _hA5_plus_vars.emplace_back((TH2D *)_AChistos_plus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA5"+sName+""+sName)));
      _hA6_plus_vars.emplace_back((TH2D *)_AChistos_plus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA6"+sName+""+sName)));
      _hA7_plus_vars.emplace_back((TH2D *)_AChistos_plus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA7"+sName+""+sName)));
      _hAUL_plus_vars.emplace_back((TH2D *)_AChistos_plus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsAUL"+sName+""+sName)));
      _htotMap_plus_vars.emplace_back((TH2D *)_AChistos_plus->Get(TString("angularCoefficients_"+_syst_kind+"/mapTot"+sName)));

      _hA0_minus_vars.emplace_back((TH2D *)_AChistos_minus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA0"+sName+""+sName)));
      _hA1_minus_vars.emplace_back((TH2D *)_AChistos_minus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA1"+sName+""+sName)));
      _hA2_minus_vars.emplace_back((TH2D *)_AChistos_minus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA2"+sName+""+sName)));
      _hA3_minus_vars.emplace_back((TH2D *)_AChistos_minus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA3"+sName+""+sName)));
      _hA4_minus_vars.emplace_back((TH2D *)_AChistos_minus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA4"+sName+""+sName)));
      _hA5_minus_vars.emplace_back((TH2D *)_AChistos_minus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA5"+sName+""+sName)));
      _hA6_minus_vars.emplace_back((TH2D *)_AChistos_minus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA6"+sName+""+sName)));
      _hA7_minus_vars.emplace_back((TH2D *)_AChistos_minus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsA7"+sName+""+sName)));
      _hAUL_minus_vars.emplace_back((TH2D *)_AChistos_minus->Get(TString("angularCoefficients_"+_syst_kind+"/harmonicsAUL"+sName+""+sName)));
      _htotMap_minus_vars.emplace_back((TH2D *)_AChistos_minus->Get(TString("angularCoefficients_"+_syst_kind+"/mapTot"+sName)));
    };
    
    
  }
  
  ~getACValues(){};

  RNode run(RNode) override;
  void getAngCoeff();
};

#endif
