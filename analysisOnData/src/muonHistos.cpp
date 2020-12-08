#include "interface/muonHistos.hpp"
#include "interface/functions.hpp"
#include "interface/TH2weightsHelper.hpp"

RNode muonHistos::run(RNode d)
{
  //for MC part this is a product of few columns// for Data this is just 1.// both are passed from python config
  auto d1 = d.Define("weight", _weight);
  if (_hcat == HistoCategory::Nominal)
    return bookNominalhistos(d1);
  else if (_hcat == HistoCategory::Corrected)
    return bookptCorrectedhistos(d1);
  else if (_hcat == HistoCategory::JME)
    return bookJMEvarhistos(d1);
  else
    std::cout << "Warning!! Histocategory undefined!!\n";
  return d1;
}

RNode muonHistos::bookNominalhistos(RNode df)
{
  TH2weightsHelper helperPt(std::string("Mu1_pt"), std::string(" ; muon p_{T} (Rochester corr.); muon charge"), _pTArr.size() - 1, _pTArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
  auto hpT = df.Filter(_filter).Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperPt), {"Mu1_pt", "Mu1_charge", "weight", _syst_weight});
  _h2Group.emplace_back(hpT);

  TH2weightsHelper helperEta(std::string("Mu1_eta"), std::string(" ; muon #{eta}; muon charge "), _etaArr.size() - 1, _etaArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
  auto heta = df.Filter(_filter).Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperEta), {"Mu1_eta", "Mu1_charge", "weight", _syst_weight});
  _h2Group.emplace_back(heta);

  TH2weightsHelper helperMT(std::string("MT"), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
  auto hMt = df.Filter(_filter).Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperMT), {"MT", "Mu1_charge", "weight", _syst_weight});
  _h2Group.emplace_back(hMt);

  return df;
}
//muon pt corrections affect both pt and MT in the filters, hence shape of all observables will vary
RNode muonHistos::bookptCorrectedhistos(RNode df)
{
  for (unsigned int i = 0; i < _colvarvec.size(); i++)
  {
    TH2weightsHelper helper_Pt(std::string("Mu1_pt" + _colvarvec[i]), std::string(" ; muon p_{T} (Rochester corr.); muon charge "), _pTArr.size() -1, _pTArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
    _h2Group.emplace_back(df.Filter(_filtervec[i]).Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper_Pt), {"Mu1_pt" + _colvarvec[i], "Mu1_charge", "weight", "Nom"}));

    //Only this is not affected//Name of the histo will change, but the nominal column will be plotted
  TH2weightsHelper helper_Eta(std::string("Mu1_eta" + _colvarvec[i]), std::string(" ; muon #{eta}; muon charge "), _etaArr.size() - 1, _etaArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
  _h2Group.emplace_back( df.Filter(_filtervec[i]).Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper_Eta), {"Mu1_eta", "Mu1_charge", "weight", "Nom"}) );

    TH2weightsHelper helper_MT(std::string("MT" + _colvarvec[i]), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
    _h2Group.emplace_back(df.Filter(_filtervec[i]).Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper_MT), {"MT" + _colvarvec[i], "Mu1_charge", "weight", "Nom"}));
  }
  return df;
}

//jme variations affect only MT in the filters, hence shape of all observables will vary
RNode muonHistos::bookJMEvarhistos(RNode df)
{
  for (unsigned int i = 0; i < _colvarvec.size(); i++)
  {
    TH2weightsHelper helper_Pt(std::string("Mu1_pt" + _colvarvec[i]), std::string(" ; muon p_{T} (Rochester corr.); muon charge "), _pTArr.size() -1, _pTArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
    _h2Group.emplace_back(df.Filter(_filtervec[i]).Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper_Pt), {"Mu1_pt", "Mu1_charge", "weight", "Nom"}));

  TH2weightsHelper helper_Eta(std::string("Mu1_eta" + _colvarvec[i]), std::string(" ; muon #{eta}; muon charge "), _etaArr.size() - 1, _etaArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
  _h2Group.emplace_back( df.Filter(_filtervec[i]).Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper_Eta), {"Mu1_eta", "Mu1_charge", "weight", "Nom"}) );

  //Only this column is affected
    TH2weightsHelper helper_MT(std::string("MT" + _colvarvec[i]), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
    _h2Group.emplace_back(df.Filter(_filtervec[i]).Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper_MT), {"MT" + _colvarvec[i], "Mu1_charge", "weight", "Nom"}));
  }
  return df;
}

void muonHistos::setAxisarrays()
{
  
  for (int i = 0; i < 61; i++)
    _pTArr[i] = 25. + i*(55.-25.)/60.;
  for (int i = 0; i < 49; i++)
    _etaArr[i] = -2.4 + i * (4.8) / 48; //eta -2.4 to 2.4
  for (int i = 0; i < 31; i++)
    _MTArr[i] = 0. + i*(150.-0.)/30;
  for (int i = 0; i < 3; i++)
    _chargeArr[i] = -2. +  i*2. ; //eta -1.5 to 1.5
}
