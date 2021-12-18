#include "interface/templates.hpp"
#include "interface/functions.hpp"
#include "interface/TH3weightsHelper.hpp"

RNode templates::run(RNode d)
{

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

RNode templates::bookNominalhistos(RNode df)
{
    TH3weightsHelper helper(std::string("templates"), std::string(" ; muon #{eta}; muon p_{T} (Rochester corr.); muon charge"), _etaArr.size() - 1, _etaArr, _pTArr.size() - 1, _pTArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
    auto h = df.Filter(_filter).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper), {"Mu1_eta", "Mu1_pt", "Mu1_charge", "weight", _syst_weight});
    _h3Group.emplace_back(h);
    
    return df;
}
//muon pt corrections affect both pt and MT
RNode templates::bookptCorrectedhistos(RNode df)
{
    for (unsigned int i = 0; i < _colvarvec.size(); i++)
    {
        TH3weightsHelper helper_Pt(std::string("templates" + _colvarvec[i]), std::string(" ; muon #{eta}; muon p_{T} (Rochester corr.); muon charge"), _etaArr.size() - 1, _etaArr, _pTArr.size() - 1, _pTArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
        _h3Group.emplace_back(df.Filter(_filtervec[i]).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper_Pt), {"Mu1_eta", "Mu1_pt" + _colvarvec[i], "Mu1_charge", "weight", "Nom"}));
    }
    return df;
}

//jme variations affect only MT
RNode templates::bookJMEvarhistos(RNode df)
{
    for (unsigned int i = 0; i < _colvarvec.size(); i++)
    {
        TH3weightsHelper helper_JME(std::string("templates" + _colvarvec[i]), std::string(" ; muon #{eta}; muon p_{T} (Rochester corr.); muon charge"), _etaArr.size() - 1, _etaArr, _pTArr.size() - 1, _pTArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
        _h3Group.emplace_back(df.Filter(_filtervec[i]).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper_JME), {"Mu1_eta", "Mu1_pt", "Mu1_charge", "weight", "Nom"}));
    }
    return df;
}

void templates::setAxisarrays()
{
  for (unsigned int i = 0; i < 61; i++){
      float binSize = (55. - 25.) / 60;
      _pTArr[i] = 25. + i*binSize;}
//   for (unsigned int i = 0; i < 601; i++){ //CHANGEBIN
//       float binSize = (55. - 25.) / 600;
//       _pTArr[i] = 25. + i*binSize;}
    for (unsigned int i = 0; i < 49; i++)
        _etaArr[i] = -2.4 + i * 4.8/48;
    for (int i = 0; i < 3; i++)
      _chargeArr[i] = -2. +  i*2. ;
}
