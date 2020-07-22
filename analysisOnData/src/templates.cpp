#include "interface/templates.hpp"
#include "interface/functions.hpp"

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
    TH3weightsHelper helper(std::string("template"), std::string(" ; muon #{eta}; muon p_{T} (Rochester corr.); muon charge"), _etaArr.size() - 1, _etaArr, _pTArr.size() - 1, _pTArr, _chargeArr.size() - 1, _chargeArr, _syst_name);
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
    for (unsigned int i = 0; i < 30; i++)
        _pTArr[i] = 26. + i;
    for (unsigned int i = 0; i < 49; i++)
        _etaArr[i] = -2.4 + i * (4.8) / 48; //eta -2.4 to 2.4
    for (int i = 0; i < 3; i++)
      _chargeArr[i] = -2. +  i*2. ;
}

std::vector<ROOT::RDF::RResultPtr<TH1D>> templates::getTH1()
{
    return _h1List;
}
std::vector<ROOT::RDF::RResultPtr<TH2D>> templates::getTH2()
{
    return _h2List;
}
std::vector<ROOT::RDF::RResultPtr<TH3D>> templates::getTH3()
{
    return _h3List;
}

std::vector<ROOT::RDF::RResultPtr<std::vector<TH1D>>> templates::getGroupTH1()
{
    return _h1Group;
}
std::vector<ROOT::RDF::RResultPtr<std::vector<TH2D>>> templates::getGroupTH2()
{
    return _h2Group;
}
std::vector<ROOT::RDF::RResultPtr<std::vector<TH3D>>> templates::getGroupTH3()
{
    return _h3Group;
}

void templates::reset()
{

    _h1List.clear();
    _h2List.clear();
    _h3List.clear();

    _h1Group.clear();
    _h2Group.clear();
    _h3Group.clear();
}
