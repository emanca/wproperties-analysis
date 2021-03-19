#include "interface/QCDhistos.hpp"
#include "interface/functions.hpp"
#include "interface/TH1weightsHelper.hpp"

RNode QCDhistos::run(RNode d)
{

    auto df = d.Define("weight", _weight);
    
    // TH1weightsHelper helperMT_pt26eta0(std::string("MT_pt26eta0"), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _syst_name);
    // auto MT_pt26eta0 = df.Filter(_filter).Filter("Mu1_charge>0 && Mu1_pt>26. && Mu1_pt<28. && Mu1_eta<0.1 && Mu1_eta>-0.1").Book<float, float,  ROOT::VecOps::RVec<float>>(std::move(helperMT_pt26eta0), {"MT", "weight", _syst_weight});
    // _h1Group.emplace_back(MT_pt26eta0);
    
    // TH1weightsHelper helperMT_pt50eta0(std::string("MT_pt50eta0"), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _syst_name);
    // auto MT_pt50eta0 = df.Filter(_filter).Filter("Mu1_charge>0 && Mu1_pt>50. && Mu1_pt<52. && Mu1_eta<0.1 && Mu1_eta>-0.1").Book<float, float,  ROOT::VecOps::RVec<float>>(std::move(helperMT_pt50eta0), {"MT", "weight", _syst_weight});
    // _h1Group.emplace_back(MT_pt50eta0);
    
    // TH1weightsHelper helperMT_pt26eta2(std::string("MT_pt26eta2"), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _syst_name);
    // auto MT_pt26eta2 = df.Filter(_filter).Filter("Mu1_charge>0 && Mu1_pt>26. && Mu1_pt<28. && Mu1_eta<2.2 && Mu1_eta>2.0").Book<float, float,  ROOT::VecOps::RVec<float>>(std::move(helperMT_pt26eta2), {"MT", "weight", _syst_weight});
    // _h1Group.emplace_back(MT_pt26eta2);
    
    // TH1weightsHelper helperMT_pt50eta2(std::string("MT_pt50eta2"), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _syst_name);
    // auto MT_pt50eta2 = df.Filter(_filter).Filter("Mu1_charge>0 && Mu1_pt>50. && Mu1_pt<52. && Mu1_eta<2.2 && Mu1_eta>2.0").Book<float, float,  ROOT::VecOps::RVec<float>>(std::move(helperMT_pt50eta2), {"MT", "weight", _syst_weight});
    // _h1Group.emplace_back(MT_pt50eta2);
    
    TH1weightsHelper helperMT_pt25(std::string("MT_pt25"), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _syst_name);
    auto MT_pt25 = df.Filter(_filter).Filter("Mu1_charge>0 && Mu1_pt>25. && Mu1_pt<35.").Book<float, float,  ROOT::VecOps::RVec<float>>(std::move(helperMT_pt25), {"MT", "weight", _syst_weight});
    _h1Group.emplace_back(MT_pt25);
    
    TH1weightsHelper helperMT_pt50(std::string("MT_pt50"), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _syst_name);
    auto MT_pt50 = df.Filter(_filter).Filter("Mu1_charge>0 && Mu1_pt>45. && Mu1_pt<55.").Book<float, float,  ROOT::VecOps::RVec<float>>(std::move(helperMT_pt50), {"MT", "weight", _syst_weight});
    _h1Group.emplace_back(MT_pt50);
    
    TH1weightsHelper helperMT_eta2(std::string("MT_eta2"), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _syst_name);
    auto MT_eta2 = df.Filter(_filter).Filter("Mu1_charge>0 && Mu1_eta<2.4 && Mu1_eta>1.6").Book<float, float,  ROOT::VecOps::RVec<float>>(std::move(helperMT_eta2), {"MT", "weight", _syst_weight});
    _h1Group.emplace_back(MT_eta2);
    
    TH1weightsHelper helperMT_eta0(std::string("MT_eta0"), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _syst_name);
    auto MT_eta0 = df.Filter(_filter).Filter("Mu1_charge>0 && Mu1_eta<0.4 && Mu1_eta>-0.4").Book<float, float,  ROOT::VecOps::RVec<float>>(std::move(helperMT_eta0), {"MT", "weight", _syst_weight});
    _h1Group.emplace_back(MT_eta0);
    
    
    
    TH1weightsHelper helperMT_int(std::string("MT_int"), std::string(" ; M_{T} (Rochester corr./smear MET); muon charge "), _MTArr.size() -1, _MTArr, _syst_name);
    auto MT_int = df.Filter(_filter).Filter("Mu1_charge>0").Book<float, float,  ROOT::VecOps::RVec<float>>(std::move(helperMT_int), {"MT", "weight", _syst_weight});
    _h1Group.emplace_back(MT_int);

    return df;
    
}
void QCDhistos::setAxisarrays()
{
  for (int i = 0; i < 21; i++)
    _MTArr[i] = 0. + i*(80.-0.)/20;
 _MTArr[21] = 90.;
 _MTArr[22] = 100.;
 _MTArr[23] = 120.;
 _MTArr[24] = 150.;
//   for (int i = 0; i < 3; i++)
//     _chargeArr[i] = -2. +  i*2. ; 
}
