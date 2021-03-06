#include "interface/zHistosROOT.hpp"
#include "interface/functions.hpp"
#include "interface/TH1weightsHelper.hpp"
#include<iostream>
//This is for validation only

RNode zHistosROOT::run(RNode d)
{

  return d;
  // auto defineNomweight = []() {
  //   ROOT::VecOps::RVec<float> One;
  //   One.emplace_back(1.);
  //   return One;
  // };

  // //for MC part this is a product of few columns// for Data this is just 1.// both are passed from python config
  // auto df = d.Define("Nom", defineNomweight).Define("weight", _weight);

  // TH1weightsHelper helperPt1(std::string("Mu1_pt"), std::string(" ; muon1 p_{T}"), _pTArr.size() - 1, _pTArr, _syst_name);
  // auto hpT1 = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperPt1), {"Mu1_pt", "weight", _syst_weight});
  // _h1Group.emplace_back(hpT1);

  // TH1weightsHelper helperPt2(std::string("Mu2_pt"), std::string(" ; muon2 p_{T}"), _pTArr.size() - 1, _pTArr, _syst_name);
  // auto hpT2 = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperPt2), {"Mu2_pt", "weight", _syst_weight});
  // _h1Group.emplace_back(hpT2);

  // TH1weightsHelper helperEta1(std::string("Mu1_eta"), std::string(" ; muon1 #{eta}; muon charge "), _etaArr.size() - 1, _etaArr, _syst_name);
  // auto heta1 = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperEta1), {"Mu1_eta", "weight", _syst_weight});
  // _h1Group.emplace_back(heta1);

  // TH1weightsHelper helperEta2(std::string("Mu2_eta"), std::string(" ; muon2 #{eta}; muon charge "), _etaArr.size() - 1, _etaArr, _syst_name);
  // auto heta2 = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperEta2), {"Mu2_eta", "weight", _syst_weight});
  // _h1Group.emplace_back(heta2);

  // TH1weightsHelper helperzMass(std::string("dimumass"), std::string(" ; diMuon mass "), _ZmassArr.size() - 1, _ZmassArr, _syst_name);
  // auto hzMass = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperzMass), {"dimuonMass", "weight", _syst_weight});
  // _h1Group.emplace_back(hzMass);

  // TH1weightsHelper helperzqt(std::string("dimuonqt"), std::string(" ; diMuon q_{T} "), _qtArr.size() - 1, _qtArr, _syst_name);
  // auto hzqt = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperzqt), {"dimuonPt", "weight", _syst_weight});
  // _h1Group.emplace_back(hzqt);

  // TH1weightsHelper helperzY(std::string("dimuonY"), std::string(" ; diMuon Y "), _etaArr.size() - 1, _etaArr, _syst_name);
  // auto hzY = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperzY), {"dimuonY", "weight", _syst_weight});
  // _h1Group.emplace_back(hzY);

  // if(_isMC) {
  //   TH1weightsHelper helperGenW(std::string("genWeight"), std::string(" ; genWeight;"), _gwArr.size() - 1, _gwArr, _syst_name);
  //   auto hGw = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperGenW), {"Generator_weight_clipped", "uno", _syst_weight});
  //   _h1Group.emplace_back(hGw);

  //   TH1weightsHelper helperGenWQt(std::string("GenV_qt"), std::string(" ; GenV #p_{T} preFSR;"), _qtArr.size() - 1, _qtArr, _syst_name);
  //   auto hGw_qt = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperGenWQt), {"Vpt_preFSR", "uno", _syst_weight});
  //   _h1Group.emplace_back(hGw_qt);
    
  //   //"Vmass_preFSR"
  //   TH1weightsHelper helperGenWY(std::string("GenV_y"), std::string(" ; GenV Y  preFSR;"), _etaArr.size() - 1, _etaArr, _syst_name);
  //   auto hGw_rap = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperGenWY), {"Vrap_preFSR", "uno", _syst_weight});
  //   _h1Group.emplace_back(hGw_rap);

  //   TH1weightsHelper helperGenCStheta(std::string("GenV_cstheta"), std::string(" ; CS #Theta;"), _cosThetaBins.size() - 1, _cosThetaBins, _syst_name);
  //   auto hGw_cstheta = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperGenCStheta), {"CStheta_preFSR", "uno", _syst_weight});
  //   _h1Group.emplace_back(hGw_cstheta);

  //   TH1weightsHelper helperGenCSphi(std::string("GenV_phi"), std::string(" ; CS #phi;"), _phiBins.size() - 1, _phiBins, _syst_name);
  //   auto hGw_csphi = df.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperGenCSphi), {"CSphi_preFSR", "uno", _syst_weight});
  //   _h1Group.emplace_back(hGw_csphi);

    
    
  // }

  // return df;
}

void zHistosROOT::setAxisarrays()
{
  
  for (int i = 0; i < 41; i++)
    _pTArr[i] = 25. + i*1.;

  for (int i = 0; i < 49; i++)
    _etaArr[i] = -2.4 + i * (4.8) / 48; //eta -2.4 to 2.4

  for (int i = 0; i < 76; i++)
    _METArr[i] = i*2.;

  for (int i = 0; i < 61; i++)
    _ZmassArr[i] = 60 + i*1.;

  for(int i = 0; i < 201; i++)
    _gwArr[i] = -12. + 200*i;
  
  auto printVec= [this](const std::vector<float>& vec) {
    std::cout << "Size=" << vec.size() << std::endl;
    for(unsigned int i = 0; i < vec.size(); i++)
      std::cout << vec[i] << ",";
    std::cout << std::endl;
  };

  printVec(_pTArr);
  printVec(_etaArr);
  printVec(_METArr);
  printVec(_ZmassArr);
  printVec(_qtArr);
  printVec(_gwArr);
}
