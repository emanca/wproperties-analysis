#include "interface/TH3weightsHelper.hpp"
#include "interface/THNweightsHelper.hpp"
#include "interface/templateBuilder.hpp"

std::vector<std::string> templateBuilder::stringMultiplication(const std::vector<std::string> &v1, const std::vector<std::string> &v2)
{
  std::vector<std::string> products;
  if (v1.size() == 0)
    return v2;
  else
  {
    products.reserve(v1.size() * v2.size());
    for (auto e1 : v1)
    {
      for (auto e2 : v2)
      {
        products.push_back(e2 + e1);
      }
    }
    return products;
  }
}

RNode templateBuilder::run(RNode d)
{
  auto d1 = d.Define("weight", _weight);
  if (_hcat == HistoCategory::Nominal)
    return bookNominalhistos(d1);
  if (_hcat == HistoCategory::Weights)
    return bookWeightVariatedhistos(d1);
  else if (_hcat == HistoCategory::Corrected)
    return bookptCorrectedhistos(d1);
  else if (_hcat == HistoCategory::JME)
    return bookJMEvarhistos(d1);
  else if (_hcat == HistoCategory::WeightsVariedCoeff)
    return bookWeightVariatedhistosVariedCoeff(d1);
  else
    std::cout << "Warning!! Histocategory undefined!!\n";
  return d1;
}

RNode templateBuilder::bookNominalhistos(RNode d)
//books nominal histos (=nominal + mass variations)
{    
  // auto cut = [](float pt, float y) { return pt < 32. && y < 2.4; };
  auto cut = [](float pt, float y) { return pt < 60. && y < 2.4; };

  // auto d1 = d.Filter(cut, {"Wpt_preFSR", "Wrap_preFSR_abs"}, "cut").Define("harmonicsWeightsMass", vecMultiplication, {"massWeights", "harmonicsWeights"});
  auto d1 = d.Filter(cut, {"Wpt_preFSR", "Wrap_preFSR_abs"}, "cut");
  // auto cutReport1 = d1.Report();
  // cutReport1->Print();

  // std::vector<std::string> mass = {"_massDown", "", "_massUp"};
  // std::vector<std::string> total = stringMultiplication(mass, helXsecs);

  // templates for the fit
  auto h = new TH2F("h", "h", nBinsY, _yArr.data(), nBinsQt, _qTArr.data());

  for (int j = 1; j < h->GetNbinsY() + 1; j++)
  { // for each W pt bin

    float lowEdgePt = h->GetYaxis()->GetBinLowEdge(j);
    float upEdgePt = h->GetYaxis()->GetBinUpEdge(j);

    auto sel = [lowEdgePt, upEdgePt](float pt) { return (pt > lowEdgePt && pt < upEdgePt); };

    TH3weightsHelper helperHelXsecs_plus(std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, helXsecs);
    // TH3weightsHelper helperHelXsecs_plus(std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
    auto htmp_plus = d1.Filter(_filter).Filter("Mu1_charge>0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_plus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeights"});
    // auto htmp_plus = d1.Filter(_filter).Filter("Mu1_charge>0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_plus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsMass"});
    _h3Group.push_back(htmp_plus);
    TH3weightsHelper helperHelXsecs_minus(std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, helXsecs);
    // TH3weightsHelper helperHelXsecs_minus(std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
    auto htmp_minus = d1.Filter(_filter).Filter("Mu1_charge<0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_minus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeights"});
    // auto htmp_minus = d1.Filter(_filter).Filter("Mu1_charge<0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_minus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsMass"});
    _h3Group.push_back(htmp_minus);
  }
  return d1;
}

RNode templateBuilder::bookWeightVariatedhistos(RNode d)
{
  // auto cut = [](float pt, float y) { return pt < 32. && y < 2.4; };
  auto cut = [](float pt, float y) { return pt < 60. && y < 2.4; };

  auto d1 = d.Filter(cut, {"Wpt_preFSR", "Wrap_preFSR_abs"}, "cut").Define("harmonicsWeightsSyst", vecMultiplication, {_syst_weight, "harmonicsWeights"});

  std::vector<std::string> total = stringMultiplication(_syst_name, helXsecs);

  // templates for the fit
  auto h = new TH2F("h", "h", nBinsY, _yArr.data(), nBinsQt, _qTArr.data());

  for (int j = 1; j < h->GetNbinsY() + 1; j++)
  { // for each W pt bin

    float lowEdgePt = h->GetYaxis()->GetBinLowEdge(j);
    float upEdgePt = h->GetYaxis()->GetBinUpEdge(j);

    auto sel = [lowEdgePt, upEdgePt](float pt) { return (pt > lowEdgePt && pt < upEdgePt); };

    TH3weightsHelper helperHelXsecs_plus(std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
    auto htmp_plus = d1.Filter(_filter).Filter("Mu1_charge>0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_plus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsSyst"});
    _h3Group.push_back(htmp_plus);
    TH3weightsHelper helperHelXsecs_minus(std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
    auto htmp_minus = d1.Filter(_filter).Filter("Mu1_charge<0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_minus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsSyst"});
    _h3Group.push_back(htmp_minus);
  }

  return d1;
}

RNode templateBuilder::bookWeightVariatedhistosVariedCoeff(RNode d)
{
  // auto cut = [](float pt, float y) { return pt < 32. && y < 2.4; };
  auto cut = [](float pt, float y) { return pt < 60. && y < 2.4; };

  auto d1 = d.Filter(cut, {"Wpt_preFSR", "Wrap_preFSR_abs"}, "cut").Define("harmonicsWeightsSyst", vecMultiplicationVariedCoeff, {_syst_weight, "harmonicsWeights"}); //this is the only line different from bookWeightVariatedhistos
  
  std::vector<std::string> total = stringMultiplication(_syst_name, helXsecs);

  // templates for the fit
  auto h = new TH2F("h", "h", nBinsY, _yArr.data(), nBinsQt, _qTArr.data());

  for (int j = 1; j < h->GetNbinsY() + 1; j++)
  { // for each W pt bin

    float lowEdgePt = h->GetYaxis()->GetBinLowEdge(j);
    float upEdgePt = h->GetYaxis()->GetBinUpEdge(j);

    auto sel = [lowEdgePt, upEdgePt](float pt) { return (pt > lowEdgePt && pt < upEdgePt); };

    TH3weightsHelper helperHelXsecs_plus(std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
    auto htmp_plus = d1.Filter(_filter).Filter("Mu1_charge>0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_plus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsSyst"});
    _h3Group.push_back(htmp_plus);
    TH3weightsHelper helperHelXsecs_minus(std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
    auto htmp_minus = d1.Filter(_filter).Filter("Mu1_charge<0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_minus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsSyst"});
    _h3Group.push_back(htmp_minus);
  }

  return d1;
}


RNode templateBuilder::bookptCorrectedhistos(RNode d)

{
  // auto cut = [](float pt, float y) { return pt < 32. && y < 2.4; };
  auto cut = [](float pt, float y) { return pt < 60. && y < 2.4; };

  auto d1 = d.Filter(cut, {"Wpt_preFSR", "Wrap_preFSR_abs"}, "cut");
  for (unsigned int i = 0; i < _colvarvec.size(); i++)
  {
    std::vector<std::string> tmp;
    tmp.emplace_back(_colvarvec[i]);
    std::vector<std::string> total = stringMultiplication(tmp, helXsecs);

    // templates for the fit
    auto h = new TH2F("h", "h", nBinsY, _yArr.data(), nBinsQt, _qTArr.data());

    for (int j = 1; j < h->GetNbinsY() + 1; j++)
    { // for each W pt bin

      float lowEdgePt = h->GetYaxis()->GetBinLowEdge(j);
      float upEdgePt = h->GetYaxis()->GetBinUpEdge(j);

      auto sel = [lowEdgePt, upEdgePt](float pt) { return (pt > lowEdgePt && pt < upEdgePt); };

      TH3weightsHelper helperHelXsecs_plus(std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
      auto htmp_plus = d1.Filter("Mu1_charge>0").Filter(_filtervec[i]).Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_plus), {"Mu1_eta", "Mu1_pt" + _colvarvec[i], "Wrap_preFSR_abs", "weight", "harmonicsWeights"});
      _h3Group.push_back(htmp_plus);
      TH3weightsHelper helperHelXsecs_minus(std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
      auto htmp_minus = d1.Filter("Mu1_charge<0").Filter(_filtervec[i]).Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_minus), {"Mu1_eta", "Mu1_pt" + _colvarvec[i], "Wrap_preFSR_abs", "weight", "harmonicsWeights"});
      _h3Group.push_back(htmp_minus);
    }
  }
  return d1;
}

RNode templateBuilder::bookJMEvarhistos(RNode d)

{
  // auto cut = [](float pt, float y) { return pt < 32. && y < 2.4; };
  auto cut = [](float pt, float y) { return pt < 60. && y < 2.4; };

  auto d1 = d.Filter(cut, {"Wpt_preFSR", "Wrap_preFSR_abs"}, "cut");
  for (unsigned int i = 0; i < _colvarvec.size(); i++)
  {
    std::vector<std::string> tmp;
    tmp.emplace_back(_colvarvec[i]);
    std::vector<std::string> total = stringMultiplication(tmp, helXsecs);

    // templates for the fit
    auto h = new TH2F("h", "h", nBinsY, _yArr.data(), nBinsQt, _qTArr.data());

    for (int j = 1; j < h->GetNbinsY() + 1; j++)
    { // for each W pt bin

      float lowEdgePt = h->GetYaxis()->GetBinLowEdge(j);
      float upEdgePt = h->GetYaxis()->GetBinUpEdge(j);

      auto sel = [lowEdgePt, upEdgePt](float pt) { return (pt > lowEdgePt && pt < upEdgePt); };

      TH3weightsHelper helperHelXsecs_plus(std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
      auto htmp_plus = d1.Filter("Mu1_charge>0").Filter(_filtervec[i]).Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_plus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeights"});
      _h3Group.push_back(htmp_plus);
      TH3weightsHelper helperHelXsecs_minus(std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
      auto htmp_minus = d1.Filter("Mu1_charge<0").Filter(_filtervec[i]).Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_minus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeights"});
      _h3Group.push_back(htmp_minus);
    }
  }
  return d1;
}

void templateBuilder::setAxisarrays()
{
  for (unsigned int i = 0; i < 61; i++)
  {
    float binSize = (55. - 25.) / 60;
    _pTArr[i] = 25. + i * binSize;
  }
  // for (unsigned int i = 0; i < 601; i++) //CHANGEBIN
  // {
  //   float binSize = (55. - 25.) / 600;
  //   _pTArr[i] = 25. + i * binSize;
  // }
  for (unsigned int i = 0; i < 49; i++)
    _etaArr[i] = -2.4 + i * 4.8 / 48;
}
