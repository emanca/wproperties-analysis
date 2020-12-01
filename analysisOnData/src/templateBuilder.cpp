#include "interface/TH3weightsHelper.hpp"
#include "interface/THNweightsHelper.hpp"
#include "interface/templateBuilder.hpp"
#include "interface/boostHweightsHelper.hpp"

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
  else
    std::cout << "Warning!! Histocategory undefined!!\n";
  return d1;
}

RNode templateBuilder::bookNominalhistos(RNode d)
//books nominal histos (=nominal + mass variations)
{
  auto cut = [](float pt, float y) { return pt < 32. && y < 2.4; };

  auto d1 = d.Filter(cut, {"Wpt_preFSR", "Wrap_preFSR_abs"}, "cut").Define("harmonicsWeightsMass", vecMultiplication, {"massWeights", "harmonicsWeights"});
  // auto cutReport1 = d1.Report();
  // cutReport1->Print();

  std::vector<std::string> mass = {"_massDown", "", "_massUp"};
  std::vector<std::string> total = stringMultiplication(mass, helXsecs);

  // // templates for the fit
  // auto h = new TH2F("h", "h", nBinsY, 0, 2.4, nBinsQt, 0, 32.);

  // for (int j = 1; j < h->GetNbinsY() + 1; j++)
  // { // for each W pt bin

  //   float lowEdgePt = h->GetYaxis()->GetBinLowEdge(j);
  //   float upEdgePt = h->GetYaxis()->GetBinUpEdge(j);

  //   auto sel = [lowEdgePt, upEdgePt](float pt) { return (pt > lowEdgePt && pt < upEdgePt); };

  //   TH3weightsHelper helperHelXsecs_plus(std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
  //   auto htmp_plus = d1.Filter("Mu1_charge>0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_plus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsMass"});
  //   _h3Group.push_back(htmp_plus);
  //   TH3weightsHelper helperHelXsecs_minus(std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
  //   auto htmp_minus = d1.Filter("Mu1_charge<0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_minus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsMass"});
  //   _h3Group.push_back(htmp_minus);
  // }

  boostHweightsHelper helper("helXsecs", helXsecs, _etaArr, _pTArr, _yArr, _qTArr, _chargeArr);
  auto vec = [](float var1, float var2, float var3, float var4, float var5) { 
        ROOT::VecOps::RVec<float> vec;
        vec.emplace_back(var1);
        vec.emplace_back(var2);
        vec.emplace_back(var3);
        vec.emplace_back(var4);
        vec.emplace_back(var5);
        return vec; };
  //auto templ = d1.Book<float, float, float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge", "weight", "harmonicsWeights"});
  auto templ = d1.Define("vec", vec, {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge"})
                   .Book<ROOT::VecOps::RVec<float>, float, ROOT::VecOps::RVec<float>>(std::move(helper), {"vec", "weight", "harmonicsWeights"});
  _hNGroup.push_back(templ);

  return d1;
}

RNode templateBuilder::bookWeightVariatedhistos(RNode d)
{
  auto cut = [](float pt, float y) { return pt < 32. && y < 2.4; };

  auto d1 = d.Filter(cut, {"Wpt_preFSR", "Wrap_preFSR_abs"}, "cut").Define("harmonicsWeightsSyst", vecMultiplication, {_syst_weight, "harmonicsWeights"});

  std::vector<std::string> total = stringMultiplication(_syst_name, helXsecs);

  // templates for the fit
  auto h = new TH2F("h", "h", nBinsY, 0, 2.4, nBinsQt, 0, 32.);

  for (int j = 1; j < h->GetNbinsY() + 1; j++)
  { // for each W pt bin

    float lowEdgePt = h->GetYaxis()->GetBinLowEdge(j);
    float upEdgePt = h->GetYaxis()->GetBinUpEdge(j);

    auto sel = [lowEdgePt, upEdgePt](float pt) { return (pt > lowEdgePt && pt < upEdgePt); };

    TH3weightsHelper helperHelXsecs_plus(std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wplus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
    auto htmp_plus = d1.Filter("Mu1_charge>0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_plus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsSyst"});
    _h3Group.push_back(htmp_plus);
    TH3weightsHelper helperHelXsecs_minus(std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), std::string("Wminus_")+std::string("qt_") + std::to_string(j) + std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
    auto htmp_minus = d1.Filter("Mu1_charge<0").Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs_minus), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsSyst"});
    _h3Group.push_back(htmp_minus);
  }

  return d1;
}

RNode templateBuilder::bookptCorrectedhistos(RNode d)

{
  auto cut = [](float pt, float y) { return pt < 32. && y < 2.4; };

  auto d1 = d.Filter(cut, {"Wpt_preFSR", "Wrap_preFSR_abs"}, "cut");
  for (unsigned int i = 0; i < _colvarvec.size(); i++)
  {
    std::vector<std::string> tmp;
    tmp.emplace_back(_colvarvec[i]);
    std::vector<std::string> total = stringMultiplication(tmp, helXsecs);

    // templates for the fit
    auto h = new TH2F("h", "h", nBinsY, 0, 2.4, nBinsQt, 0, 32.);

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
  auto cut = [](float pt, float y) { return pt < 32. && y < 2.4; };

  auto d1 = d.Filter(cut, {"Wpt_preFSR", "Wrap_preFSR_abs"}, "cut");
  for (unsigned int i = 0; i < _colvarvec.size(); i++)
  {
    std::vector<std::string> tmp;
    tmp.emplace_back(_colvarvec[i]);
    std::vector<std::string> total = stringMultiplication(tmp, helXsecs);

    // templates for the fit
    auto h = new TH2F("h", "h", nBinsY, 0, 2.4, nBinsQt, 0, 32.);

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
  for (unsigned int i = 0; i < 49; i++)
    _etaArr[i] = -2.4 + i * 4.8 / 48;
  for (unsigned int i = 0; i < 7; i++)
    _yArr[i] = 0. + i * 2.4 / 6;
  for (int i = 0; i < 3; i++)
    _chargeArr[i] = -2. + i * 2.;
  for (int i = 0; i < 9; i++)
    _qTArr[i] = 0. + i * 32./8;
}
