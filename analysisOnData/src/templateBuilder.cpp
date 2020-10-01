#include "interface/TH3weightsHelper.hpp"
#include "interface/THNweightsHelper.hpp"
#include "interface/templateBuilder.hpp"

std::vector<std::string> templateBuilder::stringMultiplication(const std::vector<std::string> &v1, const std::vector<std::string> &v2)
{
  std::vector<std::string> products;
  if (v1.size() == 0) return v2;
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

  // templates for the fit
  auto h = new TH2F("h", "h", nBinsY, 0, 2.4, nBinsQt, 0, 32.);

  for(int j=1; j<h->GetNbinsY()+1; j++){ // for each W pt bin

    float lowEdgePt = h->GetYaxis()->GetBinLowEdge(j);
    float upEdgePt = h->GetYaxis()->GetBinUpEdge(j);

    auto sel = [lowEdgePt, upEdgePt](float pt) { return (pt >lowEdgePt && pt<upEdgePt);};

    TH3weightsHelper helperHelXsecs(std::string("qt_")+std::to_string(j)+std::string("_helXsecs_"), std::string("qt_")+std::to_string(j)+std::string("_helXsecs_"), nBinsEta, _etaArr, nBinsPt, _pTArr, nBinsY, _yArr, total);
    auto htmp = d1.Filter(sel, {"Wpt_preFSR"}).Book<float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperHelXsecs), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "weight", "harmonicsWeightsMass"});
    _h3Group.push_back(htmp);

  }

  return d1;
}

RNode templateBuilder::bookWeightVariatedhistos(RNode d)
{

  auto d1 = d.Filter("Wpt_preFSR<32. && Wrap_preFSR_abs<2.4").Define("harmonicsWeightsSyst", vecMultiplication, {_syst_weight, "harmonicsWeights"});

  std::vector<std::string> total = stringMultiplication(_syst_name, helXsecs);

  THNweightsHelper helper{"helXsecs",                                        // Name
                          "helXsecs",                                        // Title
                          {nBinsEta, nBinsPt, nBinsY, nBinsQt, nBinsCharge}, // NBins
                          {-2.4, 25., 0., 0., -2},                           // Axes min values
                          {2.4, 55., 2.4, 32., 2},                           // Axes max values
                          total};

  // We book the action: it will be treated during the event loop.
  auto templ = d1.Book<float, float, float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge", "weight", "harmonicsWeightsSyst"});
  _hNGroup.push_back(templ);

  return d1;
}

RNode templateBuilder::bookptCorrectedhistos(RNode d)

{
  auto d1 = d.Filter("Wpt_preFSR<32. && Wrap_preFSR_abs<2.4");
  for (unsigned int i = 0; i < _colvarvec.size(); i++)
  {
    std::vector<std::string> tmp;
    tmp.emplace_back(_colvarvec[i]);
    std::vector<std::string> total = stringMultiplication(tmp, helXsecs);
    THNweightsHelper helper{"helXsecs",                        // Name
                            "helXsecs",                        // Title
                            {nBinsEta, nBinsPt, nBinsY, nBinsQt, nBinsCharge}, // NBins
                            {-2.4, 25., 0., 0., -2},                           // Axes min values
                            {2.4, 55., 2.4, 32., 2},                           // Axes max values
                            total};

    // We book the action: it will be treated during the event loop.
    auto templ = d1.Filter(_filtervec[i]).Book<float, float, float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper), {"Mu1_eta", "Mu1_pt" + _colvarvec[i], "Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge", "weight", "harmonicsWeights"});
    _hNGroup.push_back(templ);
  }
  return d1;
}

RNode templateBuilder::bookJMEvarhistos(RNode d)

{
  auto d1 = d.Filter("Wpt_preFSR<32. && Wrap_preFSR_abs<2.4").Define("harmonicsWeightsPt", vecMultiplication, {"massWeights", "harmonicsWeights"});

  for (unsigned int i = 0; i < _colvarvec.size(); i++)
  {
    std::vector<std::string> tmp;
    tmp.emplace_back(_colvarvec[i]);
    std::vector<std::string> total = stringMultiplication(tmp, helXsecs);

    THNweightsHelper helper{"helXsecs",                        // Name
                            "helXsecs",                        // Title
                            {nBinsEta, nBinsPt, nBinsY, nBinsQt, nBinsCharge}, // NBins
                            {-2.4, 25., 0., 0., -2},                           // Axes min values
                            {2.4, 55., 2.4, 32., 2},                           // Axes max values
                            total};

    // We book the action: it will be treated during the event loop.
    auto templ = d1.Filter(_filtervec[i]).Book<float, float, float, float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper), {"Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs", "Wpt_preFSR", "Mu1_charge", "weight", "harmonicsWeights"});
    _hNGroup.push_back(templ);
  }
  return d1;
}

void templateBuilder::setAxisarrays()
{
  for (unsigned int i = 0; i < 31; i++){
    float binSize = (55. - 25.) / 30;
    _pTArr[i] = 25. + i*binSize;}
  for (unsigned int i = 0; i < 49; i++)
    _etaArr[i] = -2.4 + i * 4.8/48;
  for (unsigned int i=0; i< 7; i++)
    _yArr[i] = 0. + i*2.4/6;
 }
