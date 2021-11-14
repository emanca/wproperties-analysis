#include "muonPrefireWeightProducer.hpp"
#include "functions.hpp"

RNode muonPrefireWeightProducer::run(RNode d) {
  auto getMuonPrefiringSF = [this](const RVec<float>& eta, const RVec<float>& pt,
				   const RVec<float>& phi) {
    float sf = 1.0;
    int nBins = _hwt->GetNbinsX();
    int prefireBin = 0;
    float prefiringProbability = 0.0;
    for (unsigned int i = 0; i < eta.size(); ++i) {
      //if (not looseId[i]) continue;
      if (eta[i] > 1.24 and eta[i] < 1.6 and phi[i] > 2.44346 and phi[i] < 2.79253) {
	prefiringProbability = _hHotspot->GetBinContent(1, 3)/(TMath::Exp( (pt[i] - _hHotspot->GetBinContent(1, 1)) / _hHotspot->GetBinContent(1, 2) ) + 1);
      } else {
	prefireBin = std::max(1, std::min(_hwt->GetXaxis()->FindFixBin(fabs(eta[i])), nBins));
	prefiringProbability = _hwt->GetBinContent(prefireBin, 3)/(TMath::Exp( (pt[i] - _hwt->GetBinContent(prefireBin, 1)) / _hwt->GetBinContent(prefireBin, 2) ) + 1);
    }
      if (prefiringProbability < 0) prefiringProbability = 0.0;
      else if (prefiringProbability > 1) prefiringProbability = 1.0;    
      sf *= (1.0 - prefiringProbability);
    }
  return sf;  
  
  };
  
  auto getMuonPrefiringSFvariationSyst = [this](const RVec<float>& eta,
						const RVec<float>& pt, const RVec<float>& phi) {
    // should be called with nSystVar = 3 (at least, can codify more variations) representing nominal, Up, Down
    RVec<float> res(2, 1.0); // initialize to 1
    int nBins = _hwt->GetNbinsX();
    int prefireBin = 0;
    double prefiringProbability = 0.0;
    for (unsigned int i = 0; i < eta.size(); ++i) {
      if (eta[i] > 1.24 and eta[i] < 1.6 and phi[i] > 2.44346 and phi[i] < 2.79253) {
	prefiringProbability = _hHotspot->GetBinContent(1, 3)/(TMath::Exp( (pt[i] - _hHotspot->GetBinContent(1, 1)) / _hHotspot->GetBinContent(1, 2) ) + 1);
      } else {
	prefireBin = std::max(1, std::min(_hwt->GetXaxis()->FindFixBin(fabs(eta[i])), nBins));
	prefiringProbability = _hwt->GetBinContent(prefireBin, 3)/(TMath::Exp( (pt[i] - _hwt->GetBinContent(prefireBin, 1)) / _hwt->GetBinContent(prefireBin, 2) ) + 1);
      }
      if (prefiringProbability < 0) prefiringProbability = 0.0;
      
      //res[0] *= (1.0 - std::min(1.0, prefiringProbability));
      res[0] *= (1.0 - std::min(1.0, 1.11*prefiringProbability));
      res[1] *= (1.0 - std::min(1.0, 0.89*prefiringProbability));
  }
    return res;  
  
  };
  
  auto d1 = d.Define("LooseMuonPt22_eta", "Muon_eta[Muon_pt > 22 && Muon_looseId]")
             .Define("LooseMuonPt22_pt", "Muon_pt[Muon_pt > 22 && Muon_looseId]")
             .Define("LooseMuonPt22_phi", "Muon_phi[Muon_pt > 22 && Muon_looseId]")
             .Define("muprefireWeight", getMuonPrefiringSF, {"LooseMuonPt22_eta", "LooseMuonPt22_pt", "LooseMuonPt22_phi"})
             .Define("muprefireWeightVars", getMuonPrefiringSFvariationSyst, {"LooseMuonPt22_eta", "LooseMuonPt22_pt", "LooseMuonPt22_phi"});
  return d1;
}
