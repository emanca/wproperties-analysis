#ifndef MUONPREFIREWEIGHTPRODUCER_H
#define MUONPREFIREWEIGHTPRODUCER_H

#include "module.hpp"
#include<string>
class muonPrefireWeightProducer : public Module {
private:
  TFile* _wtfile;
  TH2F* _hwt;
  TH2F* _hHotspot;
public:

  muonPrefireWeightProducer(TFile* wtfile, int era) {
    _wtfile = wtfile;
    if(era==1) {
      _hwt = (TH2F*)(_wtfile->Get("L1prefiring_muonparam_2016BG"));
      std::cout << "muonPrefireWeightProducer::Read histo muonPrefiring_preVFP" << std::endl;
    }
    else if(era==2) {
      _hwt = (TH2F*)(_wtfile->Get("L1prefiring_muonparam_2016postVFP"));
      std::cout << "muonPrefireWeightProducer::Read histo muonPrefiring_postVFP" << std::endl;
    }
    else 
      std::cout << "Wrong era code passed to muon prefire weight producer" << std::endl;
    _hHotspot = (TH2F*)(_wtfile->Get("L1prefiring_muonparam_2016_hotspot"));
    _hHotspot->SetDirectory(0);
  };
  ~muonPrefireWeightProducer() {};
  RNode run(RNode) override;
  
};

#endif
