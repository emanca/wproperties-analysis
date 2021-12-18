#ifndef FAKERATE_H
#define FAKERATE_H

#include "module.hpp"

class fakeRate : public Module
{

private:
  TFile *_SF;

  std::map<std::string, TH2D *> _FR;

  std::vector<std::string> _variations = {"Nominal_",
                                         "WHSF_WHSFSyst0Up", "WHSF_WHSFSyst1Up", "WHSF_WHSFSyst2Up", "WHSF_WHSFSystFlatUp", 
                                         "WHSF_WHSFSyst0Down", "WHSF_WHSFSyst1Down", "WHSF_WHSFSyst2Down", "WHSF_WHSFSystFlatDown",
                                         "LHEScaleWeight_LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_LHEScaleWeight_muR0p5_muF1p0", 
                                         "LHEScaleWeight_LHEScaleWeight_muR1p0_muF0p5","LHEScaleWeight_LHEScaleWeight_muR1p0_muF2p0", 
                                         "LHEScaleWeight_LHEScaleWeight_muR2p0_muF1p0","LHEScaleWeight_LHEScaleWeight_muR2p0_muF2p0",
                                         "ptScale_correctedUp", "ptScale_correctedDown",
					                               "jme_jesTotalUp", "jme_jesTotalDown", "jme_unclustEnUp", "jme_unclustEnDown",
                                         "PrefireWeight_PrefireWeightUp", "PrefireWeight_PrefireWeightDown",
                                         "alphaS_alphaSUp","alphaS_alphaSDown" ,
                                         "mass_massUp", "mass_massDown",
                                         "LHEPdfWeight_LHEPdfWeightHess1","LHEPdfWeight_LHEPdfWeightHess2","LHEPdfWeight_LHEPdfWeightHess3","LHEPdfWeight_LHEPdfWeightHess4","LHEPdfWeight_LHEPdfWeightHess5","LHEPdfWeight_LHEPdfWeightHess6","LHEPdfWeight_LHEPdfWeightHess7","LHEPdfWeight_LHEPdfWeightHess8","LHEPdfWeight_LHEPdfWeightHess9","LHEPdfWeight_LHEPdfWeightHess10","LHEPdfWeight_LHEPdfWeightHess11","LHEPdfWeight_LHEPdfWeightHess12","LHEPdfWeight_LHEPdfWeightHess13","LHEPdfWeight_LHEPdfWeightHess14","LHEPdfWeight_LHEPdfWeightHess15","LHEPdfWeight_LHEPdfWeightHess16","LHEPdfWeight_LHEPdfWeightHess17","LHEPdfWeight_LHEPdfWeightHess18","LHEPdfWeight_LHEPdfWeightHess19","LHEPdfWeight_LHEPdfWeightHess20","LHEPdfWeight_LHEPdfWeightHess21","LHEPdfWeight_LHEPdfWeightHess22","LHEPdfWeight_LHEPdfWeightHess23","LHEPdfWeight_LHEPdfWeightHess24","LHEPdfWeight_LHEPdfWeightHess25","LHEPdfWeight_LHEPdfWeightHess26","LHEPdfWeight_LHEPdfWeightHess27","LHEPdfWeight_LHEPdfWeightHess28","LHEPdfWeight_LHEPdfWeightHess29","LHEPdfWeight_LHEPdfWeightHess30","LHEPdfWeight_LHEPdfWeightHess31","LHEPdfWeight_LHEPdfWeightHess32","LHEPdfWeight_LHEPdfWeightHess33","LHEPdfWeight_LHEPdfWeightHess34","LHEPdfWeight_LHEPdfWeightHess35","LHEPdfWeight_LHEPdfWeightHess36","LHEPdfWeight_LHEPdfWeightHess37","LHEPdfWeight_LHEPdfWeightHess38","LHEPdfWeight_LHEPdfWeightHess39","LHEPdfWeight_LHEPdfWeightHess40","LHEPdfWeight_LHEPdfWeightHess41","LHEPdfWeight_LHEPdfWeightHess42","LHEPdfWeight_LHEPdfWeightHess43","LHEPdfWeight_LHEPdfWeightHess44","LHEPdfWeight_LHEPdfWeightHess45","LHEPdfWeight_LHEPdfWeightHess46","LHEPdfWeight_LHEPdfWeightHess47","LHEPdfWeight_LHEPdfWeightHess48","LHEPdfWeight_LHEPdfWeightHess49","LHEPdfWeight_LHEPdfWeightHess50","LHEPdfWeight_LHEPdfWeightHess51","LHEPdfWeight_LHEPdfWeightHess52","LHEPdfWeight_LHEPdfWeightHess53","LHEPdfWeight_LHEPdfWeightHess54","LHEPdfWeight_LHEPdfWeightHess55","LHEPdfWeight_LHEPdfWeightHess56","LHEPdfWeight_LHEPdfWeightHess57","LHEPdfWeight_LHEPdfWeightHess58","LHEPdfWeight_LHEPdfWeightHess59","LHEPdfWeight_LHEPdfWeightHess60",
                                         "LHEScaleWeight_WQTlow_LHEScaleWeight_muR0p5_muF0p5_WQTlow", "LHEScaleWeight_WQTlow_LHEScaleWeight_muR0p5_muF1p0_WQTlow", "LHEScaleWeight_WQTlow_LHEScaleWeight_muR1p0_muF0p5_WQTlow","LHEScaleWeight_WQTlow_LHEScaleWeight_muR1p0_muF2p0_WQTlow", "LHEScaleWeight_WQTlow_LHEScaleWeight_muR2p0_muF1p0_WQTlow","LHEScaleWeight_WQTlow_LHEScaleWeight_muR2p0_muF2p0_WQTlow",
                                         "LHEScaleWeight_WQTmid_LHEScaleWeight_muR0p5_muF0p5_WQTmid", "LHEScaleWeight_WQTmid_LHEScaleWeight_muR0p5_muF1p0_WQTmid", "LHEScaleWeight_WQTmid_LHEScaleWeight_muR1p0_muF0p5_WQTmid","LHEScaleWeight_WQTmid_LHEScaleWeight_muR1p0_muF2p0_WQTmid", "LHEScaleWeight_WQTmid_LHEScaleWeight_muR2p0_muF1p0_WQTmid","LHEScaleWeight_WQTmid_LHEScaleWeight_muR2p0_muF2p0_WQTmid",
                                         "LHEScaleWeight_WQThigh_LHEScaleWeight_muR0p5_muF0p5_WQThigh", "LHEScaleWeight_WQThigh_LHEScaleWeight_muR0p5_muF1p0_WQThigh", "LHEScaleWeight_WQThigh_LHEScaleWeight_muR1p0_muF0p5_WQThigh","LHEScaleWeight_WQThigh_LHEScaleWeight_muR1p0_muF2p0_WQThigh", "LHEScaleWeight_WQThigh_LHEScaleWeight_muR2p0_muF1p0_WQThigh","LHEScaleWeight_WQThigh_LHEScaleWeight_muR2p0_muF2p0_WQThigh",
                                         "lumi_lumiUp", "lumi_lumiDown",
                                         "topXSec_topXSecUp", "topXSec_topXSecDown",
                                         "dibosonXSec_dibosonXSecUp", "dibosonXSec_dibosonXSecDown",
                                         "tauXSec_tauXSecUp", "tauXSec_tauXSecDown",
                                         "lepVeto_lepVetoUp", "lepVeto_lepVetoDown",
                                         };
public:
  fakeRate(TFile *SF)
  {
    _SF = SF;
    for (auto s : _variations)
    {
      _FR.insert(std::pair<std::string, TH2D *>("fake_offset_" + s, (TH2D *)_SF->Get(Form("%s/fake_offset",s.c_str()))));
      _FR.insert(std::pair<std::string, TH2D *>("fake_slope_" + s, (TH2D *)_SF->Get(Form("%s/fake_slope",s.c_str()))));
      _FR.insert(std::pair<std::string, TH2D *>("prompt_offset_" + s, (TH2D *)_SF->Get(Form("%s/prompt_offset",s.c_str()))));
      _FR.insert(std::pair<std::string, TH2D *>("prompt_slope_" + s, (TH2D *)_SF->Get(Form("%s/prompt_slope",s.c_str()))));
      _FR.insert(std::pair<std::string, TH2D *>("prompt_2deg_" + s, (TH2D *)_SF->Get(Form("%s/prompt_2deg",s.c_str()))));
    }
  }
  
  // static ROOT::VecOps::RVec<float> bookFakeRate_LHEPdfWeightHess(ROOT::VecOps::RVec<float> names)
  // {
  //   ROOT::VecOps::RVec<float> fR;
  //   for(uint i=0; i<names.size(); i++){
  //     fR.emplace_back(names[i]);
  //   }
  //   return fR;
  // }

  RNode run(RNode) override;
};

#endif
