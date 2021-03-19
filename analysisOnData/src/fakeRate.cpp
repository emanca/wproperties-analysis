#include "interface/fakeRate.hpp"
#include<iostream>
RNode fakeRate::run(RNode d)
{

  RNode d1(d);
  
  // std::vector<std::string> LHEPdfWeightHess_names;
  // for(uint i=1; i<=2; i++){
    // _variations.emplace_back("LHEPdfWeight_LHEPdfWeightHess"+std::to_string(i));
    // _variations.emplace_back("LHEPdfWeight_LHEPdfWeightHess1");
    // _variations.emplace_back("LHEPdfWeight_LHEPdfWeightHess2");
  //   LHEPdfWeightHess_names.emplace_back("fakeRate_LHEPdfWeight_LHEPdfWeightHess"+std::to_string(i));
  // }
  // const ROOT::RDF::RInterface::ColumnNames_ t LHEPdfWeightHess_cols = (const ColumnNames_t)LHEPdfWeightHess_names;

  
  for (auto s : _variations)
  {
    std::cout << "Variation:" << s << std::endl;
    d1 = d1.Define(Form("fakeRate_%s", s.c_str()), [this, s](float pt, float eta, float charge, int vtype) {
      int binX = charge > 0. ? 1 : 2;
      int binY = _FR["fake_offset_" + s]->GetYaxis()->FindBin(eta);
      float f_offset = _FR["fake_offset_" + s]->GetBinContent(binX, binY);
      float f_slope = _FR["fake_slope_" + s]->GetBinContent(binX, binY);
      float p_offset = _FR["prompt_offset_" + s]->GetBinContent(binX, binY);
      float p_slope = _FR["prompt_slope_" + s]->GetBinContent(binX, binY);
      float p_2deg = _FR["prompt_2deg_" + s]->GetBinContent(binX, binY);
      float f = f_offset + (pt - 25.0) * f_slope;
      float p = p_offset * TMath::Erf(p_slope * pt + p_2deg);
      float res = vtype == 1 ? p * f / (p - f) : -(1 - p) * f / (p - f);
      
      return res;
    },
                   {"Mu1_pt", "Mu1_eta", "Mu1_charge", "Vtype"});
  }

  auto d2 = d1
                .Define("fakeRate_WHSFVars", [](float f1, float f2, float f3, float f4, float f5, float f6, float f7, float f8) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  fR.emplace_back(f3);
                  fR.emplace_back(f4);
                  fR.emplace_back(f5);
                  fR.emplace_back(f6);
                  fR.emplace_back(f7);
                  fR.emplace_back(f8);

                  return fR;
                },
                        {"fakeRate_WHSF_WHSFSyst0Up", "fakeRate_WHSF_WHSFSyst1Up", "fakeRate_WHSF_WHSFSyst2Up", "fakeRate_WHSF_WHSFSystFlatUp", "fakeRate_WHSF_WHSFSyst0Down", "fakeRate_WHSF_WHSFSyst1Down", "fakeRate_WHSF_WHSFSyst2Down", "fakeRate_WHSF_WHSFSystFlatDown"})
                .Define("fakeRate_LHEScaleWeightred", [](float f1, float f2, float f3, float f4, float f5, float f6) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  fR.emplace_back(f3);
                  fR.emplace_back(f4);
                  fR.emplace_back(f5);
                  fR.emplace_back(f6);

                  return fR;
                },
                        {"fakeRate_LHEScaleWeight_LHEScaleWeight_muR0p5_muF0p5", "fakeRate_LHEScaleWeight_LHEScaleWeight_muR0p5_muF1p0", "fakeRate_LHEScaleWeight_LHEScaleWeight_muR1p0_muF0p5", "fakeRate_LHEScaleWeight_LHEScaleWeight_muR1p0_muF2p0", "fakeRate_LHEScaleWeight_LHEScaleWeight_muR2p0_muF1p0", "fakeRate_LHEScaleWeight_LHEScaleWeight_muR2p0_muF2p0"})
                .Define("fakeRate_ptScale", [](float f1, float f2) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  return fR;
                },
                        {"fakeRate_ptScale_correctedUp", "fakeRate_ptScale_correctedDown"})
                .Define("fakeRate_jme", [](float f1, float f2, float f3, float f4) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  fR.emplace_back(f3);
                  fR.emplace_back(f4);
                  return fR;
                },
		            {"fakeRate_jme_jesTotalUp", "fakeRate_jme_jesTotalDown", "fakeRate_jme_unclustEnUp", "fakeRate_jme_unclustEnDown"})
                .Define("fakeRate_PrefireWeightVars", [](float f1, float f2) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  return fR;
		            }, {"fakeRate_PrefireWeight_PrefireWeightUp","fakeRate_PrefireWeight_PrefireWeightDown"})
               
                .Define("fakeRate_alphaSVars", [](float f1, float f2) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  return fR;
		            }, {"fakeRate_alphaS_alphaSUp","fakeRate_alphaS_alphaSDown"})
                
                .Define("fakeRate_massWeights", [](float f1, float f2) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  return fR;
		            }, {"fakeRate_mass_massUp","fakeRate_mass_massDown"})
            
                // .Define("fakeRate_LHEPdfWeight", bookFakeRate_LHEPdfWeightHess, LHEPdfWeightHess_names);
                // .Define("fakeRate_LHEPdfWeight", [](ROOT::VecOps::RVec<float> names) {
                //   ROOT::VecOps::RVec<float> fR;
                //   for(uint i=0; i<names.size(); i++){
                //     fR.emplace_back(names[i]);
                //   }
                //   return fR;
		            // },{LHEPdfWeightHess_names});
    
                .Define("fakeRate_LHEPdfWeightHess", [](float f1,float f2,float f3,float f4,float f5,float f6,float f7,float f8,float f9,float f10,float f11,float f12,float f13,float f14,float f15,float f16,float f17,float f18,float f19,float f20,float f21,float f22,float f23,float f24,float f25,float f26,float f27,float f28,float f29,float f30,float f31,float f32,float f33,float f34,float f35,float f36,float f37,float f38,float f39,float f40,float f41,float f42,float f43,float f44,float f45,float f46,float f47,float f48,float f49,float f50,float f51,float f52,float f53,float f54,float f55,float f56,float f57,float f58,float f59,float f60)
                 {                  
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);                  
                  fR.emplace_back(f2);
                  fR.emplace_back(f3);
                  fR.emplace_back(f4);
                  fR.emplace_back(f5);
                  fR.emplace_back(f6);
                  fR.emplace_back(f7);
                  fR.emplace_back(f8);
                  fR.emplace_back(f9);
                  fR.emplace_back(f10);
                  fR.emplace_back(f11);
                  fR.emplace_back(f12);
                  fR.emplace_back(f13);
                  fR.emplace_back(f14);
                  fR.emplace_back(f15);
                  fR.emplace_back(f16);
                  fR.emplace_back(f17);
                  fR.emplace_back(f18);
                  fR.emplace_back(f19);
                  fR.emplace_back(f20);
                  fR.emplace_back(f21);
                  fR.emplace_back(f22);
                  fR.emplace_back(f23);
                  fR.emplace_back(f24);
                  fR.emplace_back(f25);
                  fR.emplace_back(f26);
                  fR.emplace_back(f27);
                  fR.emplace_back(f28);
                  fR.emplace_back(f29);
                  fR.emplace_back(f30);
                  fR.emplace_back(f31);
                  fR.emplace_back(f32);
                  fR.emplace_back(f33);
                  fR.emplace_back(f34);
                  fR.emplace_back(f35);
                  fR.emplace_back(f36);
                  fR.emplace_back(f37);
                  fR.emplace_back(f38);
                  fR.emplace_back(f39);
                  fR.emplace_back(f40);
                  fR.emplace_back(f41);
                  fR.emplace_back(f42);
                  fR.emplace_back(f43);
                  fR.emplace_back(f44);
                  fR.emplace_back(f45);
                  fR.emplace_back(f46);
                  fR.emplace_back(f47);
                  fR.emplace_back(f48);
                  fR.emplace_back(f49);
                  fR.emplace_back(f50);
                  fR.emplace_back(f51);
                  fR.emplace_back(f52);
                  fR.emplace_back(f53);
                  fR.emplace_back(f54);
                  fR.emplace_back(f55);
                  fR.emplace_back(f56);
                  fR.emplace_back(f57);
                  fR.emplace_back(f58);
                  fR.emplace_back(f59);
                  fR.emplace_back(f60);
                  return fR;
		            },{"fakeRate_LHEPdfWeight_LHEPdfWeightHess1","fakeRate_LHEPdfWeight_LHEPdfWeightHess2","fakeRate_LHEPdfWeight_LHEPdfWeightHess3","fakeRate_LHEPdfWeight_LHEPdfWeightHess4","fakeRate_LHEPdfWeight_LHEPdfWeightHess5","fakeRate_LHEPdfWeight_LHEPdfWeightHess6","fakeRate_LHEPdfWeight_LHEPdfWeightHess7","fakeRate_LHEPdfWeight_LHEPdfWeightHess8","fakeRate_LHEPdfWeight_LHEPdfWeightHess9","fakeRate_LHEPdfWeight_LHEPdfWeightHess10","fakeRate_LHEPdfWeight_LHEPdfWeightHess11","fakeRate_LHEPdfWeight_LHEPdfWeightHess12","fakeRate_LHEPdfWeight_LHEPdfWeightHess13","fakeRate_LHEPdfWeight_LHEPdfWeightHess14","fakeRate_LHEPdfWeight_LHEPdfWeightHess15","fakeRate_LHEPdfWeight_LHEPdfWeightHess16","fakeRate_LHEPdfWeight_LHEPdfWeightHess17","fakeRate_LHEPdfWeight_LHEPdfWeightHess18","fakeRate_LHEPdfWeight_LHEPdfWeightHess19","fakeRate_LHEPdfWeight_LHEPdfWeightHess20","fakeRate_LHEPdfWeight_LHEPdfWeightHess21","fakeRate_LHEPdfWeight_LHEPdfWeightHess22","fakeRate_LHEPdfWeight_LHEPdfWeightHess23","fakeRate_LHEPdfWeight_LHEPdfWeightHess24","fakeRate_LHEPdfWeight_LHEPdfWeightHess25","fakeRate_LHEPdfWeight_LHEPdfWeightHess26","fakeRate_LHEPdfWeight_LHEPdfWeightHess27","fakeRate_LHEPdfWeight_LHEPdfWeightHess28","fakeRate_LHEPdfWeight_LHEPdfWeightHess29","fakeRate_LHEPdfWeight_LHEPdfWeightHess30","fakeRate_LHEPdfWeight_LHEPdfWeightHess31","fakeRate_LHEPdfWeight_LHEPdfWeightHess32","fakeRate_LHEPdfWeight_LHEPdfWeightHess33","fakeRate_LHEPdfWeight_LHEPdfWeightHess34","fakeRate_LHEPdfWeight_LHEPdfWeightHess35","fakeRate_LHEPdfWeight_LHEPdfWeightHess36","fakeRate_LHEPdfWeight_LHEPdfWeightHess37","fakeRate_LHEPdfWeight_LHEPdfWeightHess38","fakeRate_LHEPdfWeight_LHEPdfWeightHess39","fakeRate_LHEPdfWeight_LHEPdfWeightHess40","fakeRate_LHEPdfWeight_LHEPdfWeightHess41","fakeRate_LHEPdfWeight_LHEPdfWeightHess42","fakeRate_LHEPdfWeight_LHEPdfWeightHess43","fakeRate_LHEPdfWeight_LHEPdfWeightHess44","fakeRate_LHEPdfWeight_LHEPdfWeightHess45","fakeRate_LHEPdfWeight_LHEPdfWeightHess46","fakeRate_LHEPdfWeight_LHEPdfWeightHess47","fakeRate_LHEPdfWeight_LHEPdfWeightHess48","fakeRate_LHEPdfWeight_LHEPdfWeightHess49","fakeRate_LHEPdfWeight_LHEPdfWeightHess50","fakeRate_LHEPdfWeight_LHEPdfWeightHess51","fakeRate_LHEPdfWeight_LHEPdfWeightHess52","fakeRate_LHEPdfWeight_LHEPdfWeightHess53","fakeRate_LHEPdfWeight_LHEPdfWeightHess54","fakeRate_LHEPdfWeight_LHEPdfWeightHess55","fakeRate_LHEPdfWeight_LHEPdfWeightHess56","fakeRate_LHEPdfWeight_LHEPdfWeightHess57","fakeRate_LHEPdfWeight_LHEPdfWeightHess58","fakeRate_LHEPdfWeight_LHEPdfWeightHess59","fakeRate_LHEPdfWeight_LHEPdfWeightHess60"
                  })
                  
                .Define("fakeRate_LHEScaleWeightred_WQTlow", [](float f1, float f2, float f3, float f4, float f5, float f6) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  fR.emplace_back(f3);
                  fR.emplace_back(f4);
                  fR.emplace_back(f5);
                  fR.emplace_back(f6);

                  return fR;
                },
                        {"fakeRate_LHEScaleWeight_WQTlow_LHEScaleWeight_muR0p5_muF0p5_WQTlow", "fakeRate_LHEScaleWeight_WQTlow_LHEScaleWeight_muR0p5_muF1p0_WQTlow", "fakeRate_LHEScaleWeight_WQTlow_LHEScaleWeight_muR1p0_muF0p5_WQTlow", "fakeRate_LHEScaleWeight_WQTlow_LHEScaleWeight_muR1p0_muF2p0_WQTlow", "fakeRate_LHEScaleWeight_WQTlow_LHEScaleWeight_muR2p0_muF1p0_WQTlow", "fakeRate_LHEScaleWeight_WQTlow_LHEScaleWeight_muR2p0_muF2p0_WQTlow"})                  
                  
                 .Define("fakeRate_LHEScaleWeightred_WQTmid", [](float f1, float f2, float f3, float f4, float f5, float f6) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  fR.emplace_back(f3);
                  fR.emplace_back(f4);
                  fR.emplace_back(f5);
                  fR.emplace_back(f6);

                  return fR;
                },
                        {"fakeRate_LHEScaleWeight_WQTmid_LHEScaleWeight_muR0p5_muF0p5_WQTmid", "fakeRate_LHEScaleWeight_WQTmid_LHEScaleWeight_muR0p5_muF1p0_WQTmid", "fakeRate_LHEScaleWeight_WQTmid_LHEScaleWeight_muR1p0_muF0p5_WQTmid", "fakeRate_LHEScaleWeight_WQTmid_LHEScaleWeight_muR1p0_muF2p0_WQTmid", "fakeRate_LHEScaleWeight_WQTmid_LHEScaleWeight_muR2p0_muF1p0_WQTmid", "fakeRate_LHEScaleWeight_WQTmid_LHEScaleWeight_muR2p0_muF2p0_WQTmid"})                  
                      
                 .Define("fakeRate_LHEScaleWeightred_WQThigh", [](float f1, float f2, float f3, float f4, float f5, float f6) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  fR.emplace_back(f3);
                  fR.emplace_back(f4);
                  fR.emplace_back(f5);
                  fR.emplace_back(f6);

                  return fR;
                },
                        {"fakeRate_LHEScaleWeight_WQThigh_LHEScaleWeight_muR0p5_muF0p5_WQThigh", "fakeRate_LHEScaleWeight_WQThigh_LHEScaleWeight_muR0p5_muF1p0_WQThigh", "fakeRate_LHEScaleWeight_WQThigh_LHEScaleWeight_muR1p0_muF0p5_WQThigh", "fakeRate_LHEScaleWeight_WQThigh_LHEScaleWeight_muR1p0_muF2p0_WQThigh", "fakeRate_LHEScaleWeight_WQThigh_LHEScaleWeight_muR2p0_muF1p0_WQThigh", "fakeRate_LHEScaleWeight_WQThigh_LHEScaleWeight_muR2p0_muF2p0_WQThigh"})                  
                 
                .Define("fakeRate_lumi", [](float f1, float f2) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  return fR;
		            }, {"fakeRate_lumi_lumiUp", "fakeRate_lumi_lumiDown"})
                
                .Define("fakeRate_topXSec", [](float f1, float f2) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  return fR;
		            }, {"fakeRate_topXSec_topXSecUp", "fakeRate_topXSec_topXSecDown"})
                 
                .Define("fakeRate_dibosonXSec", [](float f1, float f2) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  return fR;
		            }, {"fakeRate_dibosonXSec_dibosonXSecUp", "fakeRate_dibosonXSec_dibosonXSecDown"})
                
                .Define("fakeRate_tauXSec", [](float f1, float f2) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  return fR;
		            }, {"fakeRate_tauXSec_tauXSecUp", "fakeRate_tauXSec_tauXSecDown"})
                
                .Define("fakeRate_lepVeto", [](float f1, float f2) {
                  ROOT::VecOps::RVec<float> fR;
                  fR.emplace_back(f1);
                  fR.emplace_back(f2);
                  return fR;
		            }, {"fakeRate_lepVeto_lepVetoUp", "fakeRate_lepVeto_lepVetoDown"}) 
                  ;                    
                
  
  
  return d2;
}
