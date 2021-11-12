#ifndef JMEPROD_H
#define JMEPROD_H

#include "module.hpp"
#include "TString.h"
#include "TFile.h"
#include<vector>
#include<string>
#include "JMESystematicsCalculators.h"
#include "JetCorrectorParameters.h"
#include<iostream>
class jmeProducer : public Module
{

private:
  Type1METVariationsCalculator* myJetVarCalc;
public:
  jmeProducer(int era)
    {
      std::string jecDataSrc = "/scratchnvme/wmass/CMSJMECalculators/data/";
      std::string jerDataSrc = "/scratchnvme/wmass/CMSJMECalculators/data/";
      if(era==1) {
	jecDataSrc += "Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC";
	jerDataSrc += "Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC";
      } else {
	jecDataSrc += "Summer19UL16_V7_MC/Summer19UL16_V7_MC";
	jerDataSrc += "Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC";
      }
      std::vector<JetCorrectorParameters> jecParams;
      myJetVarCalc = new Type1METVariationsCalculator();
      std::string textfilepath = jecDataSrc +  "_L1FastJet_AK4PFchs.txt";
      jecParams.push_back(JetCorrectorParameters(textfilepath));
      myJetVarCalc->setL1JEC(jecParams);
      textfilepath = jecDataSrc + "_L2Relative_AK4PFchs.txt";
      jecParams.push_back(JetCorrectorParameters(textfilepath));
      textfilepath = jecDataSrc + "_L3Absolute_AK4PFchs.txt";
      jecParams.push_back(JetCorrectorParameters(textfilepath));
      myJetVarCalc->setJEC(jecParams);
      myJetVarCalc->setL1JEC(jecParams);
      textfilepath = jecDataSrc + "_Uncertainty_AK4PFchs.txt";
      auto jecTotal = JetCorrectorParameters(textfilepath);
      myJetVarCalc->addJESUncertainty("Total", jecTotal);
      myJetVarCalc->setIsT1SmearedMET(0);
     //Now configure JER
      std::string fptResolution = jerDataSrc + "_PtResolution_AK4PFchs.txt";
      std::string fptResolutionSF = jerDataSrc + "_SF_AK4PFchs.txt";
      myJetVarCalc->setSmearing(fptResolution, fptResolutionSF, false, true, 0.2, 3.);

      for(auto& s: myJetVarCalc->available())
        std::cout << s << ",  ";
      std::cout  << std::endl;
    }

  RNode run(RNode d) override;
};

#endif
