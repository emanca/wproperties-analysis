#include "CSvariableProducer.hpp"
#include "functions.hpp"
#include "TLorentzVector.h"
#include "TMath.h"
#include<utility>
RNode CSvariableProducer::run(RNode d) {
  //Logic taken from https://github.com/WMass/nanoAOD-tools/blob/master/python/postprocessing/wmass/CSVariables.py
  
  auto getCSangles = [](const RVec<float>& pt, const RVec<float>& eta, const RVec<float>& phi, const RVec<float>& mass, const int idx1, const int idx2) {
    std::pair<float, float> cs(-99., -99.);
    if(idx1<0 || idx2<0)  return cs;

    TLorentzVector plus, minus;
    plus.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], mass[idx1]);
    minus.SetPtEtaPhiM(pt[idx2], eta[idx2], phi[idx2], mass[idx2]);
    
    TLorentzVector dilepton(plus + minus);

    int sign = dilepton.Z() ? std::abs(dilepton.Z())/dilepton.Z() : 0;

    const float ProtonMass = 0.938272;
    const float BeamEnergy = 6500.000;

    TLorentzVector p1, p2;

    p1.SetPxPyPzE(0, 0,    sign*BeamEnergy, TMath::Sqrt(BeamEnergy*BeamEnergy+ProtonMass*ProtonMass)); 
    p2.SetPxPyPzE(0, 0, -1.*sign*BeamEnergy, TMath::Sqrt(BeamEnergy*BeamEnergy+ProtonMass*ProtonMass));
		  
    p1.Boost(-dilepton.BoostVector());
    p2.Boost(-dilepton.BoostVector());

    auto CSAxis = (p1.Vect().Unit()-p2.Vect().Unit()).Unit(); //quantise along axis that bisects the boosted beams
        
    auto yAxis = (p1.Vect().Unit()).Cross((-p2.Vect().Unit())); //other axes
    yAxis = yAxis.Unit();
    auto xAxis = yAxis.Cross(CSAxis);
    xAxis = xAxis.Unit();
   
    auto boostedLep = minus;
    boostedLep.Boost(-dilepton.BoostVector());
    
    float csphi = TMath::ATan2((boostedLep.Vect()*yAxis),(boostedLep.Vect()*xAxis));
    csphi = csphi<0. ? csphi + 2.*TMath::Pi() : csphi;


    cs.first = TMath::Cos(boostedLep.Angle(CSAxis));
    cs.second = csphi;
    
    return cs;        
  };

  
  auto d1 = d.Define("CSAngles", getCSangles, {"GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass", "GenPart_preFSRLepIdx1", "GenPart_preFSRLepIdx2"})
             .Define("CStheta_preFSR", [](const std::pair<float, float>& cs){ return cs.first;}, {"CSAngles"})
             .Define("CSphi_preFSR", [](const std::pair<float, float>& cs){ return cs.second;}, {"CSAngles"});


  return d1;
}
