#include "interface/reweightFromZ.hpp"
#include "interface/functions.hpp"

RNode reweightFromZ::run(RNode d)
{

    auto getWeightPt = [this](float Pt) -> float {
        int bin = _hPt->FindBin(Pt);
        if (_hPt->IsBinOverflow(bin))
            return _hPt->GetBinContent(_hPt->GetNbinsX());
        else if (_hPt->IsBinUnderflow(bin))
            return _hPt->GetBinContent(1);
        else
            return _hPt->GetBinContent(bin);
    };
    auto getWeightY = [this](float y) -> float {
        float absy = TMath::Abs(y);
        int bin = _hY->FindBin(absy);
        if (_hY->IsBinOverflow(bin))
            return _hY->GetBinContent(_hY->GetNbinsX());
        else if (_hY->IsBinUnderflow(bin))
            return _hY->GetBinContent(1);
        else
            return _hY->GetBinContent(bin);
    };
    
    // auto getQtweights = [&](float Pt) { // used only for 10MeV study
    //     ROOT::VecOps::RVec<float> v;
    //     float vUp = 1.;//Pt;
    //     float vDown = 1.;//#Pt;
    //     if(Pt< 5) {
    //         vUp = vUp*1.04;
    //         vDown = vDown/1.04;
    //     }
    //     v.emplace_back(vUp);
    //     v.emplace_back(vDown); 
        
    //     float vUp_s1 = 1.;//Pt;
    //     float vDown_s1 = 1.;//#Pt;
    //     if(Pt< 2) {
    //         vUp_s1 = vUp_s1*1.08;
    //         vDown_s1 = vDown_s1/1.08;
    //     }
    //     v.emplace_back(vUp_s1);
    //     v.emplace_back(vDown_s1); 
        
    //     float vUp_s2 = 1.;//Pt;
    //     float vDown_s2 = 1.;//#Pt;
    //     if(Pt< 2) {
    //         vUp_s2 = vUp_s2*1.03;
    //         vDown_s2 = vDown_s2/1.03;
    //     }
    //     v.emplace_back(vUp_s2);
    //     v.emplace_back(vDown_s2);
        
    //     float vUp_s3 = 1.;//Pt;
    //     float vDown_s3 = 1.;//#Pt;
    //     if(Pt< 2) {
    //         vUp_s3 = vUp_s3*1.01;
    //         vDown_s3 = vDown_s3/1.01;
    //     }
    //     v.emplace_back(vUp_s3);
    //     v.emplace_back(vDown_s3);
        
    //     float vUp_run1 = 1.;//Pt;
    //     float vDown_run1 = 1.;//#Pt;
    //     if(Pt< 7.5) {
    //         vUp_run1 = vUp_run1*1.01;
    //         vDown_run1 = vDown_run1/1.01;
    //     }
    //     v.emplace_back(vUp_run1);
    //     v.emplace_back(vDown_run1);
        
        
         
    //     return v;
    // };

    if(!_WJets) {
        auto d1 = d.Define("Wpt_dress", "GenV_dress[0]").Define("Wrap_dress", "GenV_dress[1]")
                   .Define("weightPt", getWeightPt, {"Wpt_dress"}).Define("weightY", getWeightY, {"Wrap_dress"});
        return d1;
        }
    auto d1 = d.Define("weightPt", getWeightPt, {"Wpt_dress"}).Define("weightY", getWeightY, {"Wrap_dress"});
    // auto d1 = d.Define("weightPt", getWeightPt, {"Wpt_dress"}).Define("weightY", getWeightY, {"Wrap_dress"}).Define("Qtweights", getQtweights, {"Wpt_dress"}); // used only for 10Mev study

    return d1;
}
