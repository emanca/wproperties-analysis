#include "interface/reweightFromZ.hpp"

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

    if (!_isUL && !_isWJets)
    { //old mc case
        auto d1 = d.Define("Wpt_dress", "GenV_dress[0]").Define("Wrap_dress", "GenV_dress[1]").Define("weightPt", getWeightPt, {"Wpt_dress"}).Define("weightY", getWeightY, {"Wrap_dress"});
        return d1;
    }
    else if (!_isUL && _isWJets)
    { //WJets sample V2 already has Wpt_dress column
        //@Note for @SRC: if we use preFSR, this block can go away
        auto d1 = d.Define("Wpt_dress", "GenV_dress[0]").Define("Wrap_dress", "GenV_dress[1]").Define("weightPt", getWeightPt, {"Wpt_dress"}).Define("weightY", getWeightY, {"Wrap_dress"});
        return d1;
    }
    //
    //Note for @SRC:UL case//have to check why dress is not saved.
    auto d1 = d.Define("weightPt", getWeightPt, {"Vpt_preFSR"}).Define("weightY", getWeightY, {"Vrap_preFSR"});
    return d1;
}
