#include "lumiWeight.hpp"

RNode lumiWeight::run(RNode d)
{
    auto clipGenWeight = [](float Gen_weight)
    {
        double sign = Gen_weight / std::abs(Gen_weight);
        //return sign;
        return sign;
    };

    if(_clip){
        auto d1 = d.Define("Generator_weight_clipped", clipGenWeight, {"genWeight"}).Define("lumiweight", Form("float((%f*%f*Generator_weight_clipped)/(%f))", _targetLumi, _xsec, _genEventSumw));
        return d1;
    }
    else{
        auto d1 = d.Define("lumiweight", Form("float((%f*%f*Generator_weight)/(%f))", _targetLumi, _xsec, _genEventSumw));
        return d1;
    }
}
