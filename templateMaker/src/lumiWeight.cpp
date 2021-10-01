#include "lumiWeight.hpp"

RNode lumiWeight::run(RNode d)
{
    auto clipGenWeight = [](float Gen_weight) {
        float sign = Gen_weight / std::abs(Gen_weight);
        return sign;
    };

    //auto d1 = d.Define("Generator_weight_clipped", clipGenWeight, {"Generator_weight"}).Define("lumiweight", Form("float((%f*%f*Generator_weight_clipped)/(%f))", _targetLumi, _xsec, genEventSumw));
    auto d1 = d.Define("Generator_weight_clipped", clipGenWeight, {"genWeight"}).Define("lumiweight", Form("float((%f*%f*Generator_weight_clipped)/(%f))", _targetLumi, _xsec, _genEventSumwClipped));
    // auto d2 = d1.Display({"lumiweight"});
    // std::cout << " xsec " << _xsec << " _targetLumi" << _targetLumi << " _genEventSumwClipped " << _genEventSumwClipped << std::endl;
    // d2->Print();
    return d1;
}
