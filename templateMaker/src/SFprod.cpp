#include "SFprod.hpp"
#include "functions.hpp"

RNode SFprod::run(RNode d)
{
    auto defineSF = [this](float pt1, float eta1, float charge1, float iso1)
    {
        float sf = 1.;
        if (charge1 > 0)
        {
            if (iso1 < 0.15) // if isolated
                sf *= getValFromTH2(*_iso_plus, eta1, pt1);
            else
                sf *= getValFromTH2(*_antiiso_plus, eta1, pt1);
        }
        else
        {
            if (iso1 < 0.15) // if isolated
                sf *= getValFromTH2(*_iso_minus, eta1, pt1);
            else
                sf *= getValFromTH2(*_antiiso_minus, eta1, pt1);
        }
        return sf;
    };

    auto defineSFVar = [this](float pt1, float eta1, float charge1, float iso1)
    {
        float sf = 1.;
        if (charge1 > 0)
        {
            if (iso1 < 0.15) // if isolated
                sf *= getValFromTH2(*_iso_plus, eta1, pt1) + getErrFromTH2(*_iso_plus, eta1, pt1);
            else
                sf *= getValFromTH2(*_antiiso_plus, eta1, pt1) - getErrFromTH2(*_antiiso_plus, eta1, pt1);
        }
        else
        {
            if (iso1 < 0.15) // if isolated
                sf *= getValFromTH2(*_iso_minus, eta1, pt1) + getErrFromTH2(*_iso_minus, eta1, pt1);
            else
                sf *= getValFromTH2(*_antiiso_minus, eta1, pt1) - getErrFromTH2(*_antiiso_minus, eta1, pt1);
        }
        return sf;
    };

    auto d1 = d.Define("SF", defineSF, {"Mu1_pt", "Mu1_eta", "Mu1_charge", "Mu1_relIso"})
               .Define("SFvar", defineSFVar, {"Mu1_pt", "Mu1_eta", "Mu1_charge", "Mu1_relIso"});
    return d1;
}
