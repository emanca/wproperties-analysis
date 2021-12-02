#include "SFprod.hpp"
#include "functions.hpp"

RNode SFprod::run(RNode d)
{
    auto defineSF = [this](float pt1, float eta1, float charge1, float iso1)
    {
        float sf = 1.;
        sf *= getValFromTH2(*_reco, eta1, pt1);
        sf *= getValFromTH2(*_tracking, eta1, pt1);

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
        // std::cout<< pt1 << " " << eta1 << " " << iso1 << " "<< sf << std::endl;
        return sf;
    };

    auto defineSFStatVar = [this](float nomSF, float pt1, float eta1, float charge1, float iso1)
    {
        float sf = 0.;
        sf += TMath::Power(getErrFromTH2(*_reco, eta1, pt1) / getValFromTH2(*_reco, eta1, pt1), 2);
        sf += TMath::Power(getErrFromTH2(*_tracking, eta1, pt1) / getValFromTH2(*_tracking, eta1, pt1), 2);

        if (charge1 > 0)
        {
            if (iso1 < 0.15) // if isolated
                sf += TMath::Power(getErrFromTH2(*_iso_plus, eta1, pt1) / getValFromTH2(*_iso_plus, eta1, pt1), 2);
            else
                sf += TMath::Power(getErrFromTH2(*_antiiso_plus, eta1, pt1) / getValFromTH2(*_antiiso_plus, eta1, pt1), 2);
        }
        else
        {
            if(iso1 < 0.15) // if isolated
                sf += TMath::Power(getErrFromTH2(*_iso_minus, eta1, pt1) / getValFromTH2(*_iso_minus, eta1, pt1), 2);
            else sf += TMath::Power(getErrFromTH2(*_antiiso_minus, eta1, pt1) / getValFromTH2(*_antiiso_minus, eta1, pt1), 2);
        }
        RVec<float> var;
        var.resize(2);

        var[0] = nomSF * (1 - TMath::Sqrt(sf));
        var[1] = nomSF * (1 + TMath::Sqrt(sf));
        return var;
    };

    auto defineSFSystVar = [this](float pt1, float eta1, float charge1, float iso1)
    {
        float sf = 1.;
        sf *= getValFromTH2(*_reco_var, eta1, pt1);
        sf *= getValFromTH2(*_tracking_var, eta1, pt1);

        if (charge1 > 0)
        {
            if (iso1 < 0.15) // if isolated
                sf *= getValFromTH2(*_iso_plus_var, eta1, pt1);
            else
                sf *= getValFromTH2(*_antiiso_plus_var, eta1, pt1);
        }
        else
        {
            if (iso1 < 0.15) // if isolated
                sf *= getValFromTH2(*_iso_minus_var, eta1, pt1);
            else
                sf *= getValFromTH2(*_antiiso_minus_var, eta1, pt1);
        }
        return sf;
    };

    auto d1 = d.Define("SF", defineSF, {"Mu1_pt", "Mu1_eta", "Mu1_charge", "Mu1_relIso"})
                  .Define("SFStatvar", defineSFStatVar, {"SF", "Mu1_pt", "Mu1_eta", "Mu1_charge", "Mu1_relIso"})
                  .Define("SFSystvar", defineSFSystVar, {"Mu1_pt", "Mu1_eta", "Mu1_charge", "Mu1_relIso"});
    return d1;
}
