#include "interface/SF_ul.hpp"
#include "interface/functions.hpp"
RNode SF_ul::run(RNode d)
{

    auto defineSFZ = [this](float pt1, float eta1, int charge1, bool istrigMatched1, float pt2, float eta2, int charge2)
    {
        /*
    @SRC Note: here we choose muon of one charge(set from config) as primary or preferred muon in
    a dimuon event. This muon is trigger matched. All the SFs shall be applied for this muon.
    Whereas, for the other muon, only the tracking and ID SF is appled.
    **The default charge is +ve;
      */
        float sf = 1.;
        sf *= getValFromTH2(*_tracking, eta1, pt1);
        sf *= getValFromTH2(*_idip, eta1, pt1);
        sf *= getValFromTH2(*_iso, eta1, pt1);

        if (charge1 > 0)
            sf *= getValFromTH2(*_trigger_plus, eta1, pt1);
        else
            sf *= getValFromTH2(*_trigger_minus, eta1, pt1);

        sf *= getValFromTH2(*_tracking, eta2, pt2);
        sf *= getValFromTH2(*_idip, eta2, pt2);
        sf *= getValFromTH2(*_iso_notrig, eta2, pt2);

        return sf;
    };

    auto defineSF = [this](float pt1, float eta1, float charge1, float iso1)
    {
        float sf = 1.;
        sf *= getValFromTH2(*_reco, eta1, pt1);
        sf *= getValFromTH2(*_tracking, eta1, pt1);
        sf *= getValFromTH2(*_idip, eta1, pt1);
        if (charge1 > 0)
            sf *= getValFromTH2(*_trigger_plus, eta1, pt1);
        else
            sf *= getValFromTH2(*_trigger_minus, eta1, pt1);

        if (iso1 < 0.15)
        {
            sf *= getValFromTH2(*_iso, eta1, pt1);
        }
        else
        {
            sf *= getValFromTH2(*_antiiso, eta1, pt1);
        }
        return sf;
    };

    auto defineSFSystVar = [this](float nomSF, float pt1, float eta1, float charge1, float iso1)
    {
        float sf = 1.;
        sf *= getValFromTH2(*_reco_alt, eta1, pt1);
        sf *= getValFromTH2(*_tracking_alt, eta1, pt1);
        sf *= getValFromTH2(*_idip_alt, eta1, pt1);

        if (charge1 > 0)
            sf *= getValFromTH2(*_trigger_plus_alt, eta1, pt1);
        else
            sf *= getValFromTH2(*_trigger_minus_alt, eta1, pt1);

        if (iso1 < 0.15)
            sf *= getValFromTH2(*_iso_alt, eta1, pt1);
        else
            sf *= getValFromTH2(*_antiiso_alt, eta1, pt1);
        return sf / nomSF;
    };

    auto defineSFStatVar = [this](float nomSF, float pt1, float eta1, float charge1, float iso1)
    {
        // everything but iso
        float sf = 0.;
        sf += TMath::Power(getErrFromTH2(*_reco, eta1, pt1) / getValFromTH2(*_reco, eta1, pt1), 2);
        sf += TMath::Power(getErrFromTH2(*_tracking, eta1, pt1) / getValFromTH2(*_tracking, eta1, pt1), 2);
        sf += TMath::Power(getErrFromTH2(*_idip, eta1, pt1) / getValFromTH2(*_idip, eta1, pt1), 2);

        if (charge1 > 0)
            sf += TMath::Power(getErrFromTH2(*_trigger_plus, eta1, pt1) / getValFromTH2(*_trigger_plus, eta1, pt1), 2);
        else
            sf += TMath::Power(getErrFromTH2(*_trigger_minus, eta1, pt1) / getValFromTH2(*_trigger_minus, eta1, pt1), 2);
        RVec<float>
            var;
        var.resize(4);

        var[0] = nomSF * (1 - TMath::Sqrt(sf));
        var[1] = nomSF * (1 + TMath::Sqrt(sf));

        // now iso
        float sf_up = 1.;
        float sf_down = 1.;
        bool iso_bool=false; 

        if (iso1 < 0.15)
        {
            iso_bool = true;

            sf_up = (1 + TMath::Sqrt(TMath::Power(getErrFromTH2(*_effisodata, eta1, pt1) / getValFromTH2(*_effisodata, eta1, pt1), 2) + TMath::Power(getErrFromTH2(*_effisomc, eta1, pt1) / getValFromTH2(*_effisomc, eta1, pt1), 2)));

            sf_down =  (1 - TMath::Sqrt(TMath::Power(getErrFromTH2(*_effisodata, eta1, pt1) / getValFromTH2(*_effisodata, eta1, pt1), 2) + TMath::Power(getErrFromTH2(*_effisomc, eta1, pt1) / getValFromTH2(*_effisomc, eta1, pt1), 2)));
        }
        else
        {
            sf_up = (1 - TMath::Sqrt(TMath::Power(getErrFromTH2(*_effisodata, eta1, pt1) / (1-getValFromTH2(*_effisodata, eta1, pt1)), 2) + TMath::Power(getErrFromTH2(*_effisomc, eta1, pt1) / (1-getValFromTH2(*_effisomc, eta1, pt1)), 2)));

            sf_down = (1 + TMath::Sqrt(TMath::Power(getErrFromTH2(*_effisodata, eta1, pt1) / (1-getValFromTH2(*_effisodata, eta1, pt1)), 2) + TMath::Power(getErrFromTH2(*_effisomc, eta1, pt1) / (1-getValFromTH2(*_effisomc, eta1, pt1)), 2)));
        }

        var[2] = nomSF*sf_up;
        var[3] = nomSF*sf_down;

        // std::cout << nomSF << " " << var[2] << " " << var[3] << " " << iso_bool << " " << pt1 << " " << eta1 << std::endl;

        return var;
    };

    if (_isZ)
    {
        auto d1 = d.Define("SF", defineSFZ, {"Mu1_pt", "Mu1_eta", "Mu1_charge", "Mu1_hasTriggerMatch", "Mu2_pt", "Mu2_eta", "Mu2_charge"})
                      .Define("totalWeight", [](float lumi, float pu, float sf)
                              { return lumi * pu * sf; },
                              {"lumiweight", "puWeight", "SF"});
        return d1;
    }
    else
    {
        auto d1 = d.Define("SF", defineSF, {"Mu1_pt", "Mu1_eta", "Mu1_charge", "Mu1_relIso"}).Define("SFStatvar", defineSFStatVar, {"SF", "Mu1_pt", "Mu1_eta", "Mu1_charge", "Mu1_relIso"}).Define("SFSystvar", defineSFSystVar, {"SF", "Mu1_pt", "Mu1_eta", "Mu1_charge", "Mu1_relIso"});
        return d1;
    }
}
