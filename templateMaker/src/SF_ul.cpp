#include "SF_ul.hpp"
#include "functions.hpp"
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
        //sf *= getValFromTH2(*_tracking, eta1, pt1);
	sf *= getValFromTH1(*_pogtrk, eta1);//From POG
        sf *= getValFromTH2(*_idip, eta1, pt1);
        sf *= getValFromTH2(*_iso, eta1, pt1);

        if (charge1 > 0)
            sf *= getValFromTH2(*_trigger_plus, eta1, pt1);
        else
            sf *= getValFromTH2(*_trigger_minus, eta1, pt1);

        //sf *= getValFromTH2(*_tracking, eta2, pt2);
	sf *= getValFromTH1(*_pogtrk, eta2);//From POG
        sf *= getValFromTH2(*_idip, eta2, pt2);
        sf *= getValFromTH2(*_iso_notrig, eta2, pt2);

        return sf;
    };

    auto defineSF = [this](float pt1, float eta1, float charge1, float iso1)
    {
        float sf = 1.;
        //sf *= getValFromTH2(*_tracking, eta1, pt1);
	sf *= getValFromTH1(*_pogtrk, eta1);//From POG
        sf *= getValFromTH2(*_idip, eta1, pt1);

        if (charge1 > 0)
            sf *= getValFromTH2(*_trigger_plus, eta1, pt1);
        else
            sf *= getValFromTH2(*_trigger_minus, eta1, pt1);

        if (iso1 < 0.15)
            sf *= getValFromTH2(*_iso, eta1, pt1);
        else
            sf *= getValFromTH2(*_antiiso, eta1, pt1);

        return sf;
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
        auto d1 = d.Define("SF", defineSF, {"Mu1_pt", "Mu1_eta", "Mu1_charge", "Mu1_relIso"});
        return d1;
    }
}
