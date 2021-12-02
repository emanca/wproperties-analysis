#ifndef SFPROD_H
#define SFPROD_H

#include "module.hpp"
#include "TString.h"
#include <iostream>
class SFprod : public Module
{

private:
    TFile *_SF;
    TH2D *_iso_plus;
    TH2D *_antiiso_plus;
    TH2D *_iso_minus;
    TH2D *_antiiso_minus;
    TH2D *_iso_plus_var;
    TH2D *_antiiso_plus_var;
    TH2D *_iso_minus_var;
    TH2D *_antiiso_minus_var;
    TH2D *_reco;
    TH2D *_reco_var;
    TH2D *_tracking;
    TH2D *_tracking_var;

public:
    SFprod(TFile *SF, std::string era = "preVFP")
    {
        _SF = SF;
        std::string tag = "BtoH"; //for all 2016
        if (era == "preVFP")
            tag = "BtoF";
        else if (era == "postVFP")
            tag = "GtoH";
        std::cout << "SF tag is " << tag << std::endl;
        _iso_plus = (TH2D *)_SF->Get(Form("fullSF2D_nominal_isoTrigPlus_%s", tag.c_str()));
        _antiiso_plus = (TH2D *)_SF->Get(Form("fullSF2D_nominal_antiisoTrigPlus_%s", tag.c_str()));
        _iso_minus = (TH2D *)_SF->Get(Form("fullSF2D_nominal_isoTrigMinus_%s", tag.c_str()));
        _antiiso_minus = (TH2D *)_SF->Get(Form("fullSF2D_nominal_antiisoTrigMinus_%s", tag.c_str()));
        _reco = (TH2D *)_SF->Get(Form("SF2D_nominal_reco_%s_both", tag.c_str()));
        _tracking = (TH2D *)_SF->Get(Form("SF2D_nominal_tracking_%s_both", tag.c_str()));

        _iso_plus_var = (TH2D *)_SF->Get(Form("fullSF2D_dataAltSig_isoTrigPlus_%s", tag.c_str()));
        _antiiso_plus_var = (TH2D *)_SF->Get(Form("fullSF2D_dataAltSig_antiisoTrigPlus_%s", tag.c_str()));
        _iso_minus_var = (TH2D *)_SF->Get(Form("fullSF2D_dataAltSig_isoTrigMinus_%s", tag.c_str()));
        _antiiso_minus_var = (TH2D *)_SF->Get(Form("fullSF2D_dataAltSig_antiisoTrigMinus_%s", tag.c_str()));
        _reco_var = (TH2D *)_SF->Get(Form("SF2D_dataAltSig_reco_%s_both", tag.c_str()));
        _tracking_var = (TH2D *)_SF->Get(Form("SF2D_dataAltSig_tracking_%s_both", tag.c_str()));
    };
    ~SFprod(){};

    RNode run(RNode) override;
};

#endif
