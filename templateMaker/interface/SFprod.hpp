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

public:
    SFprod(TFile *SF, std::string era = "preVFP")
    {
        _SF = SF;
        TString tag = "BtoH"; //for all 2016
        if (era == "preVFP")
            tag = "BtoF";
        else if (era == "postVFP")
            tag = "GtoH";
        std::cout << "SF tag is " << tag << std::endl;
        _iso_plus = (TH2D *)_SF->Get("fullSF2D_isoTrigPlus_" + tag);
        _antiiso_plus = (TH2D *)_SF->Get("fullSF2D_antiisoTrigPlus_" + tag);
        _iso_minus = (TH2D *)_SF->Get("fullSF2D_isoTrigMinus_" + tag);
        _antiiso_minus = (TH2D *)_SF->Get("fullSF2D_antiisoTrigMinus_" + tag);
    };
    ~SFprod(){};

    RNode run(RNode) override;
};

#endif
