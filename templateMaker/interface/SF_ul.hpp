#ifndef SF_UL_H
#define SF_UL_H

#include "module.hpp"
#include "TString.h"
#include<iostream>
class SF_ul : public Module
{

private:
    TFile *_SF;
    TFile *_SFpogtrk;
    TH1D* _pogtrk;
    TH2D *_tracking;
    TH2D *_idip;
    TH2D *_trigger_plus;
    TH2D *_trigger_minus;
    TH2D *_iso;
    TH2D *_iso_notrig;
    TH2D *_antiiso;
    bool _isZ;
    //this is only relevant for Z studies
    int _prefCharge;

public:
  SF_ul(TFile *SF, TFile* SFpogtrk, bool isZ = false, std::string era = "preVFP", int prefCharge = 1)
    {
        _SF = SF;
	_SFpogtrk = SFpogtrk;
        _isZ = isZ;
	TString tag = "BtoH";//for all 2016
	if(era == "preVFP") {
	  tag = "BtoF";
	  _pogtrk = (TH1D*)_SFpogtrk->Get("muonPOGtrackingSF_preVFP");
	}
	else if(era == "postVFP") {
	  tag = "GtoH";
	  _pogtrk = (TH1D*)_SFpogtrk->Get("muonPOGtrackingSF_postVFP");
	}
	std::cout << "SF tag is " << tag << std::endl;
        _tracking = (TH2D *)_SF->Get("SF2D_tracking_" + tag + "_both");
        _idip = (TH2D *)_SF->Get("SF2D_idip_" + tag + "_both");
        _trigger_plus = (TH2D *)_SF->Get("SF2D_trigger_" + tag + "_plus");
        _trigger_minus = (TH2D *)_SF->Get("SF2D_trigger_" + tag + "_minus");
        _iso = (TH2D *)_SF->Get("SF2D_iso_" + tag + "_both");
        _iso_notrig = (TH2D *)_SF->Get("SF2D_isonotrig_" + tag + "_both");
        _antiiso = (TH2D *)_SF->Get("SF2D_antiiso_" + tag + "_both");
	_prefCharge = prefCharge;
    };
    ~SF_ul(){};

    RNode run(RNode) override;
};

#endif
