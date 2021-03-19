#ifndef MUONHISTOS_H
#define MUONHISTOS_H

#include "module.hpp"
#include<map>
#include<vector>


class muonHistos : public Module {

    private:

    std::vector<std::string> _syst_name;
    std::string _syst_weight;

    std::string _filter;
    std::vector<std::string> _filtervec;
    std::string _weight;
    std::string _colvar;
    std::vector<std::string> _colvarvec;
    HistoCategory _hcat;

    //Binning of WHelicity
    // std::vector<float> _pTArr = {26.0, 28.0, 30.0, 31.5, 33.0, 34.5, 36.0, 37.5, 39.0, 40.5, 42.0, 43.5, 45.0, 46.5, 48.0, 50.0, 52.0, 54.0, 56.0};
    //std::vector<float> _pTArr = {25.0,26.0, 28.0, 30.0, 31.5, 33.0, 34.5, 36.0, 37.5, 39.0, 40.5, 42.0, 43.5, 45.0, 46.5, 48.0, 50.0, 52.0, 54.0, 55.0};
    //std::vector<float> _etaArr = {-2.4, -2.1, -1.9, -1.7, -1.5, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0,
    //0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.7, 1.9, 2.1, 2.4};

    std::vector<float> _pTArr = std::vector<float>(61);
    std::vector<float> _etaArr = std::vector<float>(49);
    std::vector<float> _MTArr = std::vector<float>(31);
    std::vector<float> _chargeArr = std::vector<float>(3);
    void setAxisarrays();

    public:
    
    muonHistos(std::string filter, std::string weight, std::vector<std::string> syst_name, std::string syst_weight, HistoCategory hcat, std::string colvar = ""){
        
        _filter = filter;
        _weight = weight;
        _syst_name = syst_name;
        _syst_weight = syst_weight;
        _hcat = hcat;
        _colvar = colvar;
        setAxisarrays();

    };

    muonHistos(std::vector<std::string> filtervec, std::string weight, std::vector<std::string> syst_name, std::string syst_weight, HistoCategory hcat, std::vector<std::string> colvarvec){
        
        _filtervec = filtervec;
        _weight = weight;
        _syst_name = syst_name;
        _syst_weight = syst_weight;
        _hcat = hcat;
        _colvarvec = colvarvec;
        setAxisarrays();
    };

    ~muonHistos() {};
    RNode bookNominalhistos(RNode);
    RNode bookptCorrectedhistos(RNode);
    RNode bookJMEvarhistos(RNode);
    RNode run(RNode) override;
};

#endif
