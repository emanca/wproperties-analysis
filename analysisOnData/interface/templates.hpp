#ifndef TEMPLATES_H
#define TEMPLATES_H

#include "module.hpp"


class templates : public Module
{

private:

    std::vector<std::string> _syst_name;
    std::string _syst_weight;

    std::string _filter;
    std::vector<std::string> _filtervec;
    std::string _weight;
    std::string _colvar;
    std::vector<std::string> _colvarvec;
    HistoCategory _hcat;

    std::vector<float> _pTArr = std::vector<float>(61);
    // std::vector<float> _pTArr = std::vector<float>(601); //CHANGEBIN
    std::vector<float> _etaArr = std::vector<float>(49);
    std::vector<float> _chargeArr = std::vector<float>(3);
    void setAxisarrays();

public:
    templates(std::string filter, std::string weight, std::vector<std::string> syst_name, std::string syst_weight, HistoCategory hcat, std::string colvar = "")
    {

        _filter = filter;
        _weight = weight;
        _syst_name = syst_name;
        _syst_weight = syst_weight;
        _hcat = hcat;
        _colvar = colvar;
        setAxisarrays();
    };

    templates(std::vector<std::string> filtervec, std::string weight, std::vector<std::string> syst_name, std::string syst_weight, HistoCategory hcat, std::vector<std::string> colvarvec)
    {

        _filtervec = filtervec;
        _weight = weight;
        _syst_name = syst_name;
        _syst_weight = syst_weight;
        _hcat = hcat;
        _colvarvec = colvarvec;
        setAxisarrays();
    };

    ~templates(){};
    RNode bookNominalhistos(RNode);
    RNode bookptCorrectedhistos(RNode);
    RNode bookJMEvarhistos(RNode);

    RNode run(RNode) override;

};

#endif
