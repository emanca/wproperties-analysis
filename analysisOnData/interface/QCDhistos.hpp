#ifndef QCDHISTOS_H
#define QCDHISTOS_H

#include "module.hpp"


class QCDhistos : public Module
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

    std::vector<float> _chargeArr = std::vector<float>(3);
    std::vector<float> _MTArr = std::vector<float>(25);
    
    void setAxisarrays();

public:
    QCDhistos(std::string filter, std::string weight, std::vector<std::string> syst_name, std::string syst_weight, HistoCategory hcat, std::string colvar = "")
    {

        _filter = filter;
        _weight = weight;
        _syst_name = syst_name;
        _syst_weight = syst_weight;
        _hcat = hcat;
        _colvar = colvar;
        setAxisarrays();
    };

    ~QCDhistos(){};

    RNode run(RNode) override;

};

#endif
