#ifndef TEMPLATES_H
#define TEMPLATES_H

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TMath.h"
#include "interface/module.hpp"
#include "interface/TH1weightsHelper.hpp"
#include "interface/TH3weightsHelper.hpp"
#include <map>
#include <vector>

using RNode = ROOT::RDF::RNode;

class templates : public Module
{

private:
    std::vector<ROOT::RDF::RResultPtr<TH1D>> _h1List;
    std::vector<ROOT::RDF::RResultPtr<TH2D>> _h2List;
    std::vector<ROOT::RDF::RResultPtr<TH3D>> _h3List;

    // groups of histos
    std::vector<ROOT::RDF::RResultPtr<std::vector<TH1D>>> _h1Group;
    std::vector<ROOT::RDF::RResultPtr<std::vector<TH2D>>> _h2Group;
    std::vector<ROOT::RDF::RResultPtr<std::vector<TH3D>>> _h3Group;

    std::vector<std::string> _syst_name;
    std::string _syst_weight;

    std::string _filter;
    std::vector<std::string> _filtervec;
    std::string _weight;
    std::string _colvar;
    std::vector<std::string> _colvarvec;
    HistoCategory _hcat;

    std::vector<float> _pTArr = std::vector<float>(30);
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
    std::vector<ROOT::RDF::RResultPtr<TH1D>> getTH1() override;
    std::vector<ROOT::RDF::RResultPtr<TH2D>> getTH2() override;
    std::vector<ROOT::RDF::RResultPtr<TH3D>> getTH3() override;

    std::vector<ROOT::RDF::RResultPtr<std::vector<TH1D>>> getGroupTH1() override;
    std::vector<ROOT::RDF::RResultPtr<std::vector<TH2D>>> getGroupTH2() override;
    std::vector<ROOT::RDF::RResultPtr<std::vector<TH3D>>> getGroupTH3() override;

    void reset() override;
};

#endif
