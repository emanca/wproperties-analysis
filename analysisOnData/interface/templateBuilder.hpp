#ifndef TEMPLATEBUILDER_H
#define TEMPLATEBUILDER_H

#include "module.hpp"

class templateBuilder : public Module
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

  std::vector<std::string> helXsecs = {"L", "I", "T", "A", "P", "7", "8", "9", "UL"};
  
  const int nBinsEta = 48;
  const int nBinsPt = 60;
  // const int nBinsPt = 600; //CHANGEBIN
  const int nBinsCharge = 2;
  std::vector<float> _pTArr = std::vector<float>(61);
  // std::vector<float> _pTArr = std::vector<float>(601); //CHANGEBIN
  std::vector<float> _etaArr = std::vector<float>(49);
  // std::vector<float> _yArr = std::vector<float>(7);
  // std::vector<float> _yArr = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  // std::vector<float> _qTArr = {0., 4., 8., 12., 16., 20., 24., 28., 32.};
  // std::vector<float> _qTArr = {0.,   3.1,  5.,   7.,   9.4, 12.4, 16.5, 22.3, 32.};  //quantile
  // std::vector<float> _qTArr = {0., 2., 4., 6.,8.,   10,   12., 16., 22, 32.}; //2 GeV
  // std::vector<float> _qTArr = {0, 1., 1cdcd.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 13., 14., 15., 16.,22.,32.}; //len=32-->29 fitted bin
  
  // std::vector<float> _yArr ={0, 0.4, 0.8, 1.2 ,1.6, 2.0, 2.4, 2.8}; //Extended
  std::vector<float> _yArr ={0, 0.4, 0.8, 1.2 ,1.6, 2.0, 2.4}; //Extended
  // std::vector<float> _qTArr = {0.,  2.,  4., 6.,  8., 10., 12., 14., 16., 20., 25., 32., 45., 80};//Extended
  std::vector<float> _qTArr = {0.,  2.,  4., 6.,  8., 10., 12., 16., 20., 26., 36., 60.};//exteded-small

  
  const int nBinsY = _yArr.size()-1;
  const int nBinsQt = _qTArr.size()-1;

public:
  templateBuilder(std::string filter, std::string weight, std::vector<std::string> syst_name, std::string syst_weight, HistoCategory hcat, std::string colvar = "")
    {

        _filter = filter;
        _weight = weight;
        _syst_name = syst_name;
        _syst_weight = syst_weight;
        _hcat = hcat;
        _colvar = colvar;
        setAxisarrays();
    };

    templateBuilder(std::vector<std::string> filtervec, std::string weight, std::vector<std::string> syst_name, std::string syst_weight, HistoCategory hcat, std::vector<std::string> colvarvec)
    {

        _filtervec = filtervec;
        _weight = weight;
        _syst_name = syst_name;
        _syst_weight = syst_weight;
        _hcat = hcat;
        _colvarvec = colvarvec;
        setAxisarrays();
    };

  ~templateBuilder(){};
  RNode bookNominalhistos(RNode);
  RNode bookptCorrectedhistos(RNode);
  RNode bookJMEvarhistos(RNode);
  RNode bookWeightVariatedhistos(RNode d);
  RNode bookWeightVariatedhistosVariedCoeff(RNode d);
  RNode run(RNode) override;
  void setAxisarrays();
  std::vector<std::string> stringMultiplication(const std::vector<std::string> &v1, const std::vector<std::string> &v2);
  static ROOT::VecOps::RVec<float> vecMultiplication(const ROOT::VecOps::RVec<float> &v1, const ROOT::VecOps::RVec<float> &v2) {
    ROOT::VecOps::RVec<float> products;
    
    products.reserve(v1.size() * v2.size());
    for (auto e1 : v1)
      for (auto e2 : v2)
	products.push_back(e1 * e2);

    return products;
  
  }

static ROOT::VecOps::RVec<float> vecMultiplicationVariedCoeff(const ROOT::VecOps::RVec<float> &v1, ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> &v2) {
    ROOT::VecOps::RVec<float> products;
    
    products.reserve(v1.size() * v2[0].size());
    for (long unsigned int v=0; v!=v1.size(); v++)
      for (auto e2 : v2[v])
        products.push_back(v1[v] * e2);
   
    return products;
  }
};

#endif
