#ifndef ANGCOEFF_H
#define ANGCOEFF_H

#include "module.hpp"

class AngCoeff : public Module
{

private:
  std::vector<float> _yArr = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 6.0};
  // std::vector<float> _ptArr = {0., 4., 8., 12., 16., 20., 24., 28., 32., 40., 60., 100., 200.};
  std::vector<float> _ptArr = {0., 4.2, 7.5, 11.8, 18.5, 32., 200.};
  // std::vector<float> _yArr = {0,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.0, 6.0}; //for regularization optimization
  // std::vector<float> _ptArr = {0. ,2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 32., 40., 60., 100., 200.}; //for regularization optimization
  // std::vector<float> _ptArr;
  // std::vector<float> _yArr;

  int _nBinsY = _yArr.size()-1;
  int _nBinsPt = _ptArr.size()-1;
  std::vector<std::string> _coeff = {"A0", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "AUL"};

  std::vector<std::string> _syst_name;
  std::string _syst_weight;
  std::string _filter;

public:
  AngCoeff(std::string filter)
  {
    _syst_weight = "";
    _filter = filter;
  };

  AngCoeff(std::string filter,std::vector<std::string> syst_name, std::string syst_weight)
  {

    _syst_name = syst_name;
    _syst_weight = syst_weight;
    _filter = filter;
  };

  ~AngCoeff(){};

  RNode run(RNode) override;
  std::vector<std::string> stringMultiplication(const std::vector<std::string> &v1, const std::vector<std::string> &v2);
};

#endif
