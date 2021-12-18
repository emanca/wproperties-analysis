#ifndef ANGCOEFF_H
#define ANGCOEFF_H

#include "module.hpp"

class AngCoeff : public Module
{

private:
  // std::vector<float> _yArr = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 6.0};
  // std::vector<float> _ptArr = {0., 4., 8., 12., 16., 20., 24., 28., 32., 40., 60., 100., 200.};
  // std::vector<float> _ptArr = {0., 4.2, 7.5, 11.8, 18.5, 32., 200.}; //quantile-optimized
  // std::vector<float> _ptArr = {0., 3.1,  5.,   7.,   9.4, 12.4, 16.5, 22.3, 32.,200}; //quantile-optimized (TESTED)
  // std::vector<float> _ptArr = {0., 2., 4., 6.,8.,   10,   12., 16., 22, 32.,200}; //2 GeV
  // std::vector<float> _ptArr = {0, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 13., 14., 15., 16.,22,32,200}; //len=32-->29 fitted bin

  // std::vector<float> _yArr = {0,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.0, 6.0}; //for regularization optimization
  // std::vector<float> _ptArr = {0. ,2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 32., 40., 60., 100., 200.}; //for regularization optimization
  // std::vector<float> _ptArr;
  // std::vector<float> _yArr;
  std::vector<float> _yArr ={0, 0.4, 0.8, 1.2 ,1.6, 2.0, 2.4, 2.8, 6.}; //Extended
  // std::vector<float> _ptArr = {0.,  2.,  4., 6.,  8., 10., 12., 14., 16., 20., 25., 32., 45., 80., 200.};//Extended
  std::vector<float> _ptArr = {0.,  2.,  4., 6.,  8., 10., 12., 16., 20., 26., 36., 60., 200.};//exteded-small
  // std::vector<float> _ptArr = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80.,200};
  // std::vector<float> _Marr = std::vector<float>(1001);
  // void setAxisarrays();
  // std::vector<float> _cosThetaArr = std::vector<float>(101);
  // std::vector<float> _ptArr = {0.,2.14725709,3.75887602,5.0176509,6.35604492,7.97738581,9.8266591,12.08973559,15.1689846,19.09257703,24.54519361,33.19566976,60.,200.};
  // std::vector<float> _ptArr = {0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.,200.  };
  // std::vector<float> _yArr = {0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0};

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
    // setAxisarrays();
  };

  AngCoeff(std::string filter,std::vector<std::string> syst_name, std::string syst_weight)
  {

    _syst_name = syst_name;
    _syst_weight = syst_weight;
    _filter = filter;
    // setAxisarrays();
  };

  ~AngCoeff(){};

  RNode run(RNode) override;
  std::vector<std::string> stringMultiplication(const std::vector<std::string> &v1, const std::vector<std::string> &v2);
};

#endif
