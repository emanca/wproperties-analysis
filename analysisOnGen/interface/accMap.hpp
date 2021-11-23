#ifndef ACCMAP_H
#define ACCMAP_H

#include "module.hpp"

class accMap : public Module
{
private:
  // std::vector<float> _yArr = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 6.0};
  // std::vector<float> _ptArr = {0., 4., 8., 12., 16., 20., 24., 28., 32., 40., 60., 100., 200.};
  // std::vector<float> _ptArr = {0.,   3.1,  5.,   7.,   9.4, 12.4, 16.5, 22.3, 32.,200}; //quantile
  // std::vector<float> _ptArr = {0., 2., 4., 6.,8.,   10,   12., 16., 22, 32.,200}; //2 GeV
  // std::vector<float> _ptArr = {0, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 13., 14., 15., 16.,22,32,200}; //len=32-->29 fitted bin
  // std::vector<float> _ptArr;
  // std::vector<float> _yArr;
  
  std::vector<float> _yArr ={0, 0.4, 0.8, 1.2 ,1.6, 2.0, 2.4, 2.8, 6.}; //Extended
  // std::vector<float> _ptArr = {0.,  2.,  4., 6.,  8., 10., 12., 14., 16., 20., 25., 32., 45., 80., 200.};//Extended
  std::vector<float> _ptArr = {0.,  2.,  4., 6.,  8., 10., 12., 16., 20., 26., 36., 60., 200.};//exteded-small
  // std::vector<float> _ptArr = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80.,200};
  
  // std::vector<float> _ptArr = {0.,2.14725709,3.75887602,5.0176509,6.35604492,7.97738581,9.8266591,12.08973559,15.1689846,19.09257703,24.54519361,33.19566976,60.,200.};
  // std::vector<float> _ptArr = {0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.,200.};
  // std::vector<float> _yArr = {0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0};
  
  //acceptance study
  // std::vector<float> _ptArr = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200};
  // std::vector<float> _yArr = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2,3.6,4.,4.4,4.8,5.2,5.6,6.2};
  // std::vector<float> _ptMuArr = {25	,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
  // std::vector<float> _etaMuArr = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
  // int _nBinsPtMu = _ptMuArr.size()-1;
  // int _nBinsEtaMu = _etaMuArr.size()-1;

  int _nBinsY = _yArr.size()-1;
  int _nBinsPt = _ptArr.size()-1;
  std::string _filter;
public:
  
  accMap(std::string filter){
    _filter = filter;
  };
  ~accMap(){};

  RNode run(RNode) override;
  
};

#endif
