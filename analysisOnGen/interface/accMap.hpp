#ifndef ACCMAP_H
#define ACCMAP_H

#include "module.hpp"

class accMap : public Module
{
private:
  std::vector<float> _yArr = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 6.0};
  // std::vector<float> _ptArr = {0., 4., 8., 12., 16., 20., 24., 28., 32., 40., 60., 100., 200.};
  // std::vector<float> _ptArr = {0.,   3.1,  5.,   7.,   9.4, 12.4, 16.5, 22.3, 32.,200}; //quantile
  std::vector<float> _ptArr = {0., 2., 4., 6.,8.,   10,   12., 16., 22, 32.,200}; //2 GeV
  // std::vector<float> _ptArr;
  // std::vector<float> _yArr;

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
