#ifndef MTDEFINITIONS_H
#define MTDEFINITIONS_H

#include "module.hpp"
#include<string>
class mtDefinitions : public Module
{

private:
    bool _isMC;
    std::string _ptprefix;
    std::string _phiprefix;
public:
  mtDefinitions(bool isMC = 1, std::string ptprefix = "MET_T1_pt", std::string phiprefix = "MET_T1_phi"){
      _isMC = isMC;
      _ptprefix = ptprefix;
      _phiprefix = phiprefix;
    };
    ~mtDefinitions(){};

    RNode run(RNode) override;

};

#endif
