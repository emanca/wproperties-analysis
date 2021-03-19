#ifndef GETWEIGHTS_H
#define GETWEIGHTS_H

#include "module.hpp"

class getWeights : public Module {

    private:
        std::string _syst_kind = "";
        ROOT::VecOps::RVec<std::string> _syst_name;
        
    public:
        
        getWeights(){}
        
        getWeights( std::string syst_kind, std::vector<std::string> syst_name){
            _syst_kind = syst_kind;
            _syst_name = syst_name;
        }
    
    ~getWeights() {};
    RNode run(RNode) override;

};

#endif