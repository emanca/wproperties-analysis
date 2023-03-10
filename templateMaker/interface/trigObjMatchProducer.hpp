#ifndef TRIGOBJECTMATCHPRODUCER_H
#define TRIGOBJECTMATCHPRODUCER_H

#include "module.hpp"

class trigObjMatchProducer : public Module {

public:
  trigObjMatchProducer() {};
  ~trigObjMatchProducer() {};
  RNode run(RNode) override;
  
};

#endif
