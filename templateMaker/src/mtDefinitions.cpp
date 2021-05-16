#include "mtDefinitions.hpp"
#include "functions.hpp"

RNode mtDefinitions::run(RNode d)
{
    //define all nominal quantities // true for data and MC

   auto d1 = d.Define("MT", W_mt, {"Mu1_pt", "Mu1_phi", _ptprefix, _phiprefix});
   if(!_isMC)
     return d1;

    // //now get variations // true only for MC
   auto d1withCompvar = d1.Define("MT_jerUp", W_mt, {"Mu1_pt", "Mu1_phi", _ptprefix + "_jerUp", _phiprefix + "_jerUp"})
                          .Define("MT_jerDown", W_mt, {"Mu1_pt", "Mu1_phi", _ptprefix + "_jerDown", _phiprefix + "_jerDown"})
                          .Define("MT_jesTotalUp", W_mt, {"Mu1_pt", "Mu1_phi", _ptprefix + "_jesTotalUp", _phiprefix + "_jesTotalUp"})
                          .Define("MT_jesTotalDown", W_mt, {"Mu1_pt", "Mu1_phi", _ptprefix + "_jesTotalDown", _phiprefix + "_jesTotalDown"})
                          .Define("MT_unclustEnUp", W_mt, {"Mu1_pt", "Mu1_phi", _ptprefix + "_unclustEnUp", _phiprefix + "_unclustEnUp"})
                          .Define("MT_unclustEnDown", W_mt, {"Mu1_pt", "Mu1_phi", _ptprefix + "_unclustEnDown", _phiprefix + "_unclustEnDown"});
   return d1withCompvar;
}
