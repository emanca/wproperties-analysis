#include "interface/trigObjMatchProducer.hpp"
#include "interface/functions.hpp"

RNode trigObjMatchProducer::run(RNode d) {

  auto goodMuonTriggerCandidate = [](const RVec<int> &TrigObj_id, const RVec<float> &TrigObj_pt, const RVec<float> &TrigObj_l1pt, const RVec<float> &TrigObj_l2pt, const RVec<int> &TrigObj_filterBits) {
    RVec<bool> res(TrigObj_id.size(), false); // initialize to 0
    for (unsigned int i = 0; i < res.size(); ++i)
      {
	if (TrigObj_id[i] != 13)
	  continue;
	if (TrigObj_pt[i] < 24.)
	  continue;
	if (TrigObj_l1pt[i] < 22.)
	  continue;
	if (!((TrigObj_filterBits[i] & 8) || (TrigObj_l2pt[i] > 10. && (TrigObj_filterBits[i] & 2))))
	  continue;
	res[i] = true;
      }
    return res;
  };
  
  /*//incase we want to embed as a property of the muon obj
  auto hasTriggerMatch = [](const RVec<float> mueta, const RVec<float> muphi, const RVec<float> &TrigObj_eta, const RVec<float> &TrigObj_phi) {
    RVec<bool> muhasTrigm(mueta.size(), false);
    for(unsigned int imu = 0; imu < mueta.size(); imu++) {
      for (unsigned int jtrig = 0; jtrig < TrigObj_eta.size(); ++jtrig) {
	if (deltaR(mueta[imu], muphi[imu], TrigObj_eta[jtrig], TrigObj_phi[jtrig]) < 0.3) {
	  muhasTrigm[imu] = true;
	  break;
	}//if
      }//break from this when 1 match found
    }
    return muhasTrigm;
  };

  */
  auto d1 = d.Define("goodTrigObjs", goodMuonTriggerCandidate, {"TrigObj_id", "TrigObj_pt", "TrigObj_l1pt", "TrigObj_l2pt", "TrigObj_filterBits"})
    .Define("goodTrigObjs_eta", "TrigObj_eta[goodTrigObjs]")
    .Define("goodTrigObjs_phi", "TrigObj_phi[goodTrigObjs]");
  //.Define("Muon_hasTriggerMatch", hasTriggerMatch, {"Muon_eta", "Muon_phi", "goodTrigObjs_eta", "goodTrigObjs_phi"});

  return d1;
}
