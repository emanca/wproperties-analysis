#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"

using namespace ROOT::VecOps;

ROOT::VecOps::RVec<float> dummy(ULong64_t event);
float getFromIdx(ROOT::VecOps::RVec<float> vec, int index);
int getIntFromIdx(ROOT::VecOps::RVec<int> vec, int index);
float getCharge(ROOT::VecOps::RVec<int> vec, int idx);
float W_mt(float,float,float,float);
float Wlike_mt(float,float,float,float,float,float);
float W_hpt(float,float,float,float);
float Z_qt(float,float,float,float);
float Z_mass(float,float,float,float,float,float);
float Z_y(float,float,float,float);
ROOT::VecOps::RVec<float> float2vec(float val1, float val2, float val3, float val4, float val5);
float deltaR(float eta1, float phi1, float eta2, float phi2);
float deltaR2(float eta1, float phi1, float eta2, float phi2);
float deltaPhi(float phi1, float phi2);
bool hasTriggerMatch(const float mueta, const float muphi, const RVec<float> &TrigObj_eta, const RVec<float> &TrigObj_phi);
float getValFromTH2(const TH2 &h, const float &x, const float &y);
float getErrFromTH2(const TH2 &h, const float &x, const float &y);