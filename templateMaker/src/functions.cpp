#ifndef FUNCTIONS_H
#include "functions.hpp"
#include "TMath.h"

float getFromIdx(ROOT::VecOps::RVec<float> vec, int index){
	return vec[index];
}

int getIntFromIdx(ROOT::VecOps::RVec<int> vec, int index){
	return vec[index];
}


float getCharge(ROOT::VecOps::RVec<int> vec, int index)
{
  return float(vec[index]);
}

float W_mt(float mu_pt, float mu_phi, float met_pt, float met_phi){
  return TMath::Sqrt(2*mu_pt*met_pt*(1.0-TMath::Cos(mu_phi-met_phi)));
}

float Wlike_mt(float mu1_pt, float mu1_phi, float mu2_pt, float mu2_phi, float met_pt, float met_phi){
  float newmet_px = met_pt*TMath::Cos(met_phi) + mu2_pt*TMath::Cos(mu2_phi);
  float newmet_py = met_pt*TMath::Sin(met_phi) + mu2_pt*TMath::Sin(mu2_phi);
  float newmet_pt = TMath::Sqrt(newmet_px*newmet_px+newmet_py*newmet_py);
  float newmet_phi = TMath::ATan2(newmet_py,newmet_px);
  return TMath::Sqrt(2*mu1_pt*newmet_pt*(1.0-TMath::Cos(mu1_phi-newmet_phi)))*(80.379/91.1876);
}

float W_hpt(float mu_pt, float mu_phi, float met_pt, float met_phi){
  float px = mu_pt*TMath::Cos(mu_phi)+met_pt*TMath::Cos(met_phi);
  float py = mu_pt*TMath::Sin(mu_phi)+met_pt*TMath::Sin(met_phi);
  return TMath::Sqrt(px*px + py*py);
}

float Z_qt(float mu1_pt, float mu1_phi, float mu2_pt, float mu2_phi){
  float px = mu1_pt*TMath::Cos(mu1_phi)+mu2_pt*TMath::Cos(mu2_phi);
  float py = mu1_pt*TMath::Sin(mu1_phi)+mu2_pt*TMath::Sin(mu2_phi);
  return TMath::Sqrt(px*px + py*py);
}

float Z_mass(float mu1_pt, float mu1_eta, float mu1_phi, float mu2_pt, float mu2_eta, float mu2_phi){
  float e  = mu1_pt*TMath::CosH(mu1_eta) + mu2_pt*TMath::CosH(mu2_eta);
  float px = mu1_pt*TMath::Cos(mu1_phi)  + mu2_pt*TMath::Cos(mu2_phi);
  float py = mu1_pt*TMath::Sin(mu1_phi)  + mu2_pt*TMath::Sin(mu2_phi);
  float pz = mu1_pt*TMath::SinH(mu1_eta) + mu2_pt*TMath::SinH(mu2_eta);
  return TMath::Sqrt(e*e - px*px - py*py - pz*pz);
}

float Z_y(float mu1_pt, float mu1_eta, float mu2_pt, float mu2_eta){
  float e  = mu1_pt*TMath::CosH(mu1_eta) + mu2_pt*TMath::CosH(mu2_eta);
  float pz = mu1_pt*TMath::SinH(mu1_eta) + mu2_pt*TMath::SinH(mu2_eta);
  return 0.5*TMath::Log((e+pz)/(e-pz));
}

ROOT::VecOps::RVec<float> float2vec(float val1, float val2, float val3, float val4, float val5)
{
  ROOT::VecOps::RVec<float> vec;
  vec.push_back(val1);
  vec.push_back(val2);
  vec.push_back(val3);
  vec.push_back(val4);
  vec.push_back(val5);
  return vec;
}

float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  return std::sqrt(deltaR2(eta1, phi1, eta2, phi2));
}

float deltaR2(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = deltaPhi(phi1, phi2);
  return deta * deta + dphi * dphi;
}

float deltaPhi(float phi1, float phi2)
{
  float result = phi1 - phi2;
  while (result > float(M_PI))
    result -= float(2 * M_PI);
  while (result <= -float(M_PI))
    result += float(2 * M_PI);
  return result;
}

bool hasTriggerMatch(const float mueta, const float muphi, const RVec<float> &TrigObj_eta, const RVec<float> &TrigObj_phi) {
  bool muhasTrigm = false;
  for (unsigned int jtrig = 0; jtrig < TrigObj_eta.size(); ++jtrig)
  {
    if (deltaR(mueta, muphi, TrigObj_eta[jtrig], TrigObj_phi[jtrig]) < 0.3)
    {
      muhasTrigm = true;
      break;
    } //if
  }

  return muhasTrigm;
}

float getValFromTH2(const TH2 &h, const float &x, const float &y)
{
  //std::cout << "x,y --> " << x << "," << y << std::endl;
  int xbin = std::max(1, std::min(h.GetNbinsX(), h.GetXaxis()->FindFixBin(x)));
  int ybin = std::max(1, std::min(h.GetNbinsY(), h.GetYaxis()->FindFixBin(y)));
  //std::cout << "xbin,ybin --> " << xbin << "," << ybin << std::endl;
  return h.GetBinContent(xbin, ybin);
}

float getErrFromTH2(const TH2 &h, const float &x, const float &y)
{
  //std::cout << "x,y --> " << x << "," << y << std::endl;
  int xbin = std::max(1, std::min(h.GetNbinsX(), h.GetXaxis()->FindFixBin(x)));
  int ybin = std::max(1, std::min(h.GetNbinsY(), h.GetYaxis()->FindFixBin(y)));
  //std::cout << "xbin,ybin --> " << xbin << "," << ybin << std::endl;
  return h.GetBinError(xbin, ybin);
}

float getValFromTH1(const TH1 &h, const float &x)
{
  int xbin = std::max(1, std::min(h.GetNbinsX(), h.GetXaxis()->FindFixBin(x)));
  return h.GetBinContent(xbin);
}

ROOT::VecOps::RVec<float> vecMultiplication(const ROOT::VecOps::RVec<float> &v1, const ROOT::VecOps::RVec<float> &v2)
{
  ROOT::VecOps::RVec<float> products;

  products.reserve(v1.size() * v2.size());
  for (auto e1 : v1)
    for (auto e2 : v2)
      products.push_back(e1 * e2);

  return products;
}

#endif

