#include "genVProducer.hpp"
#include "functions.hpp"
#include "TLorentzVector.h"

RNode genVProducer::run(RNode d)
{
        //Logic equivalent to  https://github.com/WMass/nanoAOD-tools/blob/master/python/postprocessing/wmass/Vproducer.py
        //but Line 14 - how is it ensured that the idx 1 is muons?
        //to be safe I am using GenPart_mass
        auto getGenV = [](const RVec<float> &pt, const RVec<float> &eta, const RVec<float> &phi, const RVec<float> &mass, const int idx1, const int idx2)
        {
                if (idx1 < 0 || idx2 < 0)
                {
                        TLorentzVector tmp;
                        tmp.SetPtEtaPhiE(-99., -99., -99., -99.);
                        return tmp;
                }

                TLorentzVector plus, minus;
                plus.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], mass[idx1]);
                minus.SetPtEtaPhiM(pt[idx2], eta[idx2], phi[idx2], mass[idx2]);
                TLorentzVector dilepton(plus + minus);
                return dilepton;
        };

        auto d1 = d.Define("GenVP4", getGenV, {"GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass", "GenPart_preFSRLepIdx1", "GenPart_preFSRLepIdx2"})
                      .Define("Mupt_preFSR", "GenPart_pt[GenPart_preFSRLepIdx1]")
                      .Define("Mueta_preFSR", "GenPart_eta[GenPart_preFSRLepIdx1]")
                      .Define("Vpt_preFSR", [this](TLorentzVector p)
                              { return float(p.Pt()); },
                              {"GenVP4"})
                      .Define("Vrap_preFSR", [this](TLorentzVector p)
                              { return float(p.Rapidity()); },
                              {"GenVP4"})
                      .Define("Vrap_preFSR_abs", "TMath::Abs(Vrap_preFSR)")
                      .Define("Vphi_preFSR", [this](TLorentzVector p)
                              { return float(p.Phi()); },
                              {"GenVP4"})
                      .Define("Vmass_preFSR", [this](TLorentzVector p)
                              { return float(p.M()); },
                              {"GenVP4"});

        return d1;
}
