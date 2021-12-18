/// \date December 2016
/// \author Danilo Piparo
#include<iostream>
#include<fstream>
#include "TH1.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TString.h"

using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using RDf = ROOT::RDataFrame;

float W_mt(float mu_pt, float mu_phi, float met_pt, float met_phi){
  return TMath::Sqrt(2*mu_pt*met_pt*(1.0-TMath::Cos(mu_phi-met_phi)));
}


float getFromIdx(ROOT::VecOps::RVec<float> vec, int index){
  return vec[index];
}



void study_cutFlow() {
  ROOT::EnableImplicitMT(128);
  // We read the tree from the file and create a RDataFrame
  RDf d("Events", "/scratchnvme/wmass/WJetsNoCUT_v2/tree_0_1.root");//WJETS emanca
  // RDf d("Events", "/scratchnvme/wmass/NanoAOD2016-V2/SingleMuon_Run2016C/tree.root");// DATA
  // RDf d("Events", "/scratchnvme/wmass/NanoAOD2016-V2/SingleMuon_Run2016H/tree.root");
  // RDf d("Events", "/scratchnvme/wmass/NanoAOD2016-V2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/tree.root");
  // RDf d("Events", "/scratchnvme/wmass/NanoAOD2016-V2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/tree.root");

  // An optional string parameter name can be passed to the Filter method to create a named filter.
  // Named filters work as usual, but also keep track of how many entries they accept and reject.
  
  
  // auto matchOnly = d.Filter("IsTObjmatched_mu1", "trigger match only");
  // auto common  = d.Filter("HLT_SingleMu24" , "HLT cut")
  
  //  auto common  = d.Filter("HLT_SingleMu24" , "HLT cut")
                  // .Filter("IsTObjmatched_mu1", "trigger match")
  auto common = d.Filter("MET_filters==1", "MET Filter")
    // .Filter("MET_filters==1", "MET Filter")
    .Filter("nVetoElectrons==0 && 1", "veto el")
    .Filter("Idx_mu1>-1", "1 muon")
    .Define("Mu1_eta", getFromIdx, {"Muon_eta", "Idx_mu1"})
    .Define("Mu1_phi", getFromIdx, {"Muon_phi", "Idx_mu1"})
    .Define("Mu1_relIso", getFromIdx, {"Muon_pfRelIso04_all", "Idx_mu1"})
    .Define("Mu1_dz", getFromIdx, {"Muon_dz", "Idx_mu1"})
    .Define("Mu1_pt", getFromIdx, {"Muon_corrected_pt", "Idx_mu1"})
    .Define("Mu1_sip3d", getFromIdx, {"Muon_sip3d", "Idx_mu1"})
    .Define("Mu1_dxy", getFromIdx, {"Muon_dxy", "Idx_mu1"})
    .Define("MT", W_mt, {"Mu1_pt", "Mu1_phi", "MET_pt_nom", "MET_phi_nom"});
  
  // auto triggerBit = d.Filter("HLT_SingleMu24", "HLT_SingleMu24 only");
  // auto triggerBitandMatch = triggerBit.Filter("IsTObjmatched_mu1", "trigger match after trigger bit");
  
  auto signalIso = common.Filter("MT>40.", "MT>40")
                .Filter("Vtype==0", "Vtype 0")
                .Filter("HLT_SingleMu24" , "HLT sig")
                .Filter("IsTObjmatched_mu1", "trigM sig");
  auto signalAIso = common.Filter("MT>40.", "MT>40")
                .Filter("Vtype==1", "Vtype 1")
                .Filter("HLT_SingleMu24" , "HLT sigA")
                .Filter("IsTObjmatched_mu1", "trigM sigA ");
  
  // auto signalIsoHLT = signalIso.Filter("HLT_SingleMu24","111 trigger bit on signal Iso");
  // auto signalIsoMATCH = signalIsoHLT.Filter("IsTObjmatched_mu1","111 trigger matching on triggered signal iso");
  // auto signalAIsoHLT = signalAIso.Filter("HLT_SingleMu24","222 trigger bit on signal AISO");
  // auto signalAIsoMATCH = signalAIsoHLT.Filter("IsTObjmatched_mu1","222 trigger matching on triggered signal AISO");
  
  
    auto sidebandIso = common.Filter("MT<30.", "MT<30")
                    .Filter("Vtype==0", "Vtype 0")
                    .Filter("HLT_SingleMu24" , "HLT side")
                    .Filter("IsTObjmatched_mu1", "trigM side");
  auto sidebandAIso = common.Filter("MT<30.", "MT<30")
                    .Filter("Vtype==1", "Vtype 1")
                    .Filter("HLT_SingleMu24" , "HLT sideA")
                    .Filter("IsTObjmatched_mu1", "trigM sideA");
  
  // auto sidebandIsoHLT = sidebandIso.Filter("HLT_SingleMu24","333 trigger bit on sideband Iso");
  // auto sidebandIsoMATCH = sidebandIsoHLT.Filter("IsTObjmatched_mu1","333 trigger matching on triggered sideband iso");
  // auto sidebandAIsoHLT = sidebandAIso.Filter("HLT_SingleMu24","444 trigger bit on sideband AISO");
  // auto sidebandAIsoMATCH = sidebandAIsoHLT.Filter("IsTObjmatched_mu1","444 trigger matching on triggered sideband AISO");

  std::cout << "All stats:" << std::endl;
  auto allCutsReport = d.Report();
  allCutsReport->Print();


  //



}
