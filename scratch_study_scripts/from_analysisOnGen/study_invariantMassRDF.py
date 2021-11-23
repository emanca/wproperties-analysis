import ROOT
ROOT.gROOT.SetBatch()





ROOT.EnableImplicitMT(64)
df = ROOT.RDataFrame("Events","/scratchnvme/wmass/WJetsNoCUT_v2/tree_*_*.root")
# histoMass = df.Filter("Wmass_preFSR>75 && 33Wmass_preFSR<85").Histo1D("Wmass_preFSR")
histoMassPlus = df.Filter("Wmass_preFSR>75 && Wmass_preFSR<85 && GenPart_pdgId[GenPart_preFSRMuonIdx]<0").Histo1D(("Wmass_preFSR_plus","Wmass_preFSR_plus",1000,75,85),"Wmass_preFSR")
histoMassMinus = df.Filter("Wmass_preFSR>75 && Wmass_preFSR<85 && GenPart_pdgId[GenPart_preFSRMuonIdx]>0").Histo1D(("Wmass_preFSR_minus","Wmass_preFSR_minus",1000,75,85),"Wmass_preFSR")


outFile = ROOT.TFile("invariantMass.root", "recreate")
outFile.cd()
histoMassPlus.Write()
histoMassMinus.Write()

