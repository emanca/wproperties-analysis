import ROOT
import json

#fileSFul = ROOT.TFile.Open("../Common/data/allSFs.root")
##Pre-easter SF
#fileSFul = ROOT.TFile.Open("../Common/data/allSFs_eta0p1.root")
#Easter bunny SF
#fileSFul = ROOT.TFile.Open("../Common/data/2021-03-31_allSFs.root")
#SF may
fileSFul = ROOT.TFile.Open("../Common/data/2021-05-31_allSFs_nodz_dxybs.root")
# SF prod
#fileSFul = ROOT.TFile.Open("../Common/data/scaleFactorProduct_31Mar2021.root")
fileSFPogTrk=ROOT.TFile.Open("../Common/data/muonPOGtrackingSF.root")

#fileSF = ROOT.TFile.Open("../Common/data/ScaleFactors_OnTheFly.root")
filePt = ROOT.TFile.Open("../Common/data/histoUnfoldingSystPt_nsel2_dy3_rebin1_default.root")
fileY = ROOT.TFile.Open("../Common/data/histoUnfoldingSystRap_nsel2_dy3_rebin1_default.root")

pufile_mc_UL2016 = ROOT.TFile.Open("../Common/data/PileupMC_2016Legacy.root")
pufile_data_UL2016_allData = ROOT.TFile.Open("../Common/data/PileupData_2016Legacy_all2016.root")
pufile_data_UL2016_preVFP = ROOT.TFile.Open("../Common/data/PileupData_2016Legacy_upTo2016FwithHIPM.root")
pufile_data_UL2016_postVFP = ROOT.TFile.Open("../Common/data/PileupData_2016Legacy_FpostHIPMandGH.root")

datajson='../Common/data/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'


##SF file for eff correction of pre to post VFP data only
fileSFdatapre=ROOT.TFile.Open('../Common/data/scaleFactorProduct_31Mar2021.root')
##Muon prefire weights
filemuPrefire=ROOT.TFile.Open("../Common/data/muonPrefiring_prePostVFP.root")

##Common path where re-weighting files are stored
weightFoldersrc="/scratchnvme/wmass/REWEIGHT/"






