import ROOT
import math
RDF = ROOT.ROOT.RDataFrame
ROOT.EnableImplicitMT(128)


NOT COMPLETE AND NOT WORKING

samples = {
    'data' : {
        "dir": ["SingleMuon_Run2016B_ver2","SingleMuon_Run2016C","SingleMuon_Run2016D","SingleMuon_Run2016E","SingleMuon_Run2016F","SingleMuon_Run2016G","SingleMuon_Run2016H",],        
        "path" : '/scratchnvme/wmass/NanoAOD2016-V2/',
    },
    'W' :{
        "path" : '/scratchnvme/wmass/WJetsNoCUT_v2/'
    },
    }

for i in ['0','1'] :
    for j in range(1,9) :
        samples['WJetsToLNu']['dir'].append(samples['W']['path']+"/tree_"+str(i)+"_"+str(j)+".root")
samples['W']['tree'].append(samples['W']['path']+"tree_1_0.root")

for tr in samples['data']['dir']:
    samples['data']['tree'].append(path+tr+"/tree.root")

for samp in samples.item() :

# inFile = ROOT.TFile('QCD_MC_hadded.root')
# tree = inFile.Get("Events")
weight = 'float(puWeight*PrefireWeight*lumiweight)'
    
    rdfTree = RDF("Events", samples[samp]['tree'])
    
    rdfTree = rdfTree.Define("Mu1_eta", "Muon_eta[Idx_mu1]")\
        .Define("Mu1_phi", "Muon_phi[Idx_mu1]")\
        .Define("Mu1_relIso", "Muon_pfRelIso04_all[Idx_mu1]")\
        .Define("Mu1_dz", "Muon_dz[Idx_mu1]")\
        .Define("Mu1_pt", "Muon_corrected_pt[Idx_mu1]")\
        .Define("Mu1_sip3d", "Muon_sip3d[Idx_mu1]")\
        .Define("Mu1_dxy", "Muon_dxy[Idx_mu1]")\
        .Define("MT","TMath::Sqrt(2*Mu1_pt*MET_pt_nom*(1.0-TMath::Cos(Mu1_phi-MET_phi_nom)))")
    
    rdfTree = rdfTree.Filter("Idx_mu1>-1", "1 muon")\
        .Filter("MET_filters==1", "MET Filter")\
        .Filter("nVetoElectrons==0 && 1", "veto el")\
        .Filter("MT>40.", "MT>40")\
        .Filter("Vtype==0", "Vtype 0")\
        .Filter("HLT_SingleMu24" , "HLT sig")
    
    
