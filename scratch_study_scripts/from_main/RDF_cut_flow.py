import ROOT
import math

######################## DATA SELECTION EFFICIENCY ################################
print("----------------")

RDF = ROOT.ROOT.RDataFrame
ROOT.EnableImplicitMT(128)

wdict = {
        "xsec": 61526.7, 
        "dir": [],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }

pathW = '/scratchnvme/wmass/WJetsNoCUT_v2/'
for i in ['0','1'] :
    for j in range(1,9) :
        wdict['dir'].append(pathW+"/tree_"+str(i)+"_"+str(j)+".root")
wdict['dir'].append(pathW+"tree_1_0.root")


rdfTree_data = RDF("Events", wdict['dir'])

rdfTree_data = rdfTree_data.Define("Mu1_eta", "Muon_eta[Idx_mu1]")\
    .Define("Mu1_phi", "Muon_phi[Idx_mu1]")\
    .Define("Mu1_relIso", "Muon_pfRelIso04_all[Idx_mu1]")\
    .Define("Mu1_dz", "Muon_dz[Idx_mu1]")\
    .Define("Mu1_pt", "Muon_corrected_pt[Idx_mu1]")\
    .Define("Mu1_sip3d", "Muon_sip3d[Idx_mu1]")\
    .Define("Mu1_dxy", "Muon_dxy[Idx_mu1]")\
    .Define("MT","TMath::Sqrt(2*Mu1_pt*MET_pt_nom*(1.0-TMath::Cos(Mu1_phi-MET_phi_nom)))")
    
rdfTree_data = rdfTree_data.Filter("Idx_mu1>-1", "1 muon")\
    .Filter("MET_filters==1", "MET Filter")\
    .Filter("nVetoElectrons==0 && 1", "veto el")\
    .Filter("HLT_SingleMu24" , "HLT sig")\


A = rdfTree_data.Filter("MT<40. && Vtype==1", "sideband aiso")
B = rdfTree_data.Filter("MT>40. && Vtype==1", "signal aiso")
C = rdfTree_data.Filter("MT<40. && Vtype==0", "sideband")
D = rdfTree_data.Filter("MT>40. && Vtype==0", "signal")

print("...............")  
print("report A")
rep_A = A.Report()
rep_A.Print() 
print("...............")        
print("report B")
rep_B = B.Report()
rep_B.Print()
print("...............")     
print("report C")
rep_C = C.Report()
rep_C.Print() 
print("...............")  
print("report D")
rep_D = D.Report()
rep_D.Print()   