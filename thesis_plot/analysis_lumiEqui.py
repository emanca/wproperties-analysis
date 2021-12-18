import ROOT
import math


samples = {

    "TTJets_SingleLeptFromT": {
        "xsec": 182.0, 
        "dir": ["TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", "TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "ZZ": {
        "xsec": 16.523, 
        "dir": ["ZZ_TuneCUETP8M1_13TeV-pythia8"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "TTJets_DiLept": {
        "xsec": 95.02, 
        "dir": ["TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", "TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "ST_tW_antitop_5f_inclusiveDecays": {
        "xsec": 35.6, 
        "dir": ["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_ext1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "DYJetsToLL_M10to50": {
        "xsec": 1903.0, 
        "dir": ["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8","DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "ST_t-channel_antitop_4f_inclusiveDecays": {
        "xsec": 80.95, 
        "dir": ["ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "TTJets_SingleLeptFromTbar": {
        "xsec": 182.0, 
        "dir": ["TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "WZ": {
        "xsec": 47.13, 
        "dir": ["WZ_TuneCUETP8M1_13TeV-pythia8","WZ_TuneCUETP8M1_13TeV-pythia8_ext1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "ST_t-channel_top_4f_inclusiveDecays_13TeV": {
        "xsec": 136.02, 
        "dir": ["ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "ST_tW_top_5f_inclusiveDecays": {
        "xsec": 35.6, 
        "dir": ["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_ext1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "DYJetsToLL_M50": {
        "xsec": 6025.2, 
        "dir": ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8", "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext2"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "ST_s-channel_4f_leptonDecays": {
        "xsec": 3.68, 
        "dir": ["ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "WW": {
        "xsec": 115.0, 
        "dir": ["WW_TuneCUETP8M1_13TeV-pythia8", "WW_TuneCUETP8M1_13TeV-pythia8_ext1"],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }, 
    "WJetsToLNu": {
        "xsec": 61526.7, 
        "dir": [],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }
}


#---------- analysis without considering the events weight -------------------


# pathW = '/scratchnvme/wmass/WJetsNoCUT_v2/'
# for i in ['0','1'] :
#     for j in range(1,9) :
#         samples['WJetsToLNu']['dir'].append(pathW+"/tree_"+str(i)+"_"+str(j)+".root")
# samples['WJetsToLNu']['dir'].append(pathW+"tree_1_0.root")



# path = '/scratchnvme/wmass/NanoAOD2016-V2/'
# for s, samp in samples.items() :
#     print(s)
#     for tr in samp['dir'] :
#         if s!="WJetsToLNu" :
#             inFile = ROOT.TFile.Open(path+tr+"/tree.root")
#         else :
#             inFile = ROOT.TFile.Open(tr)
#         tree = inFile.Get("Events")
#         samples[s]["Nevents"].append(tree.GetEntries())
#     lumieq = float(sum(samples[s]["Nevents"]))/ samples[s]["xsec"]
#     samples[s]["lumiEqui"].append(lumieq)
    
# print("------------results------------")
# for s, samp in samples.items() :
#     print(s, "      L=",samples[s]["lumiEqui"][0]/1000, "fb^-1", "    ","Nevents=", float(sum(samples[s]["Nevents"])) ) 
    
    
########################################################################################################


###############################################  considering the event weight   ###############################################
RDF = ROOT.ROOT.RDataFrame

pathW = '/scratchnvme/wmass/WJetsNoCUT_v2/'
for i in ['0','1'] :
    for j in range(1,9) :
        samples['WJetsToLNu']['dir'].append(pathW+"/tree_"+str(i)+"_"+str(j)+".root")
samples['WJetsToLNu']['dir'].append(pathW+"tree_1_0.root")



path = '/scratchnvme/wmass/NanoAOD2016-V2/'
for s, samp in samples.items() :
    print("----------------")
    print(s)
    for tr in samp['dir'] :
        if s!="WJetsToLNu" :
            inFile = ROOT.TFile.Open(path+tr+"/tree.root")
        else :
            inFile = ROOT.TFile.Open(tr)
        tree = inFile.Get("Events")
        print(tr)
        samples[s]["Nevents"].append(tree.GetEntries())
        
        runs = RDF('Runs', inFile)
        
        if s!="WJetsToLNu" :
            runs = runs.Define("effEntries", "genEventSumw_*genEventSumw_/genEventSumw2_")
            effectiveEntries = runs.Sum("effEntries").GetValue()
        else :
            runs = runs.Define("effEntries", "genEventSumw*genEventSumw/genEventSumw2")
            effectiveEntries = runs.Sum("effEntries").GetValue()
        # print("effective entries=",effectiveEntries)
        samples[s]["Nevents_w"].append(effectiveEntries)
                
    lumieq = float(sum(samples[s]["Nevents"]))/ (samples[s]["xsec"]*1000) #sigma[fb] = 10^3sigma[pb]
    samples[s]["lumiEqui"].append(lumieq)
    
    lumieq_w = float(sum(samples[s]["Nevents_w"]))/ (samples[s]["xsec"]*1000) 
    samples[s]["lumiEqui"].append(lumieq_w)
        
print("------------results------------")
for s, samp in samples.items() :
    print(s, "      L=",samples[s]["lumiEqui"][0], "fb^-1", "    ","Nevents=", float(sum(samples[s]["Nevents"])) ) 
    print(s, "     Lw=",samples[s]["lumiEqui"][1], "fb^-1", "    ","Nevents=", float(sum(samples[s]["Nevents_w"])) ) 
    

####################### EFFICIENCY OF SELECTION #######################

ROOT.EnableImplicitMT(128)

if (1) :
    for s, samp in samples.items() :
        print("----------------")
        print(s)
        
        if s!="WJetsToLNu" : 
            treeList = []
            for tr in samp['dir'] :
                treeList.append(path+tr+"/tree.root")
        else :
            treeList =  samp['dir']
        rdfTree = RDF("Events", treeList)

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
        
        # rdfTree = defTree.Define("genW", "genEventSumw_*genEventSumw_/genEventSumw2_")
        
        rdfTree = rdfTree.Define("gen2", "Generator_weight*Generator_weight")

        selEffEntries = (rdfTree.Sum("Generator_weight").GetValue())**2/ (rdfTree.Sum("gen2").GetValue())
        selEffLumi = float(selEffEntries)/ (samples[s]["xsec"]*1000) 
        print("Selection Efficiency=", selEffLumi/samples[s]["lumiEqui"][1], ", selectedLumi=", selEffLumi, "fb^{-1}")
            



######################## DATA SELECTION EFFICIENCY ################################
print("----------------")
print('data cut flow')

dataDict = {
        "xsec": 1, 
        "dir": ["SingleMuon_Run2016B_ver2","SingleMuon_Run2016C","SingleMuon_Run2016D","SingleMuon_Run2016E","SingleMuon_Run2016F","SingleMuon_Run2016G","SingleMuon_Run2016H",],
        "lumiEqui" : [],
        "Nevents" : [],
        "Nevents_w" : []
    }

treeList = []
for tr in dataDict['dir'] :
    treeList.append(path+tr+"/tree.root")
rdfTree_data = RDF("Events", treeList)

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
    .Filter("MT>40.", "MT>40")\
    .Filter("Vtype==0", "Vtype 0")\
    .Filter("HLT_SingleMu24" , "HLT sig")

print("report")
dataCutReport = rdfTree_data.Report()
dataCutReport.Print()       