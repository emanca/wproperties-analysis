import os
import sys
import ROOT
import json
import argparse
import subprocess
sys.path.append('../RDFprocessor/framework')
from RDFtree import RDFtree
sys.path.append('python/')
sys.path.append('data/')
from systematics import systematics
from selections import selections, selectionVars, selections_bkg

from getLumiWeight import getLumiWeight

ROOT.gSystem.Load('bin/libAnalysisOnData.so')
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")

def RDFprocessWJetsMCSignalACtempl(fvec, outputDir, sample, xsec, fileSF, fileScale, ncores, pretendJob):
    ROOT.ROOT.EnableImplicitMT(ncores)
    print("running with {} cores for sample:{}".format(ncores, sample)) 
    wdecayselections = { 
        'WToMu'  : ' && (abs(genVtype) == 14)',
        'WToTau' : ' && (abs(genVtype) == 16)'
    }
    filePt = ROOT.TFile.Open("data/histoUnfoldingSystPt_nsel2_dy3_rebin1_default.root")
    fileY = ROOT.TFile.Open("data/histoUnfoldingSystRap_nsel2_dy3_rebin1_default.root")
    fileACplus = ROOT.TFile.Open("../analysisOnGen/genInput_Wplus.root")
    fileACminus = ROOT.TFile.Open("../analysisOnGen/genInput_Wminus.root")
    # fileACplus = ROOT.TFile.Open("../analysisOnGen/genInput_4GeVBinning_udShift_Wminus.root")
    # fileACminus = ROOT.TFile.Open("../analysisOnGen/genInput_4GeVBinning_udShift_Wplus.root")
    #Will only ber run on signal region
    region = 'Signal'
    weight = 'float(puWeight*PrefireWeight*lumiweight*WHSF*weightPt)'
    wtomu_cut =  selections_bkg[region] + wdecayselections['WToMu']
    # wtomu_lowAcc_cut = wtomu_cut + "&& (Wpt_preFSR>32. || Wrap_preFSR_abs>2.4)"
    wtomu_lowAcc_cut = wtomu_cut + "&& (Wpt_preFSR>60. || Wrap_preFSR_abs>2.4)"
    #print weight, "NOMINAL WEIGHT"
    nom = ROOT.vector('string')()
    nom.push_back("")

    p = RDFtree(outputDir = outputDir, inputFile = fvec, outputFile="{}_AC_plots.root".format('WToMu'), pretend=pretendJob)
    p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.reweightFromZ(filePt,fileY),ROOT.baseDefinitions(True,True),ROOT.rochesterVariations(fileScale),ROOT.weightDefinitions(fileSF),getLumiWeight(xsec=xsec, inputFile=fvec, genEvsbranch = "genEventSumw"),ROOT.Replica2Hessian()])

    steps = [ROOT.getACValues(fileACplus,fileACminus), ROOT.defineHarmonics(), ROOT.getMassWeights(), ROOT.getWeights()]
    p.branch(nodeToStart = 'defs', nodeToEnd = 'defsAC', modules = steps)
    #reco templates with AC reweighting
    p.branch(nodeToStart = 'defsAC', nodeToEnd = 'templatesAC_{}/Nominal'.format(region), modules = [ROOT.templateBuilder(wtomu_cut, weight,nom,"Nom",0)])
    #reco templates for out of acceptance events
    p.branch(nodeToStart = 'defsAC', nodeToEnd = 'templatesLowAcc_{}/Nominal'.format(region), modules = [ROOT.templates(wtomu_lowAcc_cut, weight, nom,"Nom",0)])
    # mass = ROOT.vector('string')()
    # mass.push_back("_massUp")
    # mass.push_back("_massDown")
    # p.branch(nodeToStart = 'defsAC', nodeToEnd = 'templatesLowAcc_{}/Nominal'.format(region), modules = [ROOT.templates(wtomu_lowAcc_cut, weight, mass,"massWeights",0)])
    #weight variations
    for s,variations in systematics.items():
        # if '_WQT' in s: continue  #not needed for the fit
        if 'WHSF' in s and 'aiso' in region: continue
        print("branching weight variations", s)
        if '_WQT' in s:# Wqt syst using scale
            wtomu_modcut = wtomu_cut + variations[2]
            wtomu_lowAcc_modcut = wtomu_lowAcc_cut + variations[2]
        else :
            wtomu_modcut = wtomu_cut
            wtomu_lowAcc_modcut = wtomu_lowAcc_cut
        var_weight = weight.replace(s, "1.") #do something only if s is in the weight 
        vars_vec = ROOT.vector('string')()
        for var in variations[0]:
            if '_WQT' in var :
                var = var.replace("_WQTlow","")
                var = var.replace("_WQTmid","")
                var = var.replace("_WQThigh","")
            vars_vec.push_back(var)
            # print(weight,"\t",var_weight, "MODIFIED WEIGHT")
        #reco templates with AC reweighting
        nodeToStartName = 'defsAC'
        # if s=='LHEScaleWeight' or s=='LHEPdfWeight' or s=='alphaS':
        if s=='LHEScaleWeight' or s=='LHEPdfWeight' or s=='alphaS' or s=='mass':
            steps = [ROOT.getACValues(fileACplus,fileACminus,s, variations[0]), ROOT.defineHarmonics(), ROOT.getMassWeights(), ROOT.getWeights(s, variations[0])]
            p.branch(nodeToStart = 'defs', nodeToEnd = 'defsAC{}'.format(s), modules = steps)
            nodeToStartName+=s
            p.branch(nodeToStart = nodeToStartName, nodeToEnd = 'templatesAC_{}/{}'.format(region, s), modules = [ROOT.templateBuilder(wtomu_modcut, var_weight,vars_vec,variations[1], 4)])
        else :
            p.branch(nodeToStart = nodeToStartName, nodeToEnd = 'templatesAC_{}/{}'.format(region, s), modules = [ROOT.templateBuilder(wtomu_modcut, var_weight,vars_vec,variations[1], 3)])
        #reco templates for out of acceptance events
        p.branch(nodeToStart = nodeToStartName, nodeToEnd = 'templatesLowAcc_{}/{}'.format(region,s), modules = [ROOT.templates(wtomu_lowAcc_modcut, var_weight,vars_vec,variations[1], 0)])
    #column variations#weight will be nominal, cut will vary
    for vartype, vardict in selectionVars.items():
        wtomu_cut_vec = ROOT.vector('string')()
        wtomu_var_vec = ROOT.vector('string')()
        for selvar, hcat in vardict.items() :
            wtomu_newcut = wtomu_cut.replace('MT', 'MT'+selvar)
            if 'corrected' in selvar:
                wtomu_newcut = wtomu_newcut.replace('Mu1_pt', 'Mu1_pt'+selvar)
            wtomu_cut_vec.push_back(wtomu_newcut)
            wtomu_var_vec.push_back(selvar)
        print("branching column variations:", vartype, " for region:", region) #, "\tvariations:", wtomu_var_vec
        #reco templates with AC reweighting
        p.branch(nodeToStart = 'defsAC', nodeToEnd = 'templatesAC_{}/{}'.format(region, vartype), modules = [ROOT.templateBuilder(wtomu_cut_vec, weight, nom,"Nom",hcat,wtomu_var_vec)])
        #reco templates for out of acceptance events
        print('low Acc')
        for cut in wtomu_cut_vec:
            # cut+= "&& (Wpt_preFSR>32. || Wrap_preFSR_abs>2.4)"
            cut+= "&& (Wpt_preFSR>60. || Wrap_preFSR_abs>2.4)"
            #print "Low acc cut vec vars:", wtomu_cut_vec
        p.branch(nodeToStart = 'defsAC', nodeToEnd = 'templatesLowAcc_{}/{}'.format(region,vartype), modules = [ROOT.templates(wtomu_cut_vec, weight, nom,"Nom",hcat,wtomu_var_vec)])
    p.getOutput()
    # p.saveGraph()


def RDFprocessWJetsMC(fvec, outputDir, sample, xsec, fileSF, fileScale, ncores, pretendJob, bkg,SBana=False):
    ROOT.ROOT.EnableImplicitMT(ncores)
    print("running with {} cores for sample:{}".format(ncores, sample)) 
    wdecayselections = { 
        'WToMu'  : ' && (abs(genVtype) == 14)',
        'WToTau' : ' && (abs(genVtype) == 16)'
    }
    p = RDFtree(outputDir = outputDir, inputFile = fvec, outputFile="{}_plots.root".format(sample), pretend=pretendJob)
    filePt = ROOT.TFile.Open("data/histoUnfoldingSystPt_nsel2_dy3_rebin1_default.root")
    fileY = ROOT.TFile.Open("data/histoUnfoldingSystRap_nsel2_dy3_rebin1_default.root")
    #fileAC = ROOT.TFile.Open("../analysisOnGen/genInput.root")
    p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.reweightFromZ(filePt,fileY),ROOT.baseDefinitions(True, True),ROOT.rochesterVariations(fileScale), ROOT.weightDefinitions(fileSF),getLumiWeight(xsec=xsec, inputFile=fvec, genEvsbranch = "genEventSumw"),ROOT.Replica2Hessian(),ROOT.getMassWeights()])
    for region,cut in selections_bkg.items():
        # cut = cut.replace("&& HLT_SingleMu24","&& IsTObjmatched_mu1")#do not use, IsTObjmatched_mu1 is bugged in wjets sample
        if 'aiso' in region:
            weight = 'float(puWeight*PrefireWeight*lumiweight*weightPt)'
        else:
            weight = 'float(puWeight*PrefireWeight*lumiweight*WHSF*weightPt)'
        # if 'aiso' in region:
        #     weight = 'float(puWeight*lumiweight)'
        # else:
        #     weight = 'float(puWeight*lumiweight*WHSF)'
        print(weight, "NOMINAL WEIGHT")
        nom = ROOT.vector('string')()
        nom.push_back("")
        wtomu_cut = cut + wdecayselections['WToMu']
        wtotau_cut = cut + wdecayselections['WToTau']
        print("Region=", region)
        print("WToMu cut=", wtomu_cut)
        print("WToTau cut=", wtotau_cut)
        #last argument refers to histo category - 0 = Nominal, 1 = Pt scale , 2 = MET scale
        print("branching nominal for region:", region) 
        #Nominal templates
        p.branch(nodeToStart = 'defs', nodeToEnd = '{}/templates_{}/Nominal'.format('WToMu', region), modules = [ROOT.templates(wtomu_cut, weight, nom,"Nom",0)])            
        p.branch(nodeToStart = 'defs', nodeToEnd = '{}/templates_{}/Nominal'.format('WToTau', region), modules = [ROOT.templates(wtotau_cut, weight, nom,"Nom",0)])
        if region == "Signal" or (region=='Sideband' and SBana):
            print("adding muon histo to graph for region {}".format(region))
            p.branch(nodeToStart = 'defs', nodeToEnd = '{}/prefit_{}/Nominal'.format('WToMu', region), modules = [ROOT.muonHistos(wtomu_cut, weight, nom,"Nom",0)])     
            p.branch(nodeToStart = 'defs', nodeToEnd = '{}/prefit_{}/Nominal'.format('WToTau', region), modules = [ROOT.muonHistos(wtotau_cut, weight, nom,"Nom",0)])

        #weight variations
        for s,variations in systematics.items():
            if 'WHSF' in s and 'aiso' in region: continue
            print("branching weight variations", s)
            #if "LHEScaleWeight" in s and samples[sample]['systematics'] != 2 :  continue
            if s=='LHEScaleWeight' : continue #no correlated MCscale for W
            
            if '_WQT' in s:# Wqt syst using scale
                wtomu_modcut = wtomu_cut + variations[2]
                wtotau_modcut = wtotau_cut + variations[2]
            else :
                wtomu_modcut = wtomu_cut
                wtotau_modcut = wtotau_cut
            var_weight = weight.replace(s, "1.") #do something only if s is in the weight 
            vars_vec = ROOT.vector('string')()
            for var in variations[0]:
                if '_WQT' in var :
                    var = var.replace("_WQTlow","")
                    var = var.replace("_WQTmid","")
                    var = var.replace("_WQThigh","")
                vars_vec.push_back(var)
            print(weight,"\t",var_weight, "MODIFIED WEIGHT")
            #Template vars
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = '{}/templates_{}/{}'.format('WToMu', region,s), modules = [ROOT.templates(wtomu_modcut, var_weight,vars_vec,variations[1], 0)])
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = '{}/templates_{}/{}'.format('WToTau', region,s), modules = [ROOT.templates(wtotau_modcut, var_weight,vars_vec,variations[1], 0)])
            if region == "Signal" or (region=='Sideband' and SBana):
                p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = '{}/prefit_{}/{}'.format('WToMu', region,s), modules = [ROOT.muonHistos(wtomu_modcut, var_weight,vars_vec,variations[1], 0)])
                p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = '{}/prefit_{}/{}'.format('WToTau', region,s), modules = [ROOT.muonHistos(wtotau_modcut, var_weight,vars_vec,variations[1], 0)])
                
        #column variations#weight will be nominal, cut will vary
        for vartype, vardict in selectionVars.items():
            wtomu_cut_vec = ROOT.vector('string')()
            wtomu_var_vec = ROOT.vector('string')()
            wtotau_cut_vec = ROOT.vector('string')()
            wtotau_var_vec = ROOT.vector('string')()
            for selvar, hcat in vardict.items() :
                wtomu_newcut = wtomu_cut.replace('MT', 'MT'+selvar)
                wtotau_newcut = wtotau_cut.replace('MT', 'MT'+selvar)
                if 'corrected' in selvar:
                    wtomu_newcut = wtomu_newcut.replace('Mu1_pt', 'Mu1_pt'+selvar)
                    wtotau_newcut = wtotau_newcut.replace('Mu1_pt', 'Mu1_pt'+selvar)
                wtomu_cut_vec.push_back(wtomu_newcut)
                wtomu_var_vec.push_back(selvar)
                wtotau_cut_vec.push_back(wtotau_newcut)
                wtotau_var_vec.push_back(selvar)
            print("branching column variations:", vartype, " for region:", region) #, "\tvariations:", wtomu_var_vec
            #templates (integrated over helicity xsecs)
            p.branch(nodeToStart = 'defs', nodeToEnd = '{}/templates_{}/{}'.format('WToMu', region,vartype), modules = [ROOT.templates(wtomu_cut_vec, weight, nom,"Nom",hcat,wtomu_var_vec)])  
            p.branch(nodeToStart = 'defs', nodeToEnd = '{}/templates_{}/{}'.format('WToTau', region,vartype), modules = [ROOT.templates(wtotau_cut_vec, weight, nom,"Nom",hcat,wtotau_var_vec)])  
            if region == "Signal" or (region=='Sideband' and SBana):
                p.branch(nodeToStart = 'defs', nodeToEnd = '{}/prefit_{}/{}'.format('WToMu', region,vartype), modules = [ROOT.muonHistos(wtomu_cut_vec, weight, nom,"Nom",hcat,wtomu_var_vec)])  
                p.branch(nodeToStart = 'defs', nodeToEnd = '{}/prefit_{}/{}'.format('WToTau', region,vartype), modules = [ROOT.muonHistos(wtotau_cut_vec, weight, nom,"Nom",hcat,wtotau_var_vec)])
   
    p.getOutput()
    # p.saveGraph()
    #if not pretendJob:
    #split the output file into decay modes
    os.chdir(outputDir)
    os.system("ls -ltrh")
    os.system("root -b -q ../python/splitWJets.C")
    os.chdir("../")




def main():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-p', '--pretend',type=int, default=False, help="run over a small number of event")
    parser.add_argument('-i', '--inputDir',type=str, default="/scratchnvme/wmass/NanoAOD2016-V2/", help="input dir name")
    parser.add_argument('-o', '--outputDir',type=str, default='./output/', help="output dir name")
    parser.add_argument('-c', '--ncores',type=int, default=128, help="number of cores used")
    parser.add_argument('-sb', '--SBana', type=int, default=False, help="run also on the sideband (clousure test)")
    parser.add_argument('-b', '--runBKG', type=int, default=False, help="get histograms for bkg analysis")

    args = parser.parse_args()
    pretendJob = args.pretend
    inDir = args.inputDir
    outputDir = args.outputDir
    ncores = args.ncores
    SBana = args.SBana
    bkg = args.runBKG
    print(bkg, '\t', pretendJob, '\t', SBana)
    
    if pretendJob:
        print("Running a test job over a few events")
    else:
        print("Running on full dataset")
    samples={}
    with open('data/samples_2016.json') as f:
        samples = json.load(f)
    sample='WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'
    print(sample)
    direc = samples[sample]['dir']
    xsec = samples[sample]['xsec']
    fvec=ROOT.vector('string')()
    for dirname,fname in direc.items():
        ##check if file exists or not
        #inputFile = '/scratchnvme/emanca/wproperties-analysis/analysisOnGen/test_tree_*.root'
        inputFile = '/scratchnvme/wmass/WJetsNoCUT_v2/tree_*_*.root'
        #inputFile = '{}/{}/tree.root'.format(inDir, dirname)
        isFile = os.path.isfile(inputFile)  
        if not isFile:
            print(inputFile, " does not exist")
            #continue
        fvec.push_back(inputFile)
        break
    if fvec.empty():
        print("No files found for directory:", samples[sample], " SKIPPING processing")
        sys.exit(1)
    print(fvec) 
    
    fileSF = ROOT.TFile.Open("data/ScaleFactors_OnTheFly.root")
    fileScale = ROOT.TFile.Open("data/muscales_extended.root")
    if bkg: 
        RDFprocessWJetsMC(fvec, outputDir, sample, xsec, fileSF, fileScale, ncores, pretendJob, bkg,SBana)
    else:
        RDFprocessWJetsMCSignalACtempl(fvec, outputDir, sample, xsec, fileSF, fileScale, ncores, pretendJob)

if __name__ == "__main__":
    main()
