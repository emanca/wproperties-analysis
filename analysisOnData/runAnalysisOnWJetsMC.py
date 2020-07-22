import os
import sys
import ROOT
import json

from RDFtree import RDFtree
sys.path.append('python/')
sys.path.append('data/')
from systematics import systematics
from selections import selections_whelicity, selectionVars, selections_bkg_whelicity

from getLumiWeight import getLumiWeight

ROOT.gSystem.Load('bin/libAnalysisOnData.so')

pretendJob = True if int(sys.argv[2]) == 1 else False
if pretendJob:
    print "Running a test job over a few events"
else:
    print "Running on full dataset"

runBKG = True if int(sys.argv[1]) == 1 else False 
outFtag=""
if runBKG:
    selections_whelicity = selections_bkg_whelicity
    print "Running job for preparing inputs of background study"
    outFtag="_bkgselections"

samples={}
with open('data/samples_2016.json') as f:
  samples = json.load(f)

sample='WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'
print sample
direc = samples[sample]['dir']
xsec = samples[sample]['xsec']
c = 64		
ROOT.ROOT.EnableImplicitMT(c)
print "running with {} cores".format(c)
  
fvec=ROOT.vector('string')()
for dirname,fname in direc.iteritems():
    ##check if file exists or not
    inputFile = '/scratchssd/sroychow/NanoAOD2016-V2/{}/tree.root'.format(dirname)
    isFile = os.path.isfile(inputFile)  
    if not isFile:
        print inputFile, " does not exist"
        continue
    fvec.push_back(inputFile)
    print inputFile
    WJets = False
    if "WJetsToLNu" in inputFile: WJets = True
if fvec.empty():
    print "No files found for directory:", samples[sample], " SKIPPING processing"
    sys.exit(1)
print fvec 

fileSF = ROOT.TFile.Open("data/ScaleFactors_OnTheFly.root")

wdecayselections = { 
                     'WToMu'  : ' && (abs(genVtype) == 14)',
                     'WToTau' : ' && (abs(genVtype) == 16)'
                    }

for wdecay, decaycut in wdecayselections.iteritems() :
    print "Running for Wdecay:", wdecay

    p = RDFtree(outputDir = './output/', inputFile = fvec, outputFile="{}{}_plots.root".format(wdecay, outFtag), pretend=pretendJob)
    filePt = ROOT.TFile.Open("data/histoUnfoldingSystPt_nsel2_dy3_rebin1_default.root")
    fileY = ROOT.TFile.Open("data/histoUnfoldingSystRap_nsel2_dy3_rebin1_default.root")
    p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.reweightFromZ(filePt,fileY),ROOT.baseDefinitions(),ROOT.weightDefinitions(fileSF),getLumiWeight(xsec=xsec, inputFile=fvec)])

    for region,cut in selections_whelicity.iteritems():
        
        if 'aiso' in region:
            weight = 'float(puWeight*lumiweight*weightPt*weightY)'
        else:
            weight = 'float(puWeight*lumiweight*WHSF*weightPt*weightY)'
        
        print weight, "NOMINAL WEIGHT"
        
        nom = ROOT.vector('string')()
        nom.push_back("")
        #last argument refers to histo category - 0 = Nominal, 1 = Pt scale , 2 = MET scale
        print "branching nominal for region:", region 
        cut += decaycut
        if not runBKG:             
            p.branch(nodeToStart = 'defs', nodeToEnd = 'prefit_{}/Nominal'.format(region), modules = [ROOT.muonHistos(cut, weight, nom,"Nom",0)])     
        p.branch(nodeToStart = 'defs', nodeToEnd = 'templates_{}/Nominal'.format(region), modules = [ROOT.templates(cut, weight, nom,"Nom",0)])            
   
        #weight variations
        for s,variations in systematics.iteritems():
            print "branching weight variations", s
            if "LHEScaleWeight" in s and samples[sample]['systematics'] != 2 :  continue
            if not "LHEScaleWeight" in s:
                if 'aiso' in region: continue
                var_weight = weight.replace(s, "1.")
            else: 
                var_weight = weight

            vars_vec = ROOT.vector('string')()
            for var in variations[0]:
                vars_vec.push_back(var)
                
            print weight,var_weight, "MODIFIED WEIGHT"
                
            if not runBKG: 
                p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = 'prefit_{}/{}Vars'.format(region,s), modules = [ROOT.muonHistos(cut, var_weight,vars_vec,variations[1], 0)])
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = 'templates_{}/{}Vars'.format(region,s), modules = [ROOT.templates(cut, var_weight,vars_vec,variations[1], 0)])
                    

        #column variations#weight will be nominal, cut will vary
        for vartype, vardict in selectionVars.iteritems():
            cut_vec = ROOT.vector('string')()
            var_vec = ROOT.vector('string')()
            for selvar, hcat in vardict.iteritems() :
                newcut = cut.replace('MT', 'MT'+selvar)
                    
                if 'corrected' in selvar:
                    newcut = newcut.replace('Mu1_pt', 'Mu1_pt'+selvar)
                    
                cut_vec.push_back(newcut)
                var_vec.push_back(selvar)
            print "branching column variations:", vartype, " for region:", region, "\tvariations:", var_vec
            if not runBKG: 
                p.branch(nodeToStart = 'defs', nodeToEnd = 'prefit_{}/{}Vars'.format(region,vartype), modules = [ROOT.muonHistos(cut_vec, weight, nom,"Nom",hcat,var_vec)])  
            p.branch(nodeToStart = 'defs', nodeToEnd = 'templates_{}/{}Vars'.format(region,vartype), modules = [ROOT.templates(cut_vec, weight, nom,"Nom",hcat,var_vec)])  

    p.getOutput()
    p.saveGraph()
