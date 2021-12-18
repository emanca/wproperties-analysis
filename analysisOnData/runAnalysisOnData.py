import os
import sys
import ROOT
import json
import argparse
sys.path.append('../RDFprocessor/framework')
from RDFtree import RDFtree
sys.path.append('python/')
sys.path.append('data/')
from systematics import systematics, systFlat
from selections import selections, selections_bkg, selections_fakes, selectionVars
from getLumiWeight import getLumiWeight
ROOT.gSystem.Load('bin/libAnalysisOnData.so')

#produces templates for all regions and prefit for signal
def RDFprocessData(fvec, outputDir, ncores, pretendJob=True, SBana=False, outF="SingleMuonData_plots.root"):
    ROOT.ROOT.EnableImplicitMT(ncores)
    print("running with {} cores".format(ncores))
    weight = 'float(1)'
    outF="SingleMuonData_plots.root"
    p = RDFtree(outputDir = outputDir, inputFile = fvec, outputFile=outF, pretend=pretendJob)
    p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.baseDefinitions(0)])
    for region,cut in selections_bkg.items():    
        print(region)       
        nom = ROOT.vector('string')()
        nom.push_back("")
        #last argument refers to histo category - 0 = Nominal, 1 = Pt scale , 2 = MET scale
        print("branching nominal")
        if region == "Signal" or (region=='Sideband' and SBana):
           p.branch(nodeToStart = 'defs', nodeToEnd = 'prefit_{}/Nominal'.format(region), modules = [ROOT.muonHistos(cut, weight, nom,"Nom",0)]) 
        #nominal templates
        p.branch(nodeToStart = 'defs', nodeToEnd = 'templates_{}/Nominal'.format(region), modules = [ROOT.templates(cut, weight, nom,"Nom",0)])       
    p.getOutput()
    # p.saveGraph()

#produces Fake contribution to prefit plots computed from data 
def RDFprocessfakefromData(fvec, outputDir, bkgFile, ncores, pretendJob=True, SBana=False, outF="FakeFromData_plots.root"):
    ROOT.ROOT.EnableImplicitMT(ncores)
    print("running with {} cores".format(ncores))
    weight = 'float(1)'
    #in case we want pdf variations for fakes
    #systematics.update({ "LHEPdfWeight" : ( ["_LHEPdfWeight" + str(i)  for i in range(0, 100)], "LHEPdfWeight" ) } )
    #print systematics
    p = RDFtree(outputDir = outputDir, inputFile = fvec, outputFile=outF, pretend=pretendJob)
    for region,cut in selections_fakes.items(): 
        if 'SideBand' in region and (not SBana) : continue 
        if 'SideBand' in region : bkgFile_mod = bkgFile.replace('.root','SideBand.root') 
        else :  bkgFile_mod=bkgFile
        FR=ROOT.TFile.Open(bkgFile_mod)
        p.branch(nodeToStart = 'input', nodeToEnd = 'defs', modules = [ROOT.baseDefinitions(0),ROOT.fakeRate(FR)]) 
        print(region)       
        nom = ROOT.vector('string')()
        nom.push_back("")
        weight = "float(fakeRate_Nominal_)"
        #last argument refers to histo category - 0 = Nominal, 1 = Pt scale , 2 = MET scale
        print("branching nominal")
        p.branch(nodeToStart = 'defs', nodeToEnd = 'prefit_{}/Nominal'.format(region), modules = [ROOT.muonHistos(cut, weight, nom,"Nom",0)]) 
        p.branch(nodeToStart = 'defs', nodeToEnd = 'templates_{}/Nominal'.format(region), modules = [ROOT.templates(cut, weight, nom,"Nom",0)])       

        #now add fake variations
        for s,variations in systematics.items():
            if 'Nom_WQT' in s : continue
            if s=='LHEScaleWeight' : continue
            #only required systs
            # if "LHEPdfWeight" not in s and "LHEScaleWeight" not in s: continue
            print("branching weight variations", s)
            vars_vec = ROOT.vector('string')()
            for var in variations[0]:
                vars_vec.push_back(var)
                weight="float(1)"
            print("fakeRate_"+variations[1])
            print(vars_vec)
            if "WQT" in s :
                variationMod = variations[1]+s.replace('LHEScaleWeight','')
            else :
                variationMod = variations[1]
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = 'prefit_{}/{}'.format(region,s), modules = [ROOT.muonHistos(cut,weight,vars_vec,"fakeRate_"+variationMod, 0)])
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = 'templates_{}/{}'.format(region,s), modules = [ROOT.templates(cut,weight,vars_vec,"fakeRate_"+variationMod, 0)])
        
        #fake column variations since the cut won't change in data
        for vartype, vardict in selectionVars.items():
            # if vartype != 'jme' : continue
            vars_vec = ROOT.vector('string')()
            for selvar, hcat in vardict.items() :
                vars_vec.push_back(selvar)
            print("branching fake column variations", vartype)
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = 'prefit_{}/{}'.format(region,vartype), modules = [ROOT.muonHistos(cut,weight,vars_vec,"fakeRate_"+vartype, 0)])
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = 'templates_{}/{}'.format(region,vartype), modules = [ROOT.templates(cut,weight,vars_vec,"fakeRate_"+vartype, 0)])
        
        for s,variations in systFlat.items():
            print("branching weight variations", s)
            vars_vec = ROOT.vector('string')()
            for var in variations[0]:
                vars_vec.push_back(var)
                weight="float(1)"
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = 'prefit_{}/{}'.format(region,s), modules = [ROOT.muonHistos(cut,weight,vars_vec,"fakeRate_"+variations[1], 0)])
            p.branch(nodeToStart = 'defs'.format(region), nodeToEnd = 'templates_{}/{}'.format(region,s), modules = [ROOT.templates(cut,weight,vars_vec,"fakeRate_"+variations[1], 0)]) 

    #save output
    p.getOutput()
    # p.saveGraph()

def main():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-p', '--pretend',type=int, default=False, help="run over a small number of event")
    parser.add_argument('-i', '--inputDir',type=str, default='/scratchnvme/wmass/NanoAOD2016-V2/', help="input dir name")
    parser.add_argument('-o', '--outputDir',type=str, default='./output/', help="output dir name")
    parser.add_argument('-c', '--ncores',type=int, default=128, help="number of cores used")
    parser.add_argument('-sb', '--SBana',type=int, default=False, help="run also on the sideband (clousure test)")
    parser.add_argument('-b', '--runBKG',type=int, default=False, help="prepare the input of the bkg analysis, if =false run the prefit Plots")
    parser.add_argument('-f', '--bkgFile',type=str, default='/scratch/bertacch/wmass/wproperties-analysis/bkgAnalysis/TEST_runTheMatrix/bkg_parameters_CFstatAna.root', help="bkg parameters file path/name.root")
    args = parser.parse_args()
    pretendJob = args.pretend
    inDir = args.inputDir
    outputDir = args.outputDir
    ncores = args.ncores
    SBana = args.SBana
    runBKG = args.runBKG
    bkgFile = args.bkgFile

    if pretendJob:
        print("Running a test job over a few events")
    else:
        print("Running on full dataset")
    fvec=ROOT.vector('string')()
    samples={}
    with open('data/samples_2016.json') as f:
        samples = json.load(f)
    for sample in samples:
        if not samples[sample]['datatype']=='DATA': continue
        direc = samples[sample]['dir']
        for dirname,fname in direc.items():
            ##check if file exists or not
            inputFile = '{}/{}/tree.root'.format(inDir, dirname)
            isFile = os.path.isfile(inputFile)  
            if not isFile:
                print(inputFile, " does not exist")
                continue
            fvec.push_back(inputFile)
        
    if fvec.empty():
        print("No files found for json provided\n")
        sys.exit(1)
    print(fvec)
    if runBKG : #produces templates for all regions and prefit for signal
        RDFprocessData(fvec, outputDir, ncores, pretendJob,SBana)
    else : #produces Fake contribution to prefit plots computed from data 
        RDFprocessfakefromData(fvec, outputDir, bkgFile, ncores, pretendJob,SBana)

if __name__ == "__main__":
    main()
