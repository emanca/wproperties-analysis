import os
import sys
import ROOT
import argparse
import copy
import time
from datetime import datetime
#set all paths from FWK directory
FWKBASE=os.getenv('FWK_BASE')
sys.path.append('{}/RDFprocessor/framework'.format(FWKBASE))
from RDFtree import RDFtree
sys.path.append('{}/Common/data'.format(FWKBASE))
from binning import ptBins, etaBins, mTBinsFull, mTBins,etaBins, isoBins, chargeBins, zmassBins, yBins, qtBins,metBins,pvBins, phiBins, cosThetaBins, qtBins_syst
from externals import fileSFul, fileSFPogTrk, weightFoldersrc

sys.path.append('{}/templateMaker/python'.format(FWKBASE))
from getReweightModules import *
ROOT.gSystem.Load('{}/templateMaker/bin/libAnalysisOnData.so'.format(FWKBASE))
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")


#Build the template building sequenc
def wSelectionSequence(p, systType, nodetoStart, era):
    p.EventFilter(nodeToStart=nodetoStart, nodeToEnd='defs', evfilter="1.", filtername="{:20s}".format("true"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(HLT_IsoMu24 ||  HLT_IsoTkMu24)", filtername="{:20s}".format("Pass HLT"))
    
    # note for customizeforUL(isMC=true, isWorZ=false)
    if systType == 0: #this is data
        p.branch(nodeToStart='defs', nodeToEnd='defs', modules=[ROOT.recoDefinitions(False, False), ROOT.mtDefinitions(False,ptprefix="MET_T1comp_pt", phiprefix="MET_T1comp_phi")])
        
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(vetoMuons)==1", filtername="{:20s}".format("vetomuon"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(goodMuons)==1", filtername="{:20s}".format("onemuon"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt > 25. && Mu1_pt < 55.", filtername="{:20s}".format("accept"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_mediumId[goodMuons][0] == 1 && Muon_isGlobal[goodMuons][0]", filtername="{:20s}".format("muonID"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_mediumId[goodMuons][0] == 1 && Muon_isGlobal[goodMuons][0]", filtername="{:20s}".format("muonID"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("mu1 trig matched"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_HBHENoiseIsoFilter && Flag_HBHENoiseFilter && Flag_BadPFMuonFilter", filtername="{:20s}".format("eventFilters"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(vetoElectrons) == 0;", filtername="{:20s}".format("vetoelectrons"))
        p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso"], types = ['float']*5,node='defs',histoname=ROOT.string('data_obs'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        
    elif systType < 2: #this is MC with no PDF variations
        #falling back to old lumi weight computation
        p.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.recoDefinitions(True, False), ROOT.mtDefinitions(False,ptprefix="MET_T1comp_pt", phiprefix="MET_T1comp_phi"), ROOT.SFprod(fileSFul,era=era)])
        
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(vetoMuons)==1", filtername="{:20s}".format("vetomuon"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(goodMuons)==1", filtername="{:20s}".format("onemuon"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt > 25. && Mu1_pt < 55.", filtername="{:20s}".format("accept"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_mediumId[goodMuons][0] == 1 && Muon_isGlobal[goodMuons][0]", filtername="{:20s}".format("muonID"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_mediumId[goodMuons][0] == 1 && Muon_isGlobal[goodMuons][0]", filtername="{:20s}".format("muonID"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("mu1 trig matched"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_HBHENoiseIsoFilter && Flag_HBHENoiseFilter && Flag_BadPFMuonFilter", filtername="{:20s}".format("eventFilters"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(vetoElectrons) == 0;", filtername="{:20s}".format("vetoelectrons"))

        p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])

        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SFSystvar"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_SFSystvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight"], types = ['float']*8,node='defs',histoname=ROOT.string('ewk_SFStatvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], sample=("SFStatvar",4))
        # #prefire variations
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","SF"], types = ['float']*8,node='defs',histoname=ROOT.string('ewk_prefireVar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], sample=("muprefireWeightVars",2))
        # #jec variations
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_jesTotalUp","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_jesTotalUp'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_jesTotalDown","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_jesTotalDown'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_unclustEnUp","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_unclustEnUp'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_unclustEnDown","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_unclustEnDown'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])

    else:
        p.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.defineHarmonics(),ROOT.recoDefinitions(True, False), ROOT.mtDefinitions(False,ptprefix="MET_T1comp_pt", phiprefix="MET_T1comp_phi"), ROOT.SFprod(fileSFul,era=era)])
        
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(vetoMuons)==1", filtername="{:20s}".format("vetomuon"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(goodMuons)==1", filtername="{:20s}".format("onemuon"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt > 25. && Mu1_pt < 55.", filtername="{:20s}".format("accept"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_mediumId[goodMuons][0] == 1 && Muon_isGlobal[goodMuons][0]", filtername="{:20s}".format("muonID"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_mediumId[goodMuons][0] == 1 && Muon_isGlobal[goodMuons][0]", filtername="{:20s}".format("muonID"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("mu1 trig matched"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_HBHENoiseIsoFilter && Flag_HBHENoiseFilter && Flag_BadPFMuonFilter", filtername="{:20s}".format("eventFilters"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(vetoElectrons) == 0;", filtername="{:20s}".format("vetoelectrons"))

        p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins])

        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SFSystvar"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_SFSystvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight"], types = ['float']*8,node='defs',histoname=ROOT.string('ewk_SFStatvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], sample=("SFStatvar",4))
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_LHEPdfWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], sample=("LHEPdfWeight",103))
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_LHEScaleWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], sample=("LHEScaleWeight",9))
        # #prefire variations
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","SF"], types = ['float']*8,node='defs',histoname=ROOT.string('ewk_prefireVar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], sample=("muprefireWeightVars",2))
        # # #jec variations
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_jesTotalUp","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_jesTotalUp'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_jesTotalDown","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_jesTotalDown'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_unclustEnUp","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_unclustEnUp'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_unclustEnDown","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_unclustEnDown'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])

    return p

def wSelectionHelWeightsSequence(p, nodetoStart,era):
    # here get angular coefficients
    p.branch(nodeToStart=nodetoStart, nodeToEnd='defs', modules=[ROOT.defineHarmonics()])
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight"], types = ['float']*3,node='defs',histoname=ROOT.string("xsecs"),bins=[yBins,qtBins], sample=('harmonicsVec',9))
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight"], types = ['float']*3,node='defs',histoname=ROOT.string("totxsecs"),bins=[yBins,qtBins])
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","CStheta_preFSR","lumiweight"], types = ['float']*4,node='defs',histoname=ROOT.string("qtycostheta"),bins=[yBins,qtBins,cosThetaBins])

    # p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight","yqtweight","coeffsweight"], types = ['float']*5,node='defs',histoname=ROOT.string("xsecs"),bins=[yBins,qtBins], sample=('harmonicsVec',9))
    # p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight","yqtweight","coeffsweight"], types = ['float']*5,node='defs',histoname=ROOT.string("totxsecs"),bins=[yBins,qtBins])
    
    # pdf variations
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight"], types = ['float']*3,node='defs',histoname=ROOT.string("xsecs_LHEPdfWeight"),bins=[yBins,qtBins], sample=('harmonicsVec_LHEPdfWeight',9*103))
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight"], types = ['float']*3,node='defs',histoname=ROOT.string("totxsecs_LHEPdfWeight"),bins=[yBins,qtBins], sample=("LHEPdfWeight",103))
    # scale variations
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight"], types = ['float']*3,node='defs',histoname=ROOT.string("xsecs_LHEScaleWeight"),bins=[yBins,qtBins], sample=('harmonicsVec_LHEScaleWeight',9*9))
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight"], types = ['float']*3,node='defs',histoname=ROOT.string("totxsecs_LHEScaleWeight"),bins=[yBins,qtBins], sample=("LHEScaleWeight",9))
    return p

def wSelectionDifferentialSequence(p,era,sample):
    genInfo='{}/genInfo_syst.root'.format(weightFoldersrc)
    helWeightsrc='{}/Common/data/reweight/'.format(FWKBASE)
    # helWeightsrc='{}/config/'.format(FWKBASE)
    chargeStr=''
    charge=1
    genCoeff=''
    helWeightFile=''
    mods=[]
    # here load angular coefficients and reweight
    #separate blocks are needed for Numba defines
    if 'WPlus' in sample:
        chargeStr= 'WPlus' 
        charge=1
        genCoeff='{}/genInput_v7_syst_Wplus.root'.format(weightFoldersrc)
        helWeightFile='{}/outputW_sroychow_{}/{}JetsToMuNu_helweights.hdf5'.format(weightFoldersrc, era, chargeStr)
        mods=[getHelWeightsWplus(era=era, helwtFile=helWeightFile, syst=""), \
              reweightcoeffsWplus(era=era,helWtsrcdir=weightFoldersrc, geninputF=genCoeff), \
              reweightyqtWplus(era=era, inFilehelwt=helWeightFile, genInfoFile=genInfo), \
              getMassWeightsWplus(era=era)]
    else:
        chargeStr= 'WMinus'
        charge=-1
        genCoeff='{}/genInput_v7_syst_Wminus.root'.format(weightFoldersrc)
        helWeightFile='{}/outputW_sroychow_{}/{}JetsToMuNu_helweights.hdf5'.format(weightFoldersrc, era, chargeStr)
        mods=[getHelWeightsWminus(era=era, helwtFile=helWeightFile, syst=""), \
              reweightyqtWminus(era=era, inFilehelwt=helWeightFile, genInfoFile=genInfo),
              reweightcoeffsWminus(era=era,helWtsrcdir=weightFoldersrc, geninputF=genCoeff),
              getMassWeightsWminus(era=era)]

    p.branch(nodeToStart='defs', nodeToEnd='templates', modules=mods)
    p.EventFilter(nodeToStart='templates', nodeToEnd='nominal', evfilter="Vrap_preFSR_abs<2.4 && Vpt_preFSR<60.", filtername="{:20s}".format("signal templ"))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights",9))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates_mass'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helmassweights",9*2))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SFSystvar"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates_SFSystvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights",9))
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight"], types = ['float']*10,node='nominal',histoname=ROOT.string('signalTemplates_SFStatvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("SFStatvar_helweights",9*4))
    #prefire variations
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","SF"], types = ['float']*10,node='nominal',histoname=ROOT.string('signalTemplates_prefireVars'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("muprefireWeightVars_helweights",9*2))
    # #jes
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_jesTotalUp","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates_jesTotalUp'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights",9))
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_jesTotalDown","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates_jesTotalDown'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights",9))
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_unclustEnUp","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates_unclustEnUp'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights",9))
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_unclustEnDown","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates_unclustEnDown'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights",9))
    
    # low acceptance
    p.EventFilter(nodeToStart='templates', nodeToEnd='lowacc', evfilter="(Vrap_preFSR_abs>2.4 || Vpt_preFSR>60.) && Vpt_preFSR<200.", filtername="{:20s}".format("low acc"))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*10,node='lowacc',histoname=ROOT.string('lowacc'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,qtBins_syst])
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='lowacc',histoname=ROOT.string('lowacc_mass'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins],sample=('massWeights',2))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='lowacc',histoname=ROOT.string('lowacc_LHEPdfWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins],sample=('LHEPdfWeight',103))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*10,node='lowacc',histoname=ROOT.string('lowacc_LHEScaleWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,qtBins_syst],sample=('LHEScaleWeight',9))
    
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF","yqtweight","coeffsweight"], types = ['float']*11,node='lowacc',histoname=ROOT.string('lowacc_rew'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins])
    
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SFSystvar"], types = ['float']*9,node='lowacc',histoname=ROOT.string('lowacc_SFSystvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins])
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","muprefireWeight"], types = ['float']*8,node='lowacc',histoname=ROOT.string('lowacc_SFStatvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins],sample=('SFStatvar',4))

    # #prefire
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","SF"], types = ['float']*8,node='lowacc',histoname=ROOT.string('lowacc_prefireVars'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], sample=("muprefireWeightVars",2))
    
    # #jec
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_jesTotalUp","Mu1_relIso","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*10,node='lowacc',histoname=ROOT.string('lowacc_jesTotalUp'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,qtBins_syst])
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_jesTotalDown","Mu1_relIso","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*10,node='lowacc',histoname=ROOT.string('lowacc_jesTotalDown'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,qtBins_syst])
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_unclustEnUp","Mu1_relIso","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*10,node='lowacc',histoname=ROOT.string('lowacc_unclustEnUp'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,qtBins_syst])
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT_unclustEnDown","Mu1_relIso","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*10,node='lowacc',histoname=ROOT.string('lowacc_unclustEnDown'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,qtBins_syst])

    # pdf-scale uncertainties for full templates
    if 'WPlus' in sample:
        helWeightFile='{}/powheg_acc_{}/{}JetsToMuNu_helweights.hdf5'.format(weightFoldersrc, era, chargeStr)
        mods=[getHelWeightsWplus(era=era, helwtFile=helWeightFile, syst="LHEPdfWeight")]
    else:
        helWeightFile='{}/powheg_acc_{}/{}JetsToMuNu_helweights.hdf5'.format(weightFoldersrc, era, chargeStr)
        mods=[getHelWeightsWminus(era=era, helwtFile=helWeightFile, syst="LHEPdfWeight")]
    #p.branch(nodeToStart='defs', nodeToEnd='LHEPdfWeight', modules=mods)
    #p.EventFilter(nodeToStart='LHEPdfWeight', nodeToEnd='LHEPdfWeight', evfilter="Vrap_preFSR_abs<2.4 && Vpt_preFSR<60.", filtername="{:20s}".format("signal templ"))
    #p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='LHEPdfWeight',histoname=ROOT.string('signalTemplates_LHEPdfWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights_LHEPdfWeight",9*103))

    if 'WPlus' in sample:
        helWeightFile='{}/powheg_acc_{}/{}JetsToMuNu_helweights.hdf5'.format(helWeightsrc, era, chargeStr)
        mods=[getHelWeightsWplus(era=era, helwtFile=helWeightFile, syst="LHEScaleWeight")]
    else:
        helWeightFile='{}/powheg_acc_{}/{}JetsToMuNu_helweights.hdf5'.format(helWeightsrc, era, chargeStr)
        mods=[getHelWeightsWminus(era=era, helwtFile=helWeightFile, syst="LHEScaleWeight")]
    
    #p.branch(nodeToStart='defs', nodeToEnd='LHEScaleWeight', modules=mods)
    #p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='LHEScaleWeight',histoname=ROOT.string('signalTemplates_LHEScaleWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights_LHEScaleWeight",9*9))
    
    return p
