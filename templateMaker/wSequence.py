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
from externals import fileSFul, fileSFPogTrk

sys.path.append('{}/templateMaker/python'.format(FWKBASE))
ROOT.gSystem.Load('{}/templateMaker/bin/libAnalysisOnData.so'.format(FWKBASE))
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")


#Build the template building sequenc
def wSelectionSequence(p, systType, nodetoStart, era):
    p.EventFilter(nodeToStart=nodetoStart, nodeToEnd='defs', evfilter="1.", filtername="{:20s}".format("true"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(HLT_IsoMu24 ||  HLT_IsoTkMu24)", filtername="{:20s}".format("Pass HLT"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(Muon_mediumId)", filtername="{:20s}".format("MuonID"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_dxy) < 0.05)", filtername="{:20s}".format("dxy"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_dz) < 0.2)", filtername="{:20s}".format("dz"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(Muon_pt > 25.)", filtername="{:20s}".format("Muon pt cut"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_eta) < 2.4)", filtername="{:20s}".format("Muon eta cut"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="nMuon== 1", filtername="{:20s}".format("one muon"))
    
    # note for customizeforUL(isMC=true, isWorZ=false)
    if systType == 0: #this is data
        p.branch(nodeToStart='defs', nodeToEnd='defs', modules=[ROOT.recoDefinitions(False, False), ROOT.mtDefinitions(False,ptprefix="MET_pt", phiprefix="MET_phi")])
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("mu1 trig matched"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt < 65.", filtername="{:20s}".format("mu1 pt-eta acceptance"))
        p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso"], types = ['float']*5,node='defs',histoname=ROOT.string('data_obs'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
    elif systType < 2: #this is MC with no PDF variations
        #falling back to old lumi weight computation
        p.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.recoDefinitions(True, False),ROOT.mtDefinitions(False,ptprefix="MET_pt", phiprefix="MET_phi"),ROOT.SF_ul(fileSFul, fileSFPogTrk,era=era)])
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("mu1 trig matched"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt < 65.", filtername="{:20s}".format("mu1 pt-eta acceptance"))
        p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        #p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SFvar"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_SFvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])

    else:
        p.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.defineHarmonics(),ROOT.recoDefinitions(True, False), ROOT.mtDefinitions(False,ptprefix="MET_pt", phiprefix="MET_phi"), ROOT.SF_ul(fileSFul, fileSFPogTrk,isZ=False,era=era)])
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("mu1 trig matched"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt < 65.", filtername="{:20s}".format("mu1 pt-eta acceptance"))
        p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins])
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SFvar"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_SFvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_LHEPdfWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], sample=("LHEPdfWeight",103))
        # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk_LHEScaleWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], sample=("LHEScaleWeight",9))

    return p

def wSelectionHelWeightsSequence(p, nodetoStart,era):
    from reweightyqt import reweightyqt
    from getHelWeights import getHelWeights
    from reweightyqt import reweightyqt
    from reweightcoeffs import reweightcoeffs
    # here get angular coefficients
    p.branch(nodeToStart=nodetoStart, nodeToEnd='defs', modules=[ROOT.defineHarmonics()]) #,reweightyqt(era=era),reweightycostheta(era=era),getHelWeights(era=era,syst="",pseudodata=True)
    
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
    from getHelWeights import getHelWeights
    from reweightyqt import reweightyqt
    from reweightycostheta import reweightycostheta
    from getMassWeights import getMassWeights
    # here load angular coefficients and reweight
    helWeightsrc='{}/Common/data/reweight/'.format(FWKBASE)
    helWeightFile='{}/Common/data/reweight/powheg_acc_{}/WPlusJetsToMuNu_helweights.hdf5'.format(FWKBASE, era)
    genInfo='{}/Common/data/reweight/genInfo_syst.root'.format(FWKBASE)
    genCoeff='{}/Common/data/reweight/genInput_v7_syst_Wplus.root'.format(FWKBASE)

    p.branch(nodeToStart='defs', nodeToEnd='templates', modules=[getHelWeights(era=era,helwtFile=helWeightFile,syst=""),getMassWeights(era=era)])

    p.EventFilter(nodeToStart='templates', nodeToEnd='nominal', evfilter="Vrap_preFSR_abs<2.4 && Vpt_preFSR<60.", filtername="{:20s}".format("signal templ"))

    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights",9))
    # p.branch(nodeToStart='defs', nodeToEnd='mass', modules=[ROOT.getMassWeights(),ROOT.defineHarmonics(),getHelWeights(era=era,syst="mass")])
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates_mass'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helmassweights",9*2))
    
    # low acceptance
    p.EventFilter(nodeToStart='templates', nodeToEnd='lowacc', evfilter="(Vrap_preFSR_abs>2.4 || Vpt_preFSR>60.) && Vpt_preFSR<200.", filtername="{:20s}".format("low acc"))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*10,node='lowacc',histoname=ROOT.string('lowacc'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,qtBins_syst])
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='lowacc',histoname=ROOT.string('lowacc_mass'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins],sample=('massWeights',2))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='lowacc',histoname=ROOT.string('lowacc_LHEPdfWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins],sample=('LHEPdfWeight',103))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*10,node='lowacc',histoname=ROOT.string('lowacc_LHEScaleWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,qtBins_syst],sample=('LHEScaleWeight',9))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF","yqtweight","coeffsweight"], types = ['float']*12,node='lowacc',histoname=ROOT.string('lowacc_rew'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,qtBins_syst])

    # pdf-scale uncertainties for full templates
    p.branch(nodeToStart='defs', nodeToEnd='LHEPdfWeight', modules=[getHelWeights(era=era,syst="LHEPdfWeight")])
    p.EventFilter(nodeToStart='LHEPdfWeight', nodeToEnd='LHEPdfWeight', evfilter="Vrap_preFSR_abs<2.4 && Vpt_preFSR<60.", filtername="{:20s}".format("signal templ"))
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='LHEPdfWeight',histoname=ROOT.string('signalTemplates_LHEPdfWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights_LHEPdfWeight",9*103))
    # p.branch(nodeToStart='defs', nodeToEnd='LHEScaleWeight', modules=[getHelWeights(era=era,syst="LHEScaleWeight")])
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='LHEScaleWeight',histoname=ROOT.string('signalTemplates_LHEScaleWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights_LHEScale",9*9))
    
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SFvar"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates_SFvar'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample=("helWeights",9))

    return p
