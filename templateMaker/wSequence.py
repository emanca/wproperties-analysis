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
from binning import ptBins, etaBins, mTBinsFull, mTBins,etaBins, isoBins, chargeBins, zmassBins, yBins, qtBins,metBins,pvBins, phiBins, cosThetaBins
from externals import fileSFul

sys.path.append('{}/templateMaker/python'.format(FWKBASE))
ROOT.gSystem.Load('{}/templateMaker/bin/libAnalysisOnData.so'.format(FWKBASE))
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")


#Build the template building sequenc
def wSelectionSequence(p, systType, nodetoStart, era):

    p.EventFilter(nodeToStart=nodetoStart, nodeToEnd='defs', evfilter="(HLT_IsoMu24 ||  HLT_IsoTkMu24)", filtername="{:20s}".format("Pass HLT"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(Muon_mediumId)", filtername="{:20s}".format("MuonID"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_dxy) < 0.05)", filtername="{:20s}".format("dxy"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_dz) < 0.2)", filtername="{:20s}".format("dz"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(Muon_pt > 25.)", filtername="{:20s}".format("Muon pt cut"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_eta) < 2.4)", filtername="{:20s}".format("Muon eta cut"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="nMuon== 1", filtername="{:20s}".format("one muon"))
    
    # note for customizeforUL(isMC=true, isWorZ=false)
    if systType == 0: #this is data
        p.branch(nodeToStart='defs', nodeToEnd='defs', modules=[ROOT.recoDefinitions(False, False), ROOT.mtDefinitions(False)])
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("mu1 trig matched"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt < 65.", filtername="{:20s}".format("mu1 pt-eta acceptance"))
        p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso"], types = ['float']*5,node='defs',histoname=ROOT.string('data_obs'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
    elif systType < 2: #this is MC with no PDF variations
        #falling back to old lumi weight computation
        p.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.recoDefinitions(True, False),ROOT.mtDefinitions(True,"MET_pt", "MET_phi"),ROOT.SF_ul(fileSFul,isZ=False,era=era)])
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("mu1 trig matched"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt < 65.", filtername="{:20s}".format("mu1 pt-eta acceptance"))
        p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])
    else:
        p.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.recoDefinitions(True, False), ROOT.mtDefinitions(), ROOT.SF_ul(fileSFul)])
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("mu1 trig matched"))
        p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt < 65.", filtername="{:20s}".format("mu1 pt-eta acceptance"))
        p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*9,node='defs',histoname=ROOT.string('ewk'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins], variations = [])

    return p

def wSelectionHelWeightsSequence(p, nodetoStart):
    # here get angular coefficients
    p.branch(nodeToStart=nodetoStart, nodeToEnd='defs', modules=[ROOT.defineHarmonics()])
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight"], types = ['float']*3,node='defs',histoname=ROOT.string("xsecs"),bins=[yBins,qtBins], sample='harmonicsVec')
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight"], types = ['float']*3,node='defs',histoname=ROOT.string("totxsecs"),bins=[yBins,qtBins])
    # pdf variations
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight"], types = ['float']*3,node='defs',histoname=ROOT.string("xsecs_LHEPdfWeight"),bins=[yBins,qtBins], sample='harmonicsVec_LHEPdfWeight')
    p.Histogram(columns = ["Vrap_preFSR_abs","Vpt_preFSR","lumiweight"], types = ['float']*3,node='defs',histoname=ROOT.string("totxsecs_LHEPdfWeight"),bins=[yBins,qtBins], sample='LHEPdfWeight')
    
    return p
def wSelectionDifferentialSequence(p):
    from getHelWeights import getHelWeights
    from getMassWeights import getMassWeights
    # here load angular coefficients and reweight
    p.branch(nodeToStart='defs', nodeToEnd='nominal', modules=[ROOT.defineHarmonics(),getHelWeights(),getMassWeights()])
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='nominal',histoname=ROOT.string('signalTemplates'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample="helWeights")
    # p.branch(nodeToStart='defs', nodeToEnd='LHEPdfWeight', modules=[ROOT.defineHarmonics(),getHelWeights(syst="LHEPdfWeight")])
    # p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='LHEPdfWeight',histoname=ROOT.string('signalTemplates_LHEPdfWeight'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample="helWeights_LHEPdfWeight")
    p.Histogram(columns = ["Mu1_eta","Mu1_pt","Mu1_charge","MT","Mu1_relIso", "Vrap_preFSR_abs","Vpt_preFSR","lumiweight","puWeight","muprefireWeight","SF"], types = ['float']*11,node='mass',histoname=ROOT.string('signalTemplates_mass'),bins = [etaBins,ptBins,chargeBins,mTBins,isoBins,yBins,qtBins], sample="helmassWeights")

    return p
