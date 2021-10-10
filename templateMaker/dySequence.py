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
from binning import ptBins, etaBins, mTBins, etaBins, isoBins, chargeBins, zmassBins, qtBins,metBins,pvBins,phiBins,cosThetaBins
from externals import fileSFul,filePt, fileY,fileSFPogTrk
from dataluminosity import lumi_preVFP,lumi_postVFP,lumi_total2016

sys.path.append('{}/templateMaker/python'.format(FWKBASE))
from getLumiWeight import getLumiWeight
ROOT.gSystem.Load('{}/templateMaker/bin/libAnalysisOnData.so'.format(FWKBASE))
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")


#Build the template building sequenc
def dySelectionSequence(p, xsec, systType, sumwClipped, nodetoStart, era):
    print(ptBins)
    print(zmassBins)
    luminosityN = lumi_total2016
    if era == 'preVFP' :     luminosityN = lumi_preVFP
    else: luminosityN = lumi_postVFP

    if systType not in [0,1]:#for signal MC
        p.Histogram(columns = ["CStheta_preFSR","CSphi_preFSR","Vpt_preFSR","Vrap_preFSR","Vmass_preFSR"], types = ['float']*5,node=nodetoStart,histoname=ROOT.string('genhistos'),bins = [cosThetaBins,phiBins,qtBins,etaBins,zmassBins], variations = [])
    
    p.branch(nodeToStart=nodetoStart, nodeToEnd='defs', modules=[ROOT.zVetoMuons()])
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(vetoMuons)==2 && Sum(goodMuons)==2", filtername="{:20s}".format("two muons"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(HLT_IsoMu24 ||  HLT_IsoTkMu24)", filtername="{:20s}".format("Pass HLT"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(Muon_charge[goodMuons][0] + Muon_charge[goodMuons][1]) == 0", filtername="{:20s}".format("Opposite charge"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="std::abs(Muon_eta[goodMuons][0]) < 2.4 && std::abs(Muon_eta[goodMuons][1]) < 2.4", filtername="{:20s}".format("Accept"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_mediumId[goodMuons][0] == 1 && Muon_mediumId[goodMuons][1] == 1", filtername="{:20s}".format("MuonId"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_pfRelIso04_all[goodMuons][0] < 0.15 && Muon_pfRelIso04_all[goodMuons][1] < 0.15", filtername="{:20s}".format("Isolation"))
    p.branch(nodeToStart='defs', nodeToEnd='defs', modules=[ROOT.zSelection()])
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="60. < dimuonMass && dimuonMass < 120.", filtername="{:20s}".format("mZ range"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("+ve mu trig matched"))

    if systType == 0: #this is data
        p.Histogram(columns = ["dimuonMass", "Mu1_eta", "Mu1_pt", "Mu2_eta", "Mu2_pt"], types = ['float']*5,node='defs',histoname=ROOT.string('data_muons'),bins = [zmassBins, etaBins, ptBins, etaBins, ptBins], variations = [])
        #p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu2_pt < 200. ", filtername="{:20s}".format("mu2 pt upper acceptance"))
        p.Histogram(columns = ["dimuonMass", "dimuonPt", "dimuonY", "MET_pt"], types = ['float']*4,node='defs',histoname=ROOT.string('data_dimuon'),bins = [zmassBins,qtBins, etaBins, metBins], variations = [])
        return p
    elif systType == 1:
        print("Sample will be normalized to {}/fb".format(luminosityN))
        p.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.SF_ul(fileSFul, fileSFPogTrk, isZ=True,era=era)])
        #p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu2_pt < 200. ", filtername="{:20s}".format("mu2 pt upper acceptance"))
        p.Histogram(columns = ["dimuonMass", "dimuonPt", "dimuonY", "MET_pt", "lumiweight", "puWeight", "SF"], types = ['float']*7,node='defs',histoname=ROOT.string('DY_dimuon'),bins = [zmassBins,qtBins, etaBins,metBins], variations = [])
        return p
    else:
        print("Sample will be normalized to {}/fb".format(luminosityN))
        p.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.SF_ul(fileSFul, fileSFPogTrk, isZ=True, era=era)])
        p.Histogram(columns = ["dimuonMass", "Mu1_eta", "Mu1_pt", "Mu2_eta", "Mu2_pt", "lumiweight", "puWeight", "SF"], types = ['float']*8,node='defs',histoname=ROOT.string('DY_muons'),bins = [zmassBins, etaBins, ptBins, etaBins, ptBins], variations = [])
        #p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu2_pt < 200. ", filtername="{:20s}".format("mu2 pt upper acceptance"))
        p.Histogram(columns = ["dimuonMass", "dimuonPt", "dimuonY", "MET_pt", "lumiweight", "puWeight", "SF"], types = ['float']*7,node='defs',histoname=ROOT.string('DY_dimuon'),bins = [zmassBins,qtBins, etaBins,metBins], variations = [])
        return p


    '''
    #Selection for PC 25/02/2021
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(Muon_mediumId)", filtername="{:20s}".format("MuonID"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(Muon_pfRelIso04_all < 0.15)", filtername="{:20s}".format("both mu isolated"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_dxy) < 0.05)", filtername="{:20s}".format("dxy"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_dz) < 0.2)", filtername="{:20s}".format("dz"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(Muon_pt > 25.)", filtername="{:20s}".format("Muon pt cut"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_eta) < 2.4)", filtername="{:20s}".format("Muon eta cut"))

    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="nMuon==2", filtername="{:20s}".format("twomuon"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(Muon_charge[0] + Muon_charge[1] )== 0", filtername="{:20s}".format("Opposite charge"))
    p.branch(nodeToStart='defs', nodeToEnd='defs', modules=[ROOT.zSelection()])
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="60. < dimuonMass && dimuonMass < 120.", filtername="{:20s}".format("mZ range"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("+ve mu trig matched"))
    p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt < 65.", filtername="{:20s}".format("mu1 pt-eta acceptance"))
    '''
