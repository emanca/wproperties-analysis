import os
import sys
import ROOT
import copy
import time
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from array import array
sys.path.append('../RDFprocessor/framework')
sys.path.append('../Common/data')
from RDFtree import RDFtree
sys.path.append('python/')
from getLumiWeight import getLumiWeight
from binning import qtBins, yBins, ptBins, etaBins, chargeBins
from externals import filePt, fileY, fileSF

# matplotlib stuff
plt.style.use([hep.style.ROOT, hep.style.firamath])
#hep.cms.label(loc=0)
hep.cms.text('Simulation')

ROOT.gSystem.Load('bin/libAnalysisOnData.so')
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")

ROOT.ROOT.EnableImplicitMT(128)

outputDir = 'PLOTS'
inputFile = '/scratchnvme/wmass/NanoAOD2016-UL/postNanoDec2020/WplusJetsToMuNu_preVFP_addVars/merged/*.root'

p = RDFtree(outputDir = outputDir, inputFile = inputFile, outputFile="test.root", pretend=False)
p.branch(nodeToStart='input', nodeToEnd='defs', modules=[ROOT.lumiWeight(xsec=11572.19), ROOT.customizeforUL(), ROOT.genDefinitions(), ROOT.defineHarmonics()])
harmonics = {"P0" : (20./3., 1./10),"P1": (5.,0.), "P2": (20.,0.), "P3": (4.,0.),"P4":(4.,0.),"P5":(5.,0.),"P6":(5.,0.),"P7":(4.,0.)}

for harm in harmonics:
    p.Histogram(columns = ["Wrap_preFSR_abs","Vpt_preFSR","{}".format(harm),"lumiweight","pdfWeightNNPDF0"], types = ['float']*5,node='defs',histoname=ROOT.string("xsecs_{}".format(harm)),bins = [yBins,qtBins], variations = {"pdfWeightNNPDF0":"LHEPdfWeight"})
p.Histogram(columns = ["Wrap_preFSR_abs","Vpt_preFSR","lumiweight","pdfWeightNNPDF0"], types = ['float']*4,node='defs',histoname=ROOT.string("xsecs"),bins = [yBins,qtBins], variations = {"pdfWeightNNPDF0":"LHEPdfWeight"})
p.gethdf5Output()
assert(0)
fewk = h5py.File('PLOTS/test.hdf5', mode='r+')
htot = np.array(fewk['xsecs'][:].reshape((len(yBins)-1,len(qtBins)-1),order='F'),order='C')

f = ROOT.TFile('PLOTS/angcoeff.root','recreate')
for harmonic, factors in harmonics.items():
    hewk = np.array(fewk['xsecs_{}'.format(harmonic)][:].reshape((len(yBins)-1,len(qtBins)-1),order='F'),order='C')
    hep.hist2dplot(factors[0]*(hewk[:,:]/htot[:,:]+factors[1]), yBins, qtBins)
    plt.savefig('PLOTS/xsecs_{}.png'.format(harmonic))
    plt.clf()
    #ROOT format to feed RDF
    th2 = ROOT.TH2D('xsecs_{}'.format(harmonic),'xsecs_{}'.format(harmonic), len(yBins)-1, array('d',yBins), len(qtBins)-1, array('d',qtBins))
    for i in range(1, th2.GetNbinsX()+1):
        for j in range(1, th2.GetNbinsY()+1):
            th2.SetBinContent(i,j,hewk[i-1,j-1]/htot[i-1,j-1]+factors[1])
    f.cd()
    th2.Write()

hep.hist2dplot(htot[:,:], yBins, qtBins)
plt.savefig('PLOTS/xsecs_{}.png'.format(harmonic))
plt.clf()
th2 = ROOT.TH2D('xsecs','xsecs', len(yBins)-1, array('d',yBins), len(qtBins)-1, array('d',qtBins))
for i in range(1, th2.GetNbinsX()+1):
        for j in range(1, th2.GetNbinsY()+1):
            th2.SetBinContent(i,j,htot[i-1,j-1])
f.cd()
th2.Write()
f.Close()

assert(0)
f = ROOT.TFile.Open('PLOTS/angcoeff.root')

helicities = ["L", "I", "T", "A", "P", "UL"]

p = RDFtree(outputDir = outputDir, inputFile = inputFile, outputFile="SignalTemplates.root", pretend=True)
p.branch(nodeToStart='input', nodeToEnd='defs', modules=[getLumiWeight(xsec=61526.7, inputFile=inputFile),ROOT.customizeforUL(True, True), ROOT.genDefinitions(),ROOT.defineHarmonics()])
p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="HLT_SingleMu24", filtername="{:20s}".format("Pass HLT"))
p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(Vtype==0 || Vtype==1)", filtername="{:20s}".format("Vtype selection"))
p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="HLT_SingleMu24", filtername="{:20s}".format("Pass HLT"))
p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="MET_filters==1", filtername="{:20s}".format("Pass MET filter"))
p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="nVetoElectrons==0", filtername="{:20s}".format("Electron veto"))
p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Idx_mu1>-1", filtername="{:20s}".format("Atleast 1 mu"))
p.branch(nodeToStart='defs', nodeToEnd='signalTemplates', modules=[ROOT.recoDefinitions(True,True),ROOT.recoWeightDefinitions(fileSF), ROOT.getACValues(f,f), ROOT.getMassWeights(), ROOT.getWeights()])

for helicity in helicities:
    p.Histogram(columns = ["Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs","Vpt_preFSR","weight{}".format(helicity),"lumiweight", "weightPt", "massWeights"], types = ['float']*8,node='signalTemplates',histoname=ROOT.string("xsecs_{}".format(helicity)),bins = [etaBins, ptBins, yBins,qtBins])
p.gethdf5Output()
