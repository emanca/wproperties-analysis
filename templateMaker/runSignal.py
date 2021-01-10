import os
import sys
import ROOT
import copy
import time
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from root_numpy import array2hist
from array import array
sys.path.append('../RDFprocessor/framework')
from RDFtree import RDFtree
sys.path.append('python/')
from getLumiWeight import getLumiWeight
sys.path.append('data/')
from binning import qtBins, yBins, ptBins, etaBins
from externals import filePt, fileY, fileSF

# matplotlib stuff
plt.style.use([hep.style.ROOT, hep.style.firamath])
#hep.cms.label(loc=0)
hep.cms.text('Simulation')

ROOT.gSystem.Load('bin/libAnalysisOnData.so')
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")

ROOT.ROOT.EnableImplicitMT(64)

outputDir = 'PLOTS'
inputFile = '/scratchnvme/wmass/WJetsNoCUT_v2/tree_*_*.root'

p = RDFtree(outputDir = outputDir, inputFile = inputFile, outputFile="test.root", pretend=False)
p.branch(nodeToStart='input', nodeToEnd='defs', modules=[getLumiWeight(xsec=61526.7, inputFile=inputFile),ROOT.reweightFromZ(filePt,fileY), ROOT.genDefinitions(),ROOT.defineHarmonics(), ROOT.Replica2Hessian()])
harmonics = {"P0" : (20./3., 1./10),"P1": (5.,0.), "P2": (20.,0.), "P3": (4.,0.),"P4":(4.,0.),"P5":(5.,0.),"P6":(5.,0.),"P7":(4.,0.)}

for harm in harmonics:
    if "P5" in harm or "P6" in harm or "P7" in harm: # don't make variations
        p.Histogram(columns = ["Wrap_preFSR_abs","Wpt_preFSR","{}".format(harm),"lumiweight", "weightPt"], types = ['float']*5,node='defs',histoname=ROOT.string("xsecs_{}".format(harm)),bins = [yBins,qtBins])
    else:
        p.Histogram(columns = ["Wrap_preFSR_abs","Wpt_preFSR","{}".format(harm),"lumiweight", "weightPt"], types = ['float']*5,node='defs',histoname=ROOT.string("xsecs_{}".format(harm)),bins = [yBins,qtBins])
p.Histogram(columns = ["Wrap_preFSR_abs","Wpt_preFSR","lumiweight", "weightPt"], types = ['float']*4,node='defs',histoname=ROOT.string("xsecs"),bins = [yBins,qtBins])
p.gethdf5Output()

shape = (len(yBins)-1) * (len(qtBins)-1) 

fewk = h5py.File('PLOTS/test.hdf5', mode='r+')
htot = np.array(fewk['xsecs'][:].reshape((len(yBins)-1,len(qtBins)-1),order='F'),order='C')

f = ROOT.TFile('PLOTS/angcoeff.root','recreate')
for harmonic, factors in harmonics.items():
    hewk = np.array(fewk['xsecs_{}'.format(harmonic)][:].reshape((len(yBins)-1,len(qtBins)-1),order='F'),order='C')
    hep.hist2dplot(factors[0]*(hewk/htot+factors[1]), yBins, qtBins)
    plt.savefig('PLOTS/xsecs_{}.png'.format(harmonic))
    plt.clf()
    #ROOT format to feed RDF
    th2 = ROOT.TH2D('xsecs_{}'.format(harmonic),'xsecs_{}'.format(harmonic), len(yBins)-1, array('d',yBins), len(qtBins)-1, array('d',qtBins))
    hist = array2hist(factors[0]*(hewk/htot+factors[1]), th2, errors=None)
    f.cd()
    hist.Write()

hep.hist2dplot(htot, yBins, qtBins)
plt.savefig('PLOTS/xsecs_{}.png'.format(harmonic))
plt.clf()
th2 = ROOT.TH2D('xsecs','xsecs', len(yBins)-1, array('d',yBins), len(qtBins)-1, array('d',qtBins))
hist = array2hist(htot, th2, errors=None)
f.cd()
hist.Write()
f.Close()

f = ROOT.TFile.Open('PLOTS/angcoeff.root')

helicities = ["L", "I", "T", "A", "P", "UL"]

p = RDFtree(outputDir = outputDir, inputFile = inputFile, outputFile="SignalTemplates.root", pretend=False)
p.branch(nodeToStart='input', nodeToEnd='defs', modules=[getLumiWeight(xsec=61526.7, inputFile=inputFile),ROOT.reweightFromZ(filePt,fileY), ROOT.genDefinitions(),ROOT.defineHarmonics(), ROOT.Replica2Hessian()])
p.branch(nodeToStart='defs', nodeToEnd='signalTemplates', modules=[ROOT.recoDefinitions(True,True),ROOT.recoWeightDefinitions(fileSF), ROOT.getACValues(f,f), ROOT.getMassWeights(), ROOT.getWeights()])

for helicity in helicities:
    p.Histogram(columns = ["Mu1_eta", "Mu1_pt", "Wrap_preFSR_abs","Wpt_preFSR","weight{}".format(helicity),"lumiweight", "weightPt", "massWeights", "alphaSVars", "LHEPdfWeightHess"], types = ['float']*10,node='signalTemplates',histoname=ROOT.string("xsecs_{}".format(helicity)),bins = [etaBins, ptBins, yBins,qtBins])
p.gethdf5Output()