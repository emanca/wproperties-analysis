import os
import sys
import ROOT

from RDFtree import RDFtree
sys.path.append('python/')


ROOT.gSystem.Load('bin/libSignalAnalysis.so')

c=64

ROOT.ROOT.EnableImplicitMT(c)

print("running with {} cores".format(c))

inputFile = '/scratchssd/emanca/wproperties-analysis/signalAnalysis/nanowmass_1.root'

p = RDFtree(outputDir = 'TESTPDF', inputFile = inputFile, outputFile="testpdf.root")
p.branch(nodeToStart = 'input', nodeToEnd = 'PDFreweight', modules = [ROOT.baseDefinitions(),ROOT.Replica2Hessian()])
p.getOutput()


