# bkg_systematics = {
#      #"puWeightVars"  : ["puWeightUp", "puWeightDown"],
#     "PrefireWeight"  : ["PrefireWeightUp", "PrefireWeightDown"],
#     "WHSF"  : ["WHSFSyst0Up", "WHSFSyst1Up","WHSFSyst2Up","WHSFSystFlatUp","WHSFSyst0Down", "WHSFSyst1Down","WHSFSyst2Down","WHSFSystFlatDown"],
#     "LHEScaleWeight" : ["LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0","LHEScaleWeight_muR1p0_muF0p5","LHEScaleWeight_muR1p0_muF2p0","LHEScaleWeight_muR2p0_muF1p0", "LHEScaleWeight_muR2p0_muF2p0"],
#     "ptScale" : [ "correctedUp", "correctedDown"], 
#     "jme" : ["jesTotalUp", "jesTotalDown", "unclustEnUp","unclustEnDown"],
#     "LHEPdfWeight" : ["LHEPdfWeightHess{}".format(i+1) for i in range(60)],
#     "alphaS" : ["alphaSUp", "alphaSDown"]
# }
#bkg_systematics["LHEPdfWeight"].append("alphaSUp")
#bkg_systematics["LHEPdfWeight"].append("alphaSDown")

import copy
import sys
sys.path.append('../analysisOnData/data/')
sys.path.append('../data/') #for the prefit plotter only
from systematics import systematics, systFlat

def buildBkgSyst() :
    out = {}
    for key, val in systematics.items() :
        if 'Nom_WQT' in key : continue
        out[key] = [var[1:] for var in val[0]]
    for key, val in systFlat.items() : #if doesnt work see below explicit expression
        out[key] = [var[1:] for var in val[0]]
    
    out['jme'] = ["jesTotalUp", "jesTotalDown", "unclustEnUp","unclustEnDown"]
    out['ptScale'] = [ "correctedUp", "correctedDown"]
    
    # out['lumi'] = ['lumiUp', 'lumiDown']
    # out['topXSec'] = ['topXSecUp', 'topXSecDown']
    # out['dibosonXSec'] = ['dibosonXSecUp', 'dibosonXSecDown']
    # out['tauXSec'] = ['tauXSecUp', 'tauXSecDown']
    # out['lepVeto'] = ['lepVetoUp', 'lepVetoDown']
    
    # print("dictionary=", out)
    return out 

bkg_systematics = buildBkgSyst()

# ptBinning = [ (25. + 0.5*x)  for x in range(0, 61)] 
# etaBinning = [-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4]
# mtBinning = [0,30,200]

looseCutDict = {
    '15' : [0,15],
    '30' : [15,30],
    '40' : [30,40],
}
