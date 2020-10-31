bkg_systematics = {
     #"puWeightVars"  : ["puWeightUp", "puWeightDown"],
    #"PrefireWeight"  : ["PrefireWeightUp", "PrefireWeightDown"],
    "WHSF"  : ["WHSFSyst0Up", "WHSFSyst1Up","WHSFSyst2Up","WHSFSystFlatUp","WHSFSyst0Down", "WHSFSyst1Down","WHSFSyst2Down","WHSFSystFlatDown"],
    "LHEScaleWeight" : ["LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0","LHEScaleWeight_muR1p0_muF0p5","LHEScaleWeight_muR1p0_muF2p0","LHEScaleWeight_muR2p0_muF1p0", "LHEScaleWeight_muR2p0_muF2p0"],
    "ptScale" : [ "correctedUp", "correctedDown"], 
    # "jmeVars" : [ "jerUp", "jerDown" , "jesTotalUp", "jesTotalDown", "unclustEnUp","unclustEnDown"],
    "jme" : ["jesTotalUp", "jesTotalDown", "unclustEnUp","unclustEnDown"],
    "LHEPdfWeight" : ["LHEPdfWeightHess{}".format(i+1) for i in range(60)],
    #"alphaSVars" : ["alphaSUp", "alphaSDown"]
}
#bkg_systematics["LHEPdfWeight"].append("alphaSUp")
#bkg_systematics["LHEPdfWeight"].append("alphaSDown")


#ptBinning = [25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55] 
ptBinning = [ (25. + 0.5*x)  for x in range(0, 61)] 
etaBinning = [-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4]
mtBinning = [0,30,200]

looseCutDict = {
    '15' : [0,15],
    '30' : [15,30],
    '40' : [30,40],
}
