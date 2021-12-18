systematics = {
    #####"puWeight"  : (["puWeightUp", "puWeightDown"],"puWeightVars"),
    "PrefireWeight"  : (["_PrefireWeightUp", "_PrefireWeightDown"],"PrefireWeightVars"),
    "WHSF"  : (["_WHSFSyst0Up", "_WHSFSyst1Up","_WHSFSyst2Up","_WHSFSystFlatUp","_WHSFSyst0Down", "_WHSFSyst1Down","_WHSFSyst2Down","_WHSFSystFlatDown"],"WHSFVars"),
    "LHEScaleWeight" : (["_LHEScaleWeight_muR0p5_muF0p5", "_LHEScaleWeight_muR0p5_muF1p0","_LHEScaleWeight_muR1p0_muF0p5","_LHEScaleWeight_muR1p0_muF2p0","_LHEScaleWeight_muR2p0_muF1p0","_LHEScaleWeight_muR2p0_muF2p0"], "LHEScaleWeightred"),
    "LHEPdfWeight" : (["_LHEPdfWeightHess{}".format(i+1) for i in range(60)], "LHEPdfWeightHess"),
    "alphaS"       : (["_alphaSUp", "_alphaSDown"], "alphaSVars"),
    "mass"         : (["_massUp","_massDown"], "massWeights"),
    "LHEScaleWeight_WQTlow" : (["_LHEScaleWeight_muR0p5_muF0p5_WQTlow", "_LHEScaleWeight_muR0p5_muF1p0_WQTlow","_LHEScaleWeight_muR1p0_muF0p5_WQTlow","_LHEScaleWeight_muR1p0_muF2p0_WQTlow","_LHEScaleWeight_muR2p0_muF1p0_WQTlow","_LHEScaleWeight_muR2p0_muF2p0_WQTlow"], "LHEScaleWeightred", "&& Wpt_preFSR<5"),
    "LHEScaleWeight_WQTmid" : (["_LHEScaleWeight_muR0p5_muF0p5_WQTmid", "_LHEScaleWeight_muR0p5_muF1p0_WQTmid","_LHEScaleWeight_muR1p0_muF0p5_WQTmid","_LHEScaleWeight_muR1p0_muF2p0_WQTmid","_LHEScaleWeight_muR2p0_muF1p0_WQTmid","_LHEScaleWeight_muR2p0_muF2p0_WQTmid"], "LHEScaleWeightred", "&& Wpt_preFSR>5 && Wpt_preFSR<15"),
    "LHEScaleWeight_WQThigh" : (["_LHEScaleWeight_muR0p5_muF0p5_WQThigh", "_LHEScaleWeight_muR0p5_muF1p0_WQThigh","_LHEScaleWeight_muR1p0_muF0p5_WQThigh","_LHEScaleWeight_muR1p0_muF2p0_WQThigh","_LHEScaleWeight_muR2p0_muF1p0_WQThigh","_LHEScaleWeight_muR2p0_muF2p0_WQThigh"], "LHEScaleWeightred", "&& Wpt_preFSR>15"),
    "Nom_WQTlow" : ([""], "Nom", "&& Wpt_preFSR<5"),
    "Nom_WQTmid" : ([""], "Nom", "&& Wpt_preFSR>5 && Wpt_preFSR<15"),
    "Nom_WQThigh" : ([""], "Nom", "&& Wpt_preFSR>15"),
}

systFlat = {
    'lumi' : (['_lumiUp', '_lumiDown'], 'lumi'),
    'topXSec' : (['_topXSecUp', '_topXSecDown'], 'topXSec'),
    'dibosonXSec' : (['_dibosonXSecUp', '_dibosonXSecDown'], 'dibosonXSec'),
    'tauXSec' : (['_tauXSecUp', '_tauXSecDown'], 'tauXSec'),
    'lepVeto' : (['_lepVetoUp', '_lepVetoDown'], 'lepVeto'),
}
