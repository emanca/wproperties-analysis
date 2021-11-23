systematicsDict = {
    "Nominal": {
    },
    "CMSlumi": {
       "vars":["CMSlumi"],
       "procs": ["Signal", "DYJets","DiBoson","Top","WtoTau","LowAcc"], #Fake
       "type": "lnN",
       "weight" : 1.025 
    },
    
    "Topxsec":{
       "vars":["Topxsec"],
       "procs": ["Top"],
       "type": "lnN",
       "weight" : 1.060 
    },
    "Dibosonxsec":{
       "vars":["Dibosonxsec"],
       "procs": ["DiBoson"],
       "type": "lnN",
       "weight" : 1.160 
    },
     "Tauxsec":{
       "vars":["Tauxsec"],
       "procs": ["WtoTau"],
       "type": "lnN",
       "weight" : 1.04 
    },
    
    "mass" : {
        "vars":["mass"],
        "procs": ["Signal", "LowAcc", "Fake"],
        "type": "shapeNoConstraint",
        # "type": "shape",
        "weight" : 1.
    },
    
    "WplusWHSFStat"  : {
       "vars": ["WplusWHSFSyst0Eta{}".format(i) for i in range(1, 49)]+["WplusWHSFSyst1Eta{}".format(i) for i in range(1, 49)]+["WplusWHSFSyst2Eta{}".format(i) for i in range(1, 49)],
      # "vars": ["WplusWHSFSyst0Eta{}".format(i) for i in range(1, 15)]+["WplusWHSFSyst1Eta{}".format(i) for i in range(1, 15)]+["WplusWHSFSyst2Eta{}".format(i) for i in range(1, 15)],
      #  "vars": ["WplusWHSFSyst0Eta{}".format(i) for i in range(1, 10)],
       "procs": ["WplusSignal", "WplusDYJets","WplusDiBoson","WplusTop","WplusFake","WplusWtoTau","WplusLowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    
    "WplusWHSFSyst":{
       "vars": ["WplusWHSFSystFlat"],
       "procs": ["WplusSignal","WplusDYJets","WplusDiBoson","WplusTop","WplusFake","WplusWtoTau","WplusLowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    "WminusWHSFStat"  : {
       "vars": ["WminusWHSFSyst0Eta{}".format(i) for i in range(1, 49)]+["WminusWHSFSyst1Eta{}".format(i) for i in range(1, 49)]+["WminusWHSFSyst2Eta{}".format(i) for i in range(1, 49)],
      #  "vars": ["WminusWHSFSyst0Eta{}".format(i) for i in range(1, 15)]+["WminusWHSFSyst1Eta{}".format(i) for i in range(1, 15)]+["WminusWHSFSyst2Eta{}".format(i) for i in range(1, 15)],
      #  "vars": ["WminusWHSFSyst0Eta{}".format(i) for i in range(1, 10)],
       "procs": ["WminusSignal", "WminusDYJets","WminusDiBoson","WminusTop","WminusFake","WminusWtoTau","WminusLowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    "WminusWHSFSyst":{
       "vars": ["WminusWHSFSystFlat"],
       "procs": ["WminusSignal","WminusDYJets","WminusDiBoson","WminusTop","WminusFake","WminusWtoTau","WminusLowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    
    "Wplusjme" : {
      "vars":["WplusjesTotal", "WplusunclustEn"],
      "procs": ["WplusSignal", "WplusDYJets","WplusDiBoson","WplusTop","WplusFake","WplusWtoTau","WplusLowAcc"],
      "type": "shape",
      "weight" : 1.
    },
   "Wminusjme" : {
      "vars":["WminusjesTotal", "WminusunclustEn"],
      "procs": ["WminusSignal", "WminusDYJets","WminusDiBoson","WminusTop","WminusFake","WminusWtoTau","WminusLowAcc"],
      "type": "shape",
      "weight" : 1.
    },
    
    "WplusPrefireWeight":{
       "vars":["WplusPrefireWeight"],
       "procs": ["WplusSignal", "WplusDYJets","WplusDiBoson","WplusTop","WplusFake","WplusWtoTau","WplusLowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    "WminusPrefireWeight":{
       "vars":["WminusPrefireWeight"],
       "procs": ["WminusSignal", "WminusDYJets","WminusDiBoson","WminusTop","WminusFake","WminusWtoTau","WminusLowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    
    "LHEPdfWeight" : {
       "vars":["LHEPdfWeightHess{}".format(i+1) for i in range(60)],
       "procs": ["Signal", "DYJets", "Fake", "WtoTau", "LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    "alphaS" :{
       "vars": ["alphaS"],
       "procs": ["Signal", "DYJets", "Fake", "WtoTau", "LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    
     "LHEScaleWeight" :{
       "vars": ["LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0","LHEScaleWeight_muR1p0_muF0p5","LHEScaleWeight_muR1p0_muF2p0","LHEScaleWeight_muR2p0_muF1p0","LHEScaleWeight_muR2p0_muF2p0"],
       "procs": ["DYJets"],
       "type": "shape",
       "weight" : 1.
    },
     "LHEScaleWeight_WQTlow" :{
       "vars": ["LHEScaleWeight_muR0p5_muF0p5_WQTlow", "LHEScaleWeight_muR0p5_muF1p0_WQTlow","LHEScaleWeight_muR1p0_muF0p5_WQTlow","LHEScaleWeight_muR1p0_muF2p0_WQTlow","LHEScaleWeight_muR2p0_muF1p0_WQTlow","LHEScaleWeight_muR2p0_muF2p0_WQTlow"],
       "procs": ["Fake", "WtoTau", "LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
     "LHEScaleWeight_WQTmid" :{
       "vars": ["LHEScaleWeight_muR0p5_muF0p5_WQTmid", "LHEScaleWeight_muR0p5_muF1p0_WQTmid","LHEScaleWeight_muR1p0_muF0p5_WQTmid","LHEScaleWeight_muR1p0_muF2p0_WQTmid","LHEScaleWeight_muR2p0_muF1p0_WQTmid","LHEScaleWeight_muR2p0_muF2p0_WQTmid"],
       "procs": ["Fake", "WtoTau", "LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
     "LHEScaleWeight_WQThigh" :{
       "vars": ["LHEScaleWeight_muR0p5_muF0p5_WQThigh", "LHEScaleWeight_muR0p5_muF1p0_WQThigh","LHEScaleWeight_muR1p0_muF0p5_WQThigh","LHEScaleWeight_muR1p0_muF2p0_WQThigh","LHEScaleWeight_muR2p0_muF1p0_WQThigh","LHEScaleWeight_muR2p0_muF2p0_WQThigh"],
       "procs": ["Fake", "WtoTau", "LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    
    "WplusptScale" : {
       "vars": ["Wpluscorrected"],
       "procs": ["WplusSignal", "WplusDYJets","WplusDiBoson","WplusTop","WplusFake","WplusWtoTau","WplusLowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    "WminusptScale" : {
       "vars": ["Wminuscorrected"],
       "procs": ["WminusSignal", "WminusDYJets","WminusDiBoson","WminusTop","WminusFake","WminusWtoTau","WminusLowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    
    "WplusLeptonVeto":{
       "vars":["WplusLeptonVeto"],
       "procs": ["WplusDYJets"],
       "type": "lnN",
       "weight" : 1.020 
    },
    "WminusLeptonVeto":{
       "vars":["WminusLeptonVeto"],
       "procs": ["WminusDYJets"],
       "type": "lnN",
       "weight" : 1.020 
    },
    
    "FakeLumi":{
       "vars":["lumi"],
       "procs": ["Fake"],
       "type": "shape",
       "weight" : 1.
    },
     "WplusQCDnorm":{
       "vars":["WplusQCDnorm"],
       "procs": ["WplusFake"],
       "type": "lnN",
       "weight" : 1.05
    },
     "WminusQCDnorm":{
       "vars":["WminusQCDnorm"],
       "procs": ["WminusFake"],
       "type": "lnN",
       "weight" : 1.05
    },
}