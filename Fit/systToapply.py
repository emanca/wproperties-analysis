systematicsDict = {
    "Nominal": {
    },
    "CMSlumi": {
       "vars":["CMSlumi"],
       "procs": ["Signal", "DYJets","DiBoson","Top","WtoTau","LowAcc","Fake"],
      #  "procs": [ "DYJets","DiBoson","Top","WtoTau","LowAcc","Fake"],
       "type": "lnN",
       "weight" : 1.025 
    },
    # "DYxsec":{
    #    "vars":["DYxsec"],
    #    "procs": ["DYJets"],
    #    "type": "lnN",
    #    "weight" : 1.038 
    # },
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
        "procs": ["Signal", "LowAcc"],
        # "procs": ["Signal"],
        "type": "shapeNoConstraint",
        # "type": "shape",
        "weight" : 1.
    },
    "WHSFStat"  : {
       "vars": ["WHSFSyst0Eta{}".format(i) for i in range(1, 49)]+["WHSFSyst1Eta{}".format(i) for i in range(1, 49)]+["WHSFSyst2Eta{}".format(i) for i in range(1, 49)],
       "procs": ["Signal", "DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"],
      #  "procs": [ "DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    "WHSFSyst":{
       "vars": ["WHSFSystFlat"],
       "procs": ["Signal","DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"],
      #  "procs": ["DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    "jme" : {
       "vars":["jesTotal", "unclustEn"],
       "procs": ["Signal", "DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"],
      #  "procs": ["DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    "PrefireWeight":{
       "vars":["PrefireWeight"],
       "procs": ["Signal", "DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"],
      #  "procs": ["DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    "LHEPdfWeight" : {
       "vars":["LHEPdfWeightHess{}".format(i+1) for i in range(60)],
       "procs": ["Signal", "DYJets", "Fake", "WtoTau", "LowAcc"],
      #  "procs": [ "DYJets", "Fake", "WtoTau", "LowAcc"],
      #  "procs": ["Signal", "DYJets", "Fake", "WtoTau"],
      #  "procs": ["LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
    "alphaS" :{
       "vars": ["alphaS"],
       "procs": ["Signal", "DYJets", "Fake", "WtoTau", "LowAcc"],
      #  "procs": ["DYJets", "Fake", "WtoTau", "LowAcc"],
       "type": "shape",
       "weight" : 1.
    },
     "LHEScaleWeight" :{
       "vars": ["LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0","LHEScaleWeight_muR1p0_muF0p5","LHEScaleWeight_muR1p0_muF2p0","LHEScaleWeight_muR2p0_muF1p0","LHEScaleWeight_muR2p0_muF2p0"],
       "procs": ["DYJets"],
      #  "procs": ["DYJets", "Fake", "WtoTau", "LowAcc"],
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
        #"ptScale" : {
    #    "vars": ["Eta{}zptsyst".format(j) for j in range(1, 5)] + ["Eta{}Ewksyst".format(j) for j in range(1, 5)] + ["Eta{}deltaMsyst".format(j) for j in range(1, 5)]+["Eta{}stateig{}".format(j, i) for i in range(0, 99) for j in range(1, 5)],
    #    "procs": ["Signal"],
    #    "type": "shape"
    #},
    
    "LeptonVeto":{
       "vars":["LeptonVeto"],
       "procs": ["DYJets"],
       "type": "lnN",
       "weight" : 1.020 
    },
    
    "FakeNorm":{
       "vars":["FakeNorm"],
       "procs": ["Fake"],
       "type": "lnN",
       "weight" : 1.2
       #1.05 
    },
}