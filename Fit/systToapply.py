# ["DY","Diboson","Top","fakesLowMt","fakesHighMt", "Wtau","LowAcc"]
systematicsDict = {
    "CMSlumi": {
       "vars":["CMSlumi"],
       "procs": ["Signal","LowAcc"],
       "type": "lnN",
       "weight" : 1.012
    },
    # "Topxsec":{
    #    "vars":["Topxsec"],
    #    "procs": ["Top"],
    #    "type": "lnN",
    #    "weight" : 1.060
    # },
    # "Dibosonxsec":{
    #    "vars":["Dibosonxsec"],
    #    "procs": ["Diboson"],
    #    "type": "lnN",
    #    "weight" : 1.160
    # },
    #  "Tauxsec":{
    #    "vars":["Tauxsec"],
    #    "procs": ["Wtau"],
    #    "type": "lnN",
    #    "weight" : 1.04
    # },
    "mass" : {
        "vars":["mass"],
        "procs": ["Signal","LowAcc"],
        "type": "shapeNoConstraint",
        "weight" : 1.
    },
#     "LHEScaleWeight1": {
#         "vars": ["LHEScaleWeight1_muRmuF", "LHEScaleWeight1_muR","LHEScaleWeight1_muF"],
#         "procs": ["LowAcc"],
#         "type": "shape",
#         "weight" : 1.
#     },
#     "LHEScaleWeight2": {
#         "vars": ["LHEScaleWeight2_muRmuF", "LHEScaleWeight2_muR","LHEScaleWeight2_muF"],
#         "procs": ["LowAcc"],
#         "type": "shape",
#         "weight" : 1.
#     },
#     "LHEScaleWeight3": {
#         "vars": ["LHEScaleWeight3_muRmuF", "LHEScaleWeight3_muR","LHEScaleWeight3_muF"],
#         "procs": ["LowAcc"],
#         "type": "shape",
#         "weight" : 1.
#     },
#     "LHEScaleWeight4": {
#         "vars": ["LHEScaleWeight4_muRmuF", "LHEScaleWeight4_muR","LHEScaleWeight4_muF"],
#         "procs": ["LowAcc"],
#         "type": "shape",
#         "weight" : 1.
#     },
#     "LHEScaleWeight": {
#         "vars": ["muRmuF", "muR","muF"],
#         "procs": ["DY","Wtau"],
#         "type": "shape",
#         "weight" : 1.
#     },
#     "SF"  : {
#         "vars": ["SFall{}".format(i) for i in range(624)]+["SFiso{}".format(i) for i in range(624)]+["SFSyst"],
#        "procs": ["Signal","DY","Diboson","Top","Wtau","LowAcc"],
#        "type": "shape",
#        "weight" : 1.
#     },
#     "jme" : {
#        "vars":["jes", "uncl"],
#        "procs": ["Signal","DY","Diboson","Top","Wtau","LowAcc"],
#        "type": "shape",
#        "weight" : 1.
#     },
#     "PrefireWeight":{
#        "vars":["prefire"],
#        "procs": ["Signal", "DY","Diboson","Top","Wtau","LowAcc"],
#        "type": "shape",
#        "weight" : 1.
#     },
#     "LHEPdfWeight" : {
#        "vars":["pdf{}".format(i) for i in range(1,103)],
#        "procs": ["DY","Wtau","LowAcc"],
#        "type": "shape",
#        "weight" : 1.
#     },

}
