from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt
from root_numpy import hist2array

class reweightyqtWplus(module):
   
    def __init__(self, era, inFilehelwt, genInfoFile):
        self.era=era
        self.inFilehelwt = inFilehelwt
        self.genInfoFile = genInfoFile
            
        print("reweightyqt: helWtFile=",self.inFilehelwt)
        pass
      

    def run(self,d):
        file_in = self.inFilehelwt 
        f = ROOT.TFile.Open(self.genInfoFile) 
        fin = h5py.File(file_in, mode='r+')
        
        qt_powheg = fin['qtycostheta'][:]
        qt_aMC = hist2array(f.Get('angularCoefficients_Wplus/YqTcT'))
        if self.era =="preVFP": qt_aMC*=19.514702645/35.9
        else: qt_aMC*=16.810812618/35.9

        h=np.sum(qt_aMC,axis=-1)/np.sum(qt_powheg,axis=-1)
        
        yBins = fin['edges_qtycostheta_0'][:]
        qtBins = fin['edges_qtycostheta_1'][:]

        if self.era =="preVFP":
            @ROOT.Numba.Declare(["float","float"], "float")
            def getWeightqt_preVFP_Wplus(y,pt):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                return h[biny,binpt]
            self.d = d
            self.d = self.d.Define("yqtweight", "Numba::getWeightqt_preVFP_Wplus(Vrap_preFSR_abs, Vpt_preFSR)")
        else:
            @ROOT.Numba.Declare(["float","float"], "float")
            def getWeightqt_postVFP_Wplus(y,pt):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                return h[biny,binpt]
            self.d = d
            self.d = self.d.Define("yqtweight", "Numba::getWeightqt_postVFP_Wplus(Vrap_preFSR_abs, Vpt_preFSR)")

            # data=self.d.AsNumpy(columns=["yqtweight","Vrap_preFSR_abs", "Vpt_preFSR"])
            # print(np.isinf(data["yqtweight"]).any(),np.isinf(data["Vrap_preFSR_abs"]).any(),np.isinf(data["Vpt_preFSR"]).any())
            # idx=np.where(np.isinf(data["yqtweight"]))
            # print(data["Vrap_preFSR_abs"][idx],data["Vpt_preFSR"][idx])
        return self.d

    def getTH1(self):

        return self.myTH1

    def getTH2(self):

        return self.myTH2  

    def getTH3(self):

        return self.myTH3

    def getTHN(self):

        return self.myTHN

    def getGroupTH1(self):

        return self.myTH1Group

    def getGroupTH2(self):

        return self.myTH2Group  

    def getGroupTH3(self):

        return self.myTH3Group  

    def getGroupTHN(self):

        return self.myTHNGroup

    def reset(self):

        self.myTH1 = []
        self.myTH2 = []
        self.myTH3 = [] 
        self.myTHN = [] 

        self.myTH1Group = []
        self.myTH2Group = []
        self.myTH3Group = [] 
        self.myTHNGroup = [] 
