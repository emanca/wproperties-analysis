from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt
from root_numpy import hist2array

class reweightcoeffsWplus(module):
   
    def __init__(self, era, helWtsrcdir, geninputF):
        self.era=era
        self.helwtsrcDir=helWtsrcdir
        self.getnInfputFile=geninputF
        pass
      

    def run(self,d):

        # load powheg coefficients
        file_preVFP = self.helwtsrcDir + '/outputW_sroychow_preVFP/WPlusJetsToMuNu_helweights.hdf5'
        file_postVFP = self.helwtsrcDir + '/outputW_sroychow_postVFP/WPlusJetsToMuNu_helweights.hdf5'

        f_preVFP = h5py.File(file_preVFP, mode='r+')
        f_postVFP = h5py.File(file_postVFP, mode='r+')
        htot_preVFP = f_preVFP['totxsecs'][:]
        htot_postVFP = f_postVFP['totxsecs'][:]
        h_preVFP = f_preVFP['xsecs'][:]
        h_postVFP = f_postVFP['xsecs'][:]

        htot = htot_preVFP+htot_postVFP
        h = h_preVFP+h_postVFP
        # shape h: y, qt, weights, pdf
        # shape tot: y, qt, pdf
        factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
        factors = factors[np.newaxis,np.newaxis,...]
        h = (h/htot[...,np.newaxis]+factors[...,1])*factors[...,0]
        coeffs_powheg = h[...,:-4] #remove A5,A6,A7,UL coefficient

        # load aMC@NLO coefficients
        fcoeffs = ROOT.TFile.Open(self.getnInfputFile)
        
        hists_aMC = []
        for i in range(5):
            tmp = hist2array(fcoeffs.Get('angularCoefficients/harmonicsA{}_nom_nom'.format(i)))
            hists_aMC.append(tmp)

        coeffs_aMC_stack = np.stack(hists_aMC,axis=-1)
        coeffs_aMC = np.copy(coeffs_powheg)
        coeffs_aMC[...,4]=coeffs_aMC_stack[...,4]
        # print(coeffs_powheg[...,0], 'powheg')
        # print(coeffs_aMC[...,0], 'aMC')

        yBins = f_preVFP['edges_xsecs_0'][:]
        qtBins = f_preVFP['edges_xsecs_1'][:]

        if self.era =="preVFP":
            @ROOT.Numba.Declare(["float", "float", "RVec<float>"], "float")
            def getNormRatio_Wplus(y, pt,harms):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                binpt = np.digitize(np.array([pt]), qtBins)[0]-1
    
                norm_aMC =harms[-1] #1+cos^2 theta
                norm_powheg =harms[-1] #1+cos^2 theta
                for i in range(coeffs_aMC.shape[0]):
                    norm_aMC += coeffs_aMC[biny,binpt,i]*harms[i]
                    norm_powheg += coeffs_powheg[biny,binpt,i]*harms[i]
                return norm_aMC/norm_powheg
        self.d = d
        self.d = self.d.Define("coeffsweight", "Numba::getNormRatio_Wplus(Vrap_preFSR_abs, Vpt_preFSR,harmonicsVec)")
        # data=self.d.AsNumpy(columns=["coeffsweight","Vrap_preFSR_abs", "Vpt_preFSR"])
        # print(np.isinf(data["coeffsweight"]).any(),np.isinf(data["Vrap_preFSR_abs"]).any(),np.isinf(data["Vpt_preFSR"]).any())
        # idx=np.where(np.isinf(data["coeffsweight"]))
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
