from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt

class getHelWeights(module):
   
    def __init__(self, era,syst = ""):
        self.syst = syst
        if not syst == "":
            self.syst = "_"+syst
        self.era=era
        pass
      

    def run(self,d):

        file_in = '/scratchnvme/emanca/wproperties-analysis/config/halfptBins_{}/WPlusJetsToMuNu_helweights.hdf5'.format(self.era)

        f = h5py.File(file_in, mode='r+')
        
        htot = f['totxsecs'+self.syst][:]
        h = f['xsecs'+self.syst][:]

        yBins = f['edges_totxsecs_0'][:]
        qtBins = f['edges_totxsecs_1'][:]

        # shape h: y, qt, weights, pdf
        # shape tot: y, qt, pdf
        factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
        factors = factors[np.newaxis,np.newaxis,...]
        factors_hel = np.array([2.,2*sqrt(2),4.,4.*sqrt(2),2.,2.,2.*sqrt(2),4.*sqrt(2),1.])
        factors_hel = factors_hel[np.newaxis,np.newaxis,...]
        
        if self.syst == "_LHEPdfWeight":
            h = h.reshape(len(yBins)-1, len(qtBins)-1, 9, 103)
            factors = factors[...,np.newaxis]
            factors_hel = factors_hel[...,np.newaxis]

        h = (h/htot[:,:,np.newaxis,...]+factors[:,:,:,1,...])*factors[:,:,:,0,...]
        h = h/factors_hel

        if self.syst == "":
            @ROOT.Numba.Declare(["float", "float"], "RVec<double>")
            def getCoefficients(y, pt):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                coeff = np.zeros(h.shape[-1])
                for i in range(h.shape[-1]):
                    coeff[i]=h[biny,binpt,i]
                return coeff

            @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>"], "float")
            def getNorm(y, pt, coeffs, harms):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                totMap = htot[biny,binpt]
                norm =harms[-1]*totMap
                for i in range(coeffs.shape[0]-1):
                    norm += coeffs[i]*harms[i]*totMap
                norm *=3./16./pi
                return norm

            @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>","float"], "RVec<float>")
            def getWeights(y, pt, coeffs, harms, norm):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                totMap = htot[biny,binpt]
                weights = np.zeros(h.shape[-1],dtype='float32')
                for i in range(h.shape[-1]):
                    if(norm!=0.):
                        if i!=8:
                            weights[i] = 3./16./pi * totMap * coeffs[i] *harms[i]/norm
                        else:
                            weights[i] = 3./16./pi * totMap *harms[i]/norm
                return weights
            self.d = d
                
            self.d = self.d.Define("AngCoeffVec", "Numba::getCoefficients(Vrap_preFSR_abs,Vpt_preFSR)")\
                .Define("norm", "Numba::getNorm(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec,harmonicsVec)")\
                .Define("helWeights", "Numba::getWeights(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec,harmonicsVec,norm)")\
                .Define("nhelWeights", "helWeights.size()")
        else:
            @ROOT.Numba.Declare(["float", "float"], "RVec<double>")
            def getCoefficients_LHEPdfWeight(y, pt):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                coeff = np.zeros((h.shape[-1],h.shape[2]))
                for j in range(h.shape[-1]):
                    for i in range(h.shape[2]):
                        coeff[j,i]=h[biny,binpt,i,j]
                return coeff.ravel()

            @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>"], "RVec<double>")
            def getNorm_LHEPdfWeight(y, pt, coeffs, harms):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                binpt = np.digitize(np.array([pt]), qtBins)[0]-1

                coeffs = np.ascontiguousarray(coeffs).reshape((h.shape[-1],h.shape[2]))
                harms = np.ascontiguousarray(harms).reshape((h.shape[-1],h.shape[2]))
                totMap = htot[biny,binpt,...]
                # harms (103,9) 
                # totMap (103)
                norm =harms[:,-1]*totMap
                for i in range(coeffs.shape[1]-1):
                    for j in range(coeffs.shape[0]):
                        norm[i] += coeffs[j,i]*harms[j,i]*totMap[j]
                        pass
                norm *=3./16./pi
                return norm.ravel()

            @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>","RVec<double>"], "RVec<float>")
            def getWeights_LHEPdfWeight(y, pt, coeffs, harms, norm):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                coeffs = np.ascontiguousarray(coeffs).reshape((h.shape[-1],h.shape[2]))
                harms = np.ascontiguousarray(harms).reshape((h.shape[-1],h.shape[2]))
                totMap = htot[biny,binpt,...]
                weights = np.zeros((h.shape[-1],h.shape[2]),dtype='float32')
                for j in range(h.shape[-1]):
                    for i in range(h.shape[2]):
                        if(norm[j]!=0.):
                            if i!=8:
                                weights[j,i] = 3./16./pi * totMap[j] * coeffs[j,i] *harms[j,i]/norm[j]
                            else:
                                weights[j,i] = 3./16./pi * totMap[j] *harms[j,i]/norm[j]
                return weights.ravel()
            self.d = d
            self.d = self.d.Define("AngCoeffVec_LHEPdfWeight", "Numba::getCoefficients_LHEPdfWeight(Vrap_preFSR_abs,Vpt_preFSR)")\
                .Define("norm_LHEPdfWeight", "Numba::getNorm_LHEPdfWeight(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec_LHEPdfWeight,harmonicsVec_LHEPdfWeight)")\
                .Define("helWeights_LHEPdfWeight", "Numba::getWeights_LHEPdfWeight(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec_LHEPdfWeight,harmonicsVec_LHEPdfWeight,norm_LHEPdfWeight)")\
                .Define("nhelWeights_LHEPdfWeight", "helWeights_LHEPdfWeight.size()")
        
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
