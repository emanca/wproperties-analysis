from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt

class getHelWeightsWminus(module):
   
    def __init__(self, era,helwtFile,syst = ""):
        self.syst = syst
        if not syst == "":
            self.syst = "_"+syst
        self.era=era
        self.helwtFile=helwtFile
        pass
      

    def run(self,d):
        file_in = self.helwtFile#'/scratchnvme/emanca/wproperties-analysis/config/powheg_acc_{}/WPlusJetsToMuNu_helweights.hdf5'.format(self.era)
        
        f = h5py.File(file_in, mode='r+')
        
        htot = f['totxsecs'+self.syst][:]
        h = f['xsecs'+self.syst][:]

        print(h.shape)

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
        elif self.syst == "_LHEScaleWeight":
            h = h.reshape(len(yBins)-1, len(qtBins)-1, 9, 9)
            factors = factors[...,np.newaxis]
            factors_hel = factors_hel[...,np.newaxis]

        h = (h/htot[:,:,np.newaxis,...]+factors[:,:,:,1,...])*factors[:,:,:,0,...]
        h = h/factors_hel

        if self.era =="preVFP":
            if self.syst == "":
                @ROOT.Numba.Declare(["float", "float"], "RVec<double>")
                def getCoefficients_preVFP_Wminus(y, pt):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    coeff = np.zeros(h.shape[-1])
                    for i in range(h.shape[-1]):
                        coeff[i]=h[biny,binpt,i]
                    return coeff
                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>"], "float")
                def getNorm_preVFP_Wminus(y, pt, coeffs, harms):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    totMap = htot[biny,binpt]
                    norm =harms[-1]*totMap
                    for i in range(coeffs.shape[0]-1):
                        norm += coeffs[i]*harms[i]*totMap
                    norm *=3./16./pi
                    return norm
                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>","float"], "RVec<float>")
                def getWeights_preVFP_Wminus(y, pt, coeffs, harms, norm):
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
                @ROOT.Numba.Declare(["RVec<float>","RVec<float>"], "RVec<float>")
                def multVec_preVFP_Wminus(weights, syst):
                    prod = np.zeros((weights.shape[0],syst.shape[0]),dtype='float32')
                    for i in range(weights.shape[0]):
                        for j in range(syst.shape[0]):
                            prod[i,j]=weights[i]*syst[j]
                    return prod.ravel()
                self.d = d

                self.d = self.d.Define("AngCoeffVec", "Numba::getCoefficients_preVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR)")\
                    .Define("norm", "Numba::getNorm_preVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec,harmonicsVec)")\
                    .Define("helWeights", "Numba::getWeights_preVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec,harmonicsVec,norm)")\
                    .Define("SFStatvar_helweights", "Numba::multVec_preVFP_Wminus(helWeights,SFStatvar)")
            elif self.syst == "_LHEPdfWeight":
                @ROOT.Numba.Declare(["float", "float"], "RVec<double>")
                def getCoefficients_LHEPdfWeight_preVFP_Wminus(y, pt):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    coeff = np.zeros((9,103))
                    for i in range(9):
                        for j in range(103):
                            coeff[i,j]=h[biny,binpt,i,j]
                    return coeff.ravel()
                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>"], "RVec<double>")
                def getNorm_LHEPdfWeight_preVFP_Wminus(y, pt, coeffs, harms):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1

                    coeffs = np.ascontiguousarray(coeffs).reshape((9,103))
                    harms = np.ascontiguousarray(harms).reshape((9,103))
                    totMap = htot[biny,binpt,:]

                    norm =harms[-1,:]*totMap
                    
                    for j in range(103):#pdf
                        for i in range(8):#coeff
                            norm[j] += coeffs[i,j]*harms[i,j]*totMap[j]
                            
                    norm *=3./16./pi
                    return norm.ravel()
                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>","RVec<double>"], "RVec<float>")
                def getWeights_LHEPdfWeight_preVFP_Wminus(y, pt, coeffs, harms, norm):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    coeffs = np.ascontiguousarray(coeffs).reshape((9,103))
                    harms = np.ascontiguousarray(harms).reshape((9,103))
                    totMap = htot[biny,binpt,...]
                    weights = np.zeros((9,103),dtype='float32')
                    for i in range(9):#coeff
                        for j in range(103):#pdf
                            if(norm[j]!=0.):
                                if i!=8:
                                    weights[i,j] = 3./16./pi * totMap[j] * coeffs[i,j] *harms[i,j]/norm[j]
                                else:
                                    weights[i,j] = 3./16./pi * totMap[j] *harms[i,j]/norm[j]
                    return weights.ravel()
                self.d = d
                self.d = self.d.Define("AngCoeffVec_LHEPdfWeight", "Numba::getCoefficients_LHEPdfWeight_preVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR)")\
                    .Define("norm_LHEPdfWeight", "Numba::getNorm_LHEPdfWeight_preVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec_LHEPdfWeight,harmonicsVec_LHEPdfWeight)")\
                    .Define("helWeights_LHEPdfWeight", "Numba::getWeights_LHEPdfWeight_preVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec_LHEPdfWeight,harmonicsVec_LHEPdfWeight,norm_LHEPdfWeight)")
            elif self.syst == "_LHEScaleWeight":
                @ROOT.Numba.Declare(["float", "float"], "RVec<double>")
                def getCoefficients_LHEScaleWeight_preVFP_Wminus(y, pt):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    coeff = np.zeros((9,9))
                    for i in range(9):
                        for j in range(9):
                            coeff[i,j]=h[biny,binpt,i,j]
                    return coeff.ravel()

                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>"], "RVec<double>")
                def getNorm_LHEScaleWeight_preVFP_Wminus(y, pt, coeffs, harms):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1

                    coeffs = np.ascontiguousarray(coeffs).reshape((9,9))
                    harms = np.ascontiguousarray(harms).reshape((9,9))
                    totMap = htot[biny,binpt,...]

                    norm =harms[-1,:]*totMap
                    for i in range(9):#coeff
                        for j in range(9):#pdf
                            norm[j] += coeffs[i,j]*harms[i,j]*totMap[j]
                    norm *=3./16./pi
                    return norm.ravel()

                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>","RVec<double>"], "RVec<float>")
                def getWeights_LHEScaleWeight_preVFP_Wminus(y, pt, coeffs, harms, norm):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    coeffs = np.ascontiguousarray(coeffs).reshape((9,9))
                    harms = np.ascontiguousarray(harms).reshape((9,9))
                    totMap = htot[biny,binpt,...]
                    weights = np.zeros((9,9),dtype='float32')
                    for i in range(9):#coeff
                        for j in range(9):#pdf
                            if(norm[j]!=0.):
                                if i!=8:
                                    weights[i,j] = 3./16./pi * totMap[j] * coeffs[i,j] *harms[i,j]/norm[j]
                                else:
                                    weights[i,j] = 3./16./pi * totMap[j] *harms[i,j]/norm[j]
                    return weights.ravel()
                self.d = d
                self.d = self.d.Define("AngCoeffVec_LHEScaleWeight", "Numba::getCoefficients_LHEScaleWeight_preVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR)")\
                    .Define("norm_LHEScaleWeight", "Numba::getNorm_LHEScaleWeight_preVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec_LHEScaleWeight,harmonicsVec_LHEScaleWeight)")\
                    .Define("helWeights_LHEScaleWeight", "Numba::getWeights_LHEScaleWeight_preVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec_LHEScaleWeight,harmonicsVec_LHEScaleWeight,norm_LHEScaleWeight)")\
                    .Define("nhelWeights_LHEScaleWeight", "helWeights_LHEScaleWeight.size()")
        else:
            if self.syst == "":
                @ROOT.Numba.Declare(["float", "float"], "RVec<double>")
                def getCoefficients_postVFP_Wminus(y, pt):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    coeff = np.zeros(h.shape[-1])
                    for i in range(h.shape[-1]):
                        coeff[i]=h[biny,binpt,i]
                    return coeff
                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>"], "float")
                def getNorm_postVFP_Wminus(y, pt, coeffs, harms):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    totMap = htot[biny,binpt]
                    norm =harms[-1]*totMap
                    for i in range(coeffs.shape[0]-1):
                        norm += coeffs[i]*harms[i]*totMap
                    norm *=3./16./pi
                    return norm
                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>","float"], "RVec<float>")
                def getWeights_postVFP_Wminus(y, pt, coeffs, harms, norm):
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
                @ROOT.Numba.Declare(["RVec<float>","RVec<float>"], "RVec<float>")
                def multVec_postVFP_Wminus(weights, syst):
                    prod = np.zeros((weights.shape[0],syst.shape[0]),dtype='float32')
                    for i in range(weights.shape[0]):
                        for j in range(syst.shape[0]):
                            prod[i,j]=weights[i]*syst[j]
                    return prod.ravel()
                self.d = d

                self.d = self.d.Define("AngCoeffVec", "Numba::getCoefficients_postVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR)")\
                    .Define("norm", "Numba::getNorm_postVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec,harmonicsVec)")\
                    .Define("helWeights", "Numba::getWeights_postVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec,harmonicsVec,norm)")\
                    .Define("SFStatvar_helweights", "Numba::multVec_postVFP_Wminus(helWeights,SFStatvar)")
            elif self.syst == "_LHEPdfWeight":
                @ROOT.Numba.Declare(["float", "float"], "RVec<double>")
                def getCoefficients_LHEPdfWeight_postVFP_Wminus(y, pt):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    coeff = np.zeros((9,103))
                    for i in range(9):
                        for j in range(103):
                            coeff[i,j]=h[biny,binpt,i,j]
                    return coeff.ravel()

                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>"], "RVec<double>")
                def getNorm_LHEPdfWeight_postVFP_Wminus(y, pt, coeffs, harms):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1

                    coeffs = np.ascontiguousarray(coeffs).reshape((9,103))
                    harms = np.ascontiguousarray(harms).reshape((9,103))
                    totMap = htot[biny,binpt,...]

                    norm =harms[-1,:]*totMap

                    for j in range(103):#pdf
                        for i in range(8):#coeff
                            norm[j] += coeffs[i,j]*harms[i,j]*totMap[j]
                    norm *=3./16./pi
                    return norm.ravel()

                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>","RVec<double>"], "RVec<float>")
                def getWeights_LHEPdfWeight_postVFP_Wminus(y, pt, coeffs, harms, norm):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    coeffs = np.ascontiguousarray(coeffs).reshape((9,103))
                    harms = np.ascontiguousarray(harms).reshape((9,103))
                    totMap = htot[biny,binpt,...]
                    weights = np.zeros((9,103),dtype='float32')
                    for i in range(9):#coeff
                        for j in range(103):#pdf
                            if(norm[j]!=0.):
                                if i!=8:
                                    weights[i,j] = 3./16./pi * totMap[j] * coeffs[i,j] *harms[i,j]/norm[j]
                                else:
                                    weights[i,j] = 3./16./pi * totMap[j] *harms[i,j]/norm[j]
                    return weights.ravel()
                self.d = d
                self.d = self.d.Define("AngCoeffVec_LHEPdfWeight", "Numba::getCoefficients_LHEPdfWeight_postVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR)")\
                    .Define("norm_LHEPdfWeight", "Numba::getNorm_LHEPdfWeight_postVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec_LHEPdfWeight,harmonicsVec_LHEPdfWeight)")\
                    .Define("helWeights_LHEPdfWeight", "Numba::getWeights_LHEPdfWeight_postVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec_LHEPdfWeight,harmonicsVec_LHEPdfWeight,norm_LHEPdfWeight)")\
                    .Define("nhelWeights_LHEPdfWeight", "helWeights_LHEPdfWeight.size()")

            elif self.syst == "_LHEScaleWeight":
                @ROOT.Numba.Declare(["float", "float"], "RVec<double>")
                def getCoefficients_LHEScaleWeight_postVFP_Wminus(y, pt):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    coeff = np.zeros((9,9))
                    for i in range(9):
                        for j in range(9):
                            coeff[i,j]=h[biny,binpt,i,j]
                    return coeff.ravel()

                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>"], "RVec<double>")
                def getNorm_LHEScaleWeight_postVFP_Wminus(y, pt, coeffs, harms):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1

                    coeffs = np.ascontiguousarray(coeffs).reshape((9,9))
                    harms = np.ascontiguousarray(harms).reshape((9,9))
                    totMap = htot[biny,binpt,...]

                    norm =harms[-1,:]*totMap
                    for i in range(9):#coeff
                        for j in range(9):#pdf
                            norm[j] += coeffs[i,j]*harms[i,j]*totMap[j]
                    norm *=3./16./pi
                    return norm.ravel()

                @ROOT.Numba.Declare(["float", "float", "RVec<double>", "RVec<float>","RVec<double>"], "RVec<float>")
                def getWeights_LHEScaleWeight_postVFP_Wminus(y, pt, coeffs, harms, norm):
                    biny = np.digitize(np.array([y]), yBins)[0]-1
                    binpt = np.digitize(np.array([pt]), qtBins)[0]-1
                    coeffs = np.ascontiguousarray(coeffs).reshape((9,9))
                    harms = np.ascontiguousarray(harms).reshape((9,9))
                    totMap = htot[biny,binpt,...]
                    weights = np.zeros((9,9),dtype='float32')
                    for i in range(9):#coeff
                        for j in range(9):#pdf
                            if(norm[j]!=0.):
                                if i!=8:
                                    weights[i,j] = 3./16./pi * totMap[j] * coeffs[i,j] *harms[i,j]/norm[j]
                                else:
                                    weights[i,j] = 3./16./pi * totMap[j] *harms[i,j]/norm[j]
                    return weights.ravel()
                self.d = d
                self.d = self.d.Define("AngCoeffVec_LHEScaleWeight", "Numba::getCoefficients_LHEScaleWeight_postVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR)")\
                    .Define("norm_LHEScaleWeight", "Numba::getNorm_LHEScaleWeight_postVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec_LHEScaleWeight,harmonicsVec_LHEScaleWeight)")\
                    .Define("helWeights_LHEScaleWeight", "Numba::getWeights_LHEScaleWeight_postVFP_Wminus(Vrap_preFSR_abs,Vpt_preFSR,AngCoeffVec_LHEScaleWeight,harmonicsVec_LHEScaleWeight,norm_LHEScaleWeight)")\
                    .Define("nhelWeights_LHEScaleWeight", "helWeights_LHEScaleWeight.size()")
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
