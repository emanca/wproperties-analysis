from module import *
import h5py
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

class reweightycostheta(module):
   
    def __init__(self, era):
        self.era=era
        pass
      

    def run(self,d):

        file_in = '/scratchnvme/emanca/wproperties-analysis/config/valeriosbins_{}/WPlusJetsToMuNu_costheta.hdf5'.format(self.era)
        file_in_0J = '/scratchnvme/emanca/wproperties-analysis/config/alternatesample_{}/WJetsToLNu_0J_helweights.hdf5'.format(self.era)
        file_in_1J = '/scratchnvme/emanca/wproperties-analysis/config/alternatesample_{}/WJetsToLNu_1J_helweights.hdf5'.format(self.era)
        file_in_2J = '/scratchnvme/emanca/wproperties-analysis/config/alternatesample_{}/WJetsToLNu_2J_helweights.hdf5'.format(self.era)

        f = h5py.File(file_in, mode='r+')
        f0J = h5py.File(file_in_0J, mode='r+')
        f1J = h5py.File(file_in_1J, mode='r+')
        f2J = h5py.File(file_in_2J, mode='r+')
        
        costheta_powheg = f['costheta'][:]
        costheta_aMC = f0J['costheta'][:]+f1J['costheta'][:]+f2J['costheta'][:]

        cosThetaBins = np.array([round(-1. + 2.*i/100,2) for i in range(101)])
        yBins = np.array([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0])

        fig, ax1 = plt.subplots()
        ax1.set_title("costheta powheg", fontsize=18)
        hep.hist2dplot(costheta_powheg,yBins,cosThetaBins)
        plt.tight_layout()
        plt.savefig("costheta_powheg_{}".format(self.era))
        plt.clf()

        fig, ax1 = plt.subplots()
        ax1.set_title("costheta aMC@NLO", fontsize=18)
        hep.hist2dplot(costheta_aMC,yBins,cosThetaBins)
        plt.tight_layout()
        plt.savefig("costheta_aMC{}".format(self.era))
        plt.clf()

        h = costheta_aMC/costheta_powheg
        # print(h)
        if self.era =="preVFP":
            @ROOT.Numba.Declare(["float","float"], "float")
            def getWeight_preVFP(y,costheta):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                bincostheta = np.digitize(np.array([costheta]), cosThetaBins)[0]-1
                if not (y==4.717765 and costheta==1.):
                    return h[biny,bincostheta]
                else: return 1.
            self.d = d
            self.d = self.d.Define("ycosthetaweight", "Numba::getWeight_preVFP(Vrap_preFSR_abs,CStheta_preFSR)")#.Filter("ycosthetaweight>2.")
        else:
            @ROOT.Numba.Declare(["float","float"], "float")
            def getWeight_postVFP(y,costheta):
                biny = np.digitize(np.array([y]), yBins)[0]-1
                bincostheta = np.digitize(np.array([costheta]), cosThetaBins)[0]-1
                if not (y==4.717765 and costheta==1.):
                    return h[biny,bincostheta]
                else: return 1.
            self.d = d
            self.d = self.d.Define("ycosthetaweight", "Numba::getWeight_postVFP(Vrap_preFSR_abs,CStheta_preFSR)")#.Filter("ycosthetaweight>2.")
        # node = self.d.AsNumpy()
        # print('getWeight({},{}) = {}'.format(node['Vrap_preFSR_abs'],node['CStheta_preFSR'], node['ycosthetaweight']))

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
