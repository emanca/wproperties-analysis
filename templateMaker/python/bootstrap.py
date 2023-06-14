from module import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt

ROOT.gInterpreter.Declare('#include "interface/helSystHelper.h"')

class bootstrap(module):
   
    def __init__(self):
        pass
      
    def run(self,d):

        self.d = d

        # this function creates a vector of 100 poisson random numbers per event
        ROOT.gInterpreter.Declare("""
            ROOT::VecOps::RVec<double> bootstrapping( const int & _rdfentry, const int & _rdfslot ){
            TRandom3 rnd;   
            rnd.SetSeed(_rdfentry*_rdfslot);
            ROOT::VecOps::RVec<double> rndPoisson;
            rndPoisson.reserve(100);
            for( int i =0; i < 100; ++i){
                rndPoisson.push_back((double)rnd.PoissonD(1));
            }
            
            return rndPoisson;
        };     
        """)
        self.d = self.d.Define("rndPoisson", "bootstrapping(rdfentry_,rdfslot_)")\
                    .Define("rndPoisson_tensor", f"Eigen::TensorFixedSize<double, Eigen::Sizes<100>>(wrem::vec_to_tensor_t<double, 100>(rndPoisson))")
        self.helper = ROOT.helSystHelper[100,6]()
        self.d = self.d.Define("rndPoisson_tensor_hel",self.helper,["helWeightTensor", "rndPoisson_tensor"])
        self.helper = ROOT.helMassPoissHelper[2,100]()
        self.d = self.d.Define("rndPoisson_tensor_hel_mass",self.helper,["massWeight_tensor_hel", "rndPoisson_tensor"])

        return self.d
