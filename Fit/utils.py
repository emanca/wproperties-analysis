from collections import OrderedDict
import math

helXsecs = OrderedDict()
helXsecs["L"] = "A0"
helXsecs["I"] = "A1"
helXsecs["T"] = "A2"
helXsecs["A"] = "A3"
helXsecs["P"] = "A4"
helXsecs["7"] = "A5"
helXsecs["8"] = "A6"
helXsecs["9"] = "A7"
helXsecs["UL"] = "AUL"

factors = OrderedDict()
factors["A0"] = 2.
factors["A1"] = 2.*math.sqrt(2)
factors["A2"] = 4.
factors["A3"] = 4.*math.sqrt(2)
factors["A4"] = 2.
factors["A5"] = 2.
factors["A6"] = 2.*math.sqrt(2)
factors["A7"] = 4.*math.sqrt(2)

import h5py
import numpy as np
def mergeHdf5(fsig,fdata): #then expand with other samples
    fsig_keys = list(fsig.keys())
    with h5py.File('f.hdf5', mode="w") as f:
        for key in fsig_keys:
            if 'mass' in key: continue
            if key == 'templates':
                size = fsig[key][:].shape
                dset = f.create_dataset(name='LowAcc', shape=size, dtype='float64')
                dset[...] = fsig[key][:]
            elif key == 'templates_sumw2':
                size = fsig[key][:].shape
                dset = f.create_dataset(name='LowAcc_sumw2', shape=size, dtype='float64')
                dset[...] = fsig[key][:]
            else:
                size = fsig[key][:].shape
                dset = f.create_dataset(name='{}'.format(key), shape=size, dtype='float64')
                dset[...] = fsig[key][:]
        size = fdata['templates'][:].shape
        dset = f.create_dataset(name='data_obs', shape=size, dtype='float64')
        dset[...] = fdata['templates'][:]
    return
