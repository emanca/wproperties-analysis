import os
import sys
sys.path.append('data/')
import h5py
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from binning import ptBins, etaBins, mTBins, isoBins, chargeBins
plt.style.use([hep.style.ROOT])
hep.cms.label(loc=0, year=2016, lumi=35.9, data=True)
#hep.cms.text('Simulation')

folder = "output_bkg"
os.chdir(folder)

ewkFiles = ["DYJetsToLL_M10to50.hdf5", "ST_t-channel_antitop_4f_inclusiveDecays.hdf5","ST_tW_top_5f_inclusiveDecays.hdf5","TTJets_SingleLeptFromTbar.hdf5","WZ.hdf5",\
"DYJetsToLL_M50.hdf5","ST_t-channel_top_4f_inclusiveDecays_13TeV.hdf5","TTJets_DiLept.hdf5","WJetsToLNu.hdf5","ZZ.hdf5",\
"ST_s-channel_4f_leptonDecays.hdf5","ST_tW_antitop_5f_inclusiveDecays.hdf5","TTJets_SingleLeptFromT.hdf5","WW.hdf5"]

histonames = ['ewk', 'ewk_sumw2']

def haddFiles(fileList, histonames, shape):

    f = h5py.File('ewk.hdf5', mode='w')
    for name in histonames:
        #print(name, shape)
        dset = f.create_dataset(name=name, shape=[shape], dtype='float64')
        tmp = np.zeros([shape],dtype='float64')
        for file in fileList:
            ftmp = h5py.File(file, mode='r+')
            tmp += ftmp[name][:]
        dset[...] = tmp
    return

shape = (len(etaBins)-1) * (len(ptBins)-1) * (len(chargeBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)
haddFiles(ewkFiles, histonames, shape)

fdata = h5py.File('data.hdf5', mode='r+')
fewk = h5py.File('ewk.hdf5', mode='r+')

hdata = np.array(fdata['data_obs'][:].reshape((len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1),order='F'),order='C')
hdatasumw2 = np.array(fdata['data_obs_sumw2'][:].reshape((len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1),order='F'),order='C')
hewk = np.array(fewk['ewk'][:].reshape((len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1),order='F'),order='C')
hewksumw2 = np.array(fewk['ewk_sumw2'][:].reshape((len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1),order='F'),order='C')

os.chdir('..')

fbkg = h5py.File('bkgWplus.hdf5', mode='w')

dset = fbkg.create_dataset(name='data_obs', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
dset[...] = hdata[:,:,-1,:,:].flatten() #select positive charge

dset = fbkg.create_dataset(name='ewk', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
dset[...] = hewk[:,:,-1,:,:].flatten() #select positive charge

dset = fbkg.create_dataset(name='ewk_sumw2', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
dset[...] = hewksumw2[:,:,-1,:,:].flatten() #select positive charge

hfakesLowMt = (hdata[:,:,-1,:,:]-hewk[:,:,-1,:,:])
hfakesLowMt[:,:,1,:] = np.zeros([len(etaBins)-1,len(ptBins)-1,len(isoBins)-1], dtype='float64')
dset = fbkg.create_dataset(name='fakesLowMt', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
dset[...] = hfakesLowMt.flatten()

hfakesLowMtsumw2 = (hdatasumw2[:,:,-1,:,:]+hewksumw2[:,:,-1,:,:])
hfakesLowMtsumw2[:,:,1,:] = np.zeros([len(etaBins)-1,len(ptBins)-1,len(isoBins)-1], dtype='float64')
dset = fbkg.create_dataset(name='fakesLowMt_sumw2', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
dset[...] = hfakesLowMtsumw2.flatten()

norm = np.sum(hfakesLowMt[:,:,0,1])/np.sum(hdata[:,:,-1,0,1])

hfakesHighMt = norm*(hdata[:,:,-1,:,:]-hewk[:,:,-1,:,:])
hfakesHighMt[:,:,0,:] = np.zeros([len(etaBins)-1,len(ptBins)-1,len(isoBins)-1], dtype='float64')
# keep aiso/iso ratio constant
hfakesHighMt[:,:,1,0] = norm*(hdata[:,:,-1,1,1]-hewk[:,:,-1,1,1]) * (hdata[:,:,-1,0,0]-hewk[:,:,-1,0,0])/(hdata[:,:,-1,0,1]-hewk[:,:,-1,0,1])
dset = fbkg.create_dataset(name='fakesHighMt', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
dset[...] = hfakesHighMt.flatten()

hfakesHighMtsumw2 = norm*(hdatasumw2[:,:,-1,:,:]+hewksumw2[:,:,-1,:,:])
hfakesHighMtsumw2[:,:,0,:] = np.zeros([len(etaBins)-1,len(ptBins)-1,len(isoBins)-1], dtype='float64')
dset = fbkg.create_dataset(name='fakesHighMt_sumw2', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
dset[...] = hfakesHighMtsumw2.flatten()

# plot pt, eta and mt in isolated and antiisolated region

for i in range(2):
    fig, ax = plt.subplots()
    ax.set_title("eta_iso{}_highMt_prefit".format(i), fontsize=18)
    etadata = np.sum(hdata,axis=1)[:,-1,-1,i]
    etaewk = np.sum(hewk,axis=1)[:,-1,-1,i]
    etafake = np.sum(hfakesHighMt,axis=1)[:,-1,i]

    hep.histplot([etadata],bins = etaBins, histtype = 'errorbar', color = "k", stack = False)
    hep.histplot([etafake,etaewk],bins = etaBins, histtype = 'fill',linestyle = 'solid', color =["g","r"], label=["fakes", "ewk"], stack = True)

    plt.savefig('PLOTS/eta_iso{}_highMt_prefit.png'.format(i))
    plt.cla()

    ptdata = np.sum(hdata,axis=0)[:,-1,-1,i]
    ptewk = np.sum(hewk,axis=0)[:,-1,-1,i]
    ptfake = np.sum(hfakesHighMt,axis=0)[:,-1,i]

    hep.histplot([ptdata],bins = ptBins, histtype = 'errorbar', color = "k", stack = False)
    hep.histplot([ptewk,ptfake],bins = ptBins, histtype = 'fill',linestyle = 'solid', color =["g","r"], stack = True)

    plt.savefig('PLOTS/pt_iso{}.png'.format(i))
    plt.cla()

    # mtdata = np.sum(hdata,axis=(0,1))[:,i]
    # mtewk = np.sum(hewk,axis=(0,1))[:,i]

    # hep.histplot([mtdata],bins = mTBins, histtype = 'errorbar', color = "k")
    # hep.histplot([mtewk],bins = mTBins, histtype = 'fill',linestyle = 'solid', color ="r")

    # plt.savefig('PLOTS/mt_iso{}.png'.format(i))
    # plt.clf()

for i in range(2):
    fig, ax = plt.subplots()
    ax.set_title("eta_iso{}_lowMt_prefit".format(i), fontsize=18)
    etadata = np.sum(hdata,axis=1)[:,-1,0,i]
    etaewk = np.sum(hewk,axis=1)[:,-1,0,i]
    etafake = np.sum(hfakesLowMt,axis=1)[:,0,i]

    hep.histplot([etadata],bins = etaBins, histtype = 'errorbar', color = "k", stack = False)
    hep.histplot([etafake,etaewk],bins = etaBins, histtype = 'fill',linestyle = 'solid', color =["g","r"], stack = True)

    plt.savefig('PLOTS/eta_iso{}_lowmt_prefit.png'.format(i))
    plt.cla()

    ptdata = np.sum(hdata,axis=0)[:,-1,0,i]
    ptewk = np.sum(hewk,axis=0)[:,-1,0,i]
    ptfake = np.sum(hfakesLowMt,axis=0)[:,0,i]

    hep.histplot([ptdata],bins = ptBins, histtype = 'errorbar', color = "k", stack = False)
    hep.histplot([ptewk,ptfake],bins = ptBins, histtype = 'fill',linestyle = 'solid', color =["g","r"], stack = True)

    plt.savefig('PLOTS/pt_iso{}_lowmt_prefit.png'.format(i))
    plt.cla()

    # mtdata = np.sum(hdata,axis=(0,1))[:,i]
    # mtewk = np.sum(hewk,axis=(0,1))[:,i]

    # hep.histplot([mtdata],bins = mTBins, histtype = 'errorbar', color = "k")
    # hep.histplot([mtewk],bins = mTBins, histtype = 'fill',linestyle = 'solid', color ="r")

    # plt.savefig('PLOTS/mt_iso{}.png'.format(i))
    # plt.clf()

