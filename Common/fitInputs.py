import os
import sys
sys.path.append('data/')
import h5py
from math import pi, sqrt
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from binning import ptBins, etaBins, isoBins, chargeBins, metBins, mTBins, yBins, qtBins
# from binning import mTBinsFull as mTBins
plt.style.use([hep.style.ROOT])
#hep.cms.label(loc=0, year=2016, lumi=35.9, data=True)
#hep.cms.text('Simulation')

WMuFiles = ["WPlusJetsToMuNu.hdf5"]
WTauFiles = ["WPlusJetsToTauNu.hdf5","WMinusJetsToTauNu.hdf5"]
DYFiles = ["DYJetsToMuMu_M50.hdf5","DYJetsToTauTau_M50.hdf5"]
TopFiles = ["ST_t-channel_muDecays.hdf5", "ST_t-channel_tauDecays.hdf5","ST_s-channel_4f_leptonDecays.hdf5","ST_t-channel_top_5f_InclusiveDecays.hdf5","TTToSemiLeptonic.hdf5", "TTTo2L2Nu.hdf5"]
DibosonFiles = ["WW.hdf5","WZ.hdf5"]
dataFiles = ["data.hdf5"]

def haddFiles(fileList, fname, histonames, shapes, folder, era):
    f = h5py.File('./'+fname+'_{}.hdf5'.format(era), mode='w')
    for i,name in enumerate(histonames):
        dset = f.create_dataset(name=name, shape=shapes[i], dtype='float64')
        tmp = np.zeros(shapes[i],dtype='float64')
        print(name, shape)
        for file in fileList:
            ftmp = h5py.File(folder+file, mode='r+')
            tmp += ftmp[name][:]
        dset[...] = tmp
    return

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1

eras = ["preVFP","postVFP"]
for era in eras:
    folder = "../config/rightLumi_{}/".format(era)
    shape = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1)
    haddFiles(dataFiles,"data",["data_obs","data_obs_sumw2"], [shape,shape],folder,era)
    # hadd files to bkg categories 
    histonames = ['ewk', 'ewk_sumw2']
    haddFiles(WTauFiles,"Wtau",histonames, [shape,shape],folder,era)
    haddFiles(DYFiles,"DY",histonames, [shape,shape],folder,era)
    haddFiles(TopFiles,"Top",histonames, [shape,shape],folder,era)
    haddFiles(DibosonFiles,"Diboson",histonames, [shape,shape],folder,era)

    fdata = h5py.File('./data_{}.hdf5'.format(era), mode='r+')
    fWtau = h5py.File('./Wtau_{}.hdf5'.format(era), mode='r+')
    fDY = h5py.File('./DY_{}.hdf5'.format(era), mode='r+')
    fTop = h5py.File('./Top_{}.hdf5'.format(era), mode='r+')
    fDiboson = h5py.File('./Diboson_{}.hdf5'.format(era), mode='r+')

    hdata = np.array(fdata['data_obs'][:])
    hdata_sumw2 = np.array(fdata['data_obs_sumw2'][:])
    # hWmu = np.array(fWmu['ewk'][:])
    hWtau = np.array(fWtau['ewk'][:])
    hWtau_sumw2 = np.array(fWtau['ewk_sumw2'][:])
    hDY = np.array(fDY['ewk'][:])
    hDY_sumw2 = np.array(fDY['ewk_sumw2'][:])
    hTop = np.array(fTop['ewk'][:])
    hTop_sumw2 = np.array(fTop['ewk_sumw2'][:])
    hDiboson = np.array(fDiboson['ewk'][:])
    hDiboson_sumw2 = np.array(fDiboson['ewk_sumw2'][:])

    # write down shapes as fit input
    fshapes = h5py.File('shapesWplus_{}.hdf5'.format(era), mode='w')

    # uncomment to write out real data
    # dset = fshapes.create_dataset(name='data_obs', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    # dset[...] = hdata[:,:,-1,:,:] #select positive charge

    dset = fshapes.create_dataset(name='Wtau', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hWtau[:,:,-1,:,:] #select positive charge

    dset = fshapes.create_dataset(name='Wtau_sumw2', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hWtau_sumw2[:,:,-1,:,:] #select positive charge

    dset = fshapes.create_dataset(name='DY', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hDY[:,:,-1,:,:] #select positive charge

    dset = fshapes.create_dataset(name='DY_sumw2', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hDY_sumw2[:,:,-1,:,:] #select positive charge

    dset = fshapes.create_dataset(name='Top', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hTop[:,:,-1,:,:] #select positive charge

    dset = fshapes.create_dataset(name='Top_sumw2', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hTop_sumw2[:,:,-1,:,:] #select positive charge

    dset = fshapes.create_dataset(name='Diboson', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hDiboson[:,:,-1,:,:] #select positive charge

    dset = fshapes.create_dataset(name='Diboson_sumw2', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hDiboson_sumw2[:,:,-1,:,:] #select positive charge

    # now hadd and write down W differential signal

    shape = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 9)
    shape_mass = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 9*2)
    shape_pdf = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 103)

    # histonames = ['signalTemplates', 'signalTemplates_sumw2', 'signalTemplates_mass', 'signalTemplates_LHEPdfWeight']
    histonames = ['signalTemplates', 'signalTemplates_sumw2']
    haddFiles(WMuFiles,"WmuSignal",histonames, [shape,shape,shape_mass,shape_pdf],folder,era)

    fsignal = h5py.File('./WmuSignal_{}.hdf5'.format(era), mode='r+')

    # create shape for modified pseudo-data
    qtBins = np.array(qtBins)
    qtBinsC = 0.5*(qtBins[1:]+qtBins[:-1])
    slope = -1.5*10**-3
    inter = 1.1
    qtweight_vec = (slope*qtBinsC+inter)[np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis,:,np.newaxis]
    qtBinsS = qtBins[1:]-qtBins[:-1]

    fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
    ax1.set_title("test reweight", fontsize=18)
    hep.histplot(np.sum(np.array(fsignal['signalTemplates'][:]),axis=(0,1,5,-1))[-1,-1,0],qtBins, ax=ax1)
    # hep.histplot(np.sum(np.array(fsignal['signalTemplates'][:])*qtweight_vec,axis=(0,1,5,-1))[-1,-1,0]/qtBinsS,qtBins, ax=ax1)
    ax2.set_ylim([0.8, 1.2])
    hep.histplot([np.sum(np.array(fsignal['signalTemplates'][:]),axis=(0,1,5,-1))[-1,-1,0]/np.sum(np.array(fsignal['signalTemplates'][:])*qtweight_vec,axis=(0,1,5,-1))[-1,-1,0]],bins = qtBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
    plt.tight_layout()
    plt.savefig("testprefit/test_reweight_{}.png".format(era))
    plt.clf()

    # this is pseudodata made with the sum of the templates
    dset = fshapes.create_dataset(name='data_obs', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = np.sum((np.array(fsignal['signalTemplates'][:])*qtweight_vec)[...,:threshold_y,:threshold_qt,:],axis=(5,6,-1))[:,:,-1,:,:] #select positive charge

    # signal: differential in y,pt and helicity
    hsignal = np.array(fsignal['signalTemplates'][:])[:,:,-1,...,:threshold_y,:threshold_qt,:]
    hsignal_sumw2 = np.array(fsignal['signalTemplates_sumw2'][:])[:,:,-1,...,:threshold_y,:threshold_qt,:]

    print(era, fsignal['signalTemplates'][:].shape,hsignal.shape,threshold_y)

    # # signal: mass
    # hsignal_mass = np.array(fsignal['signalTemplates_mass'][:])[:,:,-1,...,:threshold_y,:threshold_qt,:]
    # hsignal_mass=hsignal_mass.reshape(hsignal_mass.shape[:-1] + (9,2))
    # # signal: PDFs
    # # hsignal_PDF = np.array(fsignal['signalTemplates_LHEPdfWeight'][:])[:,:,-1,...,:threshold_y,:threshold_qt,:]
    # # hsignal_PDF=hsignal_PDF.reshape(hsignal_PDF.shape[:-1] + (9,103))

    # sum of all helicities = total diff xsec
    hdifftot = np.sum(np.array(fsignal['signalTemplates'][:]),axis=-1)
    hdifftot_sumw2 = np.sum(np.array(fsignal['signalTemplates_sumw2'][:]),axis=-1)

    # hdifftot_mass = np.array(fsignal['signalTemplates_mass'][:])
    # hdifftot_mass=hdifftot_mass.reshape(hdifftot_mass.shape[:-1] + (9,2))
    # hdifftot_mass =np.sum(hdifftot_mass,axis=-2)

    # hdifftot_PDF = np.array(fsignal['signalTemplates_LHEPdfWeight'][:])

    hWmu = np.sum(np.array(hdifftot[...,:threshold_y,:threshold_qt]),axis=(5,6))
    hWmu_sumw2 = np.sum(np.array(hdifftot_sumw2),axis=(5,6)) # integrated over helicity, y and qt

    # events falling out of fit range
    hlowacc = np.sum(hdifftot[...,threshold_y:,threshold_qt:],axis=(-1,-2))+np.sum(hdifftot[...,threshold_y:,:threshold_qt],axis=(-1,-2))+np.sum(hdifftot[...,:threshold_y,threshold_qt:],axis=(-1,-2))
    hlowacc_sumw2 = np.sum(hdifftot_sumw2[...,threshold_y:,threshold_qt:],axis=(-1,-2))+np.sum(hdifftot_sumw2[...,threshold_y:,:threshold_qt],axis=(-1,-2))+np.sum(hdifftot_sumw2[...,:threshold_y,threshold_qt:],axis=(-1,-2))

    # # lowacc: mass
    # hlowacc_mass = np.sum(hdifftot_mass[...,threshold_y:,threshold_qt:,:],axis=(-2,-3))+np.sum(hdifftot_mass[...,threshold_y:,:threshold_qt,:],axis=(-2,-3))+np.sum(hdifftot_mass[...,:threshold_y,threshold_qt:,:],axis=(-2,-3))
    # hlowacc_mass=hlowacc_mass.reshape(hlowacc_mass.shape[:-1] + (2,))
    # # lowacc: PDF
    # hlowacc_PDF = np.sum(hdifftot_PDF[...,threshold_y:,threshold_qt:,:],axis=(-2,-3))+np.sum(hdifftot_PDF[...,threshold_y:,:threshold_qt,:],axis=(-2,-3))+np.sum(hdifftot_PDF[...,:threshold_y,threshold_qt:,:],axis=(-2,-3))
    # hlowacc_PDF=hlowacc_PDF.reshape(hlowacc_PDF.shape[:-1] + (103,))
    # print(hlowacc_mass.shape, hlowacc_PDF.shape)

    print(hWmu.shape, hlowacc.shape, hWtau.shape)
    hewk = hWmu+hlowacc+hWtau+hDY+hTop+hDiboson

    fig, ax1 = plt.subplots()
    ax1.set_title("total xsec closure", fontsize=18)
    hep.hist2dplot(np.sum(hdifftot,axis=(5,6))[:,:,-1,1,0],etaBins,ptBins)
    plt.tight_layout()
    plt.savefig("testprefit/total_xsec_clos_{}".format(era))
    plt.clf()
    fig, ax1 = plt.subplots()
    ax1.set_title("total xsec closure", fontsize=18)
    hep.hist2dplot(hsignal[:,:,1,0,0,2,-1],etaBins,ptBins)
    plt.tight_layout()
    plt.savefig("testprefit/templ_".format(era))
    plt.clf()
    fig, ax1 = plt.subplots()
    ax1.set_title("total xsec closure", fontsize=18)
    hep.histplot(np.sum(hdifftot,axis=(1,5,6))[:,-1,1,0],etaBins)
    plt.tight_layout()
    plt.savefig("testprefit/total_xsec_clos_eta_{}".format(era))
    plt.clf()
    fig, ax1 = plt.subplots()
    ax1.set_title("total xsec closure", fontsize=18)
    hep.histplot(np.sum(hdifftot,axis=(0,5,6))[:,-1,1,0],ptBins)
    plt.tight_layout()
    plt.savefig("testprefit/total_xsec_clos_pt_{}".format(era))
    plt.clf()

    dset = fshapes.create_dataset(name='template', shape=hsignal.shape, dtype='float64')
    dset[...] = hsignal #select positive charge

    dset = fshapes.create_dataset(name='template_sumw2', shape=hsignal.shape, dtype='float64')
    dset[...] = hsignal_sumw2 #select positive charge

    dset = fshapes.create_dataset(name='lowacc', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hlowacc[:,:,-1,:,:] #select positive charge

    dset = fshapes.create_dataset(name='lowacc_sumw2', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hlowacc_sumw2[:,:,-1,:,:] #select positive charge

    # # now write shapes with systematics

    # dset = fshapes.create_dataset(name='template_mass', shape=hsignal_mass.shape, dtype='float64')
    # dset[...] = hsignal_mass #select positive charge

    # # dset = fshapes.create_dataset(name='template_LHEPdfWeight', shape=hsignal_PDF.shape, dtype='float64')
    # # dset[...] = hsignal_PDF #select positive charge

    # dset = fshapes.create_dataset(name='lowacc_mass', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1,2), dtype='float64')
    # dset[...] = hlowacc_mass[:,:,-1,:,:] #select positive charge

    # dset = fshapes.create_dataset(name='lowacc_LHEPdfWeight', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1,103), dtype='float64')
    # dset[...] = hlowacc_PDF[:,:,-1,:,:] #select positive charge

    # helicity xsecs without acceptance cuts for unfolding
    file_gen = '../config/halfptBins_{}/WPlusJetsToMuNu_helweights.hdf5'.format(era)
    f_gen = h5py.File(file_gen, mode='r+')

    # merge pre and post VFP xsecs
    htot = f_gen['totxsecs'][:]
    h = f_gen['xsecs'][:]

    yBins = f_gen['edges_totxsecs_0'][:]
    qtBins = f_gen['edges_totxsecs_1'][:]
    print(yBins)
    factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
    factors = factors[np.newaxis,np.newaxis,...]
    print(factors.shape)
    h = (h/htot[...,np.newaxis]+factors[...,1])*factors[...,0]

    factors_hel = np.array([2.,2*sqrt(2),4.,4.*sqrt(2),2.,2.,2.*sqrt(2),4.*sqrt(2),1.])
    factors_hel = factors_hel[np.newaxis,np.newaxis,...]
    h = 3./(16.*pi)*h*htot[...,np.newaxis]/factors_hel[...,:threshold_y,:threshold_qt,:]
    h[...,-1] = 3./(16.*pi)*htot
    dset = fshapes.create_dataset(name='helicity', shape=h.shape, dtype='float64')
    dset[...] = h #select positive charge

    print(np.sum(htot[...,:threshold_y,:threshold_qt],axis=0)[0])
    print(np.sum(htot,axis=0)[0])

    fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
    ax1.set_title("test qt closure", fontsize=18)
    hep.histplot(np.sum(np.array(fsignal['signalTemplates'][:])[:,:,:,...,:threshold_y,:,:],axis=(0,1,5,7))[-1,-1,0,:],qtBins, ax=ax1)
    hep.histplot(np.sum(htot,axis=0),qtBins, ax=ax1)
    ax2.set_ylim([0., 1.])
    hep.histplot(np.sum(np.array(fsignal['signalTemplates'][:]),axis=(0,1,5,7))[-1,-1,0]/(np.sum(htot,axis=0)),bins = qtBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
    plt.tight_layout()
    plt.savefig("testprefit/test_qt_clos_{}.png".format(era))
    plt.clf()

    fig, ax1 = plt.subplots()
    hep.hist2dplot(np.sum(np.array(fsignal['signalTemplates'][:]),axis=(0,1,7))[-1,-1,0]/htot,yBins,qtBins)
    plt.tight_layout()
    plt.savefig("testprefit/test_yqt_clos_{}.png".format(era))
    plt.clf()

    # get prefit estimate of shape and normalisation of QCD bkg

    hfakesLowMt = hdata[:,:,-1,:,:]-hewk[:,:,-1,:,:]
    hfakes_unc = hdata_sumw2[:,:,-1,:,:]+hWmu_sumw2[:,:,-1,:,:]
    hfakesLowMt[:,:,1,:] = np.zeros([len(etaBins)-1,len(ptBins)-1,len(isoBins)-1], dtype='float64')
    dset = fshapes.create_dataset(name='fakesLowMt', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    dset[...] = hfakesLowMt.flatten()

    dset = fshapes.create_dataset(name='fakesLowMt_sumw2', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    dset[...] = np.zeros_like(hfakesLowMt.flatten())

    # threshold = np.digitize(30.,mTBins)
    # fR is computed in low-mt region
    # fR = np.sum(hdata[:,:,-1,:threshold-1,0]-hewk[:,:,-1,:threshold-1,0], axis=2)/np.sum(hdata[:,:,-1,:threshold-1,1]-hewk[:,:,-1,:threshold-1,1],axis=2)
    fR =(hdata[:,:,-1,0,0]-hewk[:,:,-1,0,0])/(hdata[:,:,-1,0,1]-hewk[:,:,-1,0,1])

    hfakesHighMt = np.where(hdata[:,:,-1,:,:]-hewk[:,:,-1,:,:]>0, hdata[:,:,-1,:,:]-hewk[:,:,-1,:,:], 1)
    hfakesHighMt[:,:,0,:] = np.zeros([len(etaBins)-1,len(ptBins)-1,len(isoBins)-1], dtype='float64')
    # keep aiso/iso ratio constant
    hfakesHighMt[:,:,1,0] = (hdata[:,:,-1,1,1]-hewk[:,:,-1,1,1]) * fR
    dset = fshapes.create_dataset(name='fakesHighMt', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    dset[...] = hfakesHighMt.flatten()

    dset = fshapes.create_dataset(name='fakesHighMt_sumw2', shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    dset[...] = np.zeros_like(hfakesHighMt.flatten())

    # build mask
    print("building shapes")
    for i in range(hfakesLowMt.shape[0]*hfakesLowMt.shape[1]):
        mask = np.zeros(hfakesLowMt.shape[0]*hfakesLowMt.shape[1])
        mask[i,...]=1
        mask = mask.reshape((hfakesLowMt.shape[0],hfakesLowMt.shape[1]))[...,np.newaxis,np.newaxis]
        # nuisance for changing the normalisations independently

        hfakesLowMtVarUp = np.where(mask==0, hfakesLowMt, hfakesLowMt+0.5*hfakesLowMt)
        dset = fshapes.create_dataset(name='fakesLowMt_fakeNormLowMtBin{}Up'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
        dset[...] = hfakesLowMtVarUp.flatten()
        hfakesLowMtVarDown = np.where(mask==0, hfakesLowMt, hfakesLowMt-0.5*hfakesLowMt)
        dset = fshapes.create_dataset(name='fakesLowMt_fakeNormLowMtBin{}Down'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
        dset[...] = hfakesLowMtVarDown.flatten()

        hfakesHighMtVarUp = np.where(mask==0, hfakesHighMt, hfakesHighMt+0.5*hfakesHighMt)
        dset = fshapes.create_dataset(name='fakesHighMt_fakeNormHighMtBin{}Up'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
        dset[...] = hfakesHighMtVarUp.flatten()
        hfakesHighMtVarDown = np.where(mask==0, hfakesHighMt, hfakesHighMt-0.5*hfakesHighMt)
        dset = fshapes.create_dataset(name='fakesHighMt_fakeNormHighMtBin{}Down'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
        dset[...] = hfakesHighMtVarDown.flatten()

        # print('checking if any zero or negative yields for bin {}'.format(i))
        # print(np.any((hfakesLowMtVarUp+hfakesHighMtVarUp)<=0.))
        # print(np.any((hfakesLowMtVarDown+hfakesHighMtVarDown)<=0.))

        # common nuisance for changing fake shape

        norm = np.sum(hfakesLowMt[:,:,0,:],axis=2)
        ratio = hfakesLowMt[:,:,0,0]/norm #ratio iso/iso+aiso
        rate_var = 2.
        var_idx = np.nonzero(mask)
        # set to nominal
        hfakesLowMtVarUp = np.empty_like(hfakesLowMt)
        hfakesLowMtVarDown = np.empty_like(hfakesLowMt)
        np.copyto(hfakesLowMtVarUp, hfakesLowMt) # (dst, src)
        np.copyto(hfakesLowMtVarDown, hfakesLowMt) # (dst, src)
        # apply variation to isolated part
        hfakesLowMtVarUp[var_idx[0],var_idx[1],0, 0] = (rate_var*ratio*norm)[var_idx[0],var_idx[1]]
        hfakesLowMtVarDown[var_idx[0],var_idx[1],0, 0] = ((1./rate_var)*ratio*norm)[var_idx[0],var_idx[1]]
        # apply variation to anti-isolated part
        hfakesLowMtVarUp[var_idx[0],var_idx[1],0, 1] = ((1-rate_var*ratio)*norm)[var_idx[0],var_idx[1]]
        hfakesLowMtVarDown[var_idx[0],var_idx[1],0, 1] = ((1-(1./rate_var)*ratio)*norm)[var_idx[0],var_idx[1]]

        dset = fshapes.create_dataset(name='fakesLowMt_fakeShapeBin{}Up'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
        dset[...] = (hfakesLowMtVarUp/hfakes_unc).flatten()
        dset = fshapes.create_dataset(name='fakesLowMt_fakeShapeBin{}Down'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
        dset[...] = (hfakesLowMtVarDown/hfakes_unc).flatten()

        norm = np.sum(hfakesHighMt[:,:,1,:],axis=2)
        ratio = hfakesHighMt[:,:,1,0]/norm #ratio iso/iso+aiso
        rate_var = 1.2
        var_idx = np.nonzero(mask)
        # set to nominal
        hfakesHighMtVarUp = np.empty_like(hfakesHighMt)
        hfakesHighMtVarDown = np.empty_like(hfakesHighMt)
        np.copyto(hfakesHighMtVarUp, hfakesHighMt) # (dst, src)
        np.copyto(hfakesHighMtVarDown, hfakesHighMt) # (dst, src)
        # apply variation to isolated part
        hfakesHighMtVarUp[var_idx[0],var_idx[1],1, 0] = (rate_var*ratio*norm)[var_idx[0],var_idx[1]]
        hfakesHighMtVarDown[var_idx[0],var_idx[1],1, 0] = ((1./rate_var)*ratio*norm)[var_idx[0],var_idx[1]]
        # apply variation to anti-isolated part
        hfakesHighMtVarUp[var_idx[0],var_idx[1],1, 1] = ((1-rate_var*ratio)*norm)[var_idx[0],var_idx[1]]
        hfakesHighMtVarDown[var_idx[0],var_idx[1],1, 1] = ((1-(1./rate_var)*ratio)*norm)[var_idx[0],var_idx[1]]

        dset = fshapes.create_dataset(name='fakesHighMt_fakeShapeBin{}Up'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
        dset[...] = (hfakesHighMtVarUp/hfakes_unc).flatten()
        dset = fshapes.create_dataset(name='fakesHighMt_fakeShapeBin{}Down'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
        dset[...] = (hfakesHighMtVarDown/hfakes_unc).flatten()

        # print('checking if any zero or negative yields for bin {}'.format(i))
        # print(np.any((hfakesLowMtVarUp+hfakesHighMtVarUp)<=0.))
        # print(np.any((hfakesLowMtVarDown+hfakesHighMtVarDown)<=0.))

    # plot pt, eta and mt in isolated and antiisolated region
    for i in range(2):
        fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
        ax1.set_title("eta_iso{}_highMt".format(i), fontsize=18)
        ax1.set_ylabel('number of events')
        ax2.set_ylabel('data/prediction')
        ax2.set_xlabel('muon $\eta$')
        etadata = np.sum(hdata,axis=1)[:,-1,-1,i]
        etaewk = np.sum(hewk,axis=1)[:,-1,-1,i]
        etaWmu = np.sum(hWmu,axis=1)[:,-1,-1,i]
        etaWtau = np.sum(hWtau,axis=1)[:,-1,-1,i]
        etaDY = np.sum(hDY,axis=1)[:,-1,-1,i]
        etaTop = np.sum(hTop,axis=1)[:,-1,-1,i]
        etaDiboson = np.sum(hDiboson,axis=1)[:,-1,-1,i]
        etafake = np.sum(hfakesHighMt,axis=1)[:,-1,i]
        etalowacc = np.sum(hlowacc,axis=1)[:,-1,-1,i]
        hep.histplot([etadata],bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1, label = ["data"])
        hep.histplot([etaDiboson,etaTop,etaDY,etaWtau,etafake,etalowacc,etaWmu],bins = etaBins, histtype = 'fill',linestyle = 'solid', color =["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
        ax2.set_ylim([0.9, 1.1])
        hep.histplot([etadata/(etafake+etaewk)], bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
        ax1.legend(loc='upper right', frameon=True)
        plt.tight_layout()
        plt.savefig('testprefit/eta_iso{}_highMt_{}.png'.format(i,era))
        plt.cla()

        fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
        ax1.set_title("pt_iso{}_highMt".format(i), fontsize=18)
        ax1.set_ylabel('number of events')
        ax2.set_ylabel('data/prediction')
        ax2.set_xlabel('muon $p_T$')
        ptdata = np.sum(hdata,axis=0)[:,-1,-1,i]
        ptewk = np.sum(hewk,axis=0)[:,-1,-1,i]
        ptWmu = np.sum(hWmu,axis=0)[:,-1,-1,i]
        ptWtau = np.sum(hWtau,axis=0)[:,-1,-1,i]
        ptDY = np.sum(hDY,axis=0)[:,-1,-1,i]
        ptTop = np.sum(hTop,axis=0)[:,-1,-1,i]
        ptDiboson = np.sum(hDiboson,axis=0)[:,-1,-1,i]
        ptfake = np.sum(hfakesHighMt,axis=0)[:,-1,i]
        ptlowacc = np.sum(hlowacc,axis=0)[:,-1,-1,i]
        hep.histplot([ptdata],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
        hep.histplot([ptDiboson,ptTop,ptDY,ptWtau,ptfake,ptlowacc,ptWmu],bins = ptBins, histtype = 'fill',linestyle = 'solid', color =["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
        ax2.set_ylim([0.9, 1.1])
        hep.histplot([ptdata/(ptfake+ptewk)],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
        ax1.legend(loc='upper right', frameon=True)
        plt.tight_layout()
        plt.savefig('testprefit/pt_iso{}_highMt_{}.png'.format(i,era))
        plt.cla()

        fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
        ax1.set_title("eta_iso{}_lowMt".format(i), fontsize=18)
        ax1.set_ylabel('number of events')
        ax2.set_ylabel('data/prediction')
        ax2.set_xlabel('muon $\eta$')
        etadata = np.sum(hdata,axis=1)[:,-1,0,i]
        etaewk = np.sum(hewk,axis=1)[:,-1,0,i]
        etaWmu = np.sum(hWmu,axis=1)[:,-1,0,i]
        etaWtau = np.sum(hWtau,axis=1)[:,-1,0,i]
        etaDY = np.sum(hDY,axis=1)[:,-1,0,i]
        etaTop = np.sum(hTop,axis=1)[:,-1,0,i]
        etaDiboson = np.sum(hDiboson,axis=1)[:,-1,0,i]
        etafake = np.sum(hfakesLowMt,axis=1)[:,0,i]
        etalowacc = np.sum(hlowacc,axis=1)[:,-1,0,i]
        hep.histplot([etadata],bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
        hep.histplot([etaDiboson,etaTop,etaDY,etaWtau,etafake,etalowacc,etaWmu],bins = etaBins, histtype = 'fill',linestyle = 'solid', color =["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
        ax2.set_ylim([0.9, 1.1])
        hep.histplot([etadata/(etafake+etaewk)],bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
        ax1.legend(loc='upper right', frameon=True)
        plt.tight_layout()
        plt.savefig('testprefit/eta_iso{}_lowMt_{}.png'.format(i,era))
        plt.cla()

        fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
        ax1.set_title("pt_iso{}_lowMt".format(i), fontsize=18)
        ax1.set_ylabel('number of events')
        ax2.set_ylabel('data/prediction')
        ax2.set_xlabel('muon $p_T$')
        ptdata = np.sum(hdata,axis=0)[:,-1,0,i]
        ptewk = np.sum(hewk,axis=0)[:,-1,0,i]
        ptWmu = np.sum(hWmu,axis=0)[:,-1,0,i]
        ptWtau = np.sum(hWtau,axis=0)[:,-1,0,i]
        ptDY = np.sum(hDY,axis=0)[:,-1,0,i]
        ptTop = np.sum(hTop,axis=0)[:,-1,0,i]
        ptDiboson = np.sum(hDiboson,axis=0)[:,-1,0,i]
        ptfake = np.sum(hfakesLowMt,axis=0)[:,0,i]
        ptlowacc = np.sum(hlowacc,axis=0)[:,-1,0,i]
        hep.histplot([ptdata],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
        hep.histplot([ptDiboson,ptTop,ptDY,ptWtau,ptfake,ptlowacc,ptWmu],bins = ptBins, histtype = 'fill',linestyle = 'solid', color =["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
        ax2.set_ylim([0.9, 1.1])
        hep.histplot([ptdata/(ptfake+ptewk)],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
        ax1.legend(loc='upper right', frameon=True)
        plt.tight_layout()
        plt.savefig('testprefit/pt_iso{}_lowMt_{}.png'.format(i,era))
        plt.cla()

    # fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
    # ax1.set_title("mt", fontsize=18)
    # ax1.set_ylabel('number of events')
    # ax2.set_ylabel('data/prediction')
    # ax2.set_xlabel('$m_T$')
    # mtdata = np.sum(hdata[:,:,-1,:,0],axis=(0,1))
    # mtewk = np.sum(hewk[:,:,-1,:,0],axis=(0,1))
    # mtW = np.sum(hW[:,:,-1,:,0],axis=(0,1))
    # mtDY = np.sum(hDY[:,:,-1,:,0],axis=(0,1))
    # mtTop = np.sum(hTop[:,:,-1,:,0],axis=(0,1))
    # mtDiboson = np.sum(hDiboson[:,:,-1,:,0],axis=(0,1))
    # mtfake = np.einsum('kmi,km->i',hdata[:,:,-1,:,1]-hewk[:,:,-1,:,1],fR)

    # hep.histplot([mtdata],bins = mTBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
    # hep.histplot([mtDiboson,mtTop,mtDY,mtfake,mtW],bins = mTBins, histtype = 'fill',linestyle = 'solid', color =["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
    # ax2.set_ylim([0.9, 1.1])
    # hep.histplot([mtdata/(mtfake+mtewk)],bins = mTBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
    # ax1.legend(loc='upper right', frameon=True)
    # plt.tight_layout()
    # plt.savefig('testprefit/mt_{}.png'.format(era))
    # plt.cla()