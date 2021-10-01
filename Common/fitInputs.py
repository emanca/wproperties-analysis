import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import sys
sys.path.append('data/')
import h5py
from math import pi, sqrt
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import ROOT
from root_numpy import hist2array
from binning import ptBins, etaBins, isoBins, chargeBins, metBins, mTBins, yBins, qtBins, qtBins_syst
# from binning import mTBinsFull as mTBins
plt.style.use([hep.style.ROOT])
#hep.cms.label(loc=0, year=2016, lumi=35.9, data=True)
#hep.cms.text('Simulation')

WMuFiles = ["WPlusJetsToMuNu.hdf5"]
WMuFiles_rew = ["WPlusJetsToMuNu.hdf5"]
WTauFiles = ["WPlusJetsToTauNu.hdf5"]
DYFiles = ["DYJetsToMuMu_M50.hdf5","DYJetsToTauTau_M50.hdf5"]
TopFiles = ["ST_t-channel_muDecays.hdf5", "ST_t-channel_tauDecays.hdf5","ST_s-channel_4f_leptonDecays.hdf5","ST_t-channel_top_5f_InclusiveDecays.hdf5","TTToSemiLeptonic.hdf5", "TTTo2L2Nu.hdf5"]
DibosonFiles = ["WW.hdf5","WZ.hdf5"]
dataFiles = ["data.hdf5"]

def haddFiles(fileList, fname, histonames, shapes, folder, era):
    dict={}
    for i,name in enumerate(histonames):
        tmp = np.zeros(shapes[i],dtype='float64')
        print(name, shapes[i])
        for file in fileList:
            ftmp = h5py.File(folder+file, mode='r+')
            tmp += ftmp[name][:]
        dict[name]=tmp
    return dict

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1

eras = ["preVFP","postVFP"]
for era in eras:
    # folder = "../config/alternatesample_{}/".format(era)
    folder = "../config/valeriosbins_{}/".format(era)
    shape = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1)
    datadict = haddFiles(dataFiles,"data",["data_obs","data_obs_sumw2"], [shape,shape],folder,era)
    hdata = datadict['data_obs']
    hdata_sumw2 = datadict['data_obs_sumw2']

    # hadd files to bkg categories 
    histonames = ['ewk', 'ewk_sumw2']
    taudict=haddFiles(WTauFiles,"Wtau",histonames, [shape,shape],folder,era)
    dydict=haddFiles(DYFiles,"DY",histonames, [shape,shape],folder,era)
    topdict=haddFiles(TopFiles,"Top",histonames, [shape,shape],folder,era)
    dibdict=haddFiles(DibosonFiles,"Diboson",histonames, [shape,shape],folder,era)
    rewdict=haddFiles(WMuFiles_rew,"ewk",histonames, [shape,shape],"../config/powheg_qtA4_{}/".format(era),era)

    hWtau = taudict['ewk']
    hWtau_sumw2 = taudict['ewk_sumw2']
    hDY = dydict['ewk']
    hDY_sumw2 = dydict['ewk_sumw2']
    hTop = topdict['ewk']
    hTop_sumw2 = topdict['ewk_sumw2']
    hDiboson = dibdict['ewk']
    hDiboson_sumw2 = dibdict['ewk_sumw2']

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

    shape_hel = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 9)
    shape_mass = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 9*2)
    shape_pdf = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 9*103)
    shape_qcd = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 9*9)
    shape_lowacc = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1,len(qtBins_syst)-1)
    shape_mass_lowacc = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1, 2)
    shape_pdf_lowacc = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1, 103)
    shape_qcd_lowacc = (len(etaBins)-1,len(ptBins)-1,len(chargeBins)-1,len(mTBins)-1,len(isoBins)-1,len(qtBins_syst)-1,9)

    histonames = ['signalTemplates', 'signalTemplates_sumw2','signalTemplates_mass','signalTemplates_LHEPdfWeight', 'lowacc', 'lowacc_sumw2','lowacc_mass','lowacc_LHEPdfWeight','lowacc_LHEScaleWeight', 'lowacc_rew']
    wdict=haddFiles(WMuFiles,"WmuSignal",histonames, [shape_hel,shape_hel,shape_mass,shape_pdf,shape_lowacc,shape_lowacc,shape_mass_lowacc,shape_pdf_lowacc,shape_qcd_lowacc,shape_lowacc],"../config/powheg_acc_{}/".format(era),era)

    # signal: nominal
    hsignal = wdict['signalTemplates']
    hsignal_sumw2 = wdict['signalTemplates_sumw2']

    # signal: mass
    hsignal_mass = wdict['signalTemplates_mass']
    hsignal_mass=hsignal_mass.reshape(hsignal_mass.shape[:-1] + (9,2))
    
    # signal: PDFs
    hsignal_PDF = wdict['signalTemplates_LHEPdfWeight']
    hsignal_PDF=hsignal_PDF.reshape(hsignal_PDF.shape[:-1] + (9,103))

    # # # signal: QCD
    # hsignal_QCD = np.array(fsignal['signalTemplates_LHEScaleWeight'][:])[:,:,-1,...,:threshold_y,:threshold_qt,:]
    # hsignal_QCD=hsignal_QCD.reshape(hsignal_QCD.shape[:-1] + (9,9))

    hWmu = np.sum(hsignal,axis=(-1,-2,-3)) # integrated over helicity, y and qt
    hWmu_sumw2 = np.sum(hsignal_sumw2,axis=(-1,-2,-3)) # integrated over helicity, y and qt

    # helicity xsecs without acceptance cuts for unfolding
    file_gen = '../config/powheg_acc_{}/WPlusJetsToMuNu_helweights.hdf5'.format(era)
    f_gen = h5py.File(file_gen, mode='r+')

    htot = f_gen['totxsecs'][:]
    h = f_gen['xsecs'][:]

    htot_PDF = f_gen['totxsecs_LHEPdfWeight'][:]
    h_PDF = f_gen['xsecs_LHEPdfWeight'][:]

    yBins = f_gen['edges_totxsecs_0'][:]
    qtBins = f_gen['edges_totxsecs_1'][:]

    factors = np.array([[20./3., 1./10],[5.,0.],[20.,0.],[4.,0.],[4.,0.],[5.,0.],[5.,0.],[4.,0.],[1.,0.]])
    factors = factors[np.newaxis,np.newaxis,...]
    h = (h/htot[...,np.newaxis]+factors[...,1])*factors[...,0]

    fig, ax1 = plt.subplots()
    hep.hist2dplot(h[...,4],yBins,qtBins)
    plt.tight_layout()
    plt.savefig("testprefit/A0_{}.png".format(era))
    plt.clf()

    fig, ax1 = plt.subplots()
    hep.hist2dplot(htot,yBins,qtBins)
    plt.tight_layout()
    plt.savefig("testprefit/htot_{}.png".format(era))
    plt.clf()

    factors = factors[...,np.newaxis]
    h_PDF = h_PDF.reshape(len(yBins)-1, len(qtBins)-1, 9, 103)
    h_PDF = (h_PDF/htot_PDF[:,:,np.newaxis,:]+factors[:,:,:,1,...])*factors[:,:,:,0,...]

    factors_hel = np.array([2.,2*sqrt(2),4.,4.*sqrt(2),2.,2.,2.*sqrt(2),4.*sqrt(2),1.])
    factors_hel = factors_hel[np.newaxis,np.newaxis,...]
    
    h = 3./(16.*pi)*h*htot[...,np.newaxis]/factors_hel[...,:threshold_y,:threshold_qt,:] #helicity xsecs
    h[...,-1] = 3./(16.*pi)*htot

    # pseudo_data = rewdict['ewk']
    # load aMC@NLO coefficients

    f_aMC = ROOT.TFile.Open('/scratchnvme/emanca/wproperties-analysis/Common/data/genInfo_syst.root')
    qt_aMC = np.sum(hist2array(f_aMC.Get('angularCoefficients_Wplus/YqTcT')),axis=-1)
    fcoeffs = ROOT.TFile.Open('/scratchnvme/emanca/wproperties-analysis/Common/data/genInput_v7_syst_Wplus.root')
    
    hists_aMC = []
    for i in range(5):
        tmp = hist2array(fcoeffs.Get('angularCoefficients/harmonicsA{}_nom_nom'.format(i)))
        hists_aMC.append(tmp)
    coeffs_aMC = np.stack(hists_aMC,axis=-1)

    mapTot = hist2array(f_aMC.Get('angularCoefficients_Wplus/mapTot'))
    h_aMC = 3./(16.*pi)*coeffs_aMC*mapTot[...,np.newaxis]/factors_hel[...,:threshold_y,:threshold_qt,:5] #helicity xsecs aMC
    print('done')

    factors_hel = factors_hel[...,np.newaxis]
    print(h_PDF.shape,htot_PDF.shape,factors_hel[...,:threshold_y,:threshold_qt,:].shape)
    h_PDF = 3./(16.*pi)*h_PDF*htot_PDF[:,:,np.newaxis,...]/factors_hel[...,:threshold_y,:threshold_qt,:,:]
    print(h_PDF.shape)
    h_PDF[...,-1,:] = 3./(16.*pi)*htot_PDF

    pseudo_data_tot = np.copy(hsignal)
    if era =="preVFP": 
        mapTot*=19.514702645/35.9
        h_aMC*=19.514702645/35.9
    else: 
        mapTot*=16.810812618/35.9
        h_aMC*=16.810812618/35.9
    for i in range(5):
        pseudo_data_tot[...,i] *= h_aMC[...,i]/h[...,i]
    pseudo_data_tot[...,-1] *= 3./(16.*pi)*mapTot/h[...,-1]
    pseudo_data = np.sum(pseudo_data_tot,axis=(-1,-2,-3)) # integrated over helicity, y and qt

    # events falling out of fit range
    hlowacc = np.sum(wdict['lowacc'],axis=-1)
    hlowacc_sumw2 = np.sum(wdict['lowacc_sumw2'],axis=-1)
    # lowacc: mass
    hlowacc_mass = wdict['lowacc_mass']
    # lowacc: PDF
    hlowacc_PDF = wdict['lowacc_LHEPdfWeight']
    # lowacc: QCD
    hlowacc_QCD = wdict['lowacc_LHEScaleWeight']
    hlowacc_rew = np.sum(wdict['lowacc_rew'],axis=-1)

    pseudo_data+=hlowacc_rew

    # this is pseudodata made with the sum of the templates
    dset = fshapes.create_dataset(name='data_obs', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] =  pseudo_data[:,:,-1,:,:]

    hewk = hWtau+hDY+hTop+hDiboson+hWmu+hlowacc

    fig, ax1 = plt.subplots()
    ax1.set_title("total xsec closure", fontsize=9)
    hep.hist2dplot(hWmu[:,:,-1,1,0],etaBins,ptBins)
    plt.tight_layout()
    plt.savefig("testprefit/total_xsec_clos_{}".format(era))
    plt.clf()

    fig, ax1 = plt.subplots()
    ax1.set_title("pseudodata", fontsize=9)
    hep.hist2dplot(pseudo_data[:,:,-1,1,0],etaBins,ptBins)
    plt.tight_layout()
    plt.savefig("testprefit/pseudo_data_{}".format(era))
    plt.clf()

    fig, ax1 = plt.subplots()
    ax1.set_title("total xsec closure", fontsize=9)
    hep.hist2dplot(np.sum(hsignal_PDF,axis=(-2,-3,-4))[:,:,-1,1,0,1]/hWmu[:,:,-1,1,0],etaBins,ptBins)
    plt.tight_layout()
    plt.savefig("testprefit/total_xsec_clos_pdf_{}".format(era))
    plt.clf()

    dset = fshapes.create_dataset(name='template', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 9), dtype='float64')
    dset[...] = hsignal[:,:,-1,...] #select positive charge

    dset = fshapes.create_dataset(name='template_sumw2', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 9), dtype='float64')
    dset[...] = hsignal_sumw2[:,:,-1,...] #select positive charge

    dset = fshapes.create_dataset(name='lowacc', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hlowacc[:,:,-1,...] #select positive charge

    dset = fshapes.create_dataset(name='lowacc_sumw2', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1), dtype='float64')
    dset[...] = hlowacc_sumw2[:,:,-1,...] #select positive charge

    # # now write shapes with systematics

    dset = fshapes.create_dataset(name='template_mass', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 9, 2), dtype='float64')
    dset[...] = hsignal_mass[:,:,-1,...] #select positive charge

    dset = fshapes.create_dataset(name='template_LHEPdfWeight', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1, len(yBins)-1, len(qtBins)-1, 9, 103), dtype='float64')
    dset[...] = hsignal_PDF[:,:,-1,...] #select positive charge

    # dset = fshapes.create_dataset(name='template_LHEScaleWeight', shape=hsignal_QCD.shape, dtype='float64')
    # dset[...] = hsignal_QCD #select positive charge

    dset = fshapes.create_dataset(name='lowacc_mass', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1,2), dtype='float64')
    dset[...] = hlowacc_mass[:,:,-1,...]

    dset = fshapes.create_dataset(name='lowacc_LHEPdfWeight', shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1,103), dtype='float64')
    dset[...] = hlowacc_PDF[:,:,-1,...]

    # patch the qcd variations for low acc
    for i in range(len(qtBins_syst)-1):
        fig, ax1 = plt.subplots(figsize=(48, 10))
        hlowacc_QCD_patched = np.copy(wdict['lowacc'])
        histos_scale=[]
        hep.histplot(np.sum(hlowacc_QCD_patched,axis=(0,1))[-1,1,0,:],qtBins_syst)
        for j in range(9):
            hlowacc_QCD_patched[...,i] = hlowacc_QCD[...,i,j]
            histos_scale.append(np.sum(hlowacc_QCD_patched,axis=-1))
            hep.histplot(np.sum(hlowacc_QCD_patched,axis=(0,1))[-1,1,0,:],qtBins_syst)
        plt.tight_layout()
        plt.savefig("testprefit/testnuisance{}_{}.png".format(i,era))
        plt.clf()
        arrays_scale = np.stack(histos_scale,axis=-1)
        dset = fshapes.create_dataset(name='lowacc_LHEScaleWeight{}'.format(i), shape=(len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1,9), dtype='float64')
        dset[...] = arrays_scale[:,:,-1,...]

    # check LHE nuisances
    qtBins_syst=np.array(qtBins_syst)
    qtBinsS = qtBins_syst[1:]-qtBins_syst[:-1]
    qtBinsC = 0.5*(qtBins_syst[1:]+qtBins_syst[:-1])
    hlowacc_qt = np.sum(wdict['lowacc_LHEScaleWeight'][:,:,-1,-1,0,:,:],axis=(0,1))
    hlowacc_qt_rew = np.sum(wdict['lowacc_rew'][:,:,-1,-1,0,:],axis=(0,1))
    hlowacc_qt_nom = np.sum(wdict['lowacc'][:,:,-1,-1,0,:],axis=(0,1))
    hlowacc_qt_pdf = np.sum(wdict['lowacc_LHEPdfWeight'][:,:,-1,-1,0,:],axis=(0,1))
    # print(np.where(np.isnan(wdict['lowacc_rew'])), wdict['lowacc_rew'][np.where(np.isnan(wdict['lowacc_rew']))])
    fig, ax1 = plt.subplots()

    hep.histplot(hlowacc_qt_rew/qtBinsS,qtBins_syst,histtype = 'errorbar', label="pseudodata")
    print(hlowacc_qt_rew/qtBinsS,np.sum(wdict['lowacc'][:,:,-1,-1,0,:],axis=(0,1))/qtBinsS)
    qcdsyst = [0, 1, 3, 5, 7, 8]
    for iqcd in qcdsyst:
        if iqcd<4:
            if iqcd==0:
                hep.histplot(0.8*hlowacc_qt[...,iqcd]/qtBinsS,qtBins_syst, color="green", label = "qcd nuisance")
            else:
                hep.histplot(0.8*hlowacc_qt[...,iqcd]/qtBinsS,qtBins_syst, color="green")
        else:
            hep.histplot(1.2*hlowacc_qt[...,iqcd]/qtBinsS,qtBins_syst, color="red")
    # err_pdf = np.sqrt(np.sum(np.square(hlowacc_qt_nom[:,np.newaxis]-hlowacc_qt_pdf[...,:]),axis=-1))/qtBinsS
    # ax1.fill_between(qtBinsC, (hlowacc_qt_nom/qtBinsS-err_pdf).ravel(), (hlowacc_qt_nom/qtBinsS+err_pdf).ravel(), alpha=0.2, color="orange",label='pdf uncertainty')

    ax1.legend(loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig("testprefit/checknuis_{}".format(era))
    plt.clf()
    
    dset = fshapes.create_dataset(name='helicity', shape=h.shape, dtype='float64')
    dset[...] = h
    dset = fshapes.create_dataset(name='helicity_LHEPdfWeight', shape=h_PDF.shape, dtype='float64')
    dset[...] = h_PDF

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
    # for i in range(hfakesLowMt.shape[0]*hfakesLowMt.shape[1]):
    #     mask = np.zeros(hfakesLowMt.shape[0]*hfakesLowMt.shape[1])
    #     mask[i,...]=1
    #     mask = mask.reshape((hfakesLowMt.shape[0],hfakesLowMt.shape[1]))[...,np.newaxis,np.newaxis]
    #     # nuisance for changing the normalisations independently

    #     hfakesLowMtVarUp = np.where(mask==0, hfakesLowMt, hfakesLowMt+0.5*hfakesLowMt)
    #     dset = fshapes.create_dataset(name='fakesLowMt_fakeNormLowMtBin{}Up'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    #     dset[...] = hfakesLowMtVarUp.flatten()
    #     hfakesLowMtVarDown = np.where(mask==0, hfakesLowMt, hfakesLowMt-0.5*hfakesLowMt)
    #     dset = fshapes.create_dataset(name='fakesLowMt_fakeNormLowMtBin{}Down'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    #     dset[...] = hfakesLowMtVarDown.flatten()

    #     hfakesHighMtVarUp = np.where(mask==0, hfakesHighMt, hfakesHighMt+0.5*hfakesHighMt)
    #     dset = fshapes.create_dataset(name='fakesHighMt_fakeNormHighMtBin{}Up'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    #     dset[...] = hfakesHighMtVarUp.flatten()
    #     hfakesHighMtVarDown = np.where(mask==0, hfakesHighMt, hfakesHighMt-0.5*hfakesHighMt)
    #     dset = fshapes.create_dataset(name='fakesHighMt_fakeNormHighMtBin{}Down'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    #     dset[...] = hfakesHighMtVarDown.flatten()

    #     # print('checking if any zero or negative yields for bin {}'.format(i))
    #     # print(np.any((hfakesLowMtVarUp+hfakesHighMtVarUp)<=0.))
    #     # print(np.any((hfakesLowMtVarDown+hfakesHighMtVarDown)<=0.))

    #     # common nuisance for changing fake shape

    #     norm = np.sum(hfakesLowMt[:,:,0,:],axis=2)
    #     ratio = hfakesLowMt[:,:,0,0]/norm #ratio iso/iso+aiso
    #     rate_var = 2.
    #     var_idx = np.nonzero(mask)
    #     # set to nominal
    #     hfakesLowMtVarUp = np.empty_like(hfakesLowMt)
    #     hfakesLowMtVarDown = np.empty_like(hfakesLowMt)
    #     np.copyto(hfakesLowMtVarUp, hfakesLowMt) # (dst, src)
    #     np.copyto(hfakesLowMtVarDown, hfakesLowMt) # (dst, src)
    #     # apply variation to isolated part
    #     hfakesLowMtVarUp[var_idx[0],var_idx[1],0, 0] = (rate_var*ratio*norm)[var_idx[0],var_idx[1]]
    #     hfakesLowMtVarDown[var_idx[0],var_idx[1],0, 0] = ((1./rate_var)*ratio*norm)[var_idx[0],var_idx[1]]
    #     # apply variation to anti-isolated part
    #     hfakesLowMtVarUp[var_idx[0],var_idx[1],0, 1] = ((1-rate_var*ratio)*norm)[var_idx[0],var_idx[1]]
    #     hfakesLowMtVarDown[var_idx[0],var_idx[1],0, 1] = ((1-(1./rate_var)*ratio)*norm)[var_idx[0],var_idx[1]]

    #     dset = fshapes.create_dataset(name='fakesLowMt_fakeShapeBin{}Up'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    #     dset[...] = (hfakesLowMtVarUp/hfakes_unc).flatten()
    #     dset = fshapes.create_dataset(name='fakesLowMt_fakeShapeBin{}Down'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    #     dset[...] = (hfakesLowMtVarDown/hfakes_unc).flatten()

    #     norm = np.sum(hfakesHighMt[:,:,1,:],axis=2)
    #     ratio = hfakesHighMt[:,:,1,0]/norm #ratio iso/iso+aiso
    #     rate_var = 1.2
    #     var_idx = np.nonzero(mask)
    #     # set to nominal
    #     hfakesHighMtVarUp = np.empty_like(hfakesHighMt)
    #     hfakesHighMtVarDown = np.empty_like(hfakesHighMt)
    #     np.copyto(hfakesHighMtVarUp, hfakesHighMt) # (dst, src)
    #     np.copyto(hfakesHighMtVarDown, hfakesHighMt) # (dst, src)
    #     # apply variation to isolated part
    #     hfakesHighMtVarUp[var_idx[0],var_idx[1],1, 0] = (rate_var*ratio*norm)[var_idx[0],var_idx[1]]
    #     hfakesHighMtVarDown[var_idx[0],var_idx[1],1, 0] = ((1./rate_var)*ratio*norm)[var_idx[0],var_idx[1]]
    #     # apply variation to anti-isolated part
    #     hfakesHighMtVarUp[var_idx[0],var_idx[1],1, 1] = ((1-rate_var*ratio)*norm)[var_idx[0],var_idx[1]]
    #     hfakesHighMtVarDown[var_idx[0],var_idx[1],1, 1] = ((1-(1./rate_var)*ratio)*norm)[var_idx[0],var_idx[1]]

    #     dset = fshapes.create_dataset(name='fakesHighMt_fakeShapeBin{}Up'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    #     dset[...] = (hfakesHighMtVarUp/hfakes_unc).flatten()
    #     dset = fshapes.create_dataset(name='fakesHighMt_fakeShapeBin{}Down'.format(i), shape=[(len(etaBins)-1) * (len(ptBins)-1) * (len(mTBins)-1) * (len(isoBins)-1)], dtype='float64')
    #     dset[...] = (hfakesHighMtVarDown/hfakes_unc).flatten()

    #     # print('checking if any zero or negative yields for bin {}'.format(i))
    #     # print(np.any((hfakesLowMtVarUp+hfakesHighMtVarUp)<=0.))
    #     # print(np.any((hfakesLowMtVarDown+hfakesHighMtVarDown)<=0.))

    # plot pt, eta and mt in isolated and antiisolated region
    for i in range(2):
        fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
        ax1.set_title("eta_iso{}_highMt".format(i), fontsize=9)
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
        ax1.set_title("pt_iso{}_highMt".format(i), fontsize=9)
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
        ax1.set_title("eta_iso{}_lowMt".format(i), fontsize=9)
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
        ax1.set_title("pt_iso{}_lowMt".format(i), fontsize=9)
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
    # ax1.set_title("mt", fontsize=9)
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