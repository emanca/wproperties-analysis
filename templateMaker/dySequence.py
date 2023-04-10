import os
import sys
import pickle
import pathlib
import h5py
import hdf5plugin
import lz4.frame
sys.path.append('../Common/data')
sys.path.append('python')
from defineWeight import defineWeight
from getHelWeights import getHelWeights
from getMassWeights import getMassWeights
from muonSelection import muonSelection
from getSFVariations import getSFVariations
from getPrefVariations import getPrefVariations

import ROOT
# verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kDebug+10)
import argparse
import copy
import time
from datetime import datetime
import narf
import hist
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from math import pi, sqrt
from RDFtree import RDFtree
import logging
from utilities import boostHistHelpers as hh,common,output_tools,logging

import wremnants
from wremnants import muon_calibration, muon_selections

parser,initargs = common.common_parser(True)
parser = common.set_parser_default(parser, "filterProcs", common.vprocs+["dataPostVFP"])

args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
    
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles,
                                              filt=args.filterProcs,
                                              excl=args.excludeProcs, 
                                              nanoVersion="v8" if args.v8 else "v9", base_path=args.dataPath)

logger.debug(f"Will process samples {[d.name for d in datasets]}")
print([d.name for d in datasets])
print(args.filterProcs)
print(args.excludeProcs)

wremnants_dir = f"{pathlib.Path(__file__).parent}/../wremnants"
data_dir = f"{wremnants_dir}/data/"

ROOT.gSystem.Load('bin/libAnalysis.so')

# ROOT.ROOT.EnableImplicitMT()


# qts_axis = hist.axis.Variable([0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.,200.], name = "Zpt")
# ys_axis = hist.axis.Variable([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0], name = "Zrap")
# qts_axis_red = hist.axis.Variable([0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.], name = "Zpt")
# ys_axis_red = hist.axis.Variable([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4], name = "Zrap")
qts_axis = hist.axis.Variable([0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.,200.], name = "Zpt")
ys_axis = hist.axis.Variable([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0], name = "Zrap")
qts_axis_red = hist.axis.Variable([0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.], name = "Zpt")
ys_axis_red = hist.axis.Variable([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4], name = "Zrap")
costhetas_axis =  hist.axis.Regular(100, -1.,1, name = "costheta")
phis_axis =  hist.axis.Regular(100, 0.,2.*pi, name = "phi")
pts_axis =  hist.axis.Regular(60, 25.,55., name = "mupt")
etas_axis =  hist.axis.Regular(48, -2.4,2.4, name = "mueta")
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")

def build_graph(df, dataset):
    
    logger.info(f"build graph for dataset: {dataset.name}")
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs

    if isW or isZ:

        p = RDFtree(df)
        p.branch(nodeToStart='input', nodeToEnd='defs', modules=[defineWeight(dataset,True),ROOT.genLeptonSelector(), ROOT.CSvariableProducer(), ROOT.genVProducer()])

        weightsum = p.EventCount('defs', "weight")

        axes = [ys_axis,qts_axis,costhetas_axis,phis_axis]
        cols = ["Vrap_preFSR_abs","Vpt_preFSR","CStheta_preFSR","CSphi_preFSR","weight"]
        p.Histogram('defs', 'xsecs', [*cols], axes)

        results = p.getObjects()
    
    else:
        p = RDFtree(df)
        p.branch(nodeToStart='input', nodeToEnd='defs', modules=[defineWeight(dataset,True)])

        weightsum = p.EventCount('defs', "weight")
        results = []

    return results, weightsum

#give it to narf
# resultdict = narf.build_and_run(datasets, build_graph)

# print(resultdict['ZmumuPostVFP'].keys())

# outfile = "test.hdf5"
# with h5py.File(outfile, 'w') as f:
#     narf.ioutils.pickle_dump_h5py("results", resultdict, f)

infile = "test.hdf5"
with h5py.File(infile, "r") as f:
    results = narf.ioutils.pickle_load_h5py(f["results"])
    hFullAcc = results['ZmumuPostVFP']["output"]["xsecs"].get()

lumi = results['dataPostVFP']["lumi"]
xsec = results['ZmumuPostVFP']["dataset"]["xsec"]
yBins = hFullAcc.axes[0].edges
yBinsS = yBins[1:]-yBins[:-1]

# pick the right harmonics and reweight
cosThetaBins = hFullAcc.axes[2].edges
phiBins = hFullAcc.axes[-1].edges
yBins = hFullAcc.axes[0].edges
qtBins = hFullAcc.axes[1].edges

cosThetaBins = np.asarray(cosThetaBins)
cosThetaBinsC = 0.5*(cosThetaBins[1:]+cosThetaBins[:-1])
phiBins = np.asarray(phiBins)
phiBinsC = 0.5*(phiBins[1:]+phiBins[:-1])

P0w = np.outer(1. / 2. * (1. - 3. * cosThetaBinsC * cosThetaBinsC),np.ones(len(phiBinsC)))
P1w = np.outer(2. * cosThetaBinsC * np.sqrt(1. - cosThetaBinsC * cosThetaBinsC),np.cos(phiBinsC))
P2w = np.outer(1. / 2. * (1. - cosThetaBinsC * cosThetaBinsC),np.cos(2. * phiBinsC))
P3w = np.outer(np.sqrt(1. - cosThetaBinsC * cosThetaBinsC),np.cos(phiBinsC))
P4w = np.outer(cosThetaBinsC,np.ones(len(phiBinsC)))
P5w = np.outer((1. - cosThetaBinsC * cosThetaBinsC),np.sin(2. * phiBinsC))
P6w = np.outer(2. * cosThetaBinsC * np.sqrt(1. - cosThetaBinsC * cosThetaBinsC),np.sin(phiBinsC))
P7w = np.outer(np.sqrt(1. - cosThetaBinsC * cosThetaBinsC),np.sin(phiBinsC))
P8w = np.outer(1 + cosThetaBinsC * cosThetaBinsC,np.ones(len(phiBinsC)))

wharmonics = [P0w,P1w,P2w,P3w,P4w,P5w,P6w,P7w]
hharmonics = []
hFullAcc = hFullAcc*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]
totalxsec = hFullAcc.project(0,1)

factors = [(20./3., 1./10),(5.,0.),(20.,0.),(4.,0.),(4.,0.),(5.,0.),(5.,0.),(4.,0.)]
for i,hw in enumerate(wharmonics):
    htmp = np.einsum('ijkm,km->ij',hFullAcc.values(),hw)/totalxsec.values()
    htmp = factors[i][0]*(htmp+factors[i][1])
    hharmonics.append(htmp)

hharmonics.append(totalxsec.values())
hharmonics = np.stack(hharmonics,axis=-1)

hharmonics[:,:,0]*=(1./2.)
hharmonics[:,:,1]*=(1./(2.*sqrt(2)))
hharmonics[:,:,2]*=(1./4.)
hharmonics[:,:,3]*=(1./(4.*sqrt(2)))
hharmonics[:,:,4]*=(1./2.)
hharmonics[:,:,5]*=(1./2.)
hharmonics[:,:,6]*=(1./(2.*sqrt(2)))
hharmonics[:,:,7]*=(1./(4.*sqrt(2)))

#convert np array to boost histogram and write to file
#create boost histogram as placeholder using hist package
hist_coeffs = hist.Hist(ys_axis,qts_axis,hist.axis.Integer(-1, 8, underflow=False, overflow=False, name='helicity'), name = "hist_coeffs",
data = hharmonics
    )

outfile = "angCoeffZ.hdf5"
with h5py.File(outfile, 'w') as f:
    narf.ioutils.pickle_dump_h5py("angCoeffZ", hist_coeffs, f)

## helper definition - will need better organization maybe

era="2016PostVFP"
sfFile = "allSmooth_GtoH.root"
sfFile = f"{data_dir}/testMuonSF/{sfFile}"
pileup_helper = wremnants.make_pileup_helper(era = era)
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)
print(args.sfFile)
muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = pts_axis.edges[-1],
                                                                                                                                     is_w_like = True, directIsoSFsmoothing=args.directIsoSFsmoothing)
mc_jpsi_crctn_helper, data_jpsi_crctn_helper = muon_calibration.make_jpsi_crctn_helpers(args)
mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = muon_calibration.make_muon_calibration_helpers(args)
# smearing_helper = muon_calibration.make_muon_smearing_helpers() if args.smearing else None
# bias_helper = muon_calibration.make_muon_bias_helpers(args) if args.biasCalibration else None
# corr_helpers = theory_corrections.load_corr_helpers(common.vprocs, args.theoryCorr)

def massWeightNames(matches=None, proc=""):
    central=10
    nweights=21
    names = [f"massShift{int(abs(central-i)*10)}MeV{'' if i == central else ('Down' if i < central else 'Up')}" for i in range(nweights)]
    if proc and proc in common.zprocs_all:
        # This is the PDG uncertainty (turned off for now since it doesn't seem to have been read into the nano)
        names.extend(["massShift2p1MeVDown", "massShift2p1MeVUp"])
    # If name is "" it won't be stored
    return [x if not matches or any(y in x for y in matches) else "" for x in names]

def build_graph_templates(df, dataset):
    
    logger.info(f"build graph for dataset: {dataset.name}")
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs

    p = RDFtree(df)
    p.branch(nodeToStart='input', nodeToEnd='event_count', modules=[defineWeight(dataset,True)])
    p.branch(nodeToStart='event_count', nodeToEnd='calibrations', modules=[muonSelection(dataset,args,data_calibration_helper,mc_calibration_helper,data_jpsi_crctn_helper,mc_jpsi_crctn_helper),defineWeight(dataset,False,pileup_helper,muon_efficiency_helper,muon_prefiring_helper)])
    weightsum = p.EventCount('event_count', "weight")
    if isZ:
        p.branch(nodeToStart='calibrations', nodeToEnd='theoryTools', modules=[ROOT.genLeptonSelector(), ROOT.CSvariableProducer(), ROOT.genVProducer(),ROOT.defineHarmonics(),getHelWeights("angCoeffZ.hdf5"),ROOT.getWeights(),getMassWeights(isW=False)])
        
        axes = [ys_axis_red,qts_axis_red,etas_axis,pts_axis,axis_charge]
        nom_cols = ["Vrap_preFSR_abs","Vpt_preFSR", "trigMuons_eta0", "trigMuons_pt0", "trigMuons_charge0"]
        helicity_axis = hist.axis.StrCategory(["L", "I", "T", "A", "P","UL"], name = "helicities")
        
        ## signal templates decomposed by helicity

        p.EventFilter(nodeToStart='theoryTools', nodeToEnd='signalTemplates_nominal', evfilter="Vrap_preFSR_abs<2.4 && Vpt_preFSR<60.", filtername="{:20s}".format("signalTemplates_nominal"))

        # nominal histogram
        p.Histogram('signalTemplates_nominal', 'signalTemplates_nominal', [*nom_cols,'helWeightTensor'], axes, tensor_axes= [helicity_axis])
        # mass variations
        p.Histogram('signalTemplates_nominal', 'signalTemplates_mass', [*nom_cols,'massWeight_tensor_hel'], axes,tensor_axes=[helicity_axis,hist.axis.StrCategory(massWeightNames(proc=dataset.name), name="massShift")])

        # sf variations
        # p.branch(nodeToStart='signalTemplates_nominal', nodeToEnd='signalTemplates_sf', modules=[getSFVariations(True,muon_efficiency_helper_stat,muon_efficiency_helper_syst)])
        # for key,helper in muon_efficiency_helper_stat.items():
        #     p.Histogram('signalTemplates_sf', f'signalTemplates_effStatTnP_{key}', [*nom_cols,f"effStatTnP_{key}_tensor_hel"], axes, tensor_axes = [helicity_axis,*(helper.tensor_axes)])
        # p.Histogram('signalTemplates_sf', 'signalTemplates_effSystTnP', [*nom_cols,"effSystTnP_tensor_hel"], axes, tensor_axes = [helicity_axis,*(muon_efficiency_helper_syst.tensor_axes)])

        # #prefiring variations
        # p.Histogram('signalTemplates_nominal', 'signalTemplates_mass', [*nom_cols,'massWeight_tensor_hel'], axes,tensor_axes=[helicity_axis,hist.axis.StrCategory(massWeightNames(proc=dataset.name), name="massShift")])

        ## low acc
        axes = [etas_axis,pts_axis,axis_charge]
        nom_cols = ["trigMuons_eta0", "trigMuons_pt0", "trigMuons_charge0"]
        
        p.EventFilter(nodeToStart='theoryTools', nodeToEnd='lowacc', evfilter="(Vrap_preFSR_abs>2.4 || Vpt_preFSR>60.) && Vpt_preFSR<200.", filtername="{:20s}".format("low acc"))
        p.Histogram('lowacc', 'lowacc', [*nom_cols,"nominal_weight"], axes)

        # mass variation
        p.Histogram('lowacc', 'lowacc_mass', [*nom_cols,"massWeight_tensor_wnom"], axes,tensor_axes=[hist.axis.StrCategory(massWeightNames(proc=dataset.name), name="massShift")])

        # sf variations
        p.branch(nodeToStart='lowacc', nodeToEnd='lowacc_sf', modules=[getSFVariations(True,muon_efficiency_helper_stat,muon_efficiency_helper_syst)])
        # for key,helper in muon_efficiency_helper_stat.items():
        #     p.Histogram('lowacc_sf', f'lowacc_effStatTnP_{key}', [*nom_cols,f"effStatTnP_{key}_tensor"], axes, tensor_axes = helper.tensor_axes)
        # p.Histogram('lowacc_sf', 'lowacc_effSystTnP', [*nom_cols,"effSystTnP_weight"], axes,tensor_axes = muon_efficiency_helper_syst.tensor_axes)

        # prefiring variations
        p.branch(nodeToStart='lowacc', nodeToEnd='lowacc_pref', modules=[getPrefVariations(muon_prefiring_helper_stat,muon_prefiring_helper_syst)])
        p.Histogram('lowacc_pref', 'lowacc_muonL1PrefireStat', [*nom_cols,"muonL1PrefireStat_tensor"], axes,tensor_axes = muon_prefiring_helper_stat.tensor_axes)
        p.Histogram('lowacc_pref', 'lowacc_muonL1PrefireSyst', [*nom_cols,"muonL1PrefireSyst_tensor"], axes,tensor_axes = [common.down_up_axis])
        p.Histogram('lowacc_pref', 'lowacc_ecal1L1Prefire', [*nom_cols,"ecalL1Prefire_tensor"], axes,tensor_axes = [common.down_up_axis])

        results = p.getObjects()
    
    elif isW:

        results = []
    else:
        
        axes = [etas_axis,pts_axis,axis_charge]
        cols = ["trigMuons_eta0", "trigMuons_pt0", "trigMuons_charge0"]

        p.Histogram('calibrations', 'data_obs', [*cols], axes)

        results = p.getObjects()

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph_templates)
print(resultdict['ZmumuPostVFP']['output'].keys())

outfile = "templatesTest2.hdf5"
with h5py.File(outfile, 'w') as f:
    narf.ioutils.pickle_dump_h5py("results", resultdict, f)

infile = "templates.hdf5"
with h5py.File(infile, "r") as f:
    results = narf.ioutils.pickle_load_h5py(f["results"])
    print(results['ZmumuPostVFP']['output'].keys())
    h = results['ZmumuPostVFP']["output"]["signalTemplates_nominal"].get()
    print(h.axes)
    data_obs = results['dataPostVFP']["output"]["data_obs"].get()
    lowacc = results['ZmumuPostVFP']["output"]["lowacc"].get()
    hmass = results['ZmumuPostVFP']["output"]["signalTemplates_mass"].get()
    lowacc_mass = results['ZmumuPostVFP']["output"]["lowacc_mass"].get()
    signalTemplates_effStatTnP_sf_reco = results['ZmumuPostVFP']["output"]["signalTemplates_effStatTnP_sf_reco"].get()
    # full = results['ZmumuPostVFP']["output"]["full"].get()
    # acceptance = results['ZmumuPostVFP']["output"]["acceptance"].get()
    
h = h*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]
lowacc = lowacc*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]
good_idx_mass = [5,15]
hmass = hmass*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]
lowacc_mass = lowacc_mass*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]

print(hmass.values().shape)
print(hmass.values()[:,:,:,:,:,:,good_idx_mass].shape)
good_idx = [0,1,2,3,4,-1]
hhel = hharmonics[...,good_idx]

for i in range(5):
    hhel[...,i]*=hhel[...,-1]
dtype= 'float64'
with h5py.File('templatesFit.hdf5', mode="w") as f:
    dset_templ = f.create_dataset('template', h.values()[...,0,:].shape, dtype=dtype)
    dset_templ[...] = h[...,0,:].values()
    dset_sumw2 = f.create_dataset('template_sumw2', h.variances()[...,0,:].shape, dtype=dtype)
    dset_sumw2[...] = h.variances()[...,0,:]
    dset_templ = f.create_dataset('template_mass', hmass.values()[...,0,:,:][...,good_idx_mass].shape, dtype=dtype)
    dset_templ[...] = hmass.values()[...,0,:,:][...,good_idx_mass]
    dset_hel = f.create_dataset('helicity', hhel.shape, dtype=dtype)
    dset_hel[...] = hhel
    dset_lowacc = f.create_dataset('lowacc', lowacc[...,0].values().shape, dtype=dtype)
    dset_lowacc[...] = lowacc[...,0].values()
    dset_lowacc = f.create_dataset('lowacc_mass', lowacc_mass.values()[...,0,:][...,good_idx_mass].shape, dtype=dtype)
    dset_lowacc[...] = lowacc_mass.values()[...,0,:][...,good_idx_mass]
    dset_data_obs = f.create_dataset('data_obs', data_obs[...,0].values().shape, dtype=dtype)
    dset_data_obs[...] = data_obs[...,0].values()
    dset_lowacc_sumw2 = f.create_dataset('lowacc_sumw2', lowacc[...,0].variances().shape, dtype=dtype)
    dset_lowacc_sumw2[...] = lowacc[...,0].variances()
    dset_data_obs_sumw2 = f.create_dataset('data_obs_sumw2', data_obs[...,0].variances().shape, dtype=dtype)
    dset_data_obs_sumw2[...] = data_obs[...,0].variances()



# etaBins = h.axes[-3].edges
# ptBins = h.axes[-2].edges
# yBins = full.axes[0].edges
# qtBins = full.axes[1].edges

# fig, ax1 = plt.subplots(figsize=(48,48))
# ax1.set_ylabel('Events')
# ratio = acceptance.values()/full.values()
# hep.hist2dplot(np.round(ratio,1),yBins,qtBins,labels=True)
# plt.savefig('acceptanceZ.png')

