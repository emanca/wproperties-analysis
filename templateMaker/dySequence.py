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
from muonCalibration import muonCalibration
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
# custom args
parser.add_argument("-Wlike", "--Wlike", action='store_true', help="run analysis on Wlike")
parser.add_argument("-runHel", "--runHel", action='store_true', help="get helicity cross sections")

args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
    
datasets = wremnants.datasets2016.getDatasets(maxFiles=args.maxFiles,
                                              filt=args.filterProcs,
                                              excl=args.excludeProcs, 
                                              nanoVersion="v8" if args.v8 else "v9", base_path=args.dataPath)

logger.debug(f"Will process samples {[d.name for d in datasets]}")

wremnants_dir = f"{pathlib.Path(__file__).parent}/../wremnants"
data_dir = f"{wremnants_dir}/data/"

ROOT.gSystem.Load('bin/libAnalysis.so')

isWlike = args.Wlike
runHel = args.runHel
era = args.era

qts_axis_W = hist.axis.Variable([0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.,200.],underflow=False, overflow=False, name = "Vpt")
ys_axis_W = hist.axis.Variable([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0], underflow=False, overflow=False, name = "Vrap")
qts_axis_red_W = hist.axis.Variable([0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.], underflow=False, overflow=False, name = "Vpt")
ys_axis_red_W = hist.axis.Variable([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4], underflow=False, overflow=False, name = "Vrap")

qts_axis_Z = hist.axis.Variable([ 0., 3.,  4.8,  6.7, 9., 12., 16.01, 23.6,60,200.], underflow=False, overflow=False, name = "Vpt")
ys_axis_Z = hist.axis.Variable([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0], underflow=False, overflow=False, name = "Vrap")
qts_axis_red_Z = hist.axis.Variable([ 0., 3.,  4.8,  6.7, 9., 12., 16.01, 23.6,60.], underflow=False, overflow=False, name = "Vpt")
ys_axis_red_Z = hist.axis.Variable([0., 0.4, 0.8, 1.2, 1.6, 2.0], underflow=False, overflow=False, name = "Vrap")

costhetas_axis =  hist.axis.Regular(100, -1.,1, underflow=False, overflow=False, name = "costheta")
phis_axis =  hist.axis.Regular(100, 0.,2.*pi, underflow=False, overflow=False, name = "phi")
pts_axis =  hist.axis.Regular(60, 25.,55., underflow=False, overflow=False, name = "mupt")
etas_axis =  hist.axis.Regular(48, -2.4,2.4, underflow=False, overflow=False, name = "mueta")
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
axis_passIso = common.axis_passIso
axis_passMT = common.axis_passMT


def build_graph(df, dataset):
    
    logger.info(f"build graph for dataset: {dataset.name}")
    isW = dataset.name in common.wprocs_recoil #["WplusmunuPostVFP", "WminusmunuPostVFP"]
    isZ = dataset.name in common.zprocs_recoil #["ZmumuPostVFP"]

    if isW or isZ:

        p = RDFtree(df)
        p.branch(nodeToStart='input', nodeToEnd='defs', modules=[defineWeight(dataset,True),ROOT.genLeptonSelector(), ROOT.CSvariableProducer(), ROOT.genVProducer()])

        weightsum = p.EventCount('defs', "weight")

        if isW:
            ys_axis = ys_axis_W
            qts_axis = qts_axis_W
        elif isZ:
            ys_axis = ys_axis_Z
            qts_axis = qts_axis_Z

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

if runHel:
    resultdict = narf.build_and_run(datasets, build_graph)

    outfile = "momentsWZ.hdf5"
    with h5py.File(outfile, 'w') as f:
        narf.ioutils.pickle_dump_h5py("results", resultdict, f)

infile = "momentsWZ.hdf5"

def getHelXsecs(sample):

    with h5py.File(infile, "r") as f:
        results = narf.ioutils.pickle_load_h5py(f["results"])
        hFullAcc = results[sample]["output"]["xsecs"].get()

    lumi = results['dataPostVFP']["lumi"]
    xsec = results[sample]["dataset"]["xsec"]
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
    hFullAcc = hFullAcc*lumi*1000*xsec/results[sample]["weight_sum"]
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
    if sample=="ZmumuPostVFP":
        name="ZmumuPostVFP"
        ys_axis = ys_axis_Z
        qts_axis = qts_axis_Z
    elif sample == "WplusmunuPostVFP":
        name ="WplusmunuPostVFP"
        ys_axis = ys_axis_W
        qts_axis = qts_axis_W
    else:
        name ="WminusmunuPostVFP"
        ys_axis = ys_axis_W
        qts_axis = qts_axis_W
    
    hist_coeffs = hist.Hist(ys_axis,qts_axis,hist.axis.Integer(-1, 8, underflow=False, overflow=False, name='helicity'), name = f"hist_coeffs_{name}",
    data = hharmonics
        )
    return hist_coeffs

hist_coeffs_dict = {}

hist_coeffs_dict[getHelXsecs("ZmumuPostVFP").name]=getHelXsecs("ZmumuPostVFP")
hist_coeffs_dict[getHelXsecs("WplusmunuPostVFP").name]=getHelXsecs("WplusmunuPostVFP")
hist_coeffs_dict[getHelXsecs("WminusmunuPostVFP").name]=getHelXsecs("WminusmunuPostVFP")

outfile = "angCoeffWZ.hdf5"
with h5py.File(outfile, 'w') as f:
    narf.ioutils.pickle_dump_h5py("angCoeffWZ", hist_coeffs_dict, f)

with h5py.File(outfile, "r") as f:
    results = narf.ioutils.pickle_load_h5py(f["angCoeffWZ"])
    print(results)

## helper definition - will need better organization maybe

pileup_helper = wremnants.make_pileup_helper(era = era)
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)
print(args.sfFile)
muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = pts_axis.edges[-1],
                                                                                                                                     is_w_like = False, isoEfficiencySmoothing = args.isoEfficiencySmoothing)
mc_jpsi_crctn_helper, data_jpsi_crctn_helper, jpsi_crctn_MC_unc_helper, jpsi_crctn_data_unc_helper = muon_calibration.make_jpsi_crctn_helpers(args, make_uncertainty_helper=True)

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
    isW = dataset.name in common.wprocs_recoil #["WplusmunuPostVFP", "WminusmunuPostVFP"]
    isZ = dataset.name in common.zprocs_recoil #["ZmumuPostVFP"]

    if not isWlike:
        ys_axis_red = ys_axis_red_W
        qts_axis_red = qts_axis_red_W
        flag = "goodMuons"
    else:
        ys_axis_red = ys_axis_red_Z
        qts_axis_red = qts_axis_red_Z
        flag = "trigMuons"

    p = RDFtree(df)
    p.branch(nodeToStart='input', nodeToEnd='event_count', modules=[defineWeight(dataset,isWlike,True)])
    p.branch(nodeToStart='event_count', nodeToEnd='calibrations', modules=[muonSelection(dataset,args,data_calibration_helper,mc_calibration_helper,data_jpsi_crctn_helper,mc_jpsi_crctn_helper),defineWeight(dataset,isWlike,False,pileup_helper,muon_efficiency_helper,muon_prefiring_helper)])
    weightsum = p.EventCount('event_count', "weight")
    
    if (isZ and isWlike) or (isW and not isWlike):
        p.branch(nodeToStart='calibrations', nodeToEnd='theoryTools', modules=[ROOT.genLeptonSelector(), ROOT.CSvariableProducer(), ROOT.genVProducer(),ROOT.defineHarmonics(),getHelWeights(angFile="angCoeffWZ.hdf5",type=dataset.name),ROOT.getWeights(),getMassWeights(isW=False),muonCalibration(dataset,isWlike,jpsi_crctn_data_unc_helper)])

        # axes = [ys_axis_red,qts_axis_red,etas_axis,pts_axis,axis_charge]
        # nom_cols = ["Vrap_preFSR_abs","Vpt_preFSR", f"{flag}_eta0", f"{flag}_pt0", f"{flag}_charge0"]
        axes = [etas_axis,pts_axis,axis_charge,ys_axis_red,qts_axis_red]
        nom_cols = [f"{flag}_eta0", f"{flag}_pt0", f"{flag}_charge0","Vrap_preFSR_abs","Vpt_preFSR"]

        if isW:
            axes = [etas_axis,pts_axis,axis_charge,axis_passMT,axis_passIso,ys_axis_red,qts_axis_red]
            nom_cols = [f"{flag}_eta0", f"{flag}_pt0", f"{flag}_charge0","passMT","passIso","Vrap_preFSR_abs","Vpt_preFSR"]
        
        helicity_axis = hist.axis.StrCategory(["L", "I", "T", "A", "P","UL"], name = "helicities")
        
        ## signal templates decomposed by helicity
        p.EventFilter(nodeToStart='theoryTools', nodeToEnd='theoryTools', evfilter=f"{flag}_pt0> 25. && {flag}_pt0< 55. && fabs({flag}_eta0)< 2.4", filtername="{:20s}".format("accept"))
        p.EventFilter(nodeToStart='theoryTools', nodeToEnd='signal_nominal', evfilter="Vrap_preFSR_abs<2.4 && Vpt_preFSR<60.", filtername="{:20s}".format("signal_nominal"))
        # p.getCutFlowReport("signal_nominal").Print()

        # nominal histogram
        p.Histogram('signal_nominal', 'signal_nominal', [*nom_cols,'helWeightTensor'], axes, tensor_axes= [helicity_axis])

        storage = hist.storage.Double()
        # mass variations
        p.Histogram('signal_nominal', 'signal_mass_var', [*nom_cols,'massWeight_tensor_hel'], axes,tensor_axes=[helicity_axis,hist.axis.StrCategory(["mass_var"], name="mass_var"),common.down_up_axis],storage=storage)

        # # muon calibration uncertainties
        # p.Histogram('signal_nominal', 'signal_jpsi_var', [*nom_cols,'jpsiWeight_tensor_hel'], axes,tensor_axes = [helicity_axis,*(jpsi_crctn_data_unc_helper.tensor_axes)],storage=storage)
        # nom_cols_smeared = ["Vrap_preFSR_abs","Vpt_preFSR", f"{flag}_genSmearedEta", f"{flag}_genSmearedPt", f"{flag}_genSmearedCharge"]
        # if isW:
        #     nom_cols_smeared.extend(["passMT","passIso"])
        # p.Histogram('signal_nominal', 'signal_nominal_gensmear', [*nom_cols_smeared,'helWeightTensor'], axes, tensor_axes= [helicity_axis],storage=storage)

        # sf variations and prefiring variations
        p.branch(nodeToStart='signal_nominal', nodeToEnd='signal_sf', modules=[getSFVariations(isWlike=isWlike,helper_stat=muon_efficiency_helper_stat,helper_syst=muon_efficiency_helper_syst, helperPref_stat = muon_prefiring_helper_stat, helperPref_syst=muon_prefiring_helper_syst)])
        
        for key,helper in muon_efficiency_helper_stat.items():
            print(f'signal_effStatTnP_{key}',*(helper.tensor_axes))
            p.Histogram('signal_sf', f'signal_effStatTnP_{key}', [*nom_cols,f"effStatTnP_{key}_tensor_hel"], axes, tensor_axes = [helicity_axis,*(helper.tensor_axes)],storage=storage)
        p.Histogram('signal_sf', 'signal_effSystTnP', [*nom_cols,"effSystTnP_tensor_hel"], axes, tensor_axes = [helicity_axis,*(muon_efficiency_helper_syst.tensor_axes)],storage=storage)

        p.Histogram('signal_sf', 'signal_muonL1PrefireStat_tensor', [*nom_cols,'muonL1PrefireStat_tensor_hel'], axes,tensor_axes=[helicity_axis,*(muon_prefiring_helper_stat.tensor_axes)],storage=storage)
        p.Histogram('signal_sf', 'signal_muonL1PrefireSyst_tensor', [*nom_cols,'muonL1PrefireSyst_tensor_hel'], axes,tensor_axes=[helicity_axis,common.down_up_axis],storage=storage)
        p.Histogram('signal_sf', 'signal_ecalL1Prefire_tensor', [*nom_cols,'ecalL1Prefire_tensor_hel'], axes, tensor_axes=[helicity_axis,common.down_up_axis],storage=storage)

        ## low acc
        axes = [etas_axis,pts_axis,axis_charge]
        nom_cols = [f"{flag}_eta0", f"{flag}_pt0", f"{flag}_charge0"]
        
        p.EventFilter(nodeToStart='theoryTools', nodeToEnd='lowacc', evfilter="(Vrap_preFSR_abs>2.4 || Vpt_preFSR>60.) && Vpt_preFSR<200.", filtername="{:20s}".format("low acc"))
        # p.Histogram('lowacc', 'lowacc_nominal', [*nom_cols,"nominal_weight"], axes)

        # mass variation
        # p.Histogram('lowacc', 'lowacc_mass_var', [*nom_cols,"massWeight_tensor_wnom"], axes,tensor_axes=[hist.axis.StrCategory(["mass_var"], name="mass_var"),common.down_up_axis])

        # sf variations
        # p.branch(nodeToStart='lowacc', nodeToEnd='lowacc_sf', modules=[getSFVariations(isWlike=False,helper_stat=muon_efficiency_helper_stat,helper_syst=muon_efficiency_helper_syst)])
        # for key,helper in muon_efficiency_helper_stat.items():
        #     p.Histogram('lowacc_sf', f'lowacc_effStatTnP_{key}', [*nom_cols,f"effStatTnP_{key}_tensor"], axes, tensor_axes = helper.tensor_axes)
        # p.Histogram('lowacc_sf', 'lowacc_effSystTnP', [*nom_cols,"effSystTnP_weight"], axes,tensor_axes = muon_efficiency_helper_syst.tensor_axes)

        # prefiring variations
        # p.branch(nodeToStart='lowacc', nodeToEnd='lowacc_pref', modules=[getPrefVariations(muon_prefiring_helper_stat,muon_prefiring_helper_syst)])
        # p.Histogram('lowacc_pref', 'lowacc_muonL1PrefireStat', [*nom_cols,"muonL1PrefireStat_tensor"], axes,tensor_axes = muon_prefiring_helper_stat.tensor_axes)
        # p.Histogram('lowacc_pref', 'lowacc_muonL1PrefireSyst', [*nom_cols,"muonL1PrefireSyst_tensor"], axes,tensor_axes = [common.down_up_axis])
        # p.Histogram('lowacc_pref', 'lowacc_ecal1L1Prefire', [*nom_cols,"ecalL1Prefire_tensor"], axes,tensor_axes = [common.down_up_axis])

        results = p.getObjects()

    else:
        
        axes = [etas_axis,pts_axis,axis_charge]
        cols = [f"{flag}_eta0", f"{flag}_pt0", f"{flag}_charge0"]

        p.Histogram('calibrations', 'data_obs', [*cols], axes)

        results = p.getObjects()

    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph_templates)
print(resultdict['WplusmunuPostVFP']['output'].keys())

outfile = "templatesTest_testW_swapped.hdf5"
with h5py.File(outfile, 'w') as f:
    narf.ioutils.pickle_dump_h5py("results", resultdict, f)

infile = "templatesTest_testW_swapped.hdf5"
with h5py.File(infile, "r") as f:
    results = narf.ioutils.pickle_load_h5py(f["results"])
    print(results['WplusmunuPostVFP']['output'].keys())
    print(results['WplusmunuPostVFP']["output"]["signal_nominal"].get())
#     print(h.axes)
#     data_obs = results['dataPostVFP']["output"]["data_obs"].get()
#     lowacc = results['ZmumuPostVFP']["output"]["lowacc_nominal"].get()
#     hmass = results['ZmumuPostVFP']["output"]["signal_mass"].get()
#     lowacc_mass = results['ZmumuPostVFP']["output"]["lowacc_mass"].get()
#     lumi = results['dataPostVFP']["lumi"]
#     # signal_effStatTnP_sf_reco = results['ZmumuPostVFP']["output"]["signal_effStatTnP_sf_reco"].get()
#     # full = results['ZmumuPostVFP']["output"]["full"].get()
#     # acceptance = results['ZmumuPostVFP']["output"]["acceptance"].get()
    
# h = h*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]
# lowacc = lowacc*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]
# good_idx_mass = [5,15]
# hmass = hmass*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]
# lowacc_mass = lowacc_mass*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]

# print(hmass.values().shape)
# # print(hmass.values()[:,:,:,:,:,:,good_idx_mass].shape)
# good_idx = [0,1,2,3,4,-1]
# hhel = hharmonics[...,good_idx]

# for i in range(5):
#     hhel[...,i]*=hhel[...,-1]
# dtype= 'float64'
# with h5py.File('templatesFit.hdf5', mode="w") as f:
#     dset_templ = f.create_dataset('template', h.values().shape, dtype=dtype)
#     dset_templ[...] = h.values()
#     dset_sumw2 = f.create_dataset('template_sumw2', h.variances().shape, dtype=dtype)
#     dset_sumw2[...] = h.variances()
#     dset_templ = f.create_dataset('template_mass', hmass.values()[...].shape, dtype=dtype)
#     dset_templ[...] = hmass.values()[...]
#     dset_hel = f.create_dataset('helicity', hhel.shape, dtype=dtype)
#     dset_hel[...] = hhel
#     dset_lowacc = f.create_dataset('lowacc', lowacc.values().shape, dtype=dtype)
#     dset_lowacc[...] = lowacc.values()
#     dset_lowacc = f.create_dataset('lowacc_mass', lowacc_mass.values()[...].shape, dtype=dtype)
#     dset_lowacc[...] = lowacc_mass.values()[...]
#     dset_data_obs = f.create_dataset('data_obs', data_obs.values().shape, dtype=dtype)
#     dset_data_obs[...] = data_obs.values()
#     dset_lowacc_sumw2 = f.create_dataset('lowacc_sumw2', lowacc.variances().shape, dtype=dtype)
#     dset_lowacc_sumw2[...] = lowacc.variances()
#     dset_data_obs_sumw2 = f.create_dataset('data_obs_sumw2', data_obs.variances().shape, dtype=dtype)
#     dset_data_obs_sumw2[...] = data_obs.variances()

# # import pdb; pdb.set_trace()

# # etaBins = h.axes[-3].edges
# # ptBins = h.axes[-2].edges
# # yBins = full.axes[0].edges
# # qtBins = full.axes[1].edges

# # fig, ax1 = plt.subplots(figsize=(48,48))
# # ax1.set_ylabel('Events')
# # ratio = acceptance.values()/full.values()
# # hep.hist2dplot(np.round(ratio,1),yBins,qtBins,labels=True)
# # plt.savefig('acceptanceZ.png')

