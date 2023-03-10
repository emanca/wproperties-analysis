import os
import sys
import pickle
import h5py
import lz4.frame
sys.path.append('../Common/data')
import ROOT
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
from externals import fileSFNov22
from samples2016 import getDatasets

datasets = getDatasets()
ROOT.gSystem.Load('bin/libAnalysis.so')

ROOT.ROOT.EnableImplicitMT()
datasets = [dataset for dataset in datasets if dataset.name=='ZmumuPostVFP']

# p = RDFtree(outputDir = 'testZ', datasets = datasets, outputFile="test.root", pretend=False)
# p.branch(nodeToStart='input', nodeToEnd='defs', modules=[ROOT.genLeptonSelector(), ROOT.CSvariableProducer(), ROOT.genVProducer()])

qts_axis = hist.axis.Variable([0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.,200.], name = "Zpt")
ys_axis = hist.axis.Variable([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0], name = "Zrap")
costhetas_axis =  hist.axis.Regular(100, -1.,1, name = "costheta")
phis_axis =  hist.axis.Regular(100, 0.,2.*pi, name = "phi")

axes = [ys_axis,qts_axis,costhetas_axis,phis_axis]
cols = ["Vrap_preFSR_abs","Vpt_preFSR","CStheta_preFSR","CSphi_preFSR"]
# p.Histogram('defs', 'xsecs', [*cols], axes)
# p.getOutput()

fname = "testZ/test.pkl.lz4"
with (lz4.frame.open(fname, "r")) as openfile:
    resultdict_mc = pickle.load(openfile)

hFullAcc = resultdict_mc["xsecs"]
yBins = hFullAcc.axes[0].edges
yBinsS = yBins[1:]-yBins[:-1]
fig, ax1 = plt.subplots()
ax1.set_title("y", fontsize=18)
ax1.set_ylabel('Events')
ax1.set_xlabel('y')
data = np.sum(hFullAcc.values(),axis=(1,2,3))
# print(data)
hep.histplot(data/yBinsS,bins = hFullAcc.axes[0].edges, histtype = 'errorbar', color = "k", stack = False, ax=ax1)
plt.show()


# fig, ax1 = plt.subplots()
# ax1.set_title("CS theta", fontsize=18)
# ax1.set_ylabel('Events')
# ax1.set_xlabel('CS $\theta$')
# data = np.sum(hFullAcc.values(),axis=(0,1,3))
# # print(data)
# hep.histplot(data,bins = hFullAcc.axes[2].edges, histtype = 'errorbar', color = "k", stack = False, ax=ax1)
# plt.show()

# fig, ax1 = plt.subplots()
# ax1.set_title("CS phi", fontsize=18)
# ax1.set_ylabel('Events')
# ax1.set_xlabel('CS $\phi$')
# data = np.sum(hFullAcc.values(),axis=(0,1,2))
# hep.histplot([data],bins = hFullAcc.axes[-1].edges, histtype = 'errorbar', color = "k", stack = False, ax=ax1)
# plt.show()

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
totalxsec = np.sum(hFullAcc.values(),axis=(2,3))

# fig, ax1 = plt.subplots()
# ax1.set_title("total xsec", fontsize=18)
# hep.hist2dplot(totalxsec,yBins,qtBins)
# plt.savefig("total_xsec")
# plt.clf()

factors = [(20./3., 1./10),(5.,0.),(20.,0.),(4.,0.),(4.,0.),(5.,0.),(5.,0.),(4.,0.)]
for i,hw in enumerate(wharmonics):
    htmp = np.einsum('ijkm,km->ij',hFullAcc.values(),hw)/totalxsec
    htmp = factors[i][0]*(htmp+factors[i][1])
    hharmonics.append(htmp)
#     fig, ax1 = plt.subplots()
#     ax1.set_title("A{}".format(i), fontsize=18)
#     hep.hist2dplot(htmp,yBins,qtBins)
#     plt.savefig("A{}".format(i))
#     plt.cla()


# p = RDFtree(outputDir = 'testZ', datasets = datasets, outputFile="templates.root", pretend=False)
# p.branch(nodeToStart='input', nodeToEnd='defs', modules=[ROOT.genLeptonSelector(), ROOT.CSvariableProducer(), ROOT.genVProducer(),ROOT.trigObjMatchProducer()])
# p.branch(nodeToStart='defs', nodeToEnd='defs', modules=[ROOT.zVetoMuons(),ROOT.zSelection()])
# p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(vetoMuons)==2 && Sum(goodMuons)==2", filtername="{:20s}".format("two muons"))
# p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(HLT_IsoMu24 ||  HLT_IsoTkMu24)", filtername="{:20s}".format("Pass HLT"))
# p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(Muon_charge[goodMuons][0] + Muon_charge[goodMuons][1]) == 0", filtername="{:20s}".format("Opposite charge"))
# p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="std::abs(Muon_eta[goodMuons][0]) < 2.4 && std::abs(Muon_eta[goodMuons][1]) < 2.4", filtername="{:20s}".format("Accept"))
# p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_mediumId[goodMuons][0] == 1 && Muon_mediumId[goodMuons][1] == 1", filtername="{:20s}".format("MuonId"))
# p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_pfRelIso04_all[goodMuons][0] < 0.15 && Muon_pfRelIso04_all[goodMuons][1] < 0.15", filtername="{:20s}".format("Isolation"))
# p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="60. < dimuonMass && dimuonMass < 120.", filtername="{:20s}".format("mZ range"))
# p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("+ve mu trig matched"))

# pts_axis =  hist.axis.Regular(60, 25.,55., name = "mupt")
# etas_axis =  hist.axis.Regular(48, -2.4,2.4, name = "etapt")

# axes = [ys_axis,qts_axis,costhetas_axis,phis_axis,etas_axis,pts_axis]
# cols = ["Vrap_preFSR_abs","Vpt_preFSR","CStheta_preFSR","CSphi_preFSR", "Mu1_eta","Mu1_pt"]
# p.EventFilter(nodeToStart='defs', nodeToEnd='signal', evfilter="Vrap_preFSR_abs<2.4 && Vpt_preFSR<60.", filtername="{:20s}".format("signal templ"))
# p.Histogram('signal', 'signalTemplates', [*cols], axes)

# axes = [etas_axis,pts_axis]
# cols = ["Mu1_eta","Mu1_pt"]
# p.EventFilter(nodeToStart='defs', nodeToEnd='lowacc', evfilter="(Vrap_preFSR_abs>2.4 || Vpt_preFSR>60.) && Vpt_preFSR<200.", filtername="{:20s}".format("low acc"))
# p.Histogram('lowacc', 'lowacc', [*cols], axes)

# p.getOutput()

fname = "testZ/templates.pkl.lz4"
with (lz4.frame.open(fname, "r")) as openfile:
    resultdict_mc = pickle.load(openfile)
h = resultdict_mc['signalTemplates'].values()
herr = np.sqrt(h)

etaBins = resultdict_mc['signalTemplates'].axes[-2].edges
ptBins = resultdict_mc['signalTemplates'].axes[-1].edges

# fig, ax1 = plt.subplots()
# ax1.set_ylabel('Events')
# data = np.sum(h.values(),axis=(0,1,2,3,4))
# hep.histplot([data],bins = h.axes[5].edges, histtype = 'errorbar', color = "k", stack = False, ax=ax1)
# plt.show()

# now derive the weights for the templates
totalxsec = hFullAcc.values() # total differential xsec in y,qt,costheta,phi with no acceptance cuts
norm = np.zeros(totalxsec.shape, dtype='float64')
norm+=totalxsec*P8w

#rescale to get helicity xsecs
hharmonics[0]/=2.
hharmonics[1]/=(2.*sqrt(2))
hharmonics[2]/=4.
hharmonics[3]/=(4.*sqrt(2))
hharmonics[4]/=2.
hharmonics[5]/=2.
hharmonics[6]/=(2.*sqrt(2))
hharmonics[7]/=(4.*sqrt(2))

hharmonics_new = hharmonics
for i,hw in enumerate(wharmonics):
    hharmonics_new[i] = hharmonics[i][..., np.newaxis, np.newaxis]
    norm+=totalxsec * hw * hharmonics_new[i]
norm*=3./(16.*pi)

# import pdb; pdb.set_trace()
print(np.argwhere(totalxsec==0))
print(np.argwhere(norm==0))

wL = np.where(norm!=0,3./(16.*pi) * totalxsec * P0w * hharmonics_new[0]/norm,0.)
wI = np.where(norm!=0,3./(16.*pi) * totalxsec * P1w * hharmonics_new[1]/norm,0.)
wT = np.where(norm!=0,3./(16.*pi) * totalxsec * P2w * hharmonics_new[2]/norm,0.)
wA = np.where(norm!=0,3./(16.*pi) * totalxsec * P3w * hharmonics_new[3]/norm,0.)
wP = np.where(norm!=0,3./(16.*pi) * totalxsec * P4w * hharmonics_new[4]/norm,0.)
w7 = np.where(norm!=0,3./(16.*pi) * totalxsec * P5w * hharmonics_new[5]/norm,0.)
w8 = np.where(norm!=0,3./(16.*pi) * totalxsec * P6w * hharmonics_new[6]/norm,0.)
w9 = np.where(norm!=0,3./(16.*pi) * totalxsec * P7w * hharmonics_new[7]/norm,0.)
wUL = np.where(norm!=0,3./(16.*pi) * totalxsec * P8w/norm,0.)

whelicity = [wL,wI,wT,wA,wP,w7,w8,w9,wUL]
helicities = ['L', 'I', 'T', 'A', 'P', '7','8', '9', 'UL']
hhelicity = []
herrhelicity = []

templdic = {}
templw2dic = {}

# now reduce and produce the templates
for i,hw in enumerate(whelicity):
    htmp = np.einsum('ijmnkl,ijmn->ijkl',h,hw)
    hhelicity.append(htmp)
    templdic[helicities[i]]=htmp
    htmp_err = np.einsum('ijmnkl,ijmn->ijkl',herr,hw)
    herrhelicity.append(htmp_err)
    templw2dic[helicities[i]]=htmp_err
    # plot one slice of every templates
    # fig, ax1 = plt.subplots()
    # ax1.set_title("templates{}".format(i), fontsize=18)
    # hep.hist2dplot(htmp[2,2,:,:],etaBins,ptBins)
    # plt.savefig("templates{}".format(i))
    # plt.clf()

# check closure of templates sum

totalxsec_clos = np.sum(hhelicity[0]+hhelicity[1]+hhelicity[2]+hhelicity[3]+hhelicity[4]+hhelicity[5]+hhelicity[6]+hhelicity[7]+hhelicity[8], axis=(0,1))
# totalxsec = np.sum(h,axis=(0,1,2,3))
# fig, ax1 = plt.subplots()
# ax1.set_title("total xsec closure", fontsize=18)
# hep.hist2dplot(totalxsec_clos[:,:-1],etaBins,ptBins[:-1])
# plt.savefig("total_xsec_clos")
# plt.clf()
# fig, ax1 = plt.subplots()
# ax1.set_title("total xsec etapt", fontsize=18)
# ratio = totalxsec_clos/totalxsec
# hep.hist2dplot(ratio[:,:-1],etaBins,ptBins[:-1])
# plt.savefig("total_xsec_etapt")
# plt.clf()

# save templates to be read from fit and save gen quantities to unfold
dtype = 'float64'
print(hharmonics[0].shape, np.sum(totalxsec,axis=(2,3)).shape)
hharmonics[0]=np.squeeze(hharmonics[0])*np.sum(totalxsec,axis=(2,3))
hharmonics[1]=np.squeeze(hharmonics[1])*np.sum(totalxsec,axis=(2,3))
hharmonics[2]=np.squeeze(hharmonics[2])*np.sum(totalxsec,axis=(2,3))
hharmonics[3]=np.squeeze(hharmonics[3])*np.sum(totalxsec,axis=(2,3))
hharmonics[4]=np.squeeze(hharmonics[4])*np.sum(totalxsec,axis=(2,3))
hharmonics[5]=np.squeeze(hharmonics[5])*np.sum(totalxsec,axis=(2,3))
hharmonics[6]=np.squeeze(hharmonics[6])*np.sum(totalxsec,axis=(2,3))
hharmonics[7]=np.squeeze(hharmonics[7])*np.sum(totalxsec,axis=(2,3))
hharmonics.append(np.sum(totalxsec,axis=(2,3)))


threshold_y = np.digitize(2.4,yBins)
threshold_qt = np.digitize(60.,qtBins)
print(threshold_y,threshold_qt)

templates_all = np.stack(hhelicity,axis=0)
templatesw2_all = np.stack(herrhelicity,axis=0)
helicities_all = 3./(16.*pi)*np.stack(hharmonics,axis=0)
lowacc = resultdict_mc['lowacc'].values()
lowacc_err = np.sqrt(lowacc)

print(templates_all.shape, 'templates')
fig, ax1 = plt.subplots()
ax1.set_title("low acceptance template", fontsize=18)
hep.hist2dplot(lowacc,etaBins,ptBins)
plt.savefig("lowacc")
plt.clf()

data_obs = totalxsec_clos + lowacc
data_obs_err = np.sqrt(data_obs)

with h5py.File('templates.hdf5', mode="w") as f:
    dset_templ = f.create_dataset('template', templates_all.shape, dtype=dtype)
    dset_templ[...] = templates_all
    dset_sumw2 = f.create_dataset('template_sumw2', templatesw2_all.shape, dtype=dtype)
    dset_sumw2[...] = templatesw2_all
    dset_hel = f.create_dataset('helicity', helicities_all.shape, dtype=dtype)
    dset_hel[...] = helicities_all
    dset_lowacc = f.create_dataset('lowacc', lowacc.shape, dtype=dtype)
    dset_lowacc[...] = lowacc
    dset_data_obs = f.create_dataset('data_obs', data_obs.shape, dtype=dtype)
    dset_data_obs[...] = data_obs
    dset_lowacc_sumw2 = f.create_dataset('lowacc_sumw2', lowacc_err.shape, dtype=dtype)
    dset_lowacc_sumw2[...] = lowacc_err
    dset_data_obs_sumw2 = f.create_dataset('data_obs_sumw2', data_obs_err.shape, dtype=dtype)
    dset_data_obs_sumw2[...] = data_obs_err
    
# #Build the template building sequenc
# def dySelectionSequence(p, xsec, systType, sumwClipped, nodetoStart, era):
#     print(ptBins)
#     print(zmassBins)
#     luminosityN = lumi_total2016
#     if era == 'preVFP' :     luminosityN = lumi_preVFP
#     else: luminosityN = lumi_postVFP

#     if systType not in [0,1]:#for signal MC
#         p.Histogram(columns = ["CStheta_preFSR","CSphi_preFSR","Vpt_preFSR","Vrap_preFSR","Vmass_preFSR"], types = ['float']*5,node=nodetoStart,histoname=ROOT.string('genhistos'),bins = [cosThetaBins,phiBins,qtBins,etaBins,zmassBins], variations = [])
    
#     p.branch(nodeToStart=nodetoStart, nodeToEnd='defs', modules=[ROOT.zVetoMuons()])
#     p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Sum(vetoMuons)==2 && Sum(goodMuons)==2", filtername="{:20s}".format("two muons"))
#     p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(HLT_IsoMu24 ||  HLT_IsoTkMu24)", filtername="{:20s}".format("Pass HLT"))
#     p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="(Muon_charge[goodMuons][0] + Muon_charge[goodMuons][1]) == 0", filtername="{:20s}".format("Opposite charge"))
#     p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="std::abs(Muon_eta[goodMuons][0]) < 2.4 && std::abs(Muon_eta[goodMuons][1]) < 2.4", filtername="{:20s}".format("Accept"))
#     p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_mediumId[goodMuons][0] == 1 && Muon_mediumId[goodMuons][1] == 1", filtername="{:20s}".format("MuonId"))
#     p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Muon_pfRelIso04_all[goodMuons][0] < 0.15 && Muon_pfRelIso04_all[goodMuons][1] < 0.15", filtername="{:20s}".format("Isolation"))
#     p.branch(nodeToStart='defs', nodeToEnd='defs', modules=[ROOT.zSelection()])
#     p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="60. < dimuonMass && dimuonMass < 120.", filtername="{:20s}".format("mZ range"))
#     p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("+ve mu trig matched"))

#     if systType == 0: #this is data
#         p.Histogram(columns = ["dimuonMass", "Mu1_eta", "Mu1_pt", "Mu2_eta", "Mu2_pt"], types = ['float']*5,node='defs',histoname=ROOT.string('data_muons'),bins = [zmassBins, etaBins, ptBins, etaBins, ptBins], variations = [])
#         #p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu2_pt < 200. ", filtername="{:20s}".format("mu2 pt upper acceptance"))
#         p.Histogram(columns = ["dimuonMass", "dimuonPt", "dimuonY", "MET_pt"], types = ['float']*4,node='defs',histoname=ROOT.string('data_dimuon'),bins = [zmassBins,qtBins, etaBins, metBins], variations = [])
#         return p
#     elif systType == 1:
#         print("Sample will be normalized to {}/fb".format(luminosityN))
#         p.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.SF_ul(fileSFul, fileSFPogTrk, isZ=True,era=era)])
#         #p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu2_pt < 200. ", filtername="{:20s}".format("mu2 pt upper acceptance"))
#         p.Histogram(columns = ["dimuonMass", "dimuonPt", "dimuonY", "MET_pt", "lumiweight", "puWeight", "SF"], types = ['float']*7,node='defs',histoname=ROOT.string('DY_dimuon'),bins = [zmassBins,qtBins, etaBins,metBins], variations = [])
#         return p
#     else:
#         print("Sample will be normalized to {}/fb".format(luminosityN))
#         p.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.SF_ul(fileSFul, fileSFPogTrk, isZ=True, era=era)])
#         p.Histogram(columns = ["dimuonMass", "Mu1_eta", "Mu1_pt", "Mu2_eta", "Mu2_pt", "lumiweight", "puWeight", "SF"], types = ['float']*8,node='defs',histoname=ROOT.string('DY_muons'),bins = [zmassBins, etaBins, ptBins, etaBins, ptBins], variations = [])
#         #p.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu2_pt < 200. ", filtername="{:20s}".format("mu2 pt upper acceptance"))
#         p.Histogram(columns = ["dimuonMass", "dimuonPt", "dimuonY", "MET_pt", "lumiweight", "puWeight", "SF"], types = ['float']*7,node='defs',histoname=ROOT.string('DY_dimuon'),bins = [zmassBins,qtBins, etaBins,metBins], variations = [])
#         return p
