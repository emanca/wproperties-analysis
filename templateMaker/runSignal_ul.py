import os
import sys
import ROOT
import copy
import time
import h5py
from math import pi, sqrt
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from array import array
sys.path.append('../RDFprocessor/framework')
sys.path.append('../Common/data')
from RDFtree import RDFtree
sys.path.append('python/')
from getLumiWeight import getLumiWeight
from binning import ptBins, etaBins, mTBins, etaBins, isoBins, chargeBins, yBins, qtBins, cosThetaBins, phiBins
from externals import fileSFul

# matplotlib stuff
plt.style.use([hep.style.ROOT])
hep.cms.label(loc=0, year=2016, lumi=35.9, data=True)
hep.cms.text('Simulation')

ROOT.gSystem.Load('bin/libAnalysisOnData.so')

ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")

ROOT.ROOT.EnableImplicitMT(48)

outputDir = 'PLOTS'
inputFile = '/scratchnvme/wmass/NanoAOD2016-UL/postNanoDec2020/WplusJetsToMuNu_preVFP_addVars/merged/*.root'

luminosityN = 35.9
sumwclipped = 5895502305412.54

# p = RDFtree(outputDir = outputDir, inputFile = inputFile, outputFile="test.root", pretend=False)
# p.branch(nodeToStart='input', nodeToEnd='defs', modules=[ROOT.lumiWeight(xsec=11572.19, sumwclipped=sumwclipped, targetLumi = luminosityN), ROOT.customizeforUL(True,True), ROOT.genDefinitions()])
# p.Histogram(columns = ["Wrap_preFSR_abs","Vpt_preFSR","CStheta_preFSR","CSphi_preFSR","lumiweight"], types = ['float']*5,node='defs',histoname=ROOT.string("xsecs"), bins = [yBins,qtBins,cosThetaBins,phiBins])
# p.gethdf5Output()

fewk = h5py.File('PLOTS/test.hdf5', mode='r+')
hFullAcc = np.array(fewk['xsecs'][:].reshape((len(yBins)-1,len(qtBins)-1, len(cosThetaBins)-1, len(phiBins)-1),order='F'),order='C')

fig, ax1 = plt.subplots()
ax1.set_title("CS theta", fontsize=18)
ax1.set_ylabel('Events')
ax1.set_xlabel('CS $\theta$')
data = np.sum(hFullAcc,axis=(0,1,3))
hep.histplot([data],bins = cosThetaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1)
plt.savefig('CStheta.png')
plt.clf()

fig, ax1 = plt.subplots()
ax1.set_title("CS phi", fontsize=18)
ax1.set_ylabel('Events')
ax1.set_xlabel('CS $\phi$')
data = np.sum(hFullAcc,axis=(0,1,2))
hep.histplot([data],bins = phiBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1)
plt.savefig('CSphi.png')
plt.clf()

# pick the right harmonics and reweight
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
totalxsec = np.sum(hFullAcc,axis=(2,3))

fig, ax1 = plt.subplots()
ax1.set_title("total xsec", fontsize=18)
hep.hist2dplot(totalxsec,yBins,qtBins)
plt.savefig("total_xsec")
plt.clf()

factors = [(20./3., 1./10),(5.,0.),(20.,0.),(4.,0.),(4.,0.),(5.,0.),(5.,0.),(4.,0.)]
for i,hw in enumerate(wharmonics):
    htmp = np.einsum('ijkm,km->ij',hFullAcc,hw)/totalxsec
    htmp = factors[i][0]*(htmp+factors[i][1])
    hharmonics.append(htmp)
    fig, ax1 = plt.subplots()
    ax1.set_title("A{}".format(i), fontsize=18)
    hep.hist2dplot(htmp,yBins,qtBins)
    plt.savefig("A{}".format(i))
    plt.cla()

# # reprocess with acceptance cuts
# p2 = RDFtree(outputDir = outputDir, inputFile = inputFile, outputFile="templates.root", pretend=False)
# p2.EventFilter(nodeToStart='input', nodeToEnd='defs', evfilter="HLT_SingleMu24", filtername="{:20s}".format("Pass HLT"))
# p2.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(Muon_mediumId)", filtername="{:20s}".format("MuonID"))
# p2.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(Muon_pfRelIso04_all < 0.15)", filtername="{:20s}".format("both mu isolated"))
# p2.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_dxy) < 0.05)", filtername="{:20s}".format("dxy"))
# p2.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_dz) < 0.2)", filtername="{:20s}".format("dz"))
# p2.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(Muon_pt > 25.)", filtername="{:20s}".format("Muon pt cut"))
# p2.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="All(abs(Muon_eta) < 2.4)", filtername="{:20s}".format("Muon eta cut"))
# p2.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="nMuon== 1", filtername="{:20s}".format("one muon"))
# p2.branch(nodeToStart = 'defs', nodeToEnd = 'defs', modules = [ROOT.customizeforUL(True, True), ROOT.genDefinitions(),ROOT.recoDefinitions(True, False),ROOT.lumiWeight(xsec=11572.19, sumwclipped=sumwclipped, targetLumi = luminosityN), ROOT.SF_ul(fileSFul)])
# # p2.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_hasTriggerMatch", filtername="{:20s}".format("+ve mu trig matched"))
# p2.EventFilter(nodeToStart='defs', nodeToEnd='defs', evfilter="Mu1_pt < 65.", filtername="{:20s}".format("mu1 pt-eta acceptance"))
# p2.Histogram(columns = ["Wrap_preFSR_abs","Vpt_preFSR","CStheta_preFSR","CSphi_preFSR","Mu1_eta", "Mu1_pt","lumiweight"], types = ['float']*7,node='defs',histoname=ROOT.string("templates"), bins = [yBins,qtBins,cosThetaBins,phiBins,etaBins,ptBins])
# p2.gethdf5Output()

fewk = h5py.File('PLOTS/templates.hdf5', mode='r+')
h = np.array(fewk['templates'][:].reshape((len(yBins)-1,len(qtBins)-1, len(cosThetaBins)-1, len(phiBins)-1, len(etaBins)-1, len(ptBins)-1),order='F'),order='C')
herr = np.array(fewk['templates_sumw2'][:].reshape((len(yBins)-1,len(qtBins)-1, len(cosThetaBins)-1, len(phiBins)-1, len(etaBins)-1, len(ptBins)-1),order='F'),order='C')

# now derive the weights for the templates
totalxsec = hFullAcc # total differential xsec in y,qt,costheta,phi with no acceptance cuts
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
    fig, ax1 = plt.subplots()
    ax1.set_title("templates{}".format(i), fontsize=18)
    hep.hist2dplot(htmp[2,2,:,:],etaBins,ptBins)
    plt.savefig("templates{}".format(i))
    plt.clf()

# check closure of templates sum

# totalxsec_clos = np.sum(hhelicity[0]+hhelicity[1]+hhelicity[2]+hhelicity[3]+hhelicity[4]+hhelicity[5]+hhelicity[6]+hhelicity[7]+hhelicity[8], axis=(0,1))
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

hharmonics[0]=np.squeeze(hharmonics[0])*np.sum(totalxsec,axis=(2,3))
hharmonics[1]=np.squeeze(hharmonics[1])*np.sum(totalxsec,axis=(2,3))
hharmonics[2]=np.squeeze(hharmonics[2])*np.sum(totalxsec,axis=(2,3))
hharmonics[3]=np.squeeze(hharmonics[3])*np.sum(totalxsec,axis=(2,3))
hharmonics[4]=np.squeeze(hharmonics[4])*np.sum(totalxsec,axis=(2,3))
hharmonics[5]=np.squeeze(hharmonics[5])*np.sum(totalxsec,axis=(2,3))
hharmonics[6]=np.squeeze(hharmonics[6])*np.sum(totalxsec,axis=(2,3))
hharmonics[7]=np.squeeze(hharmonics[7])*np.sum(totalxsec,axis=(2,3))
hharmonics.append(np.sum(totalxsec,axis=(2,3)))
print(hharmonics[0].shape,hharmonics[8].shape, np.sum(totalxsec,axis=(2,3)).shape)

threshold_y = np.digitize(2.4,yBins)
threshold_qt = np.digitize(60.,qtBins)
print(threshold_y,threshold_qt)

print(np.stack(hhelicity,axis=0)[:threshold_y,:threshold_qt,...].shape)
templates_all = np.stack(hhelicity,axis=0)[:,:threshold_y,:threshold_qt,...]
templatesw2_all = np.stack(herrhelicity,axis=0)[:,:threshold_y,:threshold_qt,...]
helicities_all = 3./(16.*pi)*np.stack(hharmonics,axis=0)[:,:threshold_y,:threshold_qt,...]
lowacc = np.sum(h[threshold_y:,threshold_qt:,...], axis=(0,1,2,3))
lowacc_err = np.sum(herr[threshold_y:,threshold_qt:,...], axis=(0,1,2,3))

print(templates_all.shape, 'templates')
fig, ax1 = plt.subplots()
ax1.set_title("low acceptance template", fontsize=18)
hep.hist2dplot(lowacc,etaBins,ptBins)
plt.savefig("lowacc")
plt.clf()

print(templates_all.shape)

with h5py.File('templates.hdf5', mode="w") as f:
    dset_templ = f.create_dataset('templates', templates_all.shape, dtype=dtype)
    dset_templ[...] = templates_all
    dset_sumw2 = f.create_dataset('templates_sumw2', templatesw2_all.shape, dtype=dtype)
    dset_sumw2[...] = templatesw2_all
    dset_hel = f.create_dataset('helicity', helicities_all.shape, dtype=dtype)
    dset_hel[...] = helicities_all
    dset_lowacc = f.create_dataset('lowacc', lowacc.shape, dtype=dtype)
    dset_lowacc[...] = lowacc
    dset_lowacc_sumw2 = f.create_dataset('lowacc_sumw2', lowacc_err.shape, dtype=dtype)
    dset_lowacc_sumw2[...] = lowacc_err

# # rescale to make them positive defined

# hhelicity[0]+=hhelicity[-1]/6.
# hhelicity[1]+=hhelicity[-1]/6.
# hhelicity[2]+=hhelicity[-1]/6.
# hhelicity[3]+=hhelicity[-1]/6.
# hhelicity[4]+=hhelicity[-1]/6.
# hhelicity[-1]/=6.

# for i,h in enumerate(hhelicity):
#     # plot one slice of every templates
#     fig, ax1 = plt.subplots()
#     ax1.set_title("templates{}_rescaled".format(i), fontsize=18)
#     hep.hist2dplot(h[0,0,:,:-1],etaBins,ptBins[:-1])
#     plt.savefig("templates{}_rescaled".format(i))
#     plt.clf()
