import ROOT
import os
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import mplhep as hep
sys.path.append('data/')
from binning import ptBins, etaBins, mTBins, isoBins, chargeBins, yBins, qtBins

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
yBins = np.array(yBins[:threshold_y+1])
qtBins = np.array(qtBins[:threshold_qt+1])
yBinsC = 0.5*(yBins[1:]+yBins[:-1])
qtBinsC = 0.5*(qtBins[1:]+qtBins[:-1])
etaBins = np.array(etaBins)
ptBins = np.array(ptBins)
etaBinsC = 0.5*(etaBins[1:]+etaBins[:-1])
ptBinsC = 0.5*(ptBins[1:]+ptBins[:-1])

helXsecs = ['L', 'I', 'T', 'A', 'P','7','8', '9','UL']

processes = []
for hel in helXsecs:
    for i in range(len(yBinsC)):
        for j in range(len(qtBinsC)):
                proc = 'helXsecs' + hel + '_y_{}'.format(i)+'_qt_{}'.format(j)
                processes.append(proc)

plt.style.use([hep.style.ROOT])
# hep.cms.label(loc=0, year=2016, lumi=35.9, data=True)

fIn = ROOT.TFile.Open('../Fit/fit_Wplus.root')
data = fIn.Get('obs')
hdata = np.array(data)[1:-1].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))

types = ['prefit', 'postfit']
for type in types:
    hsignal = np.zeros((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
    for proc in processes:
        tmp=fIn.Get('expproc_{}_{}'.format(proc,type))
        # print(type,np.array(tmp)[1:-2].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1)))
        hsignal+= np.array(tmp)[1:-2].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
    fakesLowMt = fIn.Get('expproc_fakesLowMt_{}'.format(type))
    hfakesLowMt = np.array(fakesLowMt)[1:-2].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
    fakesHighMt = fIn.Get('expproc_fakesHighMt_{}'.format(type))
    hfakesHighMt = np.array(fakesHighMt)[1:-2].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
    DY = fIn.Get('expproc_DY_{}'.format(type))
    hDY = np.array(DY)[1:-2].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
    Diboson = fIn.Get('expproc_Diboson_{}'.format(type))
    hDiboson = np.array(Diboson)[1:-2].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
    Top = fIn.Get('expproc_Top_{}'.format(type))
    hTop = np.array(Top)[1:-2].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
    Wtau = fIn.Get('expproc_Wtau_{}'.format(type))
    hWtau = np.array(Wtau)[1:-2].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
    LowAcc = fIn.Get('expproc_LowAcc_{}'.format(type))
    hLowAcc = np.array(LowAcc)[1:-2].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))

    hewk = hsignal+hDY+hDiboson+hTop+hWtau+hLowAcc
    for i in range(2):
        fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
        ax1.set_title("eta_iso{}_highMt".format(i), fontsize=18)
        ax1.set_ylabel('number of events')
        ax2.set_ylabel('data/prediction')
        ax2.set_xlabel('muon $\eta$')
        etadata = np.sum(hdata,axis=1)[:,-1,i]
        etaewk = np.sum(hewk,axis=1)[:,-1,i]
        etaWmu = np.sum(hsignal,axis=1)[:,-1,i]
        etaWtau = np.sum(hWtau,axis=1)[:,-1,i]
        etaDY = np.sum(hDY,axis=1)[:,-1,i]
        etaTop = np.sum(hTop,axis=1)[:,-1,i]
        etaDiboson = np.sum(hDiboson,axis=1)[:,-1,i]
        etafake = np.sum(hfakesHighMt,axis=1)[:,-1,i]
        etaLowAcc = np.sum(hLowAcc,axis=1)[:,-1,i]
        hep.histplot([etadata],bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1, label = ["data"])
        hep.histplot([etaDiboson,etaTop,etaDY,etaWtau,etafake,etaWmu,etaLowAcc],bins = etaBins, histtype = 'fill',linestyle = 'solid', color =["grey","magenta","orange","green","blue","red", "aqua"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes",r'$W->\mu\nu$', "low acc"], stack = True, ax=ax1)
        ax2.set_ylim([0.9, 1.1])
        hep.histplot([etadata/(etafake+etaewk)], bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
        ax1.legend(loc='upper right', frameon=True)
        plt.tight_layout()
        plt.savefig('eta_iso{}_highMt_{}.png'.format(i,type))
        plt.cla()

        fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
        ax1.set_title("pt_iso{}_highMt".format(i), fontsize=18)
        ax1.set_ylabel('number of events')
        ax2.set_ylabel('data/prediction')
        ax2.set_xlabel('muon $p_T$')
        ptdata = np.sum(hdata,axis=0)[:,-1,i]
        ptewk = np.sum(hewk,axis=0)[:,-1,i]
        ptWmu = np.sum(hsignal,axis=0)[:,-1,i]
        ptWtau = np.sum(hWtau,axis=0)[:,-1,i]
        ptDY = np.sum(hDY,axis=0)[:,-1,i]
        ptTop = np.sum(hTop,axis=0)[:,-1,i]
        ptDiboson = np.sum(hDiboson,axis=0)[:,-1,i]
        ptfake = np.sum(hfakesHighMt,axis=0)[:,-1,i]
        ptLowAcc = np.sum(hLowAcc,axis=0)[:,-1,i]
        hep.histplot([ptdata],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
        hep.histplot([ptDiboson,ptTop,ptDY,ptWtau,ptfake,ptWmu,ptLowAcc],bins = ptBins, histtype = 'fill',linestyle = 'solid', color =["grey","magenta","orange","green","blue","red", "aqua"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes",r'$W->\mu\nu$', "low acc"], stack = True, ax=ax1)
        ax2.set_ylim([0.9, 1.1])
        hep.histplot([ptdata/(ptfake+ptewk)],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
        ax1.legend(loc='upper right', frameon=True)
        plt.tight_layout()
        plt.savefig('pt_iso{}_highMt_{}.png'.format(i,type))
        plt.cla()

        fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
        ax1.set_title("eta_iso{}_lowMt".format(i), fontsize=18)
        ax1.set_ylabel('number of events')
        ax2.set_ylabel('data/prediction')
        ax2.set_xlabel('muon $\eta$')
        etadata = np.sum(hdata,axis=1)[:,0,i]
        etaewk = np.sum(hewk,axis=1)[:,0,i]
        etaWmu = np.sum(hsignal,axis=1)[:,0,i]
        etaWtau = np.sum(hWtau,axis=1)[:,0,i]
        etaDY = np.sum(hDY,axis=1)[:,0,i]
        etaTop = np.sum(hTop,axis=1)[:,0,i]
        etaDiboson = np.sum(hDiboson,axis=1)[:,0,i]
        etafake = np.sum(hfakesLowMt,axis=1)[:,0,i]
        etaLowAcc = np.sum(hLowAcc,axis=1)[:,0,i]
        hep.histplot([etadata],bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
        hep.histplot([etaDiboson,etaTop,etaDY,etaWtau,etafake,etaWmu, etaLowAcc],bins = etaBins, histtype = 'fill',linestyle = 'solid', color =["grey","magenta","orange","green","blue","red", "aqua"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes",r'$W->\mu\nu$', "low acc"], stack = True, ax=ax1)
        ax2.set_ylim([0.9, 1.1])
        hep.histplot([etadata/(etafake+etaewk)],bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
        ax1.legend(loc='upper right', frameon=True)
        plt.tight_layout()
        plt.savefig('eta_iso{}_lowMt_{}.png'.format(i,type))
        plt.cla()

        fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
        ax1.set_title("pt_iso{}_lowMt".format(i), fontsize=18)
        ax1.set_ylabel('number of events')
        ax2.set_ylabel('data/prediction')
        ax2.set_xlabel('muon $p_T$')
        ptdata = np.sum(hdata,axis=0)[:,0,i]
        ptewk = np.sum(hewk,axis=0)[:,0,i]
        ptWmu = np.sum(hsignal,axis=0)[:,0,i]
        ptWtau = np.sum(hWtau,axis=0)[:,0,i]
        ptDY = np.sum(hDY,axis=0)[:,0,i]
        ptTop = np.sum(hTop,axis=0)[:,0,i]
        ptDiboson = np.sum(hDiboson,axis=0)[:,0,i]
        ptfake = np.sum(hfakesLowMt,axis=0)[:,0,i]
        ptLowAcc = np.sum(hLowAcc,axis=0)[:,0,i]
        hep.histplot([ptdata],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
        hep.histplot([ptDiboson,ptTop,ptDY,ptWtau,ptfake,ptWmu,ptLowAcc],bins = ptBins, histtype = 'fill',linestyle = 'solid', color =["grey","magenta","orange","green","blue","red", "aqua"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes",r'$W->\mu\nu$', "low acc"], stack = True, ax=ax1)
        ax2.set_ylim([0.9, 1.1])
        hep.histplot([ptdata/(ptfake+ptewk)],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
        ax1.legend(loc='upper right', frameon=True)
        plt.tight_layout()
        plt.savefig('pt_iso{}_lowMt_{}.png'.format(i,type))
        plt.cla()