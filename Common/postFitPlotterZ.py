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

helXsecs = ['L', 'I', 'T', 'A', 'P','UL']

processes = []
for hel in helXsecs:
    for i in range(len(yBinsC)):
        for j in range(len(qtBinsC)):
                proc = 'helXsecs' + hel + '_y_{}'.format(i)+'_qt_{}'.format(j)
                processes.append(proc)

plt.style.use([hep.style.ROOT])
# hep.cms.label(loc=0, year=2016, lumi=35.9, data=True)

fIn = ROOT.TFile.Open('../Fit/FitRes/fit_Z_toy.root')

data = fIn.Get('obs')
hdata = np.array(data)[1:-1]#.reshape((len(etaBins)-1,len(ptBins)-1))
print(hdata)

types = ['prefit', 'postfit']
for itype,type in enumerate(types):
    hsignal = np.zeros((len(etaBins)-1,len(ptBins)-1)).ravel()
    err2 = np.zeros((len(etaBins)-1,len(ptBins)-1)).flatten()

    # tmp=fIn.Get('expproc_{}_{}'.format(proc,type))
    tmp=fIn.Get('expsig_{}'.format(type))
    hsignal+= np.array(tmp)[1:-2]#.reshape((len(etaBins)-1,len(ptBins)-1))
    # err2+=np.array([tmp.GetSumw2()[i] for i in range(tmp.GetSumw2().GetSize())])[iera*11520:(iera+1)*11520]
    LowAcc = fIn.Get('expproc_LowAcc_{}'.format(type))
    hLowAcc = np.array(LowAcc)[1:-2]#.reshape((len(etaBins)-1,len(ptBins)-1))
    err2 = err2.reshape(hsignal.shape)
    # hewk = hist2array(tmp)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
    hewk = hsignal+hLowAcc

    # differential plots
    fig, (ax1, ax2) = plt.subplots(nrows=2,figsize=(48, 10),gridspec_kw={'height_ratios': [3, 1]})
    hep.cms.text('work in progress', loc=0, ax=ax1)
    # ax1.set_title("iso{}_highMt".format(i), fontsize=18)
    ax1.set_ylabel('number of events')
    ax2.set_ylabel('data/prediction')
    ax2.set_xlabel('unrolled $\eta-p_T$ bins')
    data = hdata
    print(data)
    Bins = np.linspace(0.,data.ravel().shape[0]+1, data.ravel().shape[0]+1)
    print(data.shape,Bins.shape)
    binsC = 0.5*(Bins[1:]+Bins[:-1])
    Z = hsignal
    LowAcc = hLowAcc
    hep.histplot([data],bins = Bins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
    hep.histplot([LowAcc,Z],bins = Bins, histtype = 'fill',linestyle = 'solid', color =  ["aqua","red"], label=["low acc",r'$Z->\mu\mu$'], stack = True, ax=ax1)
    ax2.set_ylim([0.9, 1.1])
    hep.histplot([data/(hewk)],bins = Bins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
    # ax2.fill_between(binsC, ((data/(hewk))-np.sqrt(err2[...,-1,i])*data/np.square(hewk)).ravel(), ((data/ (hewk))+np.sqrt(err2[...,-1,i])*data/np.square(hewk)).ravel())
    ax1.legend(loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig('prefitplots/Z_{}.png'.format(type))
    plt.cla()

    # integrated plots
    fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
    hep.cms.text('work in progress', loc=0, ax=ax1)
    # ax1.set_title("eta_iso{}_highMt".format(i), fontsize=18)
    ax1.set_ylabel('number of events')
    ax2.set_ylabel('data/prediction')
    ax2.set_xlabel('muon $\eta$')
    etadata = np.sum(hdata.reshape((len(etaBins)-1,len(ptBins)-1)),axis=1)[...]
    etaewk = np.sum(hewk.reshape((len(etaBins)-1,len(ptBins)-1)),axis=1)[...]
    etaZ = np.sum(hsignal.reshape((len(etaBins)-1,len(ptBins)-1)),axis=1)[...]
    etaLowAcc = np.sum(hLowAcc.reshape((len(etaBins)-1,len(ptBins)-1)),axis=1)[...]
    hep.histplot([etadata],bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1, label =   ["data"])
    hep.histplot([etaLowAcc,etaZ],bins = etaBins, histtype = 'fill',linestyle = 'solid', color =  ["aqua","red"], label=["low acc",r'$Z->\mu\mu$'], stack = True, ax=ax1)
    ax2.set_ylim([0.9, 1.1])
    hep.histplot([etadata/(etaewk)], bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
    ax1.legend(loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig('prefitplots/eta_Z_{}.png'.format(type))
    plt.cla()