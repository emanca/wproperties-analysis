import ROOT
import os
import sys
import h5py
import numpy as np
from root_numpy import hist2array
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

fIn = ROOT.TFile.Open('../Fit/FitRes/fit_WPlus_data.root')

eras = ["preVFP", "postVFP"]
for iera,era in enumerate(eras):

    data = fIn.Get('obs')
    hdata = hist2array(data)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))

    types = ['prefit', 'postfit']
    for itype,type in enumerate(types):
        hsignal = np.zeros((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
        err2 = np.zeros((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1)).flatten()
        for proc in processes:
            tmp=fIn.Get('expproc_{}_{}'.format(proc,type))
            # print(type,np.array(tmp)[1:-2].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1)))
            hsignal+= hist2array(tmp)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
        tmp=fIn.Get('expfull_{}'.format(type))
        # err2+=np.array([tmp.GetSumw2()[i] for i in range(tmp.GetSumw2().GetSize())])[iera*11520:(iera+1)*11520]
        fakesLowMt = fIn.Get('expproc_fakesLowMt_{}'.format(type))
        hfakesLowMt = hist2array(fakesLowMt)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)  -1))
        fakesHighMt = fIn.Get('expproc_fakesHighMt_{}'.format(type))
        hfakesHighMt = hist2array(fakesHighMt)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len (isoBins)-1))
        # err2+=hist2array([fakesLowMt.GetSumw2()[i] for i in range(fakesLowMt.GetSumw2().GetSize())])[iera*11520:(iera+1)*11520]
        # err2+=hist2array([fakesHighMt.GetSumw2()[i] for i in range(fakesHighMt.GetSumw2().GetSize())])[iera*11520:(iera+1)*11520]
        DY = fIn.Get('expproc_DY_{}'.format(type))
        hDY = hist2array(DY)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
        Diboson = fIn.Get('expproc_Diboson_{}'.format(type))
        hDiboson = hist2array(Diboson)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
        Top = fIn.Get('expproc_Top_{}'.format(type))
        hTop = hist2array(Top)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
        Wtau = fIn.Get('expproc_Wtau_{}'.format(type))
        hWtau = hist2array(Wtau)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
        LowAcc = fIn.Get('expproc_LowAcc_{}'.format(type))
        hLowAcc = hist2array(LowAcc)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
        # err2+=hist2array([LowAcc.GetSumw2()[i] for i in range(LowAcc.GetSumw2().GetSize())])[iera*11520:(iera+1)*11520]
        err2 = err2.reshape(hsignal.shape)
        print(hist2array(tmp).shape)
        # hewk = hist2array(tmp)[iera*11520:(iera+1)*11520].reshape((len(etaBins)-1,len(ptBins)-1,len(mTBins)-1,len(isoBins)-1))
        hewk = hsignal+hDY+hDiboson+hTop+hWtau+hLowAcc

        # # differential plots
        # fig, (ax1, ax2) = plt.subplots(nrows=2,figsize=(48, 10),gridspec_kw={'height_ratios': [3, 1]})
        # hep.cms.text('work in progress', loc=0, ax=ax1)
        # ax1.set_ylabel('number of events')
        # ax2.set_ylabel('data/prediction')
        # data = hdata[...,itype]
        # Bins = np.linspace(0.,data.ravel().shape[0]+1, data.ravel().shape[0]+1)
        # print(data.shape, Bins.shape)
        # binsC = 0.5*(Bins[1:]+Bins[:-1])
        # ewk = hewk[...,itype]
        # Wmu = hewk[...,itype]
        # hep.histplot([data],bins = Bins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
        # hep.histplot([Wmu],bins = Bins, histtype = 'fill',linestyle = 'solid', color =["red"], label=[r'$W->\mu\nu$'],    stack = True, ax=ax1)
        # ax2.set_ylim([0.9, 1.1])
        # hep.histplot([data/(ewk)],bins = Bins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
        # ax2.fill_between(binsC, ((data/(ewk))-np.sqrt(err2[...,-1,i])*data/np.square(ewk)).ravel(), ((data/(ewk))+np.sqrt (err2[...,-1,i])*data/np.square(ewk)).ravel())
        # ax1.legend(loc='upper right', frameon=True)
        # plt.tight_layout()
        # plt.savefig('prefitplots/testreweighting_{}_{}.png'.format(type)),era
        # plt.cla()

        for i in range(2):
             # integrated plots
            fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
            hep.cms.text('work in progress', loc=0, ax=ax1)
            # ax1.set_title("eta_iso{}_highMt".format(i), fontsize=18)
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
            hep.histplot([etadata],bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1, label =   ["data"])
            hep.histplot([etaDiboson,etaTop,etaDY,etaWtau,etafake,etaLowAcc,etaWmu],bins = etaBins, histtype = 'fill',  linestyle = 'solid', color =["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top",  "DY",r'$W->\tau\nu$',"fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
            ax2.set_ylim([0.9, 1.1])
            hep.histplot([etadata/(etafake+etaewk)], bins = etaBins, histtype = 'errorbar', color = "k", stack = False,     ax=ax2)
            ax1.legend(loc='upper right', frameon=True)
            plt.tight_layout()
            plt.savefig('prefitplots/eta_iso{}_highMt_{}_{}.png'.format(i,type,era))
            plt.cla()

            fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
            hep.cms.text('work in progress', loc=0, ax=ax1)
            # ax1.set_title("pt_iso{}_highMt".format(i), fontsize=18)
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
            hep.histplot([ptdata],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"] )
            hep.histplot([ptDiboson,ptTop,ptDY,ptWtau,ptfake,ptLowAcc,ptWmu],bins = ptBins, histtype = 'fill',linestyle =   'solid', color =["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top","DY", r'$W->\tau\nu$',"fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
            ax2.set_ylim([0.9, 1.1])
            hep.histplot([ptdata/(ptfake+ptewk)],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
            ax1.legend(loc='upper right', frameon=True)
            plt.tight_layout()
            plt.savefig('prefitplots/pt_iso{}_highMt_{}_{}.png'.format(i,type,era))
            plt.cla()
            
            fig, ax1 = plt.subplots()
            hep.cms.text('work in progress', loc=0, ax=ax1)
            ax1.set_ylabel('number of events')
            ax1.set_xlabel('muon $p_T$')
            hep.histplot([ptfake],bins = ptBins, histtype = 'fill',linestyle =   'solid', color =["blue"], label=["fakes $W^+$"], ax=ax1)
            ax1.legend(loc='upper right', frameon=True)
            plt.tight_layout()
            plt.savefig('prefitplots/fakes_pt_iso{}_highMt_{}_{}.png'.format(i,type,era))
            plt.cla()

            fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
            hep.cms.text('work in progress', loc=0, ax=ax1)
            # ax1.set_title("eta_iso{}_lowMt".format(i), fontsize=18)
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
            hep.histplot([etadata],bins = etaBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label =    ["data"])
            hep.histplot([etaDiboson,etaTop,etaDY,etaWtau,etafake,etaWmu, etaLowAcc],bins = etaBins, histtype = 'fill', linestyle = 'solid', color =["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top", "DY",r'$W->\tau\nu$',"fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
            ax2.set_ylim([0.9, 1.1])
            hep.histplot([etadata/(etafake+etaewk)],bins = etaBins, histtype = 'errorbar', color = "k", stack = False,  ax=ax2)
            ax1.legend(loc='upper right', frameon=True)
            plt.tight_layout()
            plt.savefig('prefitplots/eta_iso{}_lowMt_{}_{}.png'.format(i,type,era))
            plt.cla()

            fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
            hep.cms.text('work in progress', loc=0, ax=ax1)
            # ax1.set_title("pt_iso{}_lowMt".format(i), fontsize=18)
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
            hep.histplot([ptdata],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"] )
            hep.histplot([ptDiboson,ptTop,ptDY,ptWtau,ptfake,ptLowAcc,ptWmu],bins = ptBins, histtype = 'fill',linestyle =   'solid', color =["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top","DY", r'$W->\tau\nu$',"fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
            ax2.set_ylim([0.9, 1.1])
            hep.histplot([ptdata/(ptfake+ptewk)],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
            ax1.legend(loc='upper right', frameon=True)
            plt.tight_layout()
            plt.savefig('prefitplots/pt_iso{}_lowMt_{}_{}.png'.format(i,type,era))
            plt.cla()

            # differential plots
            fig, (ax1, ax2) = plt.subplots(nrows=2,figsize=(48, 10),gridspec_kw={'height_ratios': [3, 1]})
            hep.cms.text('work in progress', loc=0, ax=ax1)
            # ax1.set_title("iso{}_highMt".format(i), fontsize=18)
            ax1.set_ylabel('number of events')
            ax2.set_ylabel('data/prediction')
            ax2.set_xlabel('unrolled $\eta-p_T$ bins')
            data = hdata[...,-1,i]
            Bins = np.linspace(0.,data.ravel().shape[0]+1, data.ravel().shape[0]+1)
            binsC = 0.5*(Bins[1:]+Bins[:-1])
            ewk = hewk[...,-1,i]
            Wmu = hsignal[...,-1,i]
            Wtau = hWtau[...,-1,i]
            DY = hDY[...,-1,i]
            Top = hTop[...,-1,i]
            Diboson = hDiboson[...,-1,i]
            fake = hfakesHighMt[...,-1,i]
            LowAcc = hLowAcc[...,-1,i]
            hep.histplot([data],bins = Bins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
            hep.histplot([Diboson,Top,DY,Wtau,fake,LowAcc,Wmu],bins = Bins, histtype = 'fill',linestyle = 'solid', color =  ["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes", "low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
            ax2.set_ylim([0.9, 1.1])
            hep.histplot([data/(fake+ewk)],bins = Bins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
            ax2.fill_between(binsC, ((data/(fake+ewk))-np.sqrt(err2[...,-1,i])*data/np.square(fake+ewk)).ravel(), ((data/ (fake+ewk))+np.sqrt(err2[...,-1,i])*data/np.square(fake+ewk)).ravel())
            ax1.legend(loc='upper right', frameon=True)
            plt.tight_layout()
            plt.savefig('prefitplots/iso{}_highMt_{}_{}.png'.format(i,type,era))
            plt.cla()

            fig, ax1 = plt.subplots()

            hep.histplot(np.histogram([(data/(fake+ewk)).ravel()], bins=np.linspace(0.9,1.1,100)),bins = np.linspace(0.9,1.1,100), color = "b", stack = False, ax=ax1)
            plt.tight_layout()
            plt.savefig('prefitplots/iso{}_highMt_ratio_{}_{}.png'.format(i,type,era))
            plt.cla()

            # differential plots
            fig, (ax1, ax2) = plt.subplots(nrows=2,figsize=(48, 10),gridspec_kw={'height_ratios': [3, 1]})
            hep.cms.text('work in progress', loc=0, ax=ax1)
            # ax1.set_title("iso{}_lowMt".format(i), fontsize=18)
            ax1.set_ylabel('number of events')
            ax2.set_ylabel('data/prediction')
            ax2.set_xlabel('unrolled $\eta-p_T$ bins')
            # ax2.set_xlabel('muon $p_T$')
            data = hdata[...,0,i]
            Bins = np.linspace(0.,data.ravel().shape[0]+1, data.ravel().shape[0]+1)
            binsC = 0.5*(Bins[1:]+Bins[:-1])
            # bins= np.tile(etaBinsC, len(ptBinsC))
            # print(Bins,bins)
            # locator = mdates.HourLocator(interval=1)
            # locator.MAXTICKS = 100000
            # ax2.xaxis.set_minor_locator(locator)
            # ax2.set_xticklabels(bins)
            ewk = hewk[...,0,i]
            Wmu = hsignal[...,0,i]
            Wtau = hWtau[...,0,i]
            DY = hDY[...,0,i]
            Top = hTop[...,0,i]
            Diboson = hDiboson[...,0,i]
            fake = hfakesLowMt[...,0,i]
            LowAcc = hLowAcc[...,0,i]
            hep.histplot([data],bins = Bins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label = ["data"])
            hep.histplot([Diboson,Top,DY,Wtau,fake,LowAcc,Wmu],bins = Bins, histtype = 'fill',linestyle = 'solid', color =  ["grey","magenta","orange","green","blue","aqua","red"], label=["Diboson","Top","DY",r'$W->\tau\nu$',"fakes", "low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
            ax2.set_ylim([0.9, 1.1])
            hep.histplot([data/(fake+ewk)],bins = Bins, histtype = 'errorbar', color = "k", stack = False, ax=ax2)
            ax2.fill_between(binsC, ((data/(fake+ewk))-np.sqrt(err2[...,0,i])*data/np.square(fake+ewk)).ravel(), ((data/    (fake+ewk))+np.sqrt(err2[...,0,i])*data/np.square(fake+ewk)).ravel())
            ax1.legend(loc='upper right', frameon=True)
            plt.tight_layout()
            plt.savefig('prefitplots/iso{}_lowMt_{}_{}.png'.format(i,type,era))
            plt.cla()

            fig, ax1 = plt.subplots()
            hep.histplot(np.histogram([(data/(fake+ewk)).ravel()], bins=np.linspace(0.9,1.1,100)),bins = np.linspace(0.9,1.1,100), color = "b", stack = False, ax=ax1)
            plt.tight_layout()
            plt.savefig('prefitplots/iso{}_lowMt_ratio_{}_{}.png'.format(i,type,era))
            plt.cla()

            # # fake plot
            # fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
            # hep.cms.text('work in progress', loc=0, ax=ax1)
            # # ax1.set_title("pt_iso{}_highMt".format(i), fontsize=18)
            # ax1.set_ylabel('number of events')
            # ax2.set_ylabel('data/prediction')
            # ax2.set_xlabel('muon $p_T$')
            # ptdata = hdata[24,:,-1,i]
            # ptewk = hewk[24,:,-1,i]
            # ptWmu = hsignal[24,:,-1,i]
            # ptWtau = hWtau[24,:,-1,i]
            # ptDY = hDY[24,:,-1,i]
            # ptTop = hTop[24,:,-1,i]
            # ptDiboson = hDiboson[24,:,-1,i]
            # ptfake = hfakesHighMt[24,:,-1,i]
            # ptLowAcc = hLowAcc[24,:,-1,i]
            # ax1.set_yscale('log')
            # ax1.set_ylim([1.1e3, 5.e5])
            # # hep.histplot([ptdata],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label =  ["data"])
            # hep.histplot([ptfake,ptLowAcc,ptWmu],bins = ptBins, histtype = 'fill',linestyle = 'solid', color =["blue",    "aqua","red"], label=["fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
            # ax2.set_ylim([0.9, 1.1])
            # hep.histplot([ptdata/(ptfake+ptewk)],bins = ptBins, histtype = 'errorbar', color = "k", stack = False,    ax=ax2)
            # ax1.legend(loc='upper right', frameon=True)
            # plt.tight_layout()
            # plt.savefig('prefitplots/pt_iso{}_highMt_{}_onebin_{}.png'.format(i,type,era))
            # plt.cla()

            # fig, (ax1, ax2) = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]})
            # hep.cms.text('work in progress', loc=0, ax=ax1)
            # # ax1.set_title("pt_iso{}_lowMt".format(i), fontsize=18)
            # ax1.set_ylabel('number of events')
            # ax2.set_ylabel('data/prediction')
            # ax2.set_xlabel('muon $p_T$')
            # ptdata = hdata[24,:,0,i]
            # ptewk = hewk[24,:,0,i]
            # ptWmu = hsignal[24,:,0,i]
            # ptWtau = hWtau[24,:,0,i]
            # ptDY = hDY[24,:,0,i]
            # ptTop = hTop[24,:,0,i]
            # ptDiboson = hDiboson[24,:,0,i]
            # ptfake = hfakesLowMt[24,:,0,i]
            # ptLowAcc = hLowAcc[24,:,0,i]
            # ax1.set_yscale('log')
            # ax1.set_ylim([1.1e3, 5.e5])
            # # hep.histplot([ptdata],bins = ptBins, histtype = 'errorbar', color = "k", stack = False, ax=ax1,label =  ["data"])
            # hep.histplot([ptfake,ptLowAcc,ptWmu],bins = ptBins, histtype = 'fill',linestyle = 'solid', color =["blue",    "aqua","red"], label=["fakes","low acc",r'$W->\mu\nu$'], stack = True, ax=ax1)
            # ax2.set_ylim([0.9, 1.1])
            # hep.histplot([ptdata/(ptfake+ptewk)],bins = ptBins, histtype = 'errorbar', color = "k", stack = False,    ax=ax2)
            # ax1.legend(loc='upper right', frameon=True)
            # plt.tight_layout()
            # plt.savefig('prefitplots/pt_iso{}_lowMt_{}_onebin_{}.png'.format(i,type,era))
            # plt.cla()