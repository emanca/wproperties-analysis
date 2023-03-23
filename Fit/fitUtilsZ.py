import ROOT
import pickle
from termcolor import colored
import math
from HiggsAnalysis.CombinedLimit.DatacardParser import *
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import mplhep as hep
from collections import OrderedDict
import copy
from systToapply import systematicsDict
import numpy as np
from root_numpy import hist2array
import sys
sys.path.append('../Common/data')
import h5py

class fitUtilsZ:
    def __init__(self, channels =["Z_postVFP"], doSyst=False):
        
        self.doSyst = doSyst
        self.processes = []
        self.signals = []

        #combine utils
        self.channels = channels
        self.shapeMap = {}
        self.helGroups = OrderedDict()
        self.sumGroups = OrderedDict()
        self.helMetaGroups = OrderedDict()
        
        self.templSystematics = systematicsDict

        # self.helXsecs = ['L', 'I', 'T', 'A', 'P','7','8', '9','UL']
        self.helXsecs = ['L', 'I', 'T', 'A', 'P','UL']
        
        # reduce bins to acceptance
        from binning import ptBins, etaBins, mTBins, etaBins, isoBins, chargeBins, yBins, qtBins, cosThetaBins, phiBins
        self.ptBins=ptBins
        self.etaBins=etaBins
        self.threshold_y = np.digitize(2.4,yBins)-1
        self.threshold_qt = np.digitize(60.,qtBins)-1
        self.yBins = np.array(yBins[:self.threshold_y+1])
        self.qtBins = np.array(qtBins[:self.threshold_qt+1])
        self.yBinsC = 0.5*(self.yBins[1:]+self.yBins[:-1])
        self.qtBinsC = 0.5*(self.qtBins[1:]+self.qtBins[:-1])
        print(len(self.yBinsC),len(self.qtBinsC))

        # load files
        self.ftempl = {}
        self.data = {}
        self.templ = {}
        self.templw2 = {}
        self.templ_NLO = {}
        self.templw2_NLO = {}
        self.gen = {}
        self.gen_NLO = {}
        self.lowacc = {}
        self.lowaccw2 = {}
        self.Wtau = {}
        self.Wtauw2 = {}
        self.DY = {}
        self.DYw2 = {}
        self.Top = {}
        self.Topw2 = {}
        self.Diboson = {}
        self.Dibosonw2 = {}
        self.fakeslow = {}
        self.fakesloww2 = {}
        self.fakeshigh = {}
        self.fakeshighw2 = {}
        self.templSFStat ={} 
        for chan in self.channels:
            self.ftempl[chan] = h5py.File('../templateMaker/templatesFit.hdf5', mode='r+')
            self.data[chan] = self.ftempl[chan]['data_obs'][:]
            print(chan,'events in data:', np.sum(self.data[chan][...]))
            self.templ[chan] = self.ftempl[chan]['template'][:]
            print(chan,'events in signal templ:', np.sum(self.templ[chan][...]),self.templ[chan][...].shape)
            self.templw2[chan] = self.ftempl[chan]['template_sumw2'][:]
            self.gen[chan] = self.ftempl[chan]['helicity'][:self.threshold_y,:self.threshold_qt,:]
            self.lowacc[chan] = self.ftempl[chan]['lowacc'][:]
            print(chan,'events in low acc templ:', np.sum(self.lowacc[chan][...]))
            self.lowaccw2[chan] = self.ftempl[chan]['lowacc_sumw2'][:]
    
    def fillProcessList(self):
        for hel in self.helXsecs:
            for i in range(len(self.yBinsC)):
                for j in range(len(self.qtBinsC)):
                    proc = 'helXsecs' + hel + '_y_{}'.format(i)+'_qt_{}'.format(j)
                    self.processes.append(proc)
                    #if "helXsecsUL"in proc or "helXsecsP" in proc:
                    if not "helXsecs7" in proc and not "helXsecs8" in proc and not "helXsecs9" in proc:
                        print(proc)
                        self.signals.append(proc)
        # bkg_list = ["DY","Diboson","Top","fakesLowMt","fakesHighMt", "Wtau","LowAcc"]
        bkg_list = ["LowAcc"]
        self.processes.extend(bkg_list)
    
    # def mirrorShape(self, nominal,alternate):
    #     # avoid changing nominal
    #     up = np.empty_like(nominal)
    #     np.copyto(up, nominal) # (dst, src)
    #     down = np.empty_like(nominal)
    #     np.copyto(down, nominal) # (dst, src)
    #     ratio = alternate/nominal
    #     up*=ratio
    #     down*=np.reciprocal(ratio)
    #     return up,down
    
    # def decorrelateSFSyst(self, nominal, rawvarsUp, rawvarsDown):
    #     ptBins_SF = np.array([25.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 47.0, 50.0, 55.0])
    #     # ptBins_SF = np.array([25.0, 55.0])
    #     vars=np.zeros((48,60,2,2,6,8,9,624,2))
    #     isyst=0
    #     for j in range(48): # for each eta bin
    #         for i in range(len(ptBins_SF)-1): #create one histogram per variation in each macro bin of SF
    #             threshold1 = np.digitize(ptBins_SF[i],self.ptBins)-1
    #             threshold2 = np.digitize(ptBins_SF[i+1],self.ptBins)-1
    #             # print(threshold1,threshold2,self.ptBins[threshold1:threshold2+1])
    #             tmpUp = np.zeros_like(nominal)
    #             np.copyto(tmpUp,nominal) #numpy.copyto(dst, src)
    #             tmpUp[j,threshold1:threshold2+1,...]=rawvarsUp[j,threshold1:threshold2+1,...]
    #             # tmpUp[...]=rawvarsUp[...]
    #             vars[...,isyst,0]=tmpUp
    #             tmpDown = np.zeros_like(nominal)
    #             np.copyto(tmpDown,nominal) #numpy.copyto(dst, src)
    #             tmpDown[j,threshold1:threshold2+1,...]=rawvarsDown[j,threshold1:threshold2+1,...]
    #             # tmpDown[...]=rawvarsDown[...]
    #             vars[...,isyst,1]=tmpDown
    #             isyst=isyst+1
    #     print('end SF decorrelation')
    #     return vars
    
    # def decorrelateSFSystBkg(self, nominal, rawvarsUp, rawvarsDown):
    #     ptBins_SF = np.array([25.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 47.0, 50.0, 55.0])
    #     # ptBins_SF = np.array([25.0, 55.0])
    #     vars=np.zeros((48,60,2,2,624,2))
    #     isyst=0
    #     for j in range(48): # for each eta bin
    #         for i in range(len(ptBins_SF)-1): #create one histogram per variation in each macro bin of SF
    #             threshold1 = np.digitize(ptBins_SF[i],self.ptBins)-1
    #             threshold2 = np.digitize(ptBins_SF[i+1],self.ptBins)-1
    #             # print(threshold1,threshold2,self.ptBins[threshold1:threshold2+1])
    #             tmpUp = np.zeros_like(nominal)
    #             np.copyto(tmpUp,nominal) #numpy.copyto(dst, src)
    #             tmpUp[j,threshold1:threshold2+1,...]=rawvarsUp[j,threshold1:threshold2+1,...]
    #             # tmpUp[...]=rawvarsUp[...]
    #             vars[...,isyst,0]=tmpUp
    #             tmpDown = np.zeros_like(nominal)
    #             np.copyto(tmpDown,nominal) #numpy.copyto(dst, src)
    #             tmpDown[j,threshold1:threshold2+1,...]=rawvarsDown[j,threshold1:threshold2+1,...]
    #             # tmpDown[...]=rawvarsDown[...]
    #             vars[...,isyst,1]=tmpDown
    #             isyst=isyst+1
    #     print('end SF decorrelation')
    #     return vars

    def shapeFile(self):
        dtype = 'float64'
        compression = "gzip"
        for chan in self.channels:
            templ_mass = self.ftempl[chan]['template_mass'][:]
            # templ_SFSyst = self.ftempl[chan]['template_SFSystvar'][...,:self.threshold_y,:self.threshold_qt,:,:]
            # templ_prefire = self.ftempl[chan]['template_prefireVars'][...,:self.threshold_y,:self.threshold_qt,:,:]
            # templ_jesUp = self.ftempl[chan]['template_jesTotalUp'][...,:self.threshold_y,:self.threshold_qt,:,:]
            # templ_jesDown = self.ftempl[chan]['template_jesTotalDown'][...,:self.threshold_y,:self.threshold_qt,:,:]
            # templ_unclUp = self.ftempl[chan]['template_unclustEnUp'][...,:self.threshold_y,:self.threshold_qt,:,:]
            # templ_unclDown = self.ftempl[chan]['template_unclustEnDown'][...,:self.threshold_y,:self.threshold_qt,:,:]
            # shape_chan_SFStat = self.decorrelateSFSyst(self.templ[chan], self.templSFStat[chan][...,0], self.templSFStat[chan][...,1])
            # shape_chan_SFStat_iso = self.decorrelateSFSyst(self.templ[chan], self.templSFStat[chan][...,2], self.templSFStat[chan][...,3])
            with h5py.File('{}.hdf5'.format(chan), mode="w") as f:
                for proc in self.processes:
                    if "helXsecs" in proc:
                        coeff = self.helXsecs.index(proc.split('_')[0].replace('helXsecs',''))
                        iY = int(proc.split('_')[2])
                        iQt = int(proc.split('_')[4])
                        # print(chan,proc,self.templ[chan].shape)

                        dset_templ = f.create_dataset(proc, self.templ[chan][iY,iQt,...,coeff].ravel().shape, dtype=dtype,compression=compression)
                        dset_templ[...] = self.templ[chan][iY,iQt,...,coeff].ravel()
                        dset_templw2 = f.create_dataset(proc+'_sumw2', self.templw2[chan][iY,iQt,...,coeff].ravel().shape, dtype=dtype,compression=compression)
                        dset_templw2[...] = self.templw2[chan][iY,iQt,...,coeff].ravel()
                        # mass
                        dset_templ = f.create_dataset(proc+'_massUp', templ_mass[iY,iQt,...,coeff,0].ravel().shape, dtype=dtype,compression=compression)
                        dset_templ[...] = templ_mass[iY,iQt,...,coeff,0].ravel()
                        dset_templ = f.create_dataset(proc+'_massDown', templ_mass[iY,iQt,...,coeff,1].ravel().shape, dtype=dtype,compression=compression)
                        dset_templ[...] = templ_mass[iY,iQt,...,coeff,1].ravel()
                        # # # SF
                        # # for iSF in range(shape_chan_SFStat.shape[-2]):
                        # #     dset_templ = f.create_dataset(proc+'_SFall{}Up'.format(iSF), shape_chan_SFStat[...,iY,iQt,coeff,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                        # #     if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        # #         dset_templ[...] = -1*shape_chan_SFStat[...,iY,iQt,coeff,iSF,0].ravel()
                        # #     else:
                        # #         dset_templ[...] = shape_chan_SFStat[...,iY,iQt,coeff,iSF,0].ravel()
                        # #     dset_templ = f.create_dataset(proc+'_SFall{}Down'.format(iSF), shape_chan_SFStat[...,iY,iQt,coeff,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                        # #     if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        # #         dset_templ[...] = -1*shape_chan_SFStat[...,iY,iQt,coeff,iSF,1].ravel()
                        # #     else:
                        # #         dset_templ[...] = shape_chan_SFStat[...,iY,iQt,coeff,iSF,1].ravel()
                            
                        # #     dset_templ = f.create_dataset(proc+'_SFiso{}Up'.format(iSF), shape_chan_SFStat_iso[...,iY,iQt,coeff,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                        # #     if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        # #         dset_templ[...] = -1*shape_chan_SFStat_iso[...,iY,iQt,coeff,iSF,0].ravel()
                        # #     else:
                        # #         dset_templ[...] = shape_chan_SFStat_iso[...,iY,iQt,coeff,iSF,0].ravel()
                        # #     dset_templ = f.create_dataset(proc+'_SFiso{}Down'.format(iSF), shape_chan_SFStat_iso[...,iY,iQt,coeff,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                        # #     if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        # #         dset_templ[...] = -1*shape_chan_SFStat_iso[...,iY,iQt,coeff,iSF,1].ravel()
                        # #     else:
                        # #         dset_templ[...] = shape_chan_SFStat_iso[...,iY,iQt,coeff,iSF,1].ravel()
                        # # print(proc, self.templSFStat[chan].shape)
                        # dset_templ = f.create_dataset(proc+'_SFallUp', self.templSFStat[chan][...,iY,iQt,coeff,0].ravel().shape, dtype=dtype,compression=compression)
                        # dset_templ[...] = self.templSFStat[chan][...,iY,iQt,coeff,0].ravel()
                        # dset_templ = f.create_dataset(proc+'_SFallDown', self.templSFStat[chan][...,iY,iQt,coeff,1].ravel().shape, dtype=dtype,compression=compression)
                        # dset_templ[...] = self.templSFStat[chan][...,iY,iQt,coeff,1].ravel()
                            
                        # dset_templ = f.create_dataset(proc+'_SFisoUp', self.templSFStat[chan][...,iY,iQt,coeff,2].ravel().shape, dtype=dtype,compression=compression)
                        # dset_templ[...] = self.templSFStat[chan][...,iY,iQt,coeff,2].ravel()
                        # dset_templ = f.create_dataset(proc+'_SFisoDown', self.templSFStat[chan][...,iY,iQt,coeff,3].ravel().shape, dtype=dtype,compression=compression)
                        # dset_templ[...] = self.templSFStat[chan][...,iY,iQt,coeff,3].ravel()
                        # # SFSyst
                        # nominal = self.templ[chan][...,iY,iQt,coeff][:].ravel()
                        # alternate = templ_SFSyst[...,iY,iQt,coeff].ravel()
                        # up,down = self.mirrorShape(nominal,alternate)
                        # dset_templ = f.create_dataset(proc+'_SFSystUp', alternate.ravel().shape, dtype=dtype,compression=compression)
                        # if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        #     dset_templ[...] = -1*up.ravel()
                        # else:
                        #     dset_templ[...] = up.ravel()
                        # dset_templ = f.create_dataset(proc+'_SFSystDown', down.shape, dtype=dtype,compression=compression)
                        # if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        #     dset_templ[...] = -1*down
                        # else:
                        #     dset_templ[...] = down
                        # # prefire
                        # dset_templ = f.create_dataset(proc+'_prefireUp', templ_prefire[...,iY,iQt,coeff,0].ravel().shape, dtype=dtype,compression=compression)
                        # if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        #     dset_templ[...] = -1*templ_prefire[...,iY,iQt,coeff,0].ravel()
                        # else:
                        #     dset_templ[...] = templ_prefire[...,iY,iQt,coeff,0].ravel()
                        # dset_templ = f.create_dataset(proc+'_prefireDown', templ_prefire[...,iY,iQt,coeff,1].ravel().shape, dtype=dtype,compression=compression)
                        # if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        #     dset_templ[...] = -1*templ_prefire[...,iY,iQt,coeff,1].ravel()
                        # else:
                        #     dset_templ[...] = templ_prefire[...,iY,iQt,coeff,1].ravel()
                        # # JEC
                        # dset_templ = f.create_dataset(proc+'_jesUp', templ_jesUp[...,iY,iQt,coeff].ravel().shape, dtype=dtype,compression=compression)
                        # if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        #     dset_templ[...] = -1*templ_jesUp[...,iY,iQt,coeff].ravel()
                        # else:
                        #     dset_templ[...] = templ_jesUp[...,iY,iQt,coeff].ravel()
                        # dset_templ = f.create_dataset(proc+'_jesDown', templ_jesDown[...,iY,iQt,coeff].ravel().shape, dtype=dtype,compression=compression)
                        # if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        #     dset_templ[...] = -1*templ_jesDown[...,iY,iQt,coeff].ravel()
                        # else:
                        #     dset_templ[...] = templ_jesDown[...,iY,iQt,coeff].ravel()
                        
                        # dset_templ = f.create_dataset(proc+'_unclUp', templ_unclUp[...,iY,iQt,coeff].ravel().shape, dtype=dtype,compression=compression)
                        # if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        #     dset_templ[...] = -1*templ_unclUp[...,iY,iQt,coeff].ravel()
                        # else:
                        #     dset_templ[...] = templ_unclUp[...,iY,iQt,coeff].ravel()
                        # dset_templ = f.create_dataset(proc+'_unclDown', templ_unclDown[...,iY,iQt,coeff].ravel().shape, dtype=dtype,compression=compression)
                        # if proc=='helXsecsA_y_2_qt_5' and 'WPlus_preVFP' in chan:
                        #     dset_templ[...] = -1*templ_unclDown[...,iY,iQt,coeff].ravel()
                        # else:
                        #     dset_templ[...] = templ_unclDown[...,iY,iQt,coeff].ravel()

                dset_data = f.create_dataset('data_obs', self.data[chan][...].ravel().shape, dtype=dtype,compression=compression)
                dset_data[...] = self.data[chan][...].ravel()

                # dset_bkg = f.create_dataset("DY", self.DY[chan][...].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.DY[chan][...].ravel()
                # dset_bkgw2 = f.create_dataset("DY_sumw2", self.DYw2[chan][...].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkgw2[...] = self.DYw2[chan][...].ravel()
                # # pdf variations
                # for k in range(1,103):
                #     nominal = self.DY[chan][...]
                #     alternate = self.ftempl[chan]['DY_LHEPdfWeight'][...,k]
                #     up,down = self.mirrorShape(nominal,alternate)
                #     dset_templ = f.create_dataset("DY_pdf{}Up".format(k), up.ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = up.ravel()
                #     dset_templ = f.create_dataset("DY_pdf{}Down".format(k), down.ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = down.ravel()
                # #QCD
                # for j,syst in self.qcdsyst.items():
                #     dset_bkg = f.create_dataset(name='DY_{}'.format(syst), shape=self.ftempl[chan]['DY_LHEScaleWeight'][...,j].ravel().shape, dtype='float64')
                #     dset_bkg[...] = self.ftempl[chan]['DY_LHEScaleWeight'][...,j].ravel()
                # #SF
                # shape_chan_SFStat = self.decorrelateSFSystBkg(self.DY[chan], self.ftempl[chan]['DY_SFStatvar'][...,0], self.ftempl[chan]['DY_SFStatvar'][...,1])
                # shape_chan_SFStat_iso = self.decorrelateSFSystBkg(self.DY[chan], self.ftempl[chan]['DY_SFStatvar'][...,2], self.ftempl[chan]['DY_SFStatvar'][...,3])
                # for iSF in range(shape_chan_SFStat.shape[-2]):
                #     dset_templ = f.create_dataset('DY_SFall{}Up'.format(iSF), shape_chan_SFStat[...,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat[...,iSF,0].ravel()
                #     dset_templ = f.create_dataset('DY_SFall{}Down'.format(iSF), shape_chan_SFStat[...,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat[...,iSF,1].ravel()
                                                
                #     dset_templ = f.create_dataset('DY_SFiso{}Up'.format(iSF), shape_chan_SFStat_iso[...,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat_iso[...,iSF,0].ravel()
                #     dset_templ = f.create_dataset('DY_SFiso{}Down'.format(iSF), shape_chan_SFStat_iso[...,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat_iso[...,iSF,1].ravel()
                # # SFSyst
                # nominal = self.DY[chan][...].ravel()
                # alternate = self.ftempl[chan]['DY_SFSystvar'][:].ravel()
                # up,down = self.mirrorShape(nominal,alternate)
                # dset_templ = f.create_dataset('DY_SFSystUp', alternate.ravel().shape, dtype=dtype,compression=compression)
                # dset_templ[...] = up.ravel()
                # dset_templ = f.create_dataset('DY_SFSystDown', alternate.ravel().shape, dtype=dtype,compression=compression)
                # dset_templ[...] = down.ravel()
                # #prefire
                # dset_bkg = f.create_dataset("DY_prefireUp", self.ftempl[chan]['DY_prefireVar'][...,0].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['DY_prefireVar'][...,0].ravel()
                # dset_bkg = f.create_dataset("DY_prefireDown", self.ftempl[chan]['DY_prefireVar'][...,1].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['DY_prefireVar'][...,1].ravel()
                # #JEC
                # dset_bkg = f.create_dataset("DY_jesUp", self.ftempl[chan]['DY_jesTotalUp'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['DY_jesTotalUp'][:].ravel()
                # dset_bkg = f.create_dataset("DY_jesDown", self.ftempl[chan]['DY_jesTotalDown'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['DY_jesTotalDown'][:].ravel()
                # dset_bkg = f.create_dataset("DY_unclUp", self.ftempl[chan]['DY_unclustEnUp'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['DY_unclustEnUp'][:].ravel()
                # dset_bkg = f.create_dataset("DY_unclDown", self.ftempl[chan]['DY_unclustEnDown'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['DY_unclustEnDown'][:].ravel()

                # dset_bkg = f.create_dataset("Diboson", self.Diboson[chan][...].ravel().shape, dtype=dtype,        compression=compression)
                # dset_bkg[...] = self.Diboson[chan][...].ravel()
                # dset_bkgw2 = f.create_dataset("Diboson_sumw2", self.Dibosonw2[chan][...].ravel().shape, dtype=dtype,      compression=compression)
                # dset_bkgw2[...] = self.Dibosonw2[chan][...].ravel()
                # #SF
                # shape_chan_SFStat = self.decorrelateSFSystBkg(self.Diboson[chan], self.ftempl[chan]['Diboson_SFStatvar'][...,0], self.ftempl[chan]['Diboson_SFStatvar'][...,1])
                # shape_chan_SFStat_iso = self.decorrelateSFSystBkg(self.Diboson[chan], self.ftempl[chan]['Diboson_SFStatvar'][...,2], self.ftempl[chan]['Diboson_SFStatvar'][...,3])
                # for iSF in range(shape_chan_SFStat.shape[-2]):
                #     dset_templ = f.create_dataset('Diboson_SFall{}Up'.format(iSF), shape_chan_SFStat[...,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat[...,iSF,0].ravel()
                #     dset_templ = f.create_dataset('Diboson_SFall{}Down'.format(iSF), shape_chan_SFStat[...,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat[...,iSF,1].ravel()
                                                
                #     dset_templ = f.create_dataset('Diboson_SFiso{}Up'.format(iSF), shape_chan_SFStat_iso[...,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat_iso[...,iSF,0].ravel()
                #     dset_templ = f.create_dataset('Diboson_SFiso{}Down'.format(iSF), shape_chan_SFStat_iso[...,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat_iso[...,iSF,1].ravel()
                # # SFSyst
                # nominal = self.Diboson[chan][...].ravel()
                # alternate = self.ftempl[chan]['Diboson_SFSystvar'][:].ravel()
                # up,down = self.mirrorShape(nominal,alternate)
                # dset_templ = f.create_dataset('Diboson_SFSystUp', alternate.ravel().shape, dtype=dtype,compression=compression)
                # dset_templ[...] = up.ravel()
                # dset_templ = f.create_dataset('Diboson_SFSystDown', alternate.ravel().shape, dtype=dtype,compression=compression)
                # dset_templ[...] = down.ravel()
                # #prefire
                # dset_bkg = f.create_dataset("Diboson_prefireUp", self.ftempl[chan]['Diboson_prefireVar'][...,0].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Diboson_prefireVar'][...,0].ravel()
                # dset_bkg = f.create_dataset("Diboson_prefireDown", self.ftempl[chan]['Diboson_prefireVar'][...,1].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Diboson_prefireVar'][...,1].ravel()
                # #JEC
                # dset_bkg = f.create_dataset("Diboson_jesUp", self.ftempl[chan]['Diboson_jesTotalUp'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Diboson_jesTotalUp'][:].ravel()
                # dset_bkg = f.create_dataset("Diboson_jesDown", self.ftempl[chan]['Diboson_jesTotalDown'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Diboson_jesTotalDown'][:].ravel()
                # dset_bkg = f.create_dataset("Diboson_unclUp", self.ftempl[chan]['Diboson_unclustEnUp'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Diboson_unclustEnUp'][:].ravel()
                # dset_bkg = f.create_dataset("Diboson_unclDown", self.ftempl[chan]['Diboson_unclustEnDown'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Diboson_unclustEnDown'][:].ravel()

                # dset_bkg = f.create_dataset("Top", self.Top[chan][...].ravel().shape, dtype=dtype,        compression=compression)
                # dset_bkg[...] = self.Top[chan][...].ravel()
                # dset_bkgw2 = f.create_dataset("Top_sumw2", self.Topw2[chan][...].ravel().shape, dtype=dtype,      compression=compression)
                # dset_bkgw2[...] = self.Topw2[chan][...].ravel()
                # #SF
                # shape_chan_SFStat = self.decorrelateSFSystBkg(self.Top[chan], self.ftempl[chan]['Top_SFStatvar'][...,0], self.ftempl[chan]['Top_SFStatvar'][...,1])
                # shape_chan_SFStat_iso = self.decorrelateSFSystBkg(self.Top[chan], self.ftempl[chan]['Top_SFStatvar'][...,2], self.ftempl[chan]['Top_SFStatvar'][...,3])
                # for iSF in range(shape_chan_SFStat.shape[-2]):
                #     dset_templ = f.create_dataset('Top_SFall{}Up'.format(iSF), shape_chan_SFStat[...,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat[...,iSF,0].ravel()
                #     dset_templ = f.create_dataset('Top_SFall{}Down'.format(iSF), shape_chan_SFStat[...,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat[...,iSF,1].ravel()
                                                
                #     dset_templ = f.create_dataset('Top_SFiso{}Up'.format(iSF), shape_chan_SFStat_iso[...,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat_iso[...,iSF,0].ravel()
                #     dset_templ = f.create_dataset('Top_SFiso{}Down'.format(iSF), shape_chan_SFStat_iso[...,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat_iso[...,iSF,1].ravel()
                # SFSyst
                # nominal = self.Top[chan][...].ravel()
                # alternate = self.ftempl[chan]['Top_SFSystvar'][:].ravel()
                # up,down = self.mirrorShape(nominal,alternate)
                # dset_templ = f.create_dataset('Top_SFSystUp', alternate.ravel().shape, dtype=dtype,compression=compression)
                # dset_templ[...] = up.ravel()
                # dset_templ = f.create_dataset('Top_SFSystDown', alternate.ravel().shape, dtype=dtype,compression=compression)
                # dset_templ[...] = down.ravel()
                # #prefire
                # dset_bkg = f.create_dataset("Top_prefireUp", self.ftempl[chan]['Top_prefireVar'][...,0].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Top_prefireVar'][...,0].ravel()
                # dset_bkg = f.create_dataset("Top_prefireDown", self.ftempl[chan]['Top_prefireVar'][...,1].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Top_prefireVar'][...,1].ravel()
                # #JEC
                # dset_bkg = f.create_dataset("Top_jesUp", self.ftempl[chan]['Top_jesTotalUp'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Top_jesTotalUp'][:].ravel()
                # dset_bkg = f.create_dataset("Top_jesDown", self.ftempl[chan]['Top_jesTotalDown'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Top_jesTotalDown'][:].ravel()
                # dset_bkg = f.create_dataset("Top_unclUp", self.ftempl[chan]['Top_unclustEnUp'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Top_unclustEnUp'][:].ravel()
                # dset_bkg = f.create_dataset("Top_unclDown", self.ftempl[chan]['Top_unclustEnDown'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Top_unclustEnDown'][:].ravel()

                # dset_bkg = f.create_dataset("Wtau", self.Wtau[chan][...].ravel().shape, dtype=dtype,      compression=compression)
                # dset_bkg[...] = self.Wtau[chan][...].ravel()
                # dset_bkgw2 = f.create_dataset("Wtau_sumw2", self.Wtauw2[chan][...].ravel().shape, dtype=dtype,        compression=compression)
                # dset_bkgw2[...] = self.Wtauw2[chan][...].ravel()
                # # pdf variations
                # for k in range(1,103):
                #     nominal = self.Wtau[chan][...]
                #     alternate = self.ftempl[chan]['Wtau_LHEPdfWeight'][...,k]
                #     up,down = self.mirrorShape(nominal,alternate)
                #     dset_templ = f.create_dataset("Wtau_pdf{}Up".format(k), up.ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = up.ravel()
                #     dset_templ = f.create_dataset("Wtau_pdf{}Down".format(k), down.ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = down.ravel()
                # #QCD
                # for j,syst in self.qcdsyst.items():
                #     dset_bkg = f.create_dataset(name='Wtau_{}'.format(syst), shape=self.ftempl[chan]['Wtau_LHEScaleWeight'][...,j].ravel().shape, dtype='float64')
                #     dset_bkg[...] = self.ftempl[chan]['Wtau_LHEScaleWeight'][...,j].ravel()
                # #SF
                # shape_chan_SFStat = self.decorrelateSFSystBkg(self.Wtau[chan], self.ftempl[chan]['Wtau_SFStatvar'][...,0], self.ftempl[chan]['Wtau_SFStatvar'][...,1])
                # shape_chan_SFStat_iso = self.decorrelateSFSystBkg(self.Wtau[chan], self.ftempl[chan]['Wtau_SFStatvar'][...,2], self.ftempl[chan]['Wtau_SFStatvar'][...,3])
                # for iSF in range(shape_chan_SFStat.shape[-2]):
                #     dset_templ = f.create_dataset('Wtau_SFall{}Up'.format(iSF), shape_chan_SFStat[...,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat[...,iSF,0].ravel()
                #     dset_templ = f.create_dataset('Wtau_SFall{}Down'.format(iSF), shape_chan_SFStat[...,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat[...,iSF,1].ravel()
                                                
                #     dset_templ = f.create_dataset('Wtau_SFiso{}Up'.format(iSF), shape_chan_SFStat_iso[...,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat_iso[...,iSF,0].ravel()
                #     dset_templ = f.create_dataset('Wtau_SFiso{}Down'.format(iSF), shape_chan_SFStat_iso[...,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                #     dset_templ[...] = shape_chan_SFStat_iso[...,iSF,1].ravel()
                # # SFSyst
                # nominal = self.Wtau[chan][...].ravel()
                # alternate = self.ftempl[chan]['Wtau_SFSystvar'][:].ravel()
                # up,down = self.mirrorShape(nominal,alternate)
                # dset_templ = f.create_dataset('Wtau_SFSystUp', alternate.ravel().shape, dtype=dtype,compression=compression)
                # dset_templ[...] = up.ravel()
                # dset_templ = f.create_dataset('Wtau_SFSystDown', alternate.ravel().shape, dtype=dtype,compression=compression)
                # dset_templ[...] = down.ravel()
                # #prefire
                # dset_bkg = f.create_dataset("Wtau_prefireUp", self.ftempl[chan]['Wtau_prefireVar'][...,0].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Wtau_prefireVar'][...,0].ravel()
                # dset_bkg = f.create_dataset("Wtau_prefireDown", self.ftempl[chan]['Wtau_prefireVar'][...,1].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Wtau_prefireVar'][...,1].ravel()
                # #JEC
                # dset_bkg = f.create_dataset("Wtau_jesUp", self.ftempl[chan]['Wtau_jesTotalUp'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Wtau_jesTotalUp'][:].ravel()
                # dset_bkg = f.create_dataset("Wtau_jesDown", self.ftempl[chan]['Wtau_jesTotalDown'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Wtau_jesTotalDown'][:].ravel()
                # dset_bkg = f.create_dataset("Wtau_unclUp", self.ftempl[chan]['Wtau_unclustEnUp'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Wtau_unclustEnUp'][:].ravel()
                # dset_bkg = f.create_dataset("Wtau_unclDown", self.ftempl[chan]['Wtau_unclustEnDown'][:].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['Wtau_unclustEnDown'][:].ravel()

                # dset_bkg = f.create_dataset("fakesLowMt", self.fakeslow[chan].ravel().shape, dtype=dtype,     compression=compression)
                # dset_bkg[...] = self.fakeslow[chan].ravel()
                # dset_bkgw2 = f.create_dataset("fakesLowMt_sumw2", self.fakesloww2[chan].ravel().shape, dtype=dtype,       compression=compression)
                # dset_bkgw2[...] = self.fakesloww2[chan].ravel()

                # dset_bkg = f.create_dataset("fakesHighMt", self.fakeshigh[chan].ravel().shape, dtype=dtype,       compression=compression)
                # dset_bkg[...] = self.fakeshigh[chan].ravel()
                # dset_bkgw2 = f.create_dataset("fakesHighMt_sumw2", self.fakeshighw2[chan].ravel().shape, dtype=dtype,     compression=compression)
                # dset_bkgw2[...] = self.fakeshighw2[chan].ravel()

                dset_bkg = f.create_dataset("LowAcc", self.lowacc[chan][...].ravel().shape, dtype=dtype,     compression=compression)
                dset_bkg[...] = self.lowacc[chan][...].ravel()
                dset_bkgw2 = f.create_dataset("LowAcc_sumw2", self.lowacc[chan][...].ravel().shape, dtype=dtype,     compression=compression)
                dset_bkgw2[...] = self.lowaccw2[chan][...].ravel()

                low_mass = self.ftempl[chan]['lowacc_mass'][:][...]
                dset_bkg = f.create_dataset("LowAcc_massUp", low_mass[...,0].ravel().shape, dtype=dtype,      compression=compression)
                dset_bkg[...] = low_mass[...,0].ravel()
                dset_bkg = f.create_dataset("LowAcc_massDown", low_mass[...,1].ravel().shape, dtype=dtype,        compression=compression)
                dset_bkg[...] = low_mass[...,1].ravel()

                # low_pdf = self.ftempl[chan]['lowacc_LHEPdfWeight'][:][...]
                # # pdf variations
                # for k in range(1,103):
                #     nominal = self.lowacc[chan][...]
                #     alternate = low_pdf[...,k]
                #     up,down = self.mirrorShape(nominal,alternate)
                #     dset_templ = f.create_dataset("LowAcc_pdf{}Up".format(k), up.ravel().shape, dtype=dtype,      compression=compression)
                #     dset_templ[...] = up.ravel()
                #     dset_templ = f.create_dataset("LowAcc_pdf{}Down".format(k), down.ravel().shape, dtype=dtype,      compression=compression)
                #     dset_templ[...] = down.ravel()

                # # copy lowacc nuisances in shape file
                # for j,syst in self.qcdsyst.items():
                #     for i in range(4):
                #         if j<4:
                #             histo = 1./3*self.ftempl[chan]['lowacc_LHEScaleWeight{}'.format(i)][:][...,j]
                #         else:
                #             histo = 3.*self.ftempl[chan]['lowacc_LHEScaleWeight{}'.format(i)][:][...,j]
                #         dsetdown = f.create_dataset(name='LowAcc_LHEScaleWeight{}_{}'.format(i+1,syst), shape=histo.ravel().shape, dtype='float64')
                #         dsetdown[...] = histo.ravel()
                # #prefire
                # dset_bkg = f.create_dataset("LowAcc_prefireUp", self.ftempl[chan]['lowacc_prefireVars'][...,0].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['lowacc_prefireVars'][...,0].ravel()
                # dset_bkg = f.create_dataset("LowAcc_prefireDown", self.ftempl[chan]['lowacc_prefireVars'][...,1].ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = self.ftempl[chan]['lowacc_prefireVars'][...,1].ravel()
                # #JEC
                # dset_bkg = f.create_dataset("LowAcc_jesUp", np.sum(self.ftempl[chan]['lowacc_jesTotalUp'][:],axis=-1).ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = np.sum(self.ftempl[chan]['lowacc_jesTotalUp'][:],axis=-1).ravel()
                # dset_bkg = f.create_dataset("LowAcc_jesDown", np.sum(self.ftempl[chan]['lowacc_jesTotalDown'][:],axis=-1).ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = np.sum(self.ftempl[chan]['lowacc_jesTotalDown'][:],axis=-1).ravel()
                # dset_bkg = f.create_dataset("LowAcc_unclUp", np.sum(self.ftempl[chan]['lowacc_unclustEnUp'][:],axis=-1).ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = np.sum(self.ftempl[chan]['lowacc_unclustEnUp'][:],axis=-1).ravel()
                # dset_bkg = f.create_dataset("LowAcc_unclDown", np.sum(self.ftempl[chan]['lowacc_unclustEnDown'][:],axis=-1).ravel().shape, dtype=dtype,compression=compression)
                # dset_bkg[...] = np.sum(self.ftempl[chan]['lowacc_unclustEnDown'][:],axis=-1).ravel()
                # # #SF
                # # shape_chan_SFStat = self.decorrelateSFSystBkg(self.lowacc[chan], self.ftempl[chan]['lowacc_SFStatvar'][...,0], self.ftempl[chan]['lowacc_SFStatvar'][...,1])
                # # shape_chan_SFStat_iso = self.decorrelateSFSystBkg(self.lowacc[chan], self.ftempl[chan]['lowacc_SFStatvar'][...,2], self.ftempl[chan]['lowacc_SFStatvar'][...,3])
                # # for iSF in range(shape_chan_SFStat.shape[-2]):
                # #     dset_templ = f.create_dataset('LowAcc_SFall{}Up'.format(iSF), shape_chan_SFStat[...,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                # #     dset_templ[...] = shape_chan_SFStat[...,iSF,0].ravel()
                # #     dset_templ = f.create_dataset('LowAcc_SFall{}Down'.format(iSF), shape_chan_SFStat[...,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                # #     dset_templ[...] = shape_chan_SFStat[...,iSF,1].ravel()
                                                
                # #     dset_templ = f.create_dataset('LowAcc_SFiso{}Up'.format(iSF), shape_chan_SFStat_iso[...,iSF,0].ravel().shape, dtype=dtype,compression=compression)
                # #     dset_templ[...] = shape_chan_SFStat_iso[...,iSF,0].ravel()
                # #     dset_templ = f.create_dataset('LowAcc_SFiso{}Down'.format(iSF), shape_chan_SFStat_iso[...,iSF,1].ravel().shape, dtype=dtype,compression=compression)
                # #     dset_templ[...] = shape_chan_SFStat_iso[...,iSF,1].ravel()
                # # SFSyst
                # nominal = self.lowacc[chan][...].ravel()
                # alternate = self.ftempl[chan]['lowacc_SFSystvar'][:].ravel()
                # up,down = self.mirrorShape(nominal,alternate)
                # dset_templ = f.create_dataset('LowAcc_SFSystUp', alternate.ravel().shape, dtype=dtype,compression=compression)
                # dset_templ[...] = up.ravel()
                # dset_templ = f.create_dataset('LowAcc_SFSystDown', alternate.ravel().shape, dtype=dtype,compression=compression)
                # dset_templ[...] = down.ravel()
                # # copy fakes nuisances in shape file
                # for type in ["Up","Down"]:
                #     for i in range(48*30):
                #         # histo = self.ftempl['fakesLowMt_fakeNormLowMtBin{}{}'.format(i,type)][:]
                #         # dset = f.create_dataset(name='fakesLowMt_fakeNormLowMtBin{}{}'.format(i,type), shape=histo.     ravel().shape, dtype='float64')
                #         # dset[...] = histo.ravel()
                #         # histo = self.ftempl['fakesHighMt_fakeNormHighMtBin{}{}'.format(i,type)][:]
                #         # dset = f.create_dataset(name='fakesHighMt_fakeNormHighMtBin{}{}'.format(i,type,        shape=histo.ravel().shape, dtype='float64')
                #         # dset[...] = histo.ravel()
                #         histo = self.ftempl[chan]['fakesLowMt_fakeShapeBin{}{}'.format(i,type)][:]
                #         dset = f.create_dataset(name='fakesLowMt_fakeShapeBin{}{}'.format(i,type), shape=histo.ravel  ().shape, dtype=dtype,compression=compression)
                #         dset[...] = histo.ravel()

                #         histo = self.ftempl[chan]['fakesHighMt_fakeShapeBin{}{}'.format(i,type)][:]
                #         dset = f.create_dataset(name='fakesHighMt_fakeShapeBin{}{}'.format(i,type), shape=histo.ravel ().shape, dtype=dtype,compression=compression)
                #         dset[...] = histo.ravel()

    def maskedChannels(self):
        dtype = 'float64'
        compression = "gzip"
        for chan in self.channels:
            with h5py.File('{}_xsec.hdf5'.format(chan), mode="w") as f:
                for proc in self.processes:
                    if "helXsecs" in proc: #give the correct xsec to unfold
                        coeff = self.helXsecs.index(proc.split('_')[0].replace('helXsecs',''))
                        iY = int(proc.split('_')[2])
                        iQt = int(proc.split('_')[4])
                        
                        dset_masked = f.create_dataset(proc, [1], dtype=dtype,compression=compression)
                        dset_masked[...] = self.gen[chan][iY,iQt,coeff]

                    else:
                        dset_masked = f.create_dataset(proc, [1], dtype=dtype,compression=compression)
                        dset_masked[...] = 0.
                dset_masked = f.create_dataset("data_obs", [1], dtype=dtype,compression=compression)
                dset_masked[...] = 1.

    def setPreconditionVec(self):
        # f=h5py.File('fitresults_WPlus_blockAi_260292.hdf5', 'r')
        # hessian = f['hess'][:]
        # eig, U = np.linalg.eigh(hessian)
        # M1 = np.matmul(np.diag(1./np.sqrt(eig)),U.T)
        # # print(M1,np.linalg.inv(np.linalg.inv(M1)))
        # self.preconditioner = M1
        # self.invpreconditioner = np.linalg.inv(self.preconditioner)
        # print(np.matmul(self.preconditioner,np.linalg.inv(self.preconditioner)))
        # test = np.matmul(self.preconditioner,np.linalg.inv(self.preconditioner)) - np.identity(M1.shape[0])
        # print(np.max(np.abs(test)))
        self.preconditioner = np.identity(len(self.signals))
        self.invpreconditioner = np.identity(len(self.signals))

    def fillHelGroup(self):

        for i in range(len(self.yBinsC)):
            for j in range(len(self.qtBinsC)):

                s = 'y_{i}_qt_{j}'.format(i=i,j=j)
                self.helGroups[s] = []
                
                for hel in self.helXsecs:
                    if not '7' in hel and not '8' in hel and not '9' in hel:
                        self.helGroups[s].append('helXsecs'+hel+'_'+s)
                if self.helGroups[s] == []:
                    del self.helGroups[s]
    def fillHelMetaGroup(self):

        for i in range(len(self.yBinsC)):
            s = 'y_{i}'.format(i=i)
            self.helMetaGroups[s] = []
            for key in self.sumGroups:
                if s in key:
                    self.helMetaGroups[s].append(key)
            
            if self.helMetaGroups[s] == []:
                    del self.helMetaGroups[s]
        
        for j in range(len(self.qtBinsC)):
            s = 'qt_{j}'.format(j=j)
            self.helMetaGroups[s] = []
            for key in self.sumGroups:
                if 'qt' in key and key.split('_')[2]==str(j):
                    self.helMetaGroups[s].append(key)
        
            if self.helMetaGroups[s] == []:
                    del self.helMetaGroups[s]
    def fillSumGroup(self):

        for i in range(len(self.yBinsC)):
            s = 'y_{i}'.format(i=i)
            for hel in self.helXsecs:
                if '7' in hel or '8' in hel or '9' in hel: continue
                for signal in self.processes:
                    self.sumGroups['helXsecs'+hel+'_'+s] = []
                    for j in range(len(self.qtBinsC)):
                        #if 'helXsecs'+hel+'_'+'y_{i}_qt_{j}'.format(i=i,j=j) in self.processes:
                        self.sumGroups['helXsecs'+hel+'_'+s].append('helXsecs'+hel+'_'+s+'_qt_{j}'.format(j=j))
        
        for j in range(len(self.qtBinsC)):
            s = 'qt_{j}'.format(j=j)
            for hel in self.helXsecs:
                if '7' in hel or '8' in hel or '9' in hel: continue
                for signal in self.processes:
                    if signal.split('_')[0] == 'helXsecs'+hel and signal.split('_')[4] == str(j):
                        self.sumGroups['helXsecs'+hel+'_'+s] = []
                        for i in range(len(self.yBinsC)):
                            #print i, signal, 'helXsecs'+hel+'_'+'y_{i}_pt_{j}'.format(i=i,j=j)
                            #print 'append', 'helXsecs'+hel+'_y_{i}_'.format(i=i)+s, 'to', 'helXsecs'+hel+'_'+s
                            self.sumGroups['helXsecs'+hel+'_'+s].append('helXsecs'+hel+'_y_{i}_'.format(i=i)+s)
        
    def makeDatacard(self):

        self.DC = Datacard()

        ############## Setup the datacard (must be filled in) ###########################
        self.DC.bins =   copy.deepcopy(self.channels)
        self.DC.bins.extend([c+'_xsec' for c in self.channels]) # <type 'list'>
        self.DC.obs =    {} # <type 'dict'>
        self.DC.processes =  self.processes # <type 'list'>
        self.DC.signals =    self.signals # <type 'list'>
        self.DC.isSignal =   {} # <type 'dict'>
        for proc in self.processes:
            if proc in self.signals:
                self.DC.isSignal[proc] = True
            else:
                self.DC.isSignal[proc] = False
        self.DC.keyline = [] # <type 'list'> # not used by combine-tf
        self.DC.exp =    {} # <type 'dict'>
        for chan in self.channels:
            self.DC.exp[chan] = {}
            self.DC.exp[chan+'_xsec'] = {}
            for proc in self.processes:
                self.DC.exp[chan][proc] = -1.00
                self.DC.exp[chan+'_xsec'][proc] = -1.00
        self.DC.systs =  [] # <type 'list'>
        self.DC.groups = {}
        # aux = {} #each sys will have a separate aux dict
        # for i in range(48*30):
        #     for chan in self.channels:
        #         aux[chan] = {}
        #         aux[chan+'_xsec'] = {}
        #         for proc in self.processes:
        #             aux[chan][proc] = 0.
        #             aux[chan+'_xsec'][proc] = 0.
        #         aux[chan]['fakesLowMt'] = 1.
        #         aux[chan]['fakesHighMt'] = 1.
        #     self.DC.systs.append(('fakeShapeBin{}'.format(i), False, 'shapeNoConstraint', [], aux))
        # aux = {} #each sys will have a separate aux dict
        # for chan in self.channels:
        #     aux[chan] = {}
        #     aux[chan+'_xsec'] = {}
        #     for proc in self.processes:
        #             aux[chan][proc] = 0.
        #             aux[chan+'_xsec'][proc] = 0.
        #     aux[chan]['fakesLowMt'] = 1.5
        #     aux[chan]['fakesHighMt'] = 0
        # self.DC.systs.append(('fakesNormLowMt', False, 'lnNNoConstraint', [], aux))
        # aux = {} #each sys will have a separate aux dict
        # for chan in self.channels:
        #     aux[chan] = {}
        #     aux[chan+'_xsec'] = {}
        #     for proc in self.processes:
        #         aux[chan][proc] = 0.
        #         aux[chan+'_xsec'][proc] = 0.
        #     aux[chan]['fakesLowMt'] = 0.
        #     aux[chan]['fakesHighMt'] = 1.5
        # self.DC.systs.append(('fakesNormHighMt', False, 'lnNNoConstraint', [], aux))
        # aux = {} #each sys will have a separate aux dict
        # aux[chan] = {}
        # aux[chan+'_xsec'] = {}
        # for i in range(48*30):
        #     for proc in self.processes:
        #         aux[self.channel][proc] = 0.
        #         aux[self.channel+'_xsec'][proc] = 0.
        #     aux[self.channel]['fakesLowMt'] = 1.
        #     aux[self.channel]['fakesHighMt'] = 0.
        #     self.DC.systs.append(('fakeNormLowMtBin{}'.format(i), False, 'shapeNoConstraint', [], aux))
        # aux = {} #each sys will have a separate aux dict
        # aux[self.channel] = {}
        # aux[self.channel+'_xsec'] = {}
        # for i in range(48*30):
        #     for proc in self.processes:
        #         aux[self.channel][proc] = 0.
        #         aux[self.channel+'_xsec'][proc] = 0.
        #     aux[self.channel]['fakesLowMt'] = 0.
        #     aux[self.channel]['fakesHighMt'] = 1.
        #     self.DC.systs.append(('fakeNormHighMtBin{}'.format(i), False, 'shapeNoConstraint', [], aux))
        # list of [{bin : {process : [input file, path to shape, path to shape for uncertainty]}}]
        if self.doSyst:
            for syst in self.templSystematics: #loop over systematics
                for var in self.templSystematics[syst]["vars"]:
                    aux = {} #each sys will have a separate aux dict
                    for chan in self.channels:
                        aux[chan] = {}
                        aux[chan+'_xsec'] = {}
                        for proc in self.processes: 
                            if proc in self.templSystematics[syst]["procs"]:
                                aux[chan][proc] = self.templSystematics[syst]["weight"]
                                aux[chan+'_xsec'][proc] = 0.0
                            else:
                                if "Signal" in self.templSystematics[syst]["procs"] and proc in self.signals:
                                    aux[chan][proc] = self.templSystematics[syst]["weight"]
                                    aux[chan+'_xsec'][proc] = 0.0
                                else:
                                    aux[chan][proc] = 0.0
                                    aux[chan+'_xsec'][proc] = 0.0
                    self.DC.systs.append((var, False, self.templSystematics[syst]["type"], [], aux))
        self.DC.groups = {'mass': ['mass']}
        # self.DC.groups = {
        #                 'mass': ['mass'],
        #                 'pdfs': set(["pdf{}".format(i) for i in range(1,103)]),
        #                 'QCD theory': set(["LHEScaleWeight{}_muRmuF".format(i) for i in range(1,5)]+["LHEScaleWeight{}_muR".format(i) for i in range(1,5)]+["LHEScaleWeight{}_muF".format(i) for i in range(1,5)]+["muRmuF", "muR","muF"]),
        #                 'QCD background': ['fakeShapeBin{}'.format(i) for i in range(48*30)]+["fakesNormLowMt","fakesNormHighMt"],
        #                 "luminosity": ['CMSlumi'],
        #                 "SF":set(["SFall{}".format(i) for i in range(624)]+["SFiso{}".format(i) for i in range(624)]+["SFSyst"]),
        #                 'jme': set(['jes', 'uncl']),
        #                 'prefire':['prefire'],
        #                  }  # <type 'dict'>
        for chan in self.channels:
            self.DC.shapeMap = 	{chan: {'*': [chan+'.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']},\
            chan+'_xsec': {'*': [chan+'_xsec.root', '$PROCESS', '$PROCESS_$SYSTEMATIC']}} # <type 'dict'>
        self.DC.hasShapes =  True # <type 'bool'>
        self.DC.flatParamNuisances =  {} # <type 'dict'>
        self.DC.rateParams =  {} # <type 'dict'>
        self.DC.extArgs =    {} # <type 'dict'>
        self.DC.rateParamsOrder  =  set([]) # <type 'set'>
        self.DC.frozenNuisances  =  set([]) # <type 'set'>
        self.DC.systematicsShapeMap =  {} # <type 'dict'>
        self.DC.nuisanceEditLines    =  [] # <type 'list'>
        self.DC.discretes    =  [] # <type 'list'>
        self.DC.helGroups = self.helGroups
        self.DC.sumGroups = self.sumGroups
        self.DC.helMetaGroups = self.helMetaGroups
        self.DC.noiGroups = {'mass':['mass']}
        # self.DC.noiGroups = {}

        self.DC.preconditioner  = self.preconditioner 
        self.DC.invpreconditioner  = self.invpreconditioner 
        
        filehandler = open('{}.pkl'.format(chan.split("_")[0]), 'w')
        pickle.dump(self.DC, filehandler)
