import ROOT
from array import array
import math
from termcolor import colored
import copy
import argparse
import os
import sys
# sys.path.append('../../')
# sys.path.append('../')
from systToapply import systematicsDict



ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

class plotter :
    
    def __init__(self, anaKind,sign, asimov):
        
        self.anaKind = anaKind
        self.sign = sign
        self.asimov = asimov
        self.signDict = {
            'plus' :  'W^{+}',
            'minus' : 'W^{-}'
        }    
        self.yArr = [0, 0.4, 0.8, 1.2 ,1.6, 2.0, 2.4]
        # self.yArr = [0, 0.4, 0.8, 1.2 ,1.6, 2.0, 2.4, 2.8] #Exteded
        # self.qtArr = [0., 4., 8., 12., 16., 20., 24., 28., 32.] #equal size
        # self.qtArr = [0.,   3.1,  5.,   7.,   9.4, 12.4, 16.5, 22.3, 32.]#quantile
        # self.qtArr = [0.,  2.,  4., 6.,  8., 10., 12., 16., 22., 32.] #2GeV
        # self.qtArr = [0.,  2.,  4., 6.,  8., 10., 12., 14., 16., 20., 25., 32., 45., 80] #Extended
        self.qtArr =[0.,  2.,  4., 6.,  8., 10., 12., 16., 20., 26., 36., 60.,]
        # self.qtArr = [0, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 13., 14., 15., 16.,22,32] #05GeV

        self.coeffArr = [0,1,2,3,4,5,6]
        # self.coeffArr = [0,1,2,3,4,5]

        self.noiArr = ['mass']
        self.dirList = ['up','down']
        self.RatioPadYcut = 2.
        
        self.RatioPadYcutDict = {#unrolled, vs y, vs qt
            'A0' : [0.2,2,1.2] , 
            'A1' : [0.2,3,3],
            'A2' : [0.4,4,4],
            'A3' : [2.,7.,7.] ,
            'A4' : [2.,2.,2.] ,
            'unpolarizedxsec' : [2.,2.,2.]
            
        }
        
        self.lumi=35.9
        self.mass = 80.419 
        self.toyPullFit = True
        
        self.systDict = {
            "_LHEScaleWeight" : ["_LHEScaleWeight_muR0p5_muF0p5", "_LHEScaleWeight_muR0p5_muF1p0","_LHEScaleWeight_muR1p0_muF0p5","_LHEScaleWeight_muR1p0_muF2p0","_LHEScaleWeight_muR2p0_muF1p0", "_LHEScaleWeight_muR2p0_muF2p0"],
            "_LHEPdfWeight" : ["_LHEPdfWeightHess" + str(i)  for i in range(1, 61)],
        }
        
        self.vetoScaleList = []
        for scNum in self.systDict['_LHEScaleWeight'] :
            for scDen in self.systDict['_LHEScaleWeight'] :
                if ('0p5' in scNum and '2p0' in scDen) or ('2p0' in scNum and '0p5' in scDen): 
                    self.vetoScaleList.append([scNum,scDen])        
        
        if self.anaKind['angNames'] : 
            self.coeffDict = {
                'A0' : [1.,'A0', 'A_{0}', ROOT.kYellow+2],
                'A1' : [5.,'A1', 'A_{1}', ROOT.kViolet-1],
                'A2' : [20.,'A2', 'A_{2}', ROOT.kAzure+1],
                'A3' : [4.,'A3', 'A_{3}', ROOT.kGreen+2],
                'A4' : [4.,'A4', 'A_{4}', ROOT.kRed+2],
                'unpolarizedxsec' : [1,'unpolarizedxsec','#sigma^{U+L}', ROOT.kBlue-4]
            }
            self.coeffList = ['A0','A1','A2','A3','A4','unpolarizedxsec' ]
            # self.coeffList = ['A1','A2','A3','A4','unpolarizedxsec' ]
        else :
            self.coeffDict = {
                'L' : [2.,'A0', '#sigma_{L}'],
                'I' : [2*math.sqrt(2),'A1', '#sigma_{I}'],
                'T' : [4.,'A2', '#sigma_{T}'],
                'A' : [4*math.sqrt(2), 'A3', '#sigma_{A}'],
                'P' : [2.,'A4', '#sigma_{P}'],
                'UL' : [1., 'unpolarizedxsec', '#sigma_{U+L}']
            }
            self.coeffList = ['L', 'I', 'T', 'A', 'P', 'UL']
        
        #DEBUG josh
        # if not self.helXsec : 
        #     self.coeffDict = {
        #         'unpolarizedxsec' : [1,'unpolarizedxsec']
        #     }
        #     self.coeffList = ['unpolarizedxsec' ]
        # else :
        #     self.coeffDict = {
        #         'UL' : [1., 'unpolarizedxsec']
        #     }
        #     self.coeffList = ['UL']
        # self.coeffArr = [0]
        ############################################################ END OF DEBUG
        
        self.histos = {}
        self.canvas = {}
        self.leg = {}
        
        self.groupedSystColors = {
            # "WHSFVars"  : [ROOT.kGreen+1, 'Scale Factors'],
            # "LHEScaleWeightVars" : [ROOT.kViolet-2, 'MC Scale'],
            "LHEScaleWeightVars" : [ROOT.kGreen-3, 'QCD Scales'],
            # "ptScaleVars" : [ROOT.kBlue-4, 'pT Scale'],
            # "jmeVars" : [ROOT.kAzure+10, 'MET'],
            "LHEPdfWeightVars" : [ROOT.kRed+1, 'PDF+#alpha_{s}'],
            # "LHEPdfWeightVars" : [ROOT.kBlue-4, 'NNPDF3.0'],
            "Nominal" : [1, 'Fit Unc.'],
            "poiss" : [ROOT.kAzure+10,'Poiss. Unc.']
        }
        
        self.nuisanceDict = {
            "mass"      : [ROOT.kBlue-4, 'm_{W}', 29],
            "WHSFSyst"  : [ROOT.kGreen-4, 'SF syst.', 36],
            "pdfs"      : [ROOT.kRed+1, 'PDF+#alpha_{s}', 21],
            "WHSFStat"  : [ROOT.kGreen+1, 'SF stat.', 22],
            "stat"      : [ROOT.kGray+2, 'Statistical', 43],
            "PrefireWeight" : [ROOT.kSpring+10, 'Prefire weight', 28],
            "jme"       : [ROOT.kAzure+10, 'MET uncert.',45],
            "ptScale"       : [ROOT.kYellow+2, 'p_{T} Scale',43],
            # "binByBinStat" : [ROOT.kGray, 'BBB', 46],
            "QCDnorm" : [ROOT.kViolet, 'QCD norm.', 22],
            "CMSlumi" : [ROOT.kOrange-7,"Luminosity",33],
            "ewkXsec" : [ROOT.kTeal,"#sigma_{W#rightarrow#tau#nu}+#sigma_{t}+#sigma_{diboson}",3],
            "LeptonVeto" : [ROOT.kMagenta-7,"Lepton veto",23],
            "WQT": [ROOT.kViolet+7,"q_{T}^{V}",20],
            "tot" : [ROOT.kBlack, 'Total unc.',1]
            
            
            # "alphaS"    : [ROOT.kOrange-3, '#alpha_{s}', 38],
            # "LHEScaleWeight" : [ROOT.kViolet-2,"MC Scales",43],  
            # "DYxsec" : [ROOT.kCyan+2,"#sigma_{DY}",3],
            # "Topxsec" : [ROOT.kCyan-6,"#sigma_{t}",3],
            # "Dibosonxsec" : [ROOT.kCyan-1,"#sigma_{diboson}",3],
            # "Tauxsec" : [ROOT.kTeal,"#sigma_{W#rightarrow#tau#nu}",3],
        }
        
        self.NuiConstrLabels = {
            'CMSlumi' : 'Luminosity',
            'Topxsec' : '#sigma_{t}',
            'Dibosonxsec' : '#sigma_{diboson}',
            'Tauxsec' : '#sigma_{#tau}',
            'mass' : 'm_{W}',
            'WHSFStat' : 'SF stat.',
            'WHSFSystFlat' : 'SF syst.',
            'jesTotal' : 'p_{T}^{miss} (JES)',
            'unclustEn' :  'p_{T}^{miss} (E_{U})',
            'PrefireWidth' : 'Prefire weight',
            'LHEPdfWeight' : 'PDF',
            'alphaS' : '#alpha_{s}',
            'LHEScaleWeight' : 'q_{T}^{V} (MC scales)',
            'corrected' : 'p_{T}^{#mu} scale',
            'LeptonVeto' : 'Lepton veto',
            'lumi' : 'Lumi. (QCD)',
            'QCDnorm' : 'QCD norm.'
            
            
        }
        
        self.category = {
            ''      : '',
            'y'     : 'Y',
            'qt'    : 'Pt'
            }
        
        self.processDict = {
            'obs' : [ROOT.kBlack, 'Data'],
            'full' : [ROOT.kBlack,'full'],
            'sig' : [ROOT.kRed+2,       "W^{+}#rightarrow #mu^{+}#nu_{#mu}", ROOT.kRed+1, "W^{-}#rightarrow #mu^{-}#bar{#nu}_{#mu}"],
            'bkg' : [ROOT.kBlue, 'bkg'],
            'DYJets' : [ROOT.kAzure+2,     "Drell-Yan"],
            'DiBoson' : [ROOT.kViolet+2,    "Diboson"],
            'Top' : [ROOT.kGreen+3,     "Top"],
            'Fake' : [ ROOT.kGray,        "QCD"],
            'WtoTau' : [ROOT.kSpring+9,    "W^{#pm}#rightarrow #tau^{#pm}#nu_{#tau}"],
            'LowAcc' : [ROOT.kCyan-6, "Low-acceptance"],
            }
        # self.postFitVars = {
        #     ""   :   [ "unrolled #eta(p_{T})",    "dN/dp_{T}d#eta [GeV^{-1}]", " [GeV]"],
        #     "pt"   :  [ "p_{T}^{#mu}",    "dN/dp_{T} [GeV^{-1}]", " p_{T}^{#mu} [GeV]"],
        #     "eta"  :  ["#eta^{#mu}",     "dN/d#eta", "#eta^{#mu}"],
        # }
        self.postFitVars = {
            ""   :   [ "unrolled #eta(p_{T})",    "Events / ", " [GeV]", ' GeV'],
            "pt"   :  [ "p_{T}^{#mu}",    "Events / ", " p_{T}^{#mu} [GeV]", ' GeV'],
            "eta"  :  ["#eta^{#mu}",     "Events / ", "#eta^{#mu}", ''],
        }
        self.processOrder = ['DiBoson','WtoTau','Top','DYJets','Fake','LowAcc','sig']
        self.ptArr = []
        self.etaArr = []
        for ipt in range(0,61) :
            self.ptArr.append(25. + ipt*(55.-25.)/60.)
        for ieta in range(0,49) :
            self.etaArr.append(-2.4 + ieta * (4.8) / 48.)
        
        self.SPLITSIG=False  
        self.hxsecList = ['L', 'I', 'T', 'A', 'P', 'UL']
        if self.SPLITSIG: 
            self.processDict_sigSplit = copy.deepcopy(self.processDict)
            del self.processDict_sigSplit['sig'] 
            self.processDict_sigSplit['sig_L'] = [ROOT.kPink-9, "W^{+}#rightarrow #mu^{+}#nu_{#mu} #sigma_{L}",   ROOT.kPink-9, "W^{-}#rightarrow #mu^{-}#bar{#nu}_{#mu} #sigma_{L}", ]
            self.processDict_sigSplit['sig_I'] = [ROOT.kRed-8, "W^{+}#rightarrow #mu^{+}#nu_{#mu} #sigma_{I}",   ROOT.kRed-8, "W^{-}#rightarrow #mu^{-}#bar{#nu}_{#mu} #sigma_{I}", ]
            self.processDict_sigSplit['sig_T'] = [ROOT.kOrange+3, "W^{+}#rightarrow #mu^{+}#nu_{#mu} #sigma_{T}",  ROOT.kOrange+3, "W^{-}#rightarrow #mu^{-}#bar{#nu}_{#mu} #sigma_{T}", ]
            self.processDict_sigSplit['sig_A'] = [ROOT.kYellow+2, "W^{+}#rightarrow #mu^{+}#nu_{#mu} #sigma_{A}",  ROOT.kYellow+2, "W^{-}#rightarrow #mu^{-}#bar{#nu}_{#mu} #sigma_{A}", ]
            self.processDict_sigSplit['sig_P'] = [ROOT.kOrange-3, "W^{+}#rightarrow #mu^{+}#nu_{#mu} #sigma_{P}",  ROOT.kOrange-3, "W^{-}#rightarrow #mu^{-}#bar{#nu}_{#mu} #sigma_{P}", ]
            self.processDict_sigSplit['sig_UL'] = [ROOT.kRed+1, "W^{+}#rightarrow #mu^{+}#nu_{#mu} #sigma_{UL}",  ROOT.kRed+1, "W^{-}#rightarrow #mu^{-}#bar{#nu}_{#mu} #sigma_{UL}", ]
            self.processOrder_sigSplit = copy.deepcopy(self.processOrder)
            self.processOrder_sigSplit.remove('sig')
            self.processOrder_sigSplit.extend(['sig_L','sig_I','sig_T','sig_A','sig_P','sig_UL'])
            
            self.processDict = self.processDict_sigSplit
            self.processOrder  = self.processOrder_sigSplit
             
                
    
        #unrolled binning: qt large, y fine
        self.unrolledQtY= list(self.yArr)
        intervalYBin = []
        for y in self.yArr[:-1] :
            intervalYBin.append(self.yArr[self.yArr.index(y)+1]-self.yArr[self.yArr.index(y)])
        for q in range(len(self.qtArr)-2) :
            for y in intervalYBin :
                self.unrolledQtY.append(self.unrolledQtY[-1]+y)
        
        #unrolled binning: y large, qt fine        
        self.unrolledYQt= list(self.qtArr)
        intervalQtBin = []
        for q in self.qtArr[:-1] :
            intervalQtBin.append(self.qtArr[self.qtArr.index(q)+1]-self.qtArr[self.qtArr.index(q)])
        for y in range(len(self.yArr)-2) :
            for q in intervalQtBin :
                self.unrolledYQt.append(self.unrolledYQt[-1]+q)
        
        #unrolled binning: coeff large, y fine
        self.unrolledAY= list(self.yArr)
        intervalYBin = []
        for y in self.yArr[:-1] :
            intervalYBin.append(self.yArr[self.yArr.index(y)+1]-self.yArr[self.yArr.index(y)])
        for c in range(len(self.coeffArr)-2) :
            for y in intervalYBin :
                self.unrolledAY.append(self.unrolledAY[-1]+y)
        
        #unrolled binning: coeff large, qt fine
        self.unrolledAQt= list(self.qtArr)
        intervalQtBin = []
        for q in self.qtArr[:-1] :
            intervalQtBin.append(self.qtArr[self.qtArr.index(q)+1]-self.qtArr[self.qtArr.index(q)])
        for c in range(len(self.coeffArr)-2) :
            for q in intervalQtBin :
                self.unrolledAQt.append(self.unrolledAQt[-1]+q)
          
        #unrolled binning: eta large, pt fine
        self.unrolledEtaPt= list(self.ptArr)
        intervalPtBin = []
        for q in self.ptArr[:-1] :
            intervalPtBin.append(self.ptArr[self.ptArr.index(q)+1]-self.ptArr[self.ptArr.index(q)])
        for y in range(len(self.etaArr)-2) :
            for q in intervalPtBin :
                self.unrolledEtaPt.append(self.unrolledEtaPt[-1]+q)
        
        # print "unrolledQtY", self.unrolledQtY
        # print "unrolledYQt", self.unrolledYQt
        # print "unrolledAY", self.unrolledAY
        # print "unrolledAQt", self.unrolledAQt
    
    def CloneEmpty(self, histo,name) :
        hout = histo.Clone(name)
        for i in range(0, histo.GetNbinsX()+2) :
            hout.SetBinContent(i,0)
            hout.SetBinError(i,0)
        return hout     
        
        
    def getHistos(self,inFile, FitFile, acceptance, uncorrelate,suff,apoFile='',toyFile='',impact=False,postfit=False,data='') :
        
        acceptanceFile = ROOT.TFile.Open(acceptance)
        if self.sign =='plus' :
            self.histos[suff+'acceptance'] = acceptanceFile.Get('closPlus')
            self.histos[suff+'acceptance'+'sumOfTemplate'] = acceptanceFile.Get('sumOfTemplatePlus')
        elif self.sign =='minus' : 
            self.histos[suff+'acceptance'] = acceptanceFile.Get('closMinus')
            self.histos[suff+'acceptance'+'sumOfTemplate'] = acceptanceFile.Get('sumOfTemplateMinus')
        
        resFit = FitFile.fitresults
        genMod=''#_gen
                
        #MC (Pre Fit) angular coefficients and variations
        self.histos[suff+'MC'+'mapTot'] =  inFile.Get('angularCoefficients/mapTot')
        self.histos[suff+'MCy'+'mapTot'] =  inFile.Get('angularCoefficients/Y')
        self.histos[suff+'MCqt'+'mapTot'] =  inFile.Get('angularCoefficients/Pt')
        self.histos[suff+'MCy'+'acceptance'] =  self.histos[suff+'MCy'+'mapTot'].Clone("acceptance"+'MCy')#to evaluate acceptance
        self.histos[suff+'MCqt'+'acceptance'] =  self.histos[suff+'MCqt'+'mapTot'].Clone("acceptance"+'MCqt')#to evaluate acceptance
        
        for cat in self.category : #rescale maptot by lumi because mapTot=sigma_UL*Lumi
            self.histos[suff+'MC'+cat+'mapTot'].Scale(1/self.lumi)
        
        if self.anaKind['name'] == 'normXsec' :
            self.sigmaTot = {}
            self.sigmaTot[suff] = self.histos[suff+'MCy'+'mapTot'].Integral()
                    
        for coeff,div in self.coeffDict.items() :
            for cat,catName in self.category.items():
                if coeff!='UL' and coeff!='unpolarizedxsec' :
                    self.histos[suff+'MC'+cat+coeff] =  inFile.Get('angularCoefficients/harmonics'+catName+div[1]+'_nom_nom')
                    self.histos[suff+'MC'+cat+coeff].SetName(self.histos[suff+'MC'+cat+coeff].GetName().replace(div[1],coeff) )
                else :
                    self.histos[suff+'MC'+cat+coeff] = self.histos[suff+'MC'+cat+'mapTot'].Clone(self.histos[suff+'MC'+cat+'mapTot'].GetName()+coeff) 
    
                if ( (not self.anaKind['differential']) and coeff!='UL') : #sigma_hel = A_i*sigma_UL/const_i
                    self.histos[suff+'MC'+cat+coeff].Multiply(self.histos[suff+'MC'+cat+'mapTot'])
                    self.histos[suff+'MC'+cat+coeff].Scale(1/div[0])
                
                if self.anaKind['name'] == 'normXsec' :
                    self.histos[suff+'MC'+cat+coeff].Scale(1/self.sigmaTot[suff])

                
        for sKind, sList in self.systDict.items():

            sListMod = copy.deepcopy(sList)
            if sKind=='_LHEScaleWeight' and UNCORR :
                sListMod.append("_nom") #add nominal variation
            
            for sName in sListMod :
                self.histos[suff+'MC'+sName+'mapTot'] =  inFile.Get('angularCoefficients'+sKind+'/mapTot'+sName)
                self.histos[suff+'MCy'+sName+'mapTot'] =  inFile.Get('angularCoefficients'+sKind+'/Y'+sName)
                self.histos[suff+'MCqt'+sName+'mapTot'] =  inFile.Get('angularCoefficients'+sKind+'/Pt'+sName)
                
                for cat in self.category : #rescale maptot by lumi because mapTot=sigma_UL*Lumi
                    if sName=='_nom' : continue
                    self.histos[suff+'MC'+cat+sName+'mapTot'].Scale(1/self.lumi)
                    
                    if self.anaKind['name'] == 'normXsec' :
                        self.sigmaTot[suff+sName] = self.histos[suff+'MCy'+sName+'mapTot'].Integral()
                        
                for sNameDen in sListMod :
                    if sNameDen!=sName and not (sKind=='_LHEScaleWeight' and UNCORR and self.anaKind['angNames']): #PDF or correlated Scale
                        continue 
                    if sName=='_nom' and sNameDen=='_nom' : continue
                    for coeff,div in self.coeffDict.items() :
                        # if "unpol" in coeff: continue
                        for cat,catName in  self.category.items():
                            if UNCORR and sKind=='_LHEScaleWeight':
                
                                if coeff!='UL' and coeff!='unpolarizedxsec' :
                                    self.histos[suff+'MC'+cat+sName+sNameDen+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonics'+catName+div[1]+sName+sNameDen)
                                    self.histos[suff+'MC'+cat+sName+sNameDen+coeff].SetName(self.histos[suff+'MC'+cat+sName+sNameDen+coeff].GetName().replace(div[1],coeff) )
                                else :
                                    if sName=='_nom' or sNameDen=='_nom': continue
                                    self.histos[suff+'MC'+cat+sName+sNameDen+coeff] = self.histos[suff+'MC'+cat+sName+'mapTot'].Clone(self.histos[suff+'MC'+cat+sName+'mapTot'].GetName()+coeff)
                                
                                if ( (not self.anaKind['differential']) and coeff!='UL') : #sigma_hel = A_i*sigma_UL/const_i
                                    self.histos[suff+'MC'+cat+sName+sNameDen+coeff].Multiply(self.histos[suff+'MC'+cat+sName+'mapTot'])
                                    self.histos[suff+'MC'+cat+sName+sNameDen+coeff].Scale(1/div[0])
                                
                                if self.anaKind['name'] == 'normXsec' :
                                    self.histos[suff+'MC'+cat+sName+sNameDen+coeff].Scale(1/self.sigmaTot[suff+sName])

                                # self.histos[suff+'MC'+cat+sName+sNameDen+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonics'+catName+div[1]+sName+sNameDen)
                                # if self.helXsec :
                                #     if coeff!='UL' :
                                #         self.histos[suff+'MC'+cat+sName+sNameDen+coeff].Multiply(self.histos[suff+'MC'+cat+sName+'mapTot'])
                                #         self.histos[suff+'MC'+cat+sName+sNameDen+coeff].Scale(1/div[0])
                                #     else :
                                #         self.histos[suff+'MC'+cat+sName+sNameDen+coeff] = self.histos[suff+'MC'+cat+sName+'mapTot'].Clone(self.histos[suff+'MC'+cat+sName+'mapTot'].GetName()+'UL')
                            else :
                                    
                                if coeff!='UL' and coeff!='unpolarizedxsec' :
                                    self.histos[suff+'MC'+cat+sName+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonics'+catName+div[1]+sName+sNameDen)
                                    self.histos[suff+'MC'+cat+sName+coeff].SetName(self.histos[suff+'MC'+cat+sName+coeff].GetName().replace(div[1],coeff) )
                                else :
                                    if sName=='_nom' or sNameDen=='_nom': continue
                                    self.histos[suff+'MC'+cat+sName+coeff] = self.histos[suff+'MC'+cat+sName+'mapTot'].Clone(self.histos[suff+'MC'+cat+sName+'mapTot'].GetName()+coeff)
                                
                                if ( (not self.anaKind['differential']) and coeff!='UL') :#sigma_hel = A_i*sigma_UL/const_i
                                    self.histos[suff+'MC'+cat+sName+coeff].Multiply(self.histos[suff+'MC'+cat+sName+'mapTot'])
                                    self.histos[suff+'MC'+cat+sName+coeff].Scale(1/div[0])
                                
                                if self.anaKind['name'] == 'normXsec' :
                                    self.histos[suff+'MC'+cat+sName+coeff].Scale(1/self.sigmaTot[suff+sName]) 
                                    # self.histos[suff+'MC'+cat+sName+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonics'+catName+div[1]+sName+sNameDen)
                                    # if self.helXsec :
                                    #     if coeff!='UL' :
                                    #         self.histos[suff+'MC'+cat+sName+coeff].Multiply(self.histos[suff+'MC'+cat+sName+'mapTot'])
                                    #         self.histos[suff+'MC'+cat+sName+coeff].Scale(1/div[0])
                                    #     else :
                                    #         self.histos[suff+'MC'+cat+sName+coeff] = self.histos[suff+'MC'+cat+sName+'mapTot'].Clone(self.histos[suff+'MC'+cat+sName+'mapTot'].GetName()+'UL')
                            # else :
                                
                            #     self.histos[suff+'MC'+cat+sName+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonics'+catName+div[1]+sName) 
                            #     if self.helXsec :
                            #         if coeff!='UL' :
                            #             self.histos[suff+'MC'+cat+sName+coeff].Multiply(self.histos[suff+'MC'+cat+sName+'mapTot'])
                            #             self.histos[suff+'MC'+cat+sName+coeff].Scale(1/div[0])
                            #         else :
                            #             self.histos[suff+'MC'+cat+sName+coeff] = self.histos[suff+'MC'+cat+sName+'mapTot'].Clone(self.histos[suff+'MC'+cat+sName+'mapTot'].GetName()+'UL')
        
         
        #fit - helicity xsecs histos  (included in helXsec option)
        # for hel in self.hels:
        #     self.histos[suff+'FitHel'+hel] = ROOT.TH2D(suff+'FitHel{c}'.format(c=hel), suff+'FitHel{c}'.format(c=hel), len(self.yArr)-1, array('f',self.yArr), len(self.qtArr)-1, array('f',self.qtArr))
        #     self.histos[suff+'FitHel'+hel+'norm'] = ROOT.TH2D(suff+'FitHel{c}norm'.format(c=hel), suff+'FitHel{c}norm'.format(c=hel), len(self.yArr)-1, array('f',self.yArr), len(self.qtArr)-1, array('f',self.qtArr))
        #     self.histos[suff+'FitHel'+hel].GetXaxis().SetTitle('W rapidity')
        #     self.histos[suff+'FitHel'+hel].GetYaxis().SetTitle('W q_{T}')
        #     self.histos[suff+'FitHel'+hel+'norm'].GetXaxis().SetTitle('W rapidity')
        #     self.histos[suff+'FitHel'+hel+'norm'].GetYaxis().SetTitle('W q_{T}')
        #     for ev in resFit: #dummy because there's one event only
        #         for i in range(1, self.histos[suff+'FitHel'+hel].GetNbinsX()+1): #loop over rapidity bins
        #             for j in range(1, self.histos[suff+'FitHel'+hel].GetNbinsY()+1): #loop over pt bins
        #                 try:
        #                     coeff = eval('ev.helXsecs{}_y_{}_qt_{}_pmaskedexp{}'.format(hel, i, j,genMod))
        #                     coeff_err = eval('ev.helXsecs{}_y_{}_qt_{}_pmaskedexp_err'.format(hel, i, j))
        #                     self.histos[suff+'FitHel'+hel].SetBinContent(i,j,coeff)
        #                     self.histos[suff+'FitHel'+hel].SetBinError(i,j,coeff_err)

        #                     coeffnorm = eval('ev.helXsecs{}_y_{}_qt_{}_pmaskedexpnorm{}'.format(hel, i, j,genMod))
        #                     coeffnorm_err = eval('ev.helXsecs{}_y_{}_qt_{}_pmaskedexpnorm_err'.format(hel, i, j))
        #                     self.histos[suff+'FitHel'+hel+'norm'].SetBinContent(i,j,coeffnorm)
        #                     self.histos[suff+'FitHel'+hel+'norm'].SetBinError(i,j,coeffnorm_err)
                            
        #                 except AttributeError:
        #                     pass
                
        if self.anaKind['name']=='mu' : #normalize to 1 all the MC histos for mu
            for h, hval in self.histos.items() :
                try : 
                    if not (type(hval)==ROOT.TH1F or type(hval)==ROOT.TH1D or type(hval)==ROOT.TH2F or type(hval)==ROOT.TH2D) : continue
                    if type(hval)==ROOT.TH1F or type(hval)==ROOT.TH1D :
                        for bb in range(1,self.histos[h].GetNbinsX()+1) :
                            # if 'MCy' in h or 'MCqt' in h : continue 
                            self.histos[h].SetBinContent(bb,1)
                            self.histos[h].SetBinError(bb,0)
                    else :
                        for bb in range(1,self.histos[h].GetNbinsX()+1) :
                            for b2 in range(1,self.histos[h].GetNbinsY()+1) :
                                self.histos[h].SetBinContent(bb,b2,1)
                                self.histos[h].SetBinError(bb,b2,0)
                except :
                    print("histo issue", h)
                    continue 
        
        #fit - angular coefficients with band
        for c in self.coeffDict:
            self.histos[suff+'FitAC'+c] = ROOT.TH2D(suff+'FitAC{c}'.format(c=c), suff+'FitAC{c}'.format(c=c), len(self.yArr)-1, array('f',self.yArr), len(self.qtArr)-1, array('f',self.qtArr))
            self.histos[suff+'FitACqt'+c] = ROOT.TH1D(suff+'FitACqt{c}'.format(c=c), suff+'FitACqt{c}'.format(c=c), len(self.qtArr)-1, array('f',self.qtArr))
            self.histos[suff+'FitACy'+c] = ROOT.TH1D(suff+'FitACy{c}'.format(c=c), suff+'FitACy{c}'.format(c=c), len(self.yArr)-1, array('f',self.yArr))
            for ev in resFit: #dummy because there's one event only
                if eval('ev.status')!= 0 :
                        print("status 0")
                        continue #not fully-converged fit
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                    for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                        try:
                            if self.anaKind['angNames'] : coeffString = 'ev.y_{i}_qt_{j}_{c}{g}'.format(c=c, j=j, i=i,g=genMod)
                            else : coeffString = 'ev.helXsecs{c}_y_{i}_qt_{j}_{kind}{g}'.format(c=c, i=i, j=j,kind=self.anaKind['diffString'],g=genMod)
                            # else : coeffString = 'ev.helXsecs{c}_y_{i}_qt_{j}_mu{g}'.format(c=c, i=i, j=j,g=genMod)
                            coeff = eval(coeffString)
                            coeff_err = eval(coeffString+'_err') #doesn't work fit genMod!=0
                            if ('unpol' in c or (not self.anaKind['differential'])):
                                coeff = coeff/(3./16./math.pi)/self.lumi
                                coeff_err = coeff_err/(3./16./math.pi)/self.lumi
                            self.histos[suff+'FitAC'+c].SetBinContent(i,j,coeff)
                            self.histos[suff+'FitAC'+c].SetBinError(i,j,coeff_err)
                        except AttributeError: 
                            pass
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                    try:
                        if self.anaKind['angNames'] : coeffString = 'ev.qt_{j}_helmeta_{c}{g}'.format(c=c, j=j, g=genMod)
                        else : coeffString = 'ev.helXsecs{c}_qt_{j}_{kind}{g}'.format(c=c, j=j, kind=self.anaKind['intString'].replace('sumpois','sumxsec'), g=genMod)
                        coeff = eval(coeffString)
                        coeff_err = eval(coeffString+'_err') #doesn't work fit genMod!=0
                        if ('unpol' in c or (not self.anaKind['differential'])):
                                coeff = coeff/(3./16./math.pi)/self.lumi
                                coeff_err = coeff_err/(3./16./math.pi)/self.lumi
                        self.histos[suff+'FitACqt'+c].SetBinContent(j,coeff)
                        self.histos[suff+'FitACqt'+c].SetBinError(j,coeff_err)
                    except AttributeError: 
                        pass
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                    try:
                        if self.anaKind['angNames'] : coeffString = 'ev.y_{i}_helmeta_{c}{g}'.format(c=c, i=i,  g=genMod)
                        else : coeffString = 'ev.helXsecs{c}_y_{i}_{kind}{g}'.format(c=c, i=i,kind=self.anaKind['intString'].replace('sumpois','sumxsec'), g=genMod)
                        coeff = eval(coeffString)
                        coeff_err = eval(coeffString+'_err') #doesn't work fit genMod!=0
                        if ('unpol' in c or (not self.anaKind['differential'])):
                                coeff = coeff/(3./16./math.pi)/self.lumi
                                coeff_err = coeff_err/(3./16./math.pi)/self.lumi
                        self.histos[suff+'FitACy'+c].SetBinContent(i,coeff)
                        self.histos[suff+'FitACy'+c].SetBinError(i,coeff_err)
                    except AttributeError: 
                        pass
        
        if self.anaKind['angNames'] : 
            for cat in self.category :
                self.histos[suff+'FitAC'+cat+'unpolarizedxsec'+'4apo'] = self.histos[suff+'FitAC'+cat+'unpolarizedxsec'].Clone(self.histos[suff+'FitAC'+cat+'unpolarizedxsec'].GetName()+'_4apo')
                self.histos[suff+'FitAC'+cat+'unpolarizedxsec'+'4apo'].Scale((3./16./math.pi)*self.lumi)
        
        covString = [self.anaKind['diffString'],self.anaKind['intString']]
            
        #covariance and correlation matrices 
        self.histos[suff+'corrMat'] = FitFile.Get('correlation_matrix_channel'+covString[0])
        self.histos[suff+'covMat'] = FitFile.Get('covariance_matrix_channel'+covString[0])  
        self.histos[suff+'corrMat'+'Integrated'] = FitFile.Get('correlation_matrix_channel'+covString[1])
        self.histos[suff+'covMat'+'Integrated'] = FitFile.Get('covariance_matrix_channel'+covString[1])
        
        #mass
        self.histos[suff+'mass'] = ROOT.TH1F('mass'+suff,'mass'+suff,1,0,1)
        for ev in resFit: #dummy because there's one event only
            if eval('ev.status')!= 0 :
                        print("status 0")
                        continue #not fully-converged fit
            try :
                massVal = eval('ev.mass{}'.format(genMod))
                massErr = eval('ev.mass_err')
            except :
                print("missing mass in tree")
                massVal = 0
                massErr= 0
        self.histos[suff+'mass'].SetBinContent(1,massVal)
        self.histos[suff+'mass'].SetBinError(1,massErr)
        # print("WARNING: mass =0., error converted in GeV")
        
        if apoFile!='' :
            for c in self.coeffDict:   
                self.histos[suff+'apo'+c] = apoFile.Get('post-fit-regularization_'+c)
                self.histos[suff+'apo'+'qt'+c] = apoFile.Get('post-fit-regularization_'+c+'_qt')
                self.histos[suff+'apo'+'y'+c] = apoFile.Get('post-fit-regularization_'+c+'_y')
                
        if impact :
            self.histos[suff+'impact2D'+'group'+'UNR']=FitFile.Get('nuisance_group_impact_'+covString[0])
            self.histos[suff+'impact2D'+'UNR']=FitFile.Get('nuisance_impact_'+covString[0])
            self.histos[suff+'impact2D'+'group'+'int']=FitFile.Get('nuisance_group_impact_'+covString[1])
            try :
                self.histos[suff+'impact2D'+'group'+'mass']=FitFile.Get('nuisance_group_impact_nois')
            except :
                print("missing mass in nuisance")
                self.histos[suff+'impact2D'+'group'+'mass'] = 0 
        
        if postfit :
            for p,pval in self.processDict.items():
                if 'obs' in p  :
                    if data!='' :
                        fileData = ROOT.TFile.Open(data)
                        h2Data = fileData.Get("Data")
                        self.histos[suff+'prefit'+p] = FitFile.Get('obs')
                        for ee in range(1,len(self.etaArr)) : 
                            for pp in range(1,len(self.ptArr)) :
                                indexUNR = (ee-1)*(len(self.ptArr)-1)+(pp-1)
                                # print(ee,pp, indexUNR)
                                self.histos[suff+'prefit'+p].SetBinContent(indexUNR+1, h2Data.GetBinContent(ee,pp))
                                self.histos[suff+'prefit'+p].SetBinError(indexUNR+1, h2Data.GetBinError(ee,pp))
                    else :
                        self.histos[suff+'prefit'+p] = FitFile.Get('obs')
                    self.histos[suff+'postfit'+p] = FitFile.Get('obs')
                else :
                    if p in ["DYJets", "DiBoson", "Top", "Fake", "WtoTau", "LowAcc"] :
                        procS = 'proc_'
                    else : procS = ''
                    if not p.replace('sig_','') in self.hxsecList :
                        self.histos[suff+'prefit'+p] = FitFile.Get('exp'+procS+p+'_prefit')  
                        self.histos[suff+'postfit'+p] = FitFile.Get('exp'+procS+p+'_postfit')  
                    else : # SPLITSIG
                        # for hxsec in self.hxsecList :
                        hrefdy = FitFile.Get('exp'+ 'proc_'+'DYJets'+'_prefit')  
                        self.histos[suff+'prefit'+p] = self.CloneEmpty(hrefdy, 'prefit'+p)
                        self.histos[suff+'postfit'+p] = self.CloneEmpty(hrefdy, 'postfit'+p)
                        for y in range(1,len(self.yArr)) :
                            for qt in range(1,len(self.qtArr)) :
                                fullprocname = 'expproc_helXsecs'+p.replace('sig_','')+'_y_'+str(y)+'_qt_'+str(qt)
                                htemp_prefit = FitFile.Get(fullprocname+'_prefit')  
                                htemp_postfit = FitFile.Get(fullprocname+'_postfit')  
                                self.histos[suff+'prefit'+p].Add(htemp_prefit)
                                self.histos[suff+'postfit'+p].Add(htemp_postfit)
                                        
                    #removal of the last bin in the unrolled plots
                for k in ['prefit','postfit'] :
                    # print(k,p,  self.histos[suff+k+p])
                    hpostfittemp = self.histos[suff+k+p].Clone()
                    self.histos[suff+k+p] = ROOT.TH1F(self.histos[suff+k+p].GetName(), self.histos[suff+k+p].GetTitle(), len(self.unrolledEtaPt)-1, array('f',self.unrolledEtaPt))
                    for ietapt in range(1,self.histos[suff+k+p].GetNbinsX()+1) :
                        self.histos[suff+k+p].SetBinContent(ietapt, hpostfittemp.GetBinContent(ietapt))
                        self.histos[suff+k+p].SetBinError(ietapt, hpostfittemp.GetBinError(ietapt))

        
        #nuisance plots
        NbinsNui = 0
        self.NuiConstrDict = {
            'all' : [0, [],0, [0], [True,True,True]], #[Nbin,[sName],counterBin,[binwidth],[Flag_pdf,Flag_SFstat,Flag_Scale]]
            'LHEPdfWeight' : [0, [],0], 
            'WHSFStat' : [0, [],0],
            'LHEScaleWeight' : [0, [],0],
            'others' : [0, [],0]
        }
        
        for gk,group in systematicsDict.items() :
            if gk=='Nominal' : continue
            flagOthers=True
            for nuiDict_key, nuiDict_val in self.NuiConstrDict.items() :
                if nuiDict_key in gk or (nuiDict_key=='LHEPdfWeight' and gk=='alphaS'):
                    flagOthers=False
                    for sName in group['vars'] :
                        self.NuiConstrDict[nuiDict_key][0] +=1
                        self.NuiConstrDict[nuiDict_key][1].append(sName)
            for sName in group['vars'] :
                self.NuiConstrDict['all'][0] +=1
                self.NuiConstrDict['all'][1].append(sName)
                if 'LHEPdfWeight' not in gk and 'WHSFStat' not in gk and 'LHEScaleWeight' not in gk :
                    binW = 1
                else :
                    binW = 0.15
                self.NuiConstrDict['all'][3].append(self.NuiConstrDict['all'][3][-1]+binW)
            if flagOthers :
                for sName in group['vars'] :
                    self.NuiConstrDict['others'][0]+=1
                    self.NuiConstrDict['others'][1].append(sName)
                # for sName in group['vars'] :
                # NbinsNui+=1 
        for nuiDict_key, nuiDict_val in self.NuiConstrDict.items() : 
            if 'all' not in nuiDict_key :
                self.histos[suff+'NuiConstr'+nuiDict_key] = ROOT.TH1F(suff+'_NuiConstr'+nuiDict_key,suff+'_NuiConstr'+nuiDict_key,nuiDict_val[0],0,nuiDict_val[0])
            else :
                self.histos[suff+'NuiConstr'+nuiDict_key] = ROOT.TH1F(suff+'_NuiConstr'+nuiDict_key,suff+'_NuiConstr'+nuiDict_key,nuiDict_val[0],array('f',self.NuiConstrDict['all'][3]))
            self.histos[suff+'NuiConstr'+nuiDict_key].GetXaxis().SetNdivisions(-1)
        # nuiBinCount=0
        for ev in resFit: #dummy, 1ev
            if eval('ev.status')!= 0 :
                        print("status 0")
                        continue #not fully-converged fit
            for gk,group in systematicsDict.items() :
                if gk=='Nominal' : continue
                for sName in group['vars'] :
                    for nuiDict_key, nuiDict_val in self.NuiConstrDict.items() : 
                        if sName in nuiDict_val[1] :
                            self.NuiConstrDict[nuiDict_key][2]+=1
                        # nuiBinCount+=1
                            if 'all' in nuiDict_key and (('LHEPdfWeight' in gk and self.NuiConstrDict['all'][4][0]) or ('WHSFStat' in gk and self.NuiConstrDict['all'][4][1]) or ('LHEScaleWeight' in gk and self.NuiConstrDict['all'][4][2])):
                                if 'LHEPdfWeight' in gk : 
                                    self.NuiConstrDict['all'][4][0]=False
                                    shift = 0.5*float(len(self.NuiConstrDict['LHEPdfWeight'][1]))
                                if 'WHSFStat' in gk :
                                    self.NuiConstrDict['all'][4][1]=False
                                    shift = 0.5*float(len(self.NuiConstrDict['WHSFStat'][1]))
                                if 'LHEScaleWeight' in gk : 
                                    self.NuiConstrDict['all'][4][2]=False
                                    shift = 0.5*float(len(self.NuiConstrDict['LHEScaleWeight'][1]))                                
                                self.histos[suff+'NuiConstr'+nuiDict_key].GetXaxis().SetBinLabel(self.NuiConstrDict[nuiDict_key][2]+int(shift),gk)
                            elif 'all' not in nuiDict_key or ('LHEPdfWeight' not in gk and 'WHSFStat' not in gk and 'LHEScaleWeight' not in gk):
                                self.histos[suff+'NuiConstr'+nuiDict_key].GetXaxis().SetBinLabel(self.NuiConstrDict[nuiDict_key][2],sName)
                            try : 
                                nuiVal = eval('ev.{}'.format(sName))
                                nuiErr = eval('ev.{}_err'.format(sName))
                            except :
                                continue 
                            if 'mass' in sName :
                                print("w mass nuisance:", nuiVal,"+/-", nuiErr)
                            self.histos[suff+'NuiConstr'+nuiDict_key].SetBinContent(self.NuiConstrDict[nuiDict_key][2],nuiVal)
                            self.histos[suff+'NuiConstr'+nuiDict_key].SetBinError(self.NuiConstrDict[nuiDict_key][2],nuiErr)
        
        # if toyFile!= '' and self.helXsec : print("toys not implemented for helicty xsection")
        if toyFile!= '':
            print("start toy analysis")
            LoopToyTree = False 
            try :
                resFitToy = toyFile.fitresults
                LoopToyTree = True
            except :
                print('external toy histograms used (pre-looped)')
            
            
            self.histos[suff+'mass'+'toy'] = ROOT.TH1F(suff+'mass'+'toy',suff+'mass'+'toy',100,0,-1)
            self.histos[suff+'mass'+'toyPull'] = ROOT.TH1F(suff+'mass'+'toyPull',suff+'mass'+'toyPull',100,0,-1)
            self.histos[suff+'mass'+'toy'+'mean'] = ROOT.TH1F(suff+'mass'+'toy'+'mean',suff+'mass'+'toy'+'mean',1,0,1)
            self.histos[suff+'mass'+'toyPull'+'mean'] = ROOT.TH1F(suff+'mass'+'toyPull'+'mean',suff+'mass'+'toyPull'+'mean',1,0,1)
            
            for c in self.coeffDict:
                self.histos[suff+'FitAC'+c+'toy'] = ROOT.TH2D(suff+'FitAC{c}_toy'.format(c=c), suff+'FitAC{c}_toy'.format(c=c), len(self.yArr)-1, array('f',self.yArr), len(self.qtArr)-1, array('f',self.qtArr))
                self.histos[suff+'FitACqt'+c+'toy'] = ROOT.TH1D(suff+'FitACqt{c}_toy'.format(c=c), suff+'FitACqt{c}_toy'.format(c=c), len(self.qtArr)-1, array('f',self.qtArr))
                self.histos[suff+'FitACy'+c+'toy'] = ROOT.TH1D(suff+'FitACy{c}_toy'.format(c=c), suff+'FitACy{c}_toy'.format(c=c), len(self.yArr)-1, array('f',self.yArr))
                self.histos[suff+'FitAC'+c+'toyPull'] = ROOT.TH2D(suff+'FitAC{c}_toyPull'.format(c=c), suff+'FitAC{c}_toyPull'.format(c=c), len(self.yArr)-1, array('f',self.yArr), len(self.qtArr)-1, array('f',self.qtArr))
                self.histos[suff+'FitACqt'+c+'toyPull'] = ROOT.TH1D(suff+'FitACqt{c}_toyPull'.format(c=c), suff+'FitACqt{c}_toyPull'.format(c=c), len(self.qtArr)-1, array('f',self.qtArr))
                self.histos[suff+'FitACy'+c+'toyPull'] = ROOT.TH1D(suff+'FitACy{c}_toyPull'.format(c=c), suff+'FitACy{c}_toyPull'.format(c=c), len(self.yArr)-1, array('f',self.yArr))
                
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                    for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                        self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)] = ROOT.TH1D(suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j),suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j),100,0,-1)
                        self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)] = ROOT.TH1D(suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j),suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j),100,0,-1)
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                    self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)] = ROOT.TH1D(suff+'FitACqt'+c+'toy'+'qt'+str(j),suff+'FitACqt'+c+'toy'+'qt'+str(j),100,0,-1)
                    self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)] = ROOT.TH1D(suff+'FitACqt'+c+'toyPull'+'qt'+str(j),suff+'FitACqt'+c+'toyPull'+'qt'+str(j),100,0,-1)
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                    self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)] = ROOT.TH1D(suff+'FitACy'+c+'toy'+'y'+str(i),suff+'FitACy'+c+'toy'+'y'+str(i),100,0,-1)
                    self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)] = ROOT.TH1D(suff+'FitACy'+c+'toyPull'+'y'+str(i),suff+'FitACy'+c+'toyPull'+'y'+str(i),100,0,-1)
            
            if LoopToyTree :
                ievent =0     
                for ev in resFitToy:
                    ievent +=1
                    if eval('ev.status')!= 0 :
                        continue #not fully-converged fit
                    print("new event",ievent)
                    mtoy = eval('ev.mass')
                    mtoyPull = eval('(ev.mass-ev.mass_gen)/ev.mass_err')
                    self.histos[suff+'mass'+'toy'].Fill(mtoy)
                    self.histos[suff+'mass'+'toyPull'].Fill(mtoyPull)
                    
                    for c in self.coeffDict:
                        for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): 
                            for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): 
                                if self.anaKind['angNames'] : binName = 'y_'+str(i)+'_qt_'+str(j)+'_'+c
                                else : binName  = 'ev.helXsecs'+c+'_y_'+str(i)+'_qt_'+str(j)+'_'+self.anaKind['diffString']
                                unpolMult = '(3./16./3.1415926535)/'+str(self.lumi)
                                # print("double loop", c,i,j)
                                if ('unpol' in c or (not self.anaKind['differential'])) : 
                                    vtoy = eval('ev.'+binName+'/'+unpolMult)
                                else :
                                    vtoy = eval('ev.'+binName)
                                vtoyPull = eval('(ev.'+binName+'-ev.'+binName+'_gen'+')/ev.'+binName+'_err')  
                                self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)].Fill(vtoy)    
                                self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)].Fill(vtoyPull)    
                
                        for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                            if self.anaKind['angNames'] : binName = 'qt_'+str(j)+'_helmeta_'+c
                            else : binName  = 'ev.helXsecs'+c+'_qt_'+str(j)+'_'+self.anaKind['intString'].replace('sumpois','sumxsec')
                            unpolMult = '(3./16./3.1415926535)/'+str(self.lumi)
                            # print("qt loop", c,j)
                            
                            if ('unpol' in c or (not self.anaKind['differential'])) :
                                vtoy = eval('ev.'+binName+'/'+unpolMult)
                            else :
                                vtoy = eval('ev.'+binName)
                            vtoyPull = eval('(ev.'+binName+'-ev.'+binName+'_gen'+')/ev.'+binName+'_err')  
                            self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)].Fill(vtoy)    
                            self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)].Fill(vtoyPull)    
            
                        for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                            if self.anaKind['angNames'] : binName = 'y_'+str(i)+'_helmeta_'+c
                            else : binName  = 'ev.helXsecs'+c+'_y_'+str(i)+'_'+self.anaKind['diffString'].replace('sumpois','sumxsec')
                            unpolMult = '(3./16./3.1415926535)/'+str(self.lumi)
                            # print("y loop", c,i)
                            
                            if ('unpol' in c or (not self.anaKind['differential'])) : 
                                vtoy = eval('ev.'+binName+'/'+unpolMult)
                            else :
                                vtoy = eval('ev.'+binName)
                            vtoyPull = eval('(ev.'+binName+'-ev.'+binName+'_gen'+')/ev.'+binName+'_err')  
                            self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)].Fill(vtoy)    
                            self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)].Fill(vtoyPull)    
                
            else : #toy distirbutions from external file
                self.histos[suff+'mass'+'toy'] = toyFile.Get(suff+'mass'+'toy')
                self.histos[suff+'mass'+'toyPull'] = toyFile.Get(suff+'mass'+'toyPull')
                for c in self.coeffDict:
                        for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): 
                            for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): 
                                self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)] = toyFile.Get(suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j))
                                self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)] = toyFile.Get(suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j))
                        for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                            self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)]  = toyFile.Get(suff+'FitACqt'+c+'toy'+'qt'+str(j))
                            self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)] = toyFile.Get(suff+'FitACqt'+c+'toyPull'+'qt'+str(j))
                        for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1):     
                            self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)] = toyFile.Get(suff+'FitACy'+c+'toy'+'y'+str(i))
                            self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)] = toyFile.Get(suff+'FitACy'+c+'toyPull'+'y'+str(i))
                
                    
            
            self.histos[suff+'mass'+'toy'+'mean'].SetBinContent(1,self.histos[suff+'mass'+'toy'].GetMean())
            self.histos[suff+'mass'+'toy'+'mean'].SetBinError(1,self.histos[suff+'mass'+'toy'].GetStdDev())
            self.histos[suff+'mass'+'toyPull'+'mean'].SetBinContent(1,self.histos[suff+'mass'+'toyPull'].GetMean())
            self.histos[suff+'mass'+'toyPull'+'mean'].SetBinError(1,self.histos[suff+'mass'+'toyPull'].GetStdDev())
            for c in self.coeffDict:        
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                    for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                        self.histos[suff+'FitAC'+c+'toy'].SetBinContent(i,j,self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)].GetMean())
                        self.histos[suff+'FitAC'+c+'toy'].SetBinError(i,j,self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)].GetStdDev())
                        self.histos[suff+'FitAC'+c+'toyPull'].SetBinContent(i,j,self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)].GetMean())
                        self.histos[suff+'FitAC'+c+'toyPull'].SetBinError(i,j,self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)].GetStdDev())
                    
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                    self.histos[suff+'FitACqt'+c+'toy'].SetBinContent(j,self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)].GetMean())
                    self.histos[suff+'FitACqt'+c+'toy'].SetBinError(j,self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)].GetStdDev())
                    self.histos[suff+'FitACqt'+c+'toyPull'].SetBinContent(j,self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)].GetMean())
                    self.histos[suff+'FitACqt'+c+'toyPull'].SetBinError(j,self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)].GetStdDev())
                    
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                    self.histos[suff+'FitACy'+c+'toy'].SetBinContent(i,self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)].GetMean())
                    self.histos[suff+'FitACy'+c+'toy'].SetBinError(i,self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)].GetStdDev())
                    self.histos[suff+'FitACy'+c+'toyPull'].SetBinContent(i,self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)].GetMean())
                    self.histos[suff+'FitACy'+c+'toyPull'].SetBinError(i,self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)].GetStdDev())
                    
                    
                    
                    
                    
                    
                # for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                #     for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                #         binName = 'y_'+str(i)+'_qt_'+str(j)+'_'+c
                #         unpolMult = '(3./16./3.1415926535)/'+str(self.lumi)
                #         print(c,i,j)
                        
                #         self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)] = ROOT.TH1D(suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j),suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j),100,0,-1)
                #         self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)] = ROOT.TH1D(suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j),suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j),100,0,-1)
                        
                #         for ev in resFitToy:
                #             if 'unpol' in c: 
                #                 vtoy = eval('ev.'+binName+'/'+unpolMult)
                #             else :
                #                 vtoy = eval('ev.'+binName)
                #             vtoyPull = eval('(ev.'+binName+'-ev.'+binName+'_gen'+')/ev.'+binName+'_err')  
                #             self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)].Fill(vtoy)    
                #             self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)].Fill(vtoyPull)    
                        
                #         self.histos[suff+'FitAC'+c+'toy'].SetBinContent(i,j,self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)].GetMean())
                #         self.histos[suff+'FitAC'+c+'toy'].SetBinError(i,j,self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)].GetStdDev())
                #         self.histos[suff+'FitAC'+c+'toyPull'].SetBinContent(i,j,self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)].GetMean())
                #         self.histos[suff+'FitAC'+c+'toyPull'].SetBinError(i,j,self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)].GetStdDev())
                    
                            
                # for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                #         binName = 'qt_'+str(j)+'_helmeta_'+c
                #         unpolMult = '(3./16./3.1415926535)/'+str(self.lumi)
                #         print(c,j)
                        
                #         self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)] = ROOT.TH1D(suff+'FitACqt'+c+'toy'+'qt'+str(j),suff+'FitACqt'+c+'toy'+'qt'+str(j),100,0,-1)
                #         self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)] = ROOT.TH1D(suff+'FitACqt'+c+'toyPull'+'qt'+str(j),suff+'FitACqt'+c+'toyPull'+'qt'+str(j),100,0,-1)
                        
                #         for ev in resFitToy:
                #             if 'unpol' in c: 
                #                 vtoy = eval('ev.'+binName+'/'+unpolMult)
                #             else :
                #                 vtoy = eval('ev.'+binName)
                #             vtoyPull = eval('(ev.'+binName+'-ev.'+binName+'_gen'+')/ev.'+binName+'_err')  
                #             self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)].Fill(vtoy)    
                #             self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)].Fill(vtoyPull)    
                        
                #         self.histos[suff+'FitACqt'+c+'toy'].SetBinContent(j,self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)].GetMean())
                #         self.histos[suff+'FitACqt'+c+'toy'].SetBinError(j,self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)].GetStdDev())
                #         self.histos[suff+'FitACqt'+c+'toyPull'].SetBinContent(j,self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)].GetMean())
                #         self.histos[suff+'FitACqt'+c+'toyPull'].SetBinError(j,self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)].GetStdDev())
                
            
                # for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                #     binName = 'y_'+str(i)+'_helmeta_'+c
                #     unpolMult = '(3./16./3.1415926535)/'+str(self.lumi)
                #     print(c,i)
                    
                #     self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)] = ROOT.TH1D(suff+'FitACy'+c+'toy'+'y'+str(i),suff+'FitACy'+c+'toy'+'y'+str(i),100,0,-1)
                #     self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)] = ROOT.TH1D(suff+'FitACy'+c+'toyPull'+'y'+str(i),suff+'FitACy'+c+'toyPull'+'y'+str(i),100,0,-1)
                    
                #     for ev in resFitToy:
                #         if 'unpol' in c: 
                #             vtoy = eval('ev.'+binName+'/'+unpolMult)
                #         else :
                #             vtoy = eval('ev.'+binName)
                #         vtoyPull = eval('(ev.'+binName+'-ev.'+binName+'_gen'+')/ev.'+binName+'_err')  
                #         self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)].Fill(vtoy)    
                #         self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)].Fill(vtoyPull)    
                    
                #     self.histos[suff+'FitACy'+c+'toy'].SetBinContent(i,self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)].GetMean())
                #     self.histos[suff+'FitACy'+c+'toy'].SetBinError(i,self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)].GetStdDev())
                #     self.histos[suff+'FitACy'+c+'toyPull'].SetBinContent(i,self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)].GetMean())
                #     self.histos[suff+'FitACy'+c+'toyPull'].SetBinError(i,self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)].GetStdDev())
                
                
                    # toyLim = [-1,1]
                    # toyPullLim = [-1,1]
                    # for ev in resFitToy: #evaluate min e max
                    #     if 'unpol' in c: 
                    #         vtoy = eval('ev.'+binName+'/'+unpolMult)
                    #     else :
                    #         vtoy = eval('ev.'+binName)
                    #     vtoyPull = eval('(ev.'+binName+'-ev.'+binName+'_gen'+')/ev.'+binName+'_err')     
                    #     if vtoy<toyLim[0] : toyLim[0] = vtoy
                    #     if vtoy>toyLim[1] : toyLim[1] = vtoy
                    #     if vtoyPull<toyPullLim[0] : toyPullLim[0] = vtoy
                    #     if vtoyPull>toyPullLim[1] : toyPullLim[1] = vtoy
                    # toyLim[0] = 1.1*toyLim[0]
                    # toyLim[1] = 1.1*toyLim[1]
                    # toyLim[0] = 1.1*toyPullLim[0]
                    # toyPullLim[1] = 1.1*toyPullLim[1]
    
    def varBinWidth_modifier(self,suff) :
        for k,histo in self.histos.items() :
            
            if 'A0' in k or 'A1' in k  or 'A2' in k or 'A3' in k or 'A4' in k: continue 
            # if not 'unpol' in k and not 'mapTot' in k : continue #only unpol and maptot are plotted as absolute value
            # if not self.helXsec and (not 'unpol' in k and not 'mapTot' in k) : continue
            if 'corrMat' in k or 'covMat' in k: continue
            if 'impact' in k : continue 
            if 'NuiConstr' in k : continue
            if 'mass' in k : continue
            if suff+'apo' in k : continue 
            if '4apo' in k : continue 
            if 'toyPull' in k : continue 
            if 'prefit' in k or 'postfit' in k : continue 
            if 'acceptance' in k : continue 
            if self.anaKind['name'] =='mu' : continue
                        
            if not (type(histo)==ROOT.TH1F or type(histo)==ROOT.TH1D or type(histo)==ROOT.TH2F or type(histo)==ROOT.TH2D) : continue
            # print(k)
            
            #the following line rebin only if variable bin width is present in the histogram
            #########################################################################################################################
            # if 'unpol' in k :
            # if not 'MC' in k :
            #     maxBinX = histo.GetNbinsX()+1
            #     maxBinY = histo.GetNbinsY()+1
            # else : #this is needed because these histograms extend to qt=200, Y=6
            #     if 'MCqt' in k : 
            #         maxBinX = len(self.qtArr)
            #     else : 
            #         maxBinX = len(self.yArr)
            #     maxBinY = len(self.qtArr)

            # # if 'MCqt' in k : 
            # #     maxBinX = len(self.qtArr)
            # # else : 
            # #     maxBinX = len(self.yArr)
            # # maxBinY = len(self.qtArr)
            
            
            # if type(histo)==ROOT.TH1F or type(histo)==ROOT.TH1D :
            #     varWidth_X = False
            #     for i in range(1, maxBinX-1) :
            #         for j in range(i+1,maxBinX) :
            #             if abs(histo.GetXaxis().GetBinWidth(j)-histo.GetXaxis().GetBinWidth(i))>0.0001:
            #                 varWidth_X = True
            #                 break
            #     if varWidth_X :
            #         # print("1D, k=", k)
            #         for i in range(1, maxBinX):
            #             # print("pre", histo.GetBinContent(i))
            #             histo.SetBinContent(i, histo.GetBinContent(i)/histo.GetXaxis().GetBinWidth(i))
            #             histo.SetBinError(i, histo.GetBinError(i)/histo.GetXaxis().GetBinWidth(i))
            #             # print("post", histo.GetBinContent(i),histo.GetXaxis().GetBinWidth(i))
                        
            # elif type(histo)==ROOT.TH2F or type(histo)==ROOT.TH2D :
            #     varWidth_X = False
            #     for i in range(1, maxBinX-1) :
            #         for j in range(i+1,maxBinX) :
            #             if abs(histo.GetXaxis().GetBinWidth(j)-histo.GetXaxis().GetBinWidth(i))>0.0001:
            #                 varWidth_X = True
            #                 break
            #     if varWidth_X :
            #         # print("2Dx, k=", k)
            #         for y in range(1, maxBinY) :
            #             for i in range(1, maxBinX):
            #                 # print("pre", histo.GetBinContent(i,y))
            #                 histo.SetBinContent(i,y, histo.GetBinContent(i,y)/histo.GetXaxis().GetBinWidth(i))
            #                 histo.SetBinError(i,y, histo.GetBinError(i,y)/histo.GetXaxis().GetBinWidth(i))
            #                 # print("post", histo.GetBinContent(i,y), histo.GetXaxis().GetBinWidth(i))
            #     varWidth_Y = False
            #     for i in range(1, maxBinY-1) :
            #         for j in range(i+1,maxBinY) :
            #             if abs(histo.GetYaxis().GetBinWidth(j)-histo.GetYaxis().GetBinWidth(i))>0.0001:
            #                 varWidth_Y = True
            #                 break
            #     if varWidth_Y :
            #         # print("2Dy, k=", k)
            #         for j in range(1, maxBinY) :
            #             for x in range(1, maxBinX):
            #                 # print("pre", histo.GetBinContent(x,j))
            #                 histo.SetBinContent(x,j, histo.GetBinContent(x,j)/histo.GetYaxis().GetBinWidth(j))
            #                 histo.SetBinError(x,j, histo.GetBinError(x,j)/histo.GetYaxis().GetBinWidth(j))
            #                 # print("post", histo.GetBinContent(x,j), histo.GetYaxis().GetBinWidth(j))
            ####################################################################################################################
            
            #this rebin always
            if type(histo)==ROOT.TH1F or type(histo)==ROOT.TH1D :
                for i in range(1, histo.GetNbinsX()+1):
                    histo.SetBinContent(i, histo.GetBinContent(i)/histo.GetXaxis().GetBinWidth(i))
                    histo.SetBinError(i, histo.GetBinError(i)/histo.GetXaxis().GetBinWidth(i))    
            elif type(histo)==ROOT.TH2F or type(histo)==ROOT.TH2D :
                for i in range(1, histo.GetNbinsX()+1):
                    for j in range(1,  histo.GetNbinsY()+1):
                        histo.SetBinContent(i,j, histo.GetBinContent(i,j)/histo.GetXaxis().GetBinWidth(i)/histo.GetYaxis().GetBinWidth(j))
                        histo.SetBinError(i,j, histo.GetBinError(i,j)/histo.GetXaxis().GetBinWidth(i)/histo.GetYaxis().GetBinWidth(j)) 
                
                 
                
                
            
    def AngCoeffPlots(self,inputFile, fitFile, uncorrelate,suff,aposteriori,toy,impact,postfit, cleanNuisance, massImp,data,acceptance) :
        
        FitFile = ROOT.TFile.Open(fitFile)
        inFile = ROOT.TFile.Open(inputFile)
        # if aposteriori!='' :
        #     apoFile = ROOT.TFile.Open(aposteriori)
        #     self.getHistos(inFile=inFile, FitFile=FitFile, uncorrelate=uncorrelate,suff=suff,apoFile=apoFile)
        # else :
        #     self.getHistos(inFile=inFile, FitFile=FitFile, uncorrelate=uncorrelate,suff=suff)
        # print "WARNING: syst.band done with the fit central value (ok if asimov only)"
        
        if aposteriori!='' :
            apoFile = ROOT.TFile.Open(aposteriori) 
        else :
            apoFile = ''
        if toy!='' :
            toyFile = ROOT.TFile.Open(toy)
        else :
            toyFile = ''
        self.getHistos(inFile=inFile, FitFile=FitFile, uncorrelate=uncorrelate,suff=suff,apoFile=apoFile,toyFile=toyFile,impact=impact, postfit=postfit,data=data,acceptance=acceptance)
        self.varBinWidth_modifier(suff=suff)
        
        if cleanNuisance:
            del self.nuisanceDict['WHSFSyst']
            del self.nuisanceDict['PrefireWeight']
            del self.nuisanceDict['jme']
            del self.nuisanceDict['ptScale']
            del self.nuisanceDict['QCDnorm']
            del self.nuisanceDict['ewkXsec']
            del self.nuisanceDict['LeptonVeto']
            del self.nuisanceDict['mass']
            if massImp!='' :
                self.nuisanceDict["otherNoMass"] = [ROOT.kGreen+10, "other",23] 
            else : 
                self.nuisanceDict["other"] = [ROOT.kGreen+10, "other",23] 
             
        


        
            
        # ------------- build the bands for the angular coefficient ---------------------------# 
        for c in self.coeffDict:
            for cat in self.category : 
                if c!='UL' and c!='unpolarizedxsec' : 
                    self.histos[suff+'FitBand'+cat+c] = self.histos[suff+'MC'+cat+c].Clone('FitBand'+cat+c)
                    self.histos[suff+'FitBandPDF'+cat+c] = self.histos[suff+'MC'+cat+c].Clone('FitBandPDF'+cat+c) 
                    self.histos[suff+'FitBandScale'+cat+c] = self.histos[suff+'MC'+cat+c].Clone('FitBandScale'+cat+c) 
                else :
                    self.histos[suff+'FitBand'+cat+c] = self.histos[suff+'MC'+cat+'mapTot'].Clone('FitBand'+cat+c)
                    self.histos[suff+'FitBandPDF'+cat+c] = self.histos[suff+'MC'+cat+'mapTot'].Clone('FitBandPDF'+cat+c) 
                    self.histos[suff+'FitBandScale'+cat+c] = self.histos[suff+'MC'+cat+'mapTot'].Clone('FitBandScale'+cat+c) 
                    
                    if self.anaKind['name'] == 'normXsec' :
                        # if 'UL' in c :
                        self.histos[suff+'FitBand'+cat+c].Scale(1/self.sigmaTot[suff])
                        self.histos[suff+'FitBandPDF'+cat+c].Scale(1/self.sigmaTot[suff])
                        self.histos[suff+'FitBandScale'+cat+c].Scale(1/self.sigmaTot[suff])
                                

            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                    # print("not closure of", c, i, j , ",   (mc-fit)/mc=", (self.histos[suff+'FitBand'+c].GetBinContent(i,j)-self.histos[suff+'FitAC'+c].GetBinContent(i,j))/self.histos[suff+'FitBand'+c].GetBinContent(i,j), 'fit/mc=', self.histos[suff+'FitAC'+c].GetBinContent(i,j)/self.histos[suff+'FitBand'+c].GetBinContent(i,j))

                    #unclosure check
                    if abs(self.histos[suff+'FitBand'+c].GetBinContent(i,j))>0.0000001 :
                        if abs((self.histos[suff+'FitBand'+c].GetBinContent(i,j)-self.histos[suff+'FitAC'+c].GetBinContent(i,j))/self.histos[suff+'FitBand'+c].GetBinContent(i,j))>0.0000001 :
                                print("not closure of", c, i, j , ",   (mc-fit)/mc=", (self.histos[suff+'FitBand'+c].GetBinContent(i,j)-self.histos[suff+'FitAC'+c].GetBinContent(i,j))/self.histos[suff+'FitBand'+c].GetBinContent(i,j), 'fit/mc=', self.histos[suff+'FitAC'+c].GetBinContent(i,j)/self.histos[suff+'FitBand'+c].GetBinContent(i,j))
                    if i==1 :
                        if abs(self.histos[suff+'FitBandqt'+c].GetBinContent(j))>0.0000001 :
                            if abs((self.histos[suff+'FitBandqt'+c].GetBinContent(j)-self.histos[suff+'FitACqt'+c].GetBinContent(j))/self.histos[suff+'FitBandqt'+c].GetBinContent(j))>0.0000001 :
                                print("not closure of", c, j , "(qT),   (mc-fit)/mc=", (self.histos[suff+'FitBandqt'+c].GetBinContent(j)-self.histos[suff+'FitACqt'+c].GetBinContent(j))/self.histos[suff+'FitBandqt'+c].GetBinContent(j), 'fit/mc=', self.histos[suff+'FitACqt'+c].GetBinContent(j)/self.histos[suff+'FitBandqt'+c].GetBinContent(j))
                    if j==1 : 
                        if abs(self.histos[suff+'FitBandy'+c].GetBinContent(i))>0.0000001 :
                            if abs((self.histos[suff+'FitBandy'+c].GetBinContent(i)-self.histos[suff+'FitACy'+c].GetBinContent(i))/self.histos[suff+'FitBandy'+c].GetBinContent(i))>0.0000001 :
                                print("not closure of", c, i , "(y),   (mc-fit)/mc=", (self.histos[suff+'FitBandy'+c].GetBinContent(i)-self.histos[suff+'FitACy'+c].GetBinContent(i))/self.histos[suff+'FitBandy'+c].GetBinContent(i), 'fit/mc=', self.histos[suff+'FitACy'+c].GetBinContent(i)/self.histos[suff+'FitBandy'+c].GetBinContent(i))
                    
                    errPDF = 0.
                    if i==1 : errPDFqt = 0.
                    if j==1 : errPDFy = 0.
                    MCVal = self.histos[suff+'FitBand'+c].GetBinContent(i,j) #like MC, but already lumi scaled
                    if i==1 : MCValqt = self.histos[suff+'FitBandqt'+c].GetBinContent(j)
                    if j==1 : MCValy = self.histos[suff+'FitBandy'+c].GetBinContent(i)
                    
                    for sName in self.systDict['_LHEPdfWeight']:
                        # if 'unpol' in c:
                        #     systVal=self.histos[suff+'MC'+sName+'mapTot'].GetBinContent(i,j)
                        #     if i==1 : systValqt=self.histos[suff+'MCqt'+sName+'mapTot'].GetBinContent(j)
                        #     if j==1 : systValy=self.histos[suff+'MCy'+sName+'mapTot'].GetBinContent(i)
                        # else:
                        systVal = self.histos[suff+'MC'+sName+c].GetBinContent(i,j)
                        if i==1 : systValqt = self.histos[suff+'MCqt'+sName+c].GetBinContent(j)
                        if j==1 : systValy = self.histos[suff+'MCy'+sName+c].GetBinContent(i)
                        errPDF+= (MCVal - systVal)**2  
                        if i==1 : errPDFqt+= (MCValqt - systValqt)**2  
                        if j==1 : errPDFy+= (MCValy - systValy)**2                  
                    self.histos[suff+'FitBandPDF'+c].SetBinError(i,j,math.sqrt(errPDF)) 
                    if i==1 : self.histos[suff+'FitBandPDFqt'+c].SetBinError(j,math.sqrt(errPDFqt)) 
                    if j==1 : self.histos[suff+'FitBandPDFy'+c].SetBinError(i,math.sqrt(errPDFy)) 
                            
                    sListMod = copy.deepcopy(self.systDict['_LHEScaleWeight'])
                    if UNCORR :
                        sListMod.append("_nom") #add nominal variation
                    errScale = 0.
                    if i==1 : errScaleqt = 0.
                    if j==1 : errScaley = 0.
                    for sName in sListMod:
                        for sNameDen in sListMod :
                            if sNameDen!=sName and not (UNCORR and self.anaKind['angNames']): continue
                            if sNameDen!=sName and 'unpol' in c : continue
                            if sName=='_nom' and sNameDen=='_nom' : continue
                            if ([sName,sNameDen] in self.vetoScaleList) : continue  #extremal cases
                            # if 'unpol' in c:
                            #     systVal=self.histos[suff+'MC'+sName+'mapTot'].GetBinContent(i,j)
                            #     if i==1 : systValqt =self.histos[suff+'MCqt'+sName+'mapTot'].GetBinContent(j)
                            #     if j==1 : systValy =self.histos[suff+'MCy'+sName+'mapTot'].GetBinContent(i)
                            # else:
                            if UNCORR :
                                systVal = self.histos[suff+'MC'+sName+sNameDen+c].GetBinContent(i,j)
                                if i==1 : systValqt = self.histos[suff+'MCqt'+sName+sNameDen+c].GetBinContent(j)
                                if j==1 : systValy = self.histos[suff+'MCy'+sName+sNameDen+c].GetBinContent(i)
                            else : 
                                systVal = self.histos[suff+'MC'+sName+c].GetBinContent(i,j)
                                if i==1 : systValqt = self.histos[suff+'MCqt'+sName+c].GetBinContent(j)
                                if j==1 : systValy = self.histos[suff+'MCy'+sName+c].GetBinContent(i)
                            err_temp= (MCVal - systVal)**2
                            if i==1 : err_tempqt= (MCValqt - systValqt)**2
                            if j==1 : err_tempy= (MCValy - systValy)**2
                            if err_temp>errScale : 
                                errScale = err_temp
                            if i==1 and err_tempqt>errScaleqt : errScaleqt = err_tempqt
                            if j==1 and err_tempy>errScaley : errScaley = err_tempy
                                
                    self.histos[suff+'FitBandScale'+c].SetBinError(i,j,math.sqrt(errScale)) 
                    if i==1 : self.histos[suff+'FitBandScaleqt'+c].SetBinError(j,math.sqrt(errScaleqt)) 
                    if j==1 : self.histos[suff+'FitBandScaley'+c].SetBinError(i,math.sqrt(errScaley)) 

                    err=errPDF+errScale
                    if i==1 : errqt=errPDFqt+errScaleqt
                    if j==1 : erry=errPDFy+errScaley
                    
                    self.histos[suff+'FitBand'+c].SetBinError(i,j,math.sqrt(err)) 
                    if i==1 : self.histos[suff+'FitBandqt'+c].SetBinError(j,math.sqrt(errqt)) 
                    if j==1 : self.histos[suff+'FitBandy'+c].SetBinError(i,math.sqrt(erry)) 
        
        #--------------- build the relative uncertainity breakdown plots -----------------------#
        for c in self.coeffDict:
            self.histos[suff+'FitErr'+c] = self.histos[suff+'FitAC'+c].Clone('FitErr'+c) 
            self.histos[suff+'FitErrPDF'+c] = self.histos[suff+'FitBandPDF'+c].Clone('FitErrPDF'+c) 
            self.histos[suff+'FitErrScale'+c] = self.histos[suff+'FitBandScale'+c].Clone('FitErrScale'+c) 
            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): 
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1):
                    valCentral = self.histos[suff+'FitAC'+c].GetBinContent(i,j)
                    if valCentral==0.0 : 
                        valCentral =1
                        print("WARNING: null fit prediction")
                    self.histos[suff+'FitErr'+c].SetBinContent(i,j, self.histos[suff+'FitErr'+c].GetBinError(i,j)/abs(valCentral))
                    self.histos[suff+'FitErrPDF'+c].SetBinContent(i,j,self.histos[suff+'FitErrPDF'+c].GetBinError(i,j)/abs(valCentral))
                    self.histos[suff+'FitErrScale'+c].SetBinContent(i,j,self.histos[suff+'FitErrScale'+c].GetBinError(i,j)/abs(valCentral))
                    self.histos[suff+'FitErr'+c].SetBinError(i,j,0)
                    self.histos[suff+'FitErrPDF'+c].SetBinError(i,j,0)
                    self.histos[suff+'FitErrScale'+c].SetBinError(i,j,0)   
                    # print suff, c, i, j, "central=",valCentral, "main=", self.histos[suff+'FitErr'+c].GetBinContent(i,j), "pdf=",self.histos[suff+'FitErrPDF'+c].GetBinContent(i,j),"scale=",self.histos[suff+'FitErrScale'+c].GetBinContent(i,j)
        
        #--------------- build the relative uncertainity breakdown plots - Y TREND -----------------------#
        for c in self.coeffDict:
            self.histos[suff+'FitErry'+c] = self.histos[suff+'FitACy'+c].Clone('FitErry'+c) 
            self.histos[suff+'FitErrPDFy'+c] = self.histos[suff+'FitBandPDFy'+c].Clone('FitErrPDFy'+c) 
            self.histos[suff+'FitErrScaley'+c] = self.histos[suff+'FitBandScaley'+c].Clone('FitErrScaley'+c) 
            for i in range(1, self.histos[suff+'FitACy'+c].GetNbinsX()+1): 
                valCentral = self.histos[suff+'FitACy'+c].GetBinContent(i)
                if valCentral==0.0 : 
                        valCentral =1
                        print("WARNING: null fit prediction")
                self.histos[suff+'FitErry'+c].SetBinContent(i, self.histos[suff+'FitErry'+c].GetBinError(i)/abs(valCentral))
                self.histos[suff+'FitErrPDFy'+c].SetBinContent(i,self.histos[suff+'FitErrPDFy'+c].GetBinError(i)/abs(valCentral))
                self.histos[suff+'FitErrScaley'+c].SetBinContent(i,self.histos[suff+'FitErrScaley'+c].GetBinError(i)/abs(valCentral))
                self.histos[suff+'FitErry'+c].SetBinError(i,0)
                self.histos[suff+'FitErrPDFy'+c].SetBinError(i,0)
                self.histos[suff+'FitErrScaley'+c].SetBinError(i,0)     
        
        #--------------- build the relative uncertainity breakdown plots - qt TREND -----------------------#
        for c in self.coeffDict:
            self.histos[suff+'FitErrqt'+c] = self.histos[suff+'FitACqt'+c].Clone('FitErrqt'+c) 
            self.histos[suff+'FitErrPDFqt'+c] = self.histos[suff+'FitBandPDFqt'+c].Clone('FitErrPDFqt'+c) 
            self.histos[suff+'FitErrScaleqt'+c] = self.histos[suff+'FitBandScaleqt'+c].Clone('FitErrScaleqt'+c) 
            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1):
                valCentral = self.histos[suff+'FitACqt'+c].GetBinContent(i)
                if valCentral==0.0 : 
                        valCentral =1
                        print("WARNING: null fit prediction")
                self.histos[suff+'FitErrqt'+c].SetBinContent(i, self.histos[suff+'FitErrqt'+c].GetBinError(i)/abs(valCentral))
                self.histos[suff+'FitErrPDFqt'+c].SetBinContent(i,self.histos[suff+'FitErrPDFqt'+c].GetBinError(i)/abs(valCentral))
                self.histos[suff+'FitErrScaleqt'+c].SetBinContent(i,self.histos[suff+'FitErrScaleqt'+c].GetBinError(i)/abs(valCentral))
                self.histos[suff+'FitErrqt'+c].SetBinError(i,0)
                self.histos[suff+'FitErrPDFqt'+c].SetBinError(i,0)
                self.histos[suff+'FitErrScaleqt'+c].SetBinError(i,0)  
        
        #--------------- build the ratio plots -----------------------#
        for c in self.coeffDict :
            self.histos[suff+'FitRatioAC'+c] = self.histos[suff+'FitAC'+c].Clone('FitRatioAC'+c)
            self.histos[suff+'FitRatio'+c] = self.histos[suff+'FitBand'+c].Clone('FitRatio'+c)
            self.histos[suff+'FitRatioPDF'+c] = self.histos[suff+'FitBandPDF'+c].Clone('FitRatioPDF'+c) 
            self.histos[suff+'FitRatioScale'+c] = self.histos[suff+'FitBandScale'+c].Clone('FitRatioScale'+c) 
            for d in self.dirList : 
                self.histos[suff+'FitRatioPDF'+d+c] = self.histos[suff+'FitBandPDF'+c].Clone('FitRatioPDF'+c) 
                self.histos[suff+'FitRatioScale'+d+c] = self.histos[suff+'FitBandScale'+c].Clone('FitRatioScale'+c) 
            if apoFile!='' :
                self.histos[suff+'apo'+'Ratio'+c] = self.histos[suff+'apo'+c].Clone("apoRatio"+c)
            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): 
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1):
                    # if not "unpol" in c:   
                    valCentral = self.histos[suff+'MC'+c].GetBinContent(i,j)
                    if ('unpol' in c  or (not self.anaKind['differential'])): 
                        if valCentral==0.0 : 
                            valCentral =1
                            print("WARNING: null fit prediction")                    
                        self.histos[suff+'FitRatioAC'+c].SetBinContent(i,j, self.histos[suff+'FitRatioAC'+c].GetBinContent(i,j)/valCentral)
                        self.histos[suff+'FitRatio'+c].SetBinContent(i,j, self.histos[suff+'FitRatio'+c].GetBinContent(i,j)/valCentral)
                        self.histos[suff+'FitRatioPDF'+c].SetBinContent(i,j,self.histos[suff+'FitRatioPDF'+c].GetBinContent(i,j)/valCentral)
                        self.histos[suff+'FitRatioScale'+c].SetBinContent(i,j,self.histos[suff+'FitRatioScale'+c].GetBinContent(i,j)/valCentral)
                        self.histos[suff+'FitRatioAC'+c].SetBinError(i,j,self.histos[suff+'FitRatioAC'+c].GetBinError(i,j)/valCentral)
                        self.histos[suff+'FitRatio'+c].SetBinError(i,j,self.histos[suff+'FitRatio'+c].GetBinError(i,j)/valCentral)
                        self.histos[suff+'FitRatioPDF'+c].SetBinError(i,j,self.histos[suff+'FitRatioPDF'+c].GetBinError(i,j)/valCentral)
                        self.histos[suff+'FitRatioScale'+c].SetBinError(i,j,self.histos[suff+'FitRatioScale'+c].GetBinError(i,j)/valCentral) 
                        for d in self.dirList :
                            for var in ['PDF','Scale'] :
                                if d=='up' :
                                    self.histos[suff+'FitRatio'+var+d+c].SetBinContent(i,j,self.histos[suff+'FitBand'+var+c].GetBinContent(i,j)/valCentral+abs(self.histos[suff+'FitBand'+var+c].GetBinError(i,j)/valCentral))
                                else :
                                    self.histos[suff+'FitRatio'+var+d+c].SetBinContent(i,j,self.histos[suff+'FitBand'+var+c].GetBinContent(i,j)/valCentral-abs(self.histos[suff+'FitBand'+var+c].GetBinError(i,j)/valCentral))
                                self.histos[suff+'FitRatio'+var+c].SetBinError(i,j,0)
                        if apoFile!='' :
                            self.histos[suff+'apo'+'Ratio'+c].SetBinContent(i,j, self.histos[suff+'apo'+c].GetBinContent(i,j)/valCentral)
                            self.histos[suff+'apo'+'Ratio'+c].SetBinError(i,j, self.histos[suff+'apo'+c].GetBinError(i,j)/valCentral)
                    else : #no ratio but difference for Ai
                        self.histos[suff+'FitRatioAC'+c].SetBinContent(i,j, self.histos[suff+'FitRatioAC'+c].GetBinContent(i,j)-valCentral)
                        self.histos[suff+'FitRatio'+c].SetBinContent(i,j, self.histos[suff+'FitRatio'+c].GetBinContent(i,j)-valCentral)
                        self.histos[suff+'FitRatioPDF'+c].SetBinContent(i,j,self.histos[suff+'FitRatioPDF'+c].GetBinContent(i,j)-valCentral)
                        self.histos[suff+'FitRatioScale'+c].SetBinContent(i,j,self.histos[suff+'FitRatioScale'+c].GetBinContent(i,j)-valCentral)
                        for d in self.dirList :
                            for var in ['PDF','Scale'] :
                                if d=='up' :
                                    self.histos[suff+'FitRatio'+var+d+c].SetBinContent(i,j,self.histos[suff+'FitBand'+var+c].GetBinContent(i,j)-valCentral+abs(self.histos[suff+'FitBand'+var+c].GetBinError(i,j)))
                                else :
                                    self.histos[suff+'FitRatio'+var+d+c].SetBinContent(i,j,self.histos[suff+'FitBand'+var+c].GetBinContent(i,j)-valCentral-abs(self.histos[suff+'FitBand'+var+c].GetBinError(i,j)))
                                self.histos[suff+'FitRatio'+var+c].SetBinError(i,j,0)
                        if apoFile!='' :
                            self.histos[suff+'apo'+'Ratio'+c].SetBinContent(i,j, self.histos[suff+'apo'+c].GetBinContent(i,j)-valCentral)
                        
                    
                
            
            #--------------- build the ratio plots  Y trends -----------------------#
            self.histos[suff+'FitRatioACy'+c] = self.histos[suff+'FitAC'+'y'+c].Clone('FitRatioACy'+c)
            self.histos[suff+'FitRatioy'+c] = self.histos[suff+'FitBand'+'y'+c].Clone('FitRatioy'+c)
            self.histos[suff+'FitRatioPDFy'+c] = self.histos[suff+'FitBandPDF'+'y'+c].Clone('FitRatioPDFy'+c)
            self.histos[suff+'FitRatioScaley'+c] = self.histos[suff+'FitBandScale'+'y'+c].Clone('FitRatioScaley'+c)
            for d in self.dirList : 
                self.histos[suff+'FitRatioPDF'+'y'+d+c] = self.histos[suff+'FitBandPDF'+'y'+c].Clone('FitRatioPDF'+c) 
                self.histos[suff+'FitRatioScale'+'y'+d+c] = self.histos[suff+'FitBandScale'+'y'+c].Clone('FitRatioScale'+c) 
            if apoFile!='' :
                self.histos[suff+'apo'+'Ratio'+'y'+c] = self.histos[suff+'apo'+'y'+c].Clone("apoRatio"+'y'+c)
            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): 
                # if not "unpol" in c:   
                valCentral = self.histos[suff+'MC'+'y'+c].GetBinContent(i)
                if ('unpol' in c  or (not self.anaKind['differential'])): 
                    if valCentral==0.0 : 
                            valCentral =1
                            print("WARNING: null fit prediction")
                    self.histos[suff+'FitRatioAC'+'y'+c].SetBinContent(i, self.histos[suff+'FitRatioAC'+'y'+c].GetBinContent(i)/valCentral)
                    self.histos[suff+'FitRatio'+'y'+c].SetBinContent(i, self.histos[suff+'FitRatio'+'y'+c].GetBinContent(i)/valCentral)
                    self.histos[suff+'FitRatioPDF'+'y'+c].SetBinContent(i,self.histos[suff+'FitRatioPDF'+'y'+c].GetBinContent(i)/valCentral)
                    self.histos[suff+'FitRatioScale'+'y'+c].SetBinContent(i,self.histos[suff+'FitRatioScale'+'y'+c].GetBinContent(i)/valCentral)
                    self.histos[suff+'FitRatioAC'+'y'+c].SetBinError(i,self.histos[suff+'FitRatioAC'+'y'+c].GetBinError(i)/valCentral)
                    self.histos[suff+'FitRatio'+'y'+c].SetBinError(i,self.histos[suff+'FitRatio'+'y'+c].GetBinError(i)/valCentral)
                    self.histos[suff+'FitRatioPDF'+'y'+c].SetBinError(i,self.histos[suff+'FitRatioPDF'+'y'+c].GetBinError(i)/valCentral)
                    self.histos[suff+'FitRatioScale'+'y'+c].SetBinError(i,self.histos[suff+'FitRatioScale'+'y'+c].GetBinError(i)/valCentral) 
                    for d in self.dirList :
                        for var in ['PDF','Scale'] :
                            if d=='up' :
                                self.histos[suff+'FitRatio'+var+'y'+d+c].SetBinContent(i,self.histos[suff+'FitBand'+var+'y'+c].GetBinContent(i)/valCentral+abs(self.histos[suff+'FitBand'+var+'y'+c].GetBinError(i)/valCentral))
                            else :
                                self.histos[suff+'FitRatio'+var+'y'+d+c].SetBinContent(i,self.histos[suff+'FitBand'+var+'y'+c].GetBinContent(i)/valCentral-abs(self.histos[suff+'FitBand'+var+'y'+c].GetBinError(i)/valCentral))
                            self.histos[suff+'FitRatio'+var+'y'+c].SetBinError(i,0)
                    if apoFile!='' :
                        self.histos[suff+'apo'+'Ratio'+'y'+c].SetBinContent(i, self.histos[suff+'apo'+'y'+c].GetBinContent(i)/valCentral)
                        self.histos[suff+'apo'+'Ratio'+'y'+c].SetBinError(i, self.histos[suff+'apo'+'y'+c].GetBinError(i)/valCentral)
                else : #no ratio but difference for Ai
                    self.histos[suff+'FitRatioAC'+'y'+c].SetBinContent(i, self.histos[suff+'FitRatioAC'+'y'+c].GetBinContent(i)-valCentral)
                    self.histos[suff+'FitRatio'+'y'+c].SetBinContent(i, self.histos[suff+'FitRatio'+'y'+c].GetBinContent(i)-valCentral)
                    self.histos[suff+'FitRatioPDF'+'y'+c].SetBinContent(i,self.histos[suff+'FitRatioPDF'+'y'+c].GetBinContent(i)-valCentral)
                    self.histos[suff+'FitRatioScale'+'y'+c].SetBinContent(i,self.histos[suff+'FitRatioScale'+'y'+c].GetBinContent(i)-valCentral)
                    for d in self.dirList :
                        for var in ['PDF','Scale'] :
                            if d=='up' :
                                self.histos[suff+'FitRatio'+var+'y'+d+c].SetBinContent(i,self.histos[suff+'FitBand'+var+'y'+c].GetBinContent(i)-valCentral+abs(self.histos[suff+'FitBand'+var+'y'+c].GetBinError(i)))
                            else :
                                self.histos[suff+'FitRatio'+var+'y'+d+c].SetBinContent(i,self.histos[suff+'FitBand'+var+'y'+c].GetBinContent(i)-valCentral-abs(self.histos[suff+'FitBand'+var+'y'+c].GetBinError(i)))
                            self.histos[suff+'FitRatio'+var+'y'+c].SetBinError(i,0)
                    if apoFile!='' :
                        self.histos[suff+'apo'+'Ratio'+'y'+c].SetBinContent(i, self.histos[suff+'apo'+'y'+c].GetBinContent(i)-valCentral)
            
            
            #--------------- build the ratio plots  qt trends -----------------------#
            self.histos[suff+'FitRatioACqt'+c] = self.histos[suff+'FitAC'+'qt'+c].Clone('FitRatioACqt'+c)
            self.histos[suff+'FitRatioqt'+c] = self.histos[suff+'FitBand'+'qt'+c].Clone('FitRatioqt'+c)
            self.histos[suff+'FitRatioPDFqt'+c] = self.histos[suff+'FitBandPDF'+'qt'+c].Clone('FitRatioPDFqt'+c)
            self.histos[suff+'FitRatioScaleqt'+c] = self.histos[suff+'FitBandScale'+'qt'+c].Clone('FitBandScaleqt'+c)
            for d in self.dirList : 
                self.histos[suff+'FitRatioPDF'+'qt'+d+c] = self.histos[suff+'FitBandPDF'+'qt'+c].Clone('FitRatioPDF'+c) 
                self.histos[suff+'FitRatioScale'+'qt'+d+c] = self.histos[suff+'FitBandScale'+'qt'+c].Clone('FitRatioScale'+c) 
            if apoFile!='' :
                self.histos[suff+'apo'+'Ratio'+'qt'+c] = self.histos[suff+'apo'+'qt'+c].Clone("apoRatio"+'qt'+c)
            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): 
                # if not "unpol" in c:   
                valCentral = self.histos[suff+'MC'+'qt'+c].GetBinContent(i)
                if ('unpol' in c  or (not self.anaKind['differential'])): 
                    if valCentral==0.0 : 
                        valCentral =1
                        print("WARNING: null fit prediction")
                    self.histos[suff+'FitRatioAC'+'qt'+c].SetBinContent(i, self.histos[suff+'FitRatioAC'+'qt'+c].GetBinContent(i)/valCentral)
                    self.histos[suff+'FitRatio'+'qt'+c].SetBinContent(i, self.histos[suff+'FitRatio'+'qt'+c].GetBinContent(i)/valCentral)
                    self.histos[suff+'FitRatioPDF'+'qt'+c].SetBinContent(i,self.histos[suff+'FitRatioPDF'+'qt'+c].GetBinContent(i)/valCentral)
                    self.histos[suff+'FitRatioScale'+'qt'+c].SetBinContent(i,self.histos[suff+'FitRatioScale'+'qt'+c].GetBinContent(i)/valCentral)
                    self.histos[suff+'FitRatioAC'+'qt'+c].SetBinError(i,self.histos[suff+'FitRatioAC'+'qt'+c].GetBinError(i)/valCentral)
                    self.histos[suff+'FitRatio'+'qt'+c].SetBinError(i,self.histos[suff+'FitRatio'+'qt'+c].GetBinError(i)/valCentral)
                    self.histos[suff+'FitRatioPDF'+'qt'+c].SetBinError(i,self.histos[suff+'FitRatioPDF'+'qt'+c].GetBinError(i)/valCentral)
                    self.histos[suff+'FitRatioScale'+'qt'+c].SetBinError(i,self.histos[suff+'FitRatioScale'+'qt'+c].GetBinError(i)/valCentral) 
                    for d in self.dirList :
                        for var in ['PDF','Scale'] :
                            if d=='up' :
                                self.histos[suff+'FitRatio'+var+'qt'+d+c].SetBinContent(i,self.histos[suff+'FitBand'+var+'qt'+c].GetBinContent(i)/valCentral+abs(self.histos[suff+'FitBand'+var+'qt'+c].GetBinError(i)/valCentral))
                            else :
                                self.histos[suff+'FitRatio'+var+'qt'+d+c].SetBinContent(i,self.histos[suff+'FitBand'+var+'qt'+c].GetBinContent(i)/valCentral-abs(self.histos[suff+'FitBand'+var+'qt'+c].GetBinError(i)/valCentral))
                            self.histos[suff+'FitRatio'+var+'qt'+c].SetBinError(i,0)
                    if apoFile!='' :
                        self.histos[suff+'apo'+'Ratio'+'qt'+c].SetBinContent(i, self.histos[suff+'apo'+'qt'+c].GetBinContent(i)/valCentral)
                        self.histos[suff+'apo'+'Ratio'+'qt'+c].SetBinError(i,self.histos[suff+'apo'+'qt'+c].GetBinError(i)/valCentral)
                else : #no ratio but difference for Ai
                    self.histos[suff+'FitRatioAC'+'qt'+c].SetBinContent(i, self.histos[suff+'FitRatioAC'+'qt'+c].GetBinContent(i)-valCentral)
                    self.histos[suff+'FitRatio'+'qt'+c].SetBinContent(i, self.histos[suff+'FitRatio'+'qt'+c].GetBinContent(i)-valCentral)
                    self.histos[suff+'FitRatioPDF'+'qt'+c].SetBinContent(i,self.histos[suff+'FitRatioPDF'+'qt'+c].GetBinContent(i)-valCentral)
                    self.histos[suff+'FitRatioScale'+'qt'+c].SetBinContent(i,self.histos[suff+'FitRatioScale'+'qt'+c].GetBinContent(i)-valCentral)
                    for d in self.dirList :
                        for var in ['PDF','Scale'] :
                            if d=='up' :
                                self.histos[suff+'FitRatio'+var+'qt'+d+c].SetBinContent(i,self.histos[suff+'FitBand'+var+'qt'+c].GetBinContent(i)-valCentral+abs(self.histos[suff+'FitBand'+var+'qt'+c].GetBinError(i)))
                            else :
                                self.histos[suff+'FitRatio'+var+'qt'+d+c].SetBinContent(i,self.histos[suff+'FitBand'+var+'qt'+c].GetBinContent(i)-valCentral-abs(self.histos[suff+'FitBand'+var+'qt'+c].GetBinError(i)))
                            self.histos[suff+'FitRatio'+var+'qt'+c].SetBinError(i,0)
                    if apoFile!='' :
                        self.histos[suff+'apo'+'Ratio'+'qt'+c].SetBinContent(i, self.histos[suff+'apo'+'qt'+c].GetBinContent(i)-valCentral)

        
                
                    
        #---------------------- build the covariance and correlation matrices --------- #
        for c in self.coeffDict:# y large, qt fine        
            self.histos[suff+'corrMat'+c] = ROOT.TH2D(suff+'corrMat_{c}'.format(c=c), suff+'corrMat_{c}'.format(c=c), len(self.unrolledYQt)-1, array('f',self.unrolledYQt), len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            self.histos[suff+'covMat'+c] = ROOT.TH2D(suff+'covMat_{c}'.format(c=c), suff+'covMat_{c}'.format(c=c), len(self.unrolledYQt)-1, array('f',self.unrolledYQt), len(self.unrolledYQt)-1, array('f',self.unrolledYQt)) 
            for q1 in self.qtArr[:-1]: 
                for y1 in self.yArr[:-1] :
                    for q2 in self.qtArr[:-1]: 
                        for y2 in self.yArr[:-1] : 
                            indexUNR_X = self.yArr.index(float(y1))*(len(self.qtArr)-1)+self.qtArr.index(float(q1))
                            indexUNR_Y = self.yArr.index(float(y2))*(len(self.qtArr)-1)+self.qtArr.index(float(q2))
                            if self.anaKind['angNames'] :
                                nameX = 'y_'+str(self.yArr.index(y1)+1)+'_qt_'+str(self.qtArr.index(q1)+1)+'_'+c
                                nameY = 'y_'+str(self.yArr.index(y2)+1)+'_qt_'+str(self.qtArr.index(q2)+1)+'_'+c
                            else :
                                nameX = 'helXsecs'+c+'_y_'+str(self.yArr.index(y1)+1)+'_qt_'+str(self.qtArr.index(q1)+1)+'_'+self.anaKind['diffString']
                                nameY = 'helXsecs'+c+'_y_'+str(self.yArr.index(y2)+1)+'_qt_'+str(self.qtArr.index(q2)+1)+'_'+self.anaKind['diffString']
                            try :
                                corrVal = self.histos[suff+'corrMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                                covVal = self.histos[suff+'covMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                            except :
                                print("WARNING: missing matrices",nameX,nameY)
                                corrVal =0
                                covVal= 0
                            self.histos[suff+'corrMat'+c].SetBinContent(indexUNR_X+1,indexUNR_Y+1,corrVal)
                            self.histos[suff+'covMat'+c].SetBinContent(indexUNR_X+1,indexUNR_Y+1,covVal)
                            
        for i in range(1, self.histos[suff+'FitAC'+self.coeffList[0]].GetNbinsX()+1): #loop over y bins
            self.histos[suff+'corrMat'+'y'+str(i)] = ROOT.TH2D(suff+'corrMat_y{c}'.format(c=i), suff+'corrMat_y{c}'.format(c=i), len(self.unrolledAQt)-1, array('f',self.unrolledAQt), len(self.unrolledAQt)-1, array('f',self.unrolledAQt))
            self.histos[suff+'covMat'+'y'+str(i)] = ROOT.TH2D(suff+'covMat_y{c}'.format(c=i), suff+'covMat_y{c}'.format(c=i), len(self.unrolledAQt)-1, array('f',self.unrolledAQt), len(self.unrolledAQt)-1, array('f',self.unrolledAQt))
            for q1 in self.qtArr[:-1]: 
                for c1 in self.coeffArr[:-1] :
                    for q2 in self.qtArr[:-1]: 
                        for c2 in self.coeffArr[:-1] : 
                            indexUNR_X = self.coeffArr.index(float(c1))*(len(self.qtArr)-1)+self.qtArr.index(float(q1))
                            indexUNR_Y = self.coeffArr.index(float(c2))*(len(self.qtArr)-1)+self.qtArr.index(float(q2))
                            if self.anaKind['angNames'] : 
                                nameX = 'y_'+str(i)+'_qt_'+str(self.qtArr.index(q1)+1)+'_'+self.coeffList[self.coeffArr.index(c1)]
                                nameY = 'y_'+str(i)+'_qt_'+str(self.qtArr.index(q2)+1)+'_'+self.coeffList[self.coeffArr.index(c2)]
                            else :
                                nameX = 'helXsecs'+self.coeffList[self.coeffArr.index(c1)]+'_y_'+str(i)+'_qt_'+str(self.qtArr.index(q1)+1)+'_'+self.anaKind['diffString']
                                nameY = 'helXsecs'+self.coeffList[self.coeffArr.index(c2)]+'_y_'+str(i)+'_qt_'+str(self.qtArr.index(q2)+1)+'_'+self.anaKind['diffString']
                            try :
                                corrVal = self.histos[suff+'corrMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                                covVal = self.histos[suff+'covMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                            except :
                                print("WARNING: missing matrices",nameX,nameY)
                                corrVal =0
                                covVal= 0
                            self.histos[suff+'corrMat'+'y'+str(i)].SetBinContent(indexUNR_X+1,indexUNR_Y+1,corrVal)
                            self.histos[suff+'covMat'+'y'+str(i)].SetBinContent(indexUNR_X+1,indexUNR_Y+1,covVal)
        
        for j in range(1, self.histos[suff+'FitAC'+self.coeffList[0]].GetNbinsY()+1): #loop over qt bins
            self.histos[suff+'corrMat'+'qt'+str(j)] = ROOT.TH2D(suff+'corrMat_qt{c}'.format(c=j), suff+'corrMat_qt{c}'.format(c=j), len(self.unrolledAY)-1, array('f',self.unrolledAY), len(self.unrolledAY)-1, array('f',self.unrolledAY))
            self.histos[suff+'covMat'+'qt'+str(j)] = ROOT.TH2D(suff+'covMat_qt{c}'.format(c=j), suff+'covMat_qt{c}'.format(c=j), len(self.unrolledAY)-1, array('f',self.unrolledAY), len(self.unrolledAY)-1, array('f',self.unrolledAY))
            for y1 in self.yArr[:-1]: 
                for c1 in self.coeffArr[:-1] :
                    for y2 in self.yArr[:-1]: 
                        for c2 in self.coeffArr[:-1] : 
                            indexUNR_X = self.coeffArr.index(float(c1))*(len(self.yArr)-1)+self.yArr.index(float(y1))
                            indexUNR_Y = self.coeffArr.index(float(c2))*(len(self.yArr)-1)+self.yArr.index(float(y2))
                            if self.anaKind['angNames'] : 
                                nameX = 'y_'+str(self.yArr.index(y1)+1)+'_qt_'+str(j)+'_'+self.coeffList[self.coeffArr.index(c1)]
                                nameY = 'y_'+str(self.yArr.index(y2)+1)+'_qt_'+str(j)+'_'+self.coeffList[self.coeffArr.index(c2)]
                            else :
                                nameX = 'helXsecs'+self.coeffList[self.coeffArr.index(c1)]+'_y_'+str(self.yArr.index(y1)+1)+'_qt_'+str(j)+'_'+self.anaKind['diffString']
                                nameY = 'helXsecs'+self.coeffList[self.coeffArr.index(c2)]+'_y_'+str(self.yArr.index(y2)+1)+'_qt_'+str(j)+'_'+self.anaKind['diffString']
                            try :
                                corrVal = self.histos[suff+'corrMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                                covVal = self.histos[suff+'covMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                            except :
                                print("WARNING: missing matrices",nameX,nameY)
                                corrVal =0
                                covVal= 0
                            self.histos[suff+'corrMat'+'qt'+str(j)].SetBinContent(indexUNR_X+1,indexUNR_Y+1,corrVal)
                            self.histos[suff+'covMat'+'qt'+str(j)].SetBinContent(indexUNR_X+1,indexUNR_Y+1,covVal)        
        
        #---------------------- build the covariance and correlation matrices Y and Qt trend --------- #
        # for c in self.coeffDict:# y large, qt fine        
        #     self.histos[suff+'corrMat'+'Integrated'+c] = ROOT.TH2D(suff+'corrMat_Integrated_{c}'.format(c=c), suff+'corrMat_Integrated_{c}'.format(c=c), len(self.unrolledYQt)-1, array('f',self.unrolledYQt), len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
        #     self.histos[suff+'covMat'+'Integrated'+c] = ROOT.TH2D(suff+'covMat_Integrated_{c}'.format(c=c), suff+'covMat_Integrated_{c}'.format(c=c), len(self.unrolledYQt)-1, array('f',self.unrolledYQt), len(self.unrolledYQt)-1, array('f',self.unrolledYQt)) 
        #     for q1 in self.qtArr[:-1]: 
        #         for y1 in self.yArr[:-1] :
        #             for q2 in self.qtArr[:-1]: 
        #                 for y2 in self.yArr[:-1] : 
        #                     indexUNR_X = self.yArr.index(float(y1))*(len(self.qtArr)-1)+self.qtArr.index(float(q1))
        #                     indexUNR_Y = self.yArr.index(float(y2))*(len(self.qtArr)-1)+self.qtArr.index(float(q2))
        #                     nameX = 'y_'+str(self.yArr.index(y1)+1)+'_qt_'+str(self.qtArr.index(q1)+1)+'_'+c
        #                     nameY = 'y_'+str(self.yArr.index(y2)+1)+'_qt_'+str(self.qtArr.index(q2)+1)+'_'+c
        #                     corrVal = self.histos[suff+'corrMat'+'Integrated'].GetBinContent(self.histos[suff+'corrMat'+'Integrated'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'+'Integrated'].GetYaxis().FindBin(nameY))
        #                     covVal = self.histos[suff+'covMat'+'Integrated'].GetBinContent(self.histos[suff+'corrMat'+'Integrated'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'+'Integrated'].GetYaxis().FindBin(nameY))
        #                     self.histos[suff+'corrMat'+'Integrated'+c].SetBinContent(indexUNR_X,indexUNR_Y,corrVal)
        #                     self.histos[suff+'covMat'+'Integrated'+c].SetBinContent(indexUNR_X,indexUNR_Y,covVal)
                            
        
        self.histos[suff+'corrMat'+'Integrated'+'y'] = ROOT.TH2D(suff+'corrMat_Integrated_y', suff+'corrMat_Integrated_y', len(self.unrolledAY)-1, array('f',self.unrolledAY), len(self.unrolledAY)-1, array('f',self.unrolledAY))
        self.histos[suff+'covMat'+'Integrated'+'y'] = ROOT.TH2D(suff+'covMat_Integrated_y', suff+'covMat_Integrated_y', len(self.unrolledAY)-1, array('f',self.unrolledAY), len(self.unrolledAY)-1, array('f',self.unrolledAY))
        for y1 in self.yArr[:-1]: 
            for c1 in self.coeffArr[:-1] :
                for y2 in self.yArr[:-1]: 
                    for c2 in self.coeffArr[:-1] : 
                        indexUNR_X = self.coeffArr.index(float(c1))*(len(self.yArr)-1)+self.yArr.index(float(y1))
                        indexUNR_Y = self.coeffArr.index(float(c2))*(len(self.yArr)-1)+self.yArr.index(float(y2))
                        if self.anaKind['angNames'] :  
                            nameX = 'y_'+str(self.yArr.index(y1)+1)+'_helmeta_'+self.coeffList[self.coeffArr.index(c1)]
                            nameY = 'y_'+str(self.yArr.index(y2)+1)+'_helmeta_'+self.coeffList[self.coeffArr.index(c2)]
                        else :
                            nameX = 'helXsecs'+self.coeffList[self.coeffArr.index(c1)]+'_y_'+str(self.yArr.index(y1)+1)+'_'+self.anaKind['intString'].replace('sumpois','sumxsec')
                            nameY = 'helXsecs'+self.coeffList[self.coeffArr.index(c2)]+'_y_'+str(self.yArr.index(y2)+1)+'_'+self.anaKind['intString'].replace('sumpois','sumxsec')
                        try :
                            corrVal = self.histos[suff+'corrMat'+'Integrated'].GetBinContent(self.histos[suff+'corrMat'+'Integrated'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'+'Integrated'].GetYaxis().FindBin(nameY))
                            covVal = self.histos[suff+'covMat'+'Integrated'].GetBinContent(self.histos[suff+'corrMat'+'Integrated'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'+'Integrated'].GetYaxis().FindBin(nameY))
                        except :
                                print("WARNING: missing matrices",nameX,nameY)
                                corrVal =0
                                covVal= 0
                        self.histos[suff+'corrMat'+'Integrated'+'y'].SetBinContent(indexUNR_X+1,indexUNR_Y+1,corrVal)
                        self.histos[suff+'covMat'+'Integrated'+'y'].SetBinContent(indexUNR_X+1,indexUNR_Y+1,covVal)
        
        
        self.histos[suff+'corrMat'+'Integrated'+'qt'] = ROOT.TH2D(suff+'corrMat_Integrated_qt', suff+'corrMat_Integrated_qt', len(self.unrolledAQt)-1, array('f',self.unrolledAQt), len(self.unrolledAQt)-1, array('f',self.unrolledAQt))
        self.histos[suff+'covMat'+'Integrated'+'qt'] = ROOT.TH2D(suff+'covMat_Integrated_qt', suff+'covMat_Integrated_qt', len(self.unrolledAQt)-1, array('f',self.unrolledAQt), len(self.unrolledAQt)-1, array('f',self.unrolledAQt))
        for q1 in self.qtArr[:-1]: 
            for c1 in self.coeffArr[:-1] :
                for q2 in self.qtArr[:-1]: 
                    for c2 in self.coeffArr[:-1] : 
                        indexUNR_X = self.coeffArr.index(float(c1))*(len(self.qtArr)-1)+self.qtArr.index(float(q1))
                        indexUNR_Y = self.coeffArr.index(float(c2))*(len(self.qtArr)-1)+self.qtArr.index(float(q2))
                        if self.anaKind['angNames'] : 
                            nameX = 'qt_'+str(self.qtArr.index(q1)+1)+'_helmeta_'+self.coeffList[self.coeffArr.index(c1)]
                            nameY = 'qt_'+str(self.qtArr.index(q2)+1)+'_helmeta_'+self.coeffList[self.coeffArr.index(c2)]
                        else :
                            nameX = 'helXsecs'+self.coeffList[self.coeffArr.index(c1)]+'_qt_'+str(self.qtArr.index(q1)+1)+'_'+self.anaKind['intString'].replace('sumpois','sumxsec')
                            nameY = 'helXsecs'+self.coeffList[self.coeffArr.index(c2)]+'_qt_'+str(self.qtArr.index(q2)+1)+'_'+self.anaKind['intString'].replace('sumpois','sumxsec')
                        try :
                            corrVal = self.histos[suff+'corrMat'+'Integrated'].GetBinContent(self.histos[suff+'corrMat'+'Integrated'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'+'Integrated'].GetYaxis().FindBin(nameY))
                            covVal = self.histos[suff+'covMat'+'Integrated'].GetBinContent(self.histos[suff+'corrMat'+'Integrated'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'+'Integrated'].GetYaxis().FindBin(nameY))
                        except :
                                print("WARNING: missing matrices",nameX,nameY)
                                corrVal =0
                                covVal= 0
                        self.histos[suff+'corrMat'+'Integrated'+'qt'].SetBinContent(indexUNR_X+1,indexUNR_Y+1,corrVal)
                        self.histos[suff+'covMat'+'Integrated'+'qt'].SetBinContent(indexUNR_X+1,indexUNR_Y+1,covVal)
        
        
        #------------------------ build the impact plots ----------------------------------- #
        
        if impact :
            impactVals = {}
            
            skipDiffImpact = False
            try : 
                for qty in range(1, self.histos[suff+'impact2D'+'group'+'UNR'].GetNbinsX()+1):    
                    impBin = self.histos[suff+'impact2D'+'group'+'UNR'].GetXaxis().GetBinLabel(qty)          
                    for g in range(1, self.histos[suff+'impact2D'+'group'+'UNR'].GetNbinsY()+1):
                        nuiGroup=self.histos[suff+'impact2D'+'group'+'UNR'].GetYaxis().GetBinLabel(g)
                        impactVals[impBin+nuiGroup] = self.histos[suff+'impact2D'+'group'+'UNR'].GetBinContent(qty, g)
            except :
                skipDiffImpact = True
                print('missing differential impact')
                
            skipIntImpact = False
            try : 
                for qty in range(1, self.histos[suff+'impact2D'+'group'+'int'].GetNbinsX()+1):    
                    impBin = self.histos[suff+'impact2D'+'group'+'int'].GetXaxis().GetBinLabel(qty)          
                    for g in range(1, self.histos[suff+'impact2D'+'group'+'int'].GetNbinsY()+1):
                        nuiGroup=self.histos[suff+'impact2D'+'group'+'int'].GetYaxis().GetBinLabel(g)
                        impactVals[impBin+nuiGroup] = self.histos[suff+'impact2D'+'group'+'int'].GetBinContent(qty, g)
            except :
                skipIntImpact = True
                print('missing integrated impact')
            
            skipMassImpact = False
            try :    
                for g in range(1,self.histos[suff+'impact2D'+'group'+'mass'].GetNbinsY()+1):
                    nuiGroup=self.histos[suff+'impact2D'+'group'+'mass'].GetYaxis().GetBinLabel(g)
                    impactVals['mass'+nuiGroup] = self.histos[suff+'impact2D'+'group'+'mass'].GetBinContent(1, g)
            except :
                skipMassImpact = True
                print('missing mass impact')
            
            ###### This for rebin only if variable width histogram
            # unp = self.coeffList[-1]
            # varWidth_X = False 
            # for i in range(1, len(self.yArr)-1) :
            #         for j in range(i+1,len(self.yArr)) :
            #             if abs(self.histos[suff+'FitAC'+unp].GetXaxis().GetBinWidth(j)-self.histos[suff+'FitAC'+unp].GetXaxis().GetBinWidth(i))>0.0001:
            #                 varWidth_X = True
            #                 break
            # varWidth_Y = False
            # for i in range(1, len(self.qtArr)-1) :
            #     for j in range(i+1,len(self.qtArr)) :
            #         if abs(self.histos[suff+'FitAC'+unp].GetYaxis().GetBinWidth(j)-self.histos[suff+'FitAC'+unp].GetYaxis().GetBinWidth(i))>0.0001:
            #             varWidth_Y = True
            #             break
            # varWidth_qt = False
            # for i in range(1, len(self.qtArr)-1) :
            #     for j in range(i+1,len(self.qtArr)) :
            #         if abs(self.histos[suff+'FitACqt'+unp].GetXaxis().GetBinWidth(j)-self.histos[suff+'FitACqt'+unp].GetXaxis().GetBinWidth(i))>0.0001:
            #             varWidth_qt = True
            #             break
            # varWidth_y = False
            # for i in range(1, len(self.yArr)-1) :
            #     for j in range(i+1,len(self.yArr)) :
            #         if abs(self.histos[suff+'FitACy'+unp].GetXaxis().GetBinWidth(j)-self.histos[suff+'FitACy'+unp].GetXaxis().GetBinWidth(i))>0.0001:
            #             varWidth_y = True
            #             break
            
                     
            for c in self.coeffDict:
                for nui in self.nuisanceDict:
                    
                    self.histos[suff+'impact'+'UNRqty'+c+nui] = ROOT.TH1D('impact_'+c+'_'+nui+'_UNRqty', 'impact_'+c+'_'+nui+'_UNRqty', len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
                    for q in range(1, len(self.qtArr)) :
                        if skipDiffImpact : continue 
                        for y in range(1, len(self.yArr)): 
                            indexUNRqty = (q-1)*(len(self.yArr)-1)+y
                            if self.anaKind['angNames'] : impBin = 'y_'+str(y)+'_qt_'+str(q)+'_'+c
                            else : impBin = 'helXsecs'+c+'_'+'y_'+str(y)+'_qt_'+str(q)+'_'+self.anaKind['diffString']  
                            if nui!='tot' :                          
                                relImp =   abs(impactVals[impBin+nui])
                            else :
                                relImp = self.histos[suff+'FitAC'+c].GetBinError(y, q)
                            if ('unpol' in c  or (not self.anaKind['differential'])) :
                                relImp = relImp/abs(self.histos[suff+'FitAC'+c].GetBinContent(y, q))                   
                            if (('unpol' in c  or (not self.anaKind['differential'])) and nui!='tot') :
                                relImp = relImp/(self.lumi*3./16./math.pi)
                                relImp = relImp/self.histos[suff+'FitAC'+c].GetYaxis().GetBinWidth(q)/self.histos[suff+'FitAC'+c].GetXaxis().GetBinWidth(y)   
                            self.histos[suff+'impact'+'UNRqty'+c+nui].SetBinContent(indexUNRqty,relImp)
                    
                    self.histos[suff+'impact'+'UNRyqt'+c+nui] = ROOT.TH1D('impact_'+c+'_'+nui+'_UNRyqt', 'impact_'+c+'_'+nui+'_UNRyqt', len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
                    for y in range(1, len(self.yArr)): 
                        if skipDiffImpact : continue 
                        for q in range(1, len(self.qtArr)) :
                            indexUNRyqt = (y-1)*(len(self.qtArr)-1)+q
                            if self.anaKind['angNames'] : impBin = 'y_'+str(y)+'_qt_'+str(q)+'_'+c
                            else : impBin = 'helXsecs'+c+'_'+'y_'+str(y)+'_qt_'+str(q)+'_'+self.anaKind['diffString']          
                            if nui!='tot' :                   
                                relImp =   abs(impactVals[impBin+nui])
                            else :
                                relImp = self.histos[suff+'FitAC'+c].GetBinError(y, q)
                            if ('unpol' in c  or (not self.anaKind['differential'])) :
                                relImp = relImp/abs(self.histos[suff+'FitAC'+c].GetBinContent(y, q))                    
                            if (('unpol' in c  or (not self.anaKind['differential'])) and nui!='tot') :
                                relImp = relImp/(self.lumi*3./16./math.pi)
                                relImp = relImp/self.histos[suff+'FitAC'+c].GetYaxis().GetBinWidth(q)/self.histos[suff+'FitAC'+c].GetXaxis().GetBinWidth(y)   
                            self.histos[suff+'impact'+'UNRyqt'+c+nui].SetBinContent(indexUNRyqt,relImp)
                    
            
                    self.histos[suff+'impact'+'y'+c+nui] = ROOT.TH1D('impact_'+c+'_'+nui+'_y', 'impact_'+c+'_'+nui+'_y', len(self.yArr)-1, array('f',self.yArr))
                    for y in range(1, len(self.yArr)):   
                        if skipIntImpact : continue  
                        if self.anaKind['angNames'] :  impBin = 'y_'+str(y)+'_helmeta_'+c
                        else : impBin = 'helXsecs'+c+'_'+'y_'+str(y)+'_'+self.anaKind['intString'].replace('sumpois','sumxsec')   
                        if nui!='tot' :
                            relImp = abs(impactVals[impBin+nui])
                        else :
                            relImp = abs(self.histos[suff+'FitACy'+c].GetBinError(y))
                        if ('unpol' in c  or (not self.anaKind['differential'])) :
                              relImp = relImp/abs(self.histos[suff+'FitACy'+c].GetBinContent(y))
                        if (('unpol' in c  or (not self.anaKind['differential'])) and nui!='tot') :
                            relImp = relImp/(self.lumi*3./16./math.pi) 
                            relImp = relImp/self.histos[suff+'FitACy'+c].GetXaxis().GetBinWidth(y)  
                        self.histos[suff+'impact'+'y'+c+nui].SetBinContent(y,relImp)
                    
                    self.histos[suff+'impact'+'qt'+c+nui] = ROOT.TH1D('impact_'+c+'_'+nui+'_qt', 'impact_'+c+'_'+nui+'_qt', len(self.qtArr)-1, array('f',self.qtArr))
                    for qt in range(1, len(self.qtArr)) :
                        if skipIntImpact : continue  
                        if self.anaKind['angNames'] :  impBin = 'qt_'+str(qt)+'_helmeta_'+c
                        else : impBin = 'helXsecs'+c+'_'+'qt_'+str(qt)+'_'+self.anaKind['intString'].replace('sumpois','sumxsec')   
                        if nui!='tot' :
                            relImp = abs(impactVals[impBin+nui])
                        else :
                            relImp = abs(self.histos[suff+'FitACqt'+c].GetBinError(qt))
                        if ('unpol' in c  or (not self.anaKind['differential'])) :
                            relImp = relImp/abs(self.histos[suff+'FitACqt'+c].GetBinContent(qt))
                        if (('unpol' in c  or (not self.anaKind['differential'])) and nui!='tot') :
                            relImp = relImp/(self.lumi*3./16./math.pi) 
                            relImp = relImp/self.histos[suff+'FitACqt'+c].GetXaxis().GetBinWidth(qt) 
                        self.histos[suff+'impact'+'qt'+c+nui].SetBinContent(qt,relImp)
            
            #print("WARNING: used hardcoded value of the Wmass for the impact, mW=",self.mass," GeV")    
            print("WARNING: used hardcoded value of the up/down mass variation (50 MeV)")    
            for nui in self.nuisanceDict:
                self.histos[suff+'impact'+'mass'+nui] = ROOT.TH1D('impact_mass_'+nui, 'impact_mass_'+nui, 1, 0,1)
                if skipMassImpact : continue
                if nui=='tot' : continue     
                if nui=='mass ': continue 
                # relImp = abs(impactVals['mass'+nui]/self.histos[suff+'mass'].GetBinContent(1))
                # relImp = abs(impactVals['mass'+nui]/self.mass)
                relImp = abs(impactVals['mass'+nui]*50) #in MeV, since the weight is +/- is 50MeV.
                self.histos[suff+'impact'+'mass'+nui].SetBinContent(1,relImp) 
                self.histos[suff+'impact'+'mass'+nui].SetBinError(1,0.0000000001) 
            
            self.histos[suff+'impact'+'mass'+'tot'].SetBinContent(1,self.histos[suff+'mass'].GetBinError(1)*50)
            self.histos[suff+'impact'+'mass'+'tot'].SetBinError(1,0.0000000001)
            
            if massImp!='' :
                fileMassStat = ROOT.TFile.Open(massImp)
                treeMassStat = fileMassStat.Get('fitresults')
                for ev in treeMassStat:
                    mWstat = eval('ev.mass_err')
                self.histos[suff+'impact'+'mass'+'stat'].SetBinContent(1,mWstat*50)    
            
            
            
            
        #########sqrt(N) plot as additional nuisances to extra comparison
        #building of acceptance maps integrated
        # self.histos[suff+'acceptance'+'sumOfTemplate'+'y'] = self.histos[suff+'acceptance'+'sumOfTemplate'].ProjectionX('acceptance'+'sumOfTemplate'+'y')
        # self.histos[suff+'acceptance'+'sumOfTemplate'+'qt'] = self.histos[suff+'acceptance'+'sumOfTemplate'].ProjectionY('acceptance'+'sumOfTemplate'+'qt')
        self.histos[suff+'acceptance'+'y'] = self.histos[suff+'acceptance'+'sumOfTemplate'].ProjectionX('acceptance'+'y')
        self.histos[suff+'acceptance'+'qt'] = self.histos[suff+'acceptance'+'sumOfTemplate'].ProjectionY('acceptance'+'qt')
        self.histos[suff+'acceptance'+'y'].Divide(self.histos[suff+'MCy'+'acceptance'])
        self.histos[suff+'acceptance'+'qt'].Divide(self.histos[suff+'MCqt'+'acceptance'])
        
        c='unpolarizedxsec'
        self.histos[suff+'impact'+'UNRqty'+c+'poiss'] = ROOT.TH1D('impact_'+c+'_'+'poiss'+'_UNRqty', 'impact_'+c+'_'+'poiss'+'_UNRqty', len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
        for q in range(1, len(self.qtArr)) :
            if impact and skipDiffImpact : continue 
            if not self.anaKind['angNames'] : continue 
            for y in range(1, len(self.yArr)): 
                indexUNRqty = (q-1)*(len(self.yArr)-1)+y
                # nev = self.histos[suff+'FitAC'+c].GetBinContent(y, q)* (self.lumi*3./16./math.pi) *self.histos[suff+'FitAC'+c].GetYaxis().GetBinWidth(q)*self.histos[suff+'FitAC'+c].GetXaxis().GetBinWidth(y)  
                nev = self.histos[suff+'FitAC'+c].GetBinContent(y, q)* self.lumi *self.histos[suff+'FitAC'+c].GetYaxis().GetBinWidth(q)*self.histos[suff+'FitAC'+c].GetXaxis().GetBinWidth(y)*self.histos[suff+'acceptance'].GetBinContent(y,q)
                if nev<0 :
                    print("WARNING: negative unpolarized cross section") 
                    nev = 1
                if nev!=0 : poisserr = math.sqrt(nev)/nev 
                else : poisserr = 0 
                self.histos[suff+'impact'+'UNRqty'+c+'poiss'].SetBinContent(indexUNRqty,poisserr)
                    
        
        self.histos[suff+'impact'+'UNRyqt'+c+'poiss'] = ROOT.TH1D('impact_'+c+'_'+'poiss'+'_UNRyqt', 'impact_'+c+'_'+'poiss'+'_UNRyqt', len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
        for y in range(1, len(self.yArr)): 
                if impact and skipDiffImpact : continue 
                if not self.anaKind['angNames'] : continue
                for q in range(1, len(self.qtArr)) :
                    indexUNRyqt = (y-1)*(len(self.qtArr)-1)+q
                    # nev = self.histos[suff+'FitAC'+c].GetBinContent(y, q)* (self.lumi*3./16./math.pi) *self.histos[suff+'FitAC'+c].GetYaxis().GetBinWidth(q)*self.histos[suff+'FitAC'+c].GetXaxis().GetBinWidth(y)  
                    nev = self.histos[suff+'FitAC'+c].GetBinContent(y, q)*self.lumi*self.histos[suff+'FitAC'+c].GetYaxis().GetBinWidth(q)*self.histos[suff+'FitAC'+c].GetXaxis().GetBinWidth(y)*self.histos[suff+'acceptance'].GetBinContent(y,q)
                    if nev<0 : 
                        print("WARNING: negative unpolarized cross section")
                        nev = 1
                    if nev!=0 : poisserr = math.sqrt(nev)/nev 
                    else : poisserr = 0 
                    self.histos[suff+'impact'+'UNRyqt'+c+'poiss'].SetBinContent(indexUNRyqt,poisserr)
        
        self.histos[suff+'impact'+'y'+c+'poiss'] = ROOT.TH1D('impact_'+c+'_'+'poiss'+'_y', 'impact_'+c+'_'+'poiss'+'_y', len(self.yArr)-1, array('f',self.yArr))
        for y in range(1, len(self.yArr)):   
            if impact and skipIntImpact : continue
            if not self.anaKind['angNames'] : continue
            # nev = self.histos[suff+'FitACy'+c].GetBinContent(y)* (self.lumi*3./16./math.pi) *self.histos[suff+'FitACy'+c].GetXaxis().GetBinWidth(y)  
            nev = self.histos[suff+'FitACy'+c].GetBinContent(y) *self.lumi*self.histos[suff+'FitACy'+c].GetXaxis().GetBinWidth(y)*self.histos[suff+'acceptance'+'y'].GetBinContent(y)
            if nev<0 :
                print("WARNING: negative unpolarized cross section") 
                nev = 1
            if nev!=0 : poisserr = math.sqrt(nev)/nev 
            else : poisserr = 0 
            self.histos[suff+'impact'+'y'+c+'poiss'].SetBinContent(y,poisserr)
        
        self.histos[suff+'impact'+'qt'+c+'poiss'] = ROOT.TH1D('impact_'+c+'_'+'poiss'+'_qt', 'impact_'+c+'_'+'poiss'+'_qt', len(self.qtArr)-1, array('f',self.qtArr))
        for qt in range(1, len(self.qtArr)) :
            if impact and skipIntImpact : continue  
            if not self.anaKind['angNames'] : continue
            # nev = self.histos[suff+'FitACqt'+c].GetBinContent(qt)* (self.lumi*3./16./math.pi) *self.histos[suff+'FitACqt'+c].GetXaxis().GetBinWidth(qt)  
            nev = self.histos[suff+'FitACqt'+c].GetBinContent(qt)*self.lumi*self.histos[suff+'FitACqt'+c].GetXaxis().GetBinWidth(qt) *self.histos[suff+'acceptance'+'qt'].GetBinContent(qt) 
            if nev<0 : 
                print("WARNING: negative unpolarized cross section")
                nev = 1
            if nev!=0 : poisserr = math.sqrt(nev)/nev 
            else : poisserr = 0 
            self.histos[suff+'impact'+'qt'+c+'poiss'].SetBinContent(qt,poisserr)        
                
            # print("DEBUG impact (sum in quadrature")
            # sumOfPdf = 0.
            # labelX = self.histos[suff+'impact2D'+'UNR'].GetXaxis().GetBinLabel(2)
            # unpol_divider = self.histos[suff+'FitAC'+'unpolarizedxsec'].GetBinContent(1, 2)/self.histos[suff+'FitAC'+'unpolarizedxsec'].GetYaxis().GetBinWidth(2)/(self.lumi*3./16./math.pi)
            # # unpol_divider=1
            # for y in range(0,self.histos[suff+'impact2D'+'UNR'].GetNbinsY()+1) :
            #     labelY = self.histos[suff+'impact2D'+'UNR'].GetYaxis().GetBinLabel(y)
            #     val = self.histos[suff+'impact2D'+'UNR'].GetBinContent(2,y)
            #     val = val/unpol_divider
            #     # if "unclustEn" in labelY or "jes" in labelY:
            #     if "Pdf" in labelY:
            #         sumOfPdf = sumOfPdf+val**2
            #         # sumOfPdf = sumOfPdf+abs(val)
            #         # print(y, labelY,labelX,val)
            # sumOfPdf = math.sqrt(sumOfPdf)
            # valGroupPdf = self.histos[suff+'impact2D'+'group'+'UNR'].GetBinContent(self.histos[suff+'impact2D'+'group'+'UNR'].GetXaxis().FindBin(labelX),self.histos[suff+'impact2D'+'group'+'UNR'].GetYaxis().FindBin("pdfs"))
            # # print("sum of PDF=",sumOfPdf) 
            # # print("grouped=",valGroupPdf/unpol_divider)
            # # print("value",  self.histos[suff+'impact'+'UNRqty'+'unpolarizedxsec'+'pdfs'].GetBinContent(7))
            
            # print("DEBUG impact: max pdf") #josh debug
            # for x in range(1,self.histos[suff+'impact2D'+'UNR'].GetNbinsX()+1) :
            #     labelX = self.histos[suff+'impact2D'+'UNR'].GetXaxis().GetBinLabel(x)
            #     if not 'A4' in labelX : continue
            #     maxCounter=0 
            #     maxBinY = 0
            #     maxBinqT=0 
            #     maxName = ''
            #     quadSum = 0
            #     for x2 in range(1,self.histos[suff+'impact2D'+'UNR'].GetNbinsY()+1) :
            #         labelY = self.histos[suff+'impact2D'+'UNR'].GetYaxis().GetBinLabel(x2)
            #         if "Pdf" in labelY: 
            #             impVal = abs(self.histos[suff+'impact2D'+'UNR'].GetBinContent(x,x2)  )
            #             y = int(labelX.split('_')[1])
            #             qt = int(labelX.split('_')[3])
            #             impVal = impVal/abs(self.histos[suff+'FitAC'+'A4'].GetBinContent(y, qt))  
            #             quadSum+=impVal**2
            #             # if labelY=='LHEPdfWeightHess3'  : print("inside", impVal, maxCounter)
            #             if abs(impVal) > maxCounter :
            #                 maxCounter = abs(impVal)
            #                 maxBinqtY = labelX
            #                 maxName = labelY
            #             refVal = abs(self.histos[suff+'impact2D'+'UNR'].GetBinContent(x,self.histos[suff+'impact2D'+'UNR'].GetYaxis().FindBin("LHEPdfWeightHess3")))
            #             refVal = abs(refVal/self.histos[suff+'FitAC'+'A4'].GetBinContent(y, qt))   
            #             # refVal = (MCVal - refVal)**2/(MCVal**2)   
            #             if labelX=='y_1_qt_1_A4' :
            #                 print(labelX, labelY, impVal,impVal*abs(self.histos[suff+'FitAC'+'A4'].GetBinContent(y, qt)), x,x2, y, qt, self.histos[suff+'FitAC'+'A4'].GetBinContent(y, qt))
            #     print("Impact: name=",maxName, "value=", maxCounter, "bin qty=", maxBinqtY, "hess3=", refVal, "quadSum=",math.sqrt(quadSum))#"val3/val=",refVal/maxCounter)
        
        #------------- prefit and postfit folding -------------#
        if postfit : 
            for p,pval in self.processDict.items():
                for k in ['prefit','postfit'] :
                    
                    
                    self.histos[suff+k+p+'2D'] = ROOT.TH2F(self.histos[suff+k+p].GetName()+'_2D',self.histos[suff+k+p].GetName()+'_2D', len(self.etaArr)-1, array('f',self.etaArr), len(self.ptArr)-1, array('f',self.ptArr))
                    for ieta in range(1,self.histos[suff+k+p+'2D'].GetNbinsX()+1) :
                        for ipt in range(1,self.histos[suff+k+p+'2D'].GetNbinsY()+1) :
                            # indexUNRetpt = self.etaArr.index(float(ieta))*(len(self.ptArr)-1)+self.ptArr.index(float(ipt))
                            indexUNRetpt = (ieta-1)*(len(self.ptArr)-1)+(ipt-1)
                            self.histos[suff+k+p+'2D'].SetBinContent(ieta,ipt, self.histos[suff+k+p].GetBinContent(indexUNRetpt+1))
                            self.histos[suff+k+p+'2D'].SetBinError(ieta,ipt, self.histos[suff+k+p].GetBinError(indexUNRetpt+1))
                    
                    #projections
                    self.histos[suff+k+p+'eta'] = self.histos[suff+k+p+'2D'].ProjectionX(self.histos[suff+k+p].GetName()+'_eta',0,-1)
                    self.histos[suff+k+p+'eta'].SetTitle(self.histos[suff+k+p].GetName()+'_eta')
                    # self.histos[suff+k+p+'eta'].Scale(1/self.histos[suff+k+p+'eta'].GetBinWidth(1))
                    self.histos[suff+k+p+'pt'] = self.histos[suff+k+p+'2D'].ProjectionY(self.histos[suff+k+p].GetName()+'_pt',0,-1)
                    self.histos[suff+k+p+'pt'].SetTitle(self.histos[suff+k+p].GetName()+'_pt')  
                    # self.histos[suff+k+p+'pt'].Scale(1/self.histos[suff+k+p+'pt'].GetBinWidth(1))   
                    # self.histos[suff+k+p].Scale(1/(self.histos[suff+k+p].GetXaxis().GetBinWidth(1)*self.histos[suff+k+p].GetYaxis().GetBinWidth(1)))       
            
            #------------ prefit and postfit stacked and ratio band ---------------#
            for k in ['prefit','postfit'] :
                for var,varInfo in self.postFitVars.items() :
                            
                    self.leg[suff+k+var] = ROOT.TLegend(0.5, 0.55, 0.90, 0.85)            
                    self.histos[suff+k+var+'stack'] = ROOT.THStack('stack_'+suff+k+var, 'stack_'+suff+k+var)
                    # hSum = self.histos[suff+k+'full'+var] #COMMENT ME!!!
                    # hData = self.histos[suff+k+'obs'+var] #COMMENT ME!!!

                    #build stacked plot
                    for p in self.processOrder :
                        self.histos[suff+k+p+var].SetFillStyle(1001)
                        self.histos[suff+k+p+var].SetLineWidth(1)
                        self.histos[suff+k+p+var].SetFillColor(self.processDict[p][0])
                        # self.histos[suff+k+p+var].SetLineColor(self.processDict[p][0])
                        self.histos[suff+k+var+'stack'].Add(self.histos[suff+k+p+var])
                        # hSum.Add(self.histos[suff+k+p+var])
                    self.leg[suff+k+var].AddEntry(self.histos[suff+k+'obs'+var], 'Data', "PE1")
                    for p in self.processOrder[::-1] :
                        if p =='sig' and 'minus' in inputFile: #signal has different color and legend entry
                            self.histos[suff+k+p+var].SetFillColor(self.processDict[p][2])
                            legTag = self.processDict[p][3]
                        else :
                            legTag = self.processDict[p][1]
                        self.leg[suff+k+var].AddEntry(self.histos[suff+k+p+var],legTag,  "f")                    
                    
                    #build ratio plot Data/pred
                    self.histos[suff+k+var+'ratio'] = self.histos[suff+k+'obs'+var].Clone('hRatio_'+var)
                    self.histos[suff+k+var+'ratioBand'] = self.histos[suff+k+'full'+var].Clone('hRatio_'+var)

                    if k=='prefit' and data!='' :
                        for xx in range(1,self.histos[suff+k+var+'ratio'].GetNbinsX()+1) :
                            self.histos[suff+k+var+'ratio'].SetBinContent(xx, self.histos[suff+k+var+'ratio'].GetBinContent(xx)/self.histos[suff+k+'full'+var].GetBinContent(xx) )
                            self.histos[suff+k+var+'ratio'].SetBinError(xx, self.histos[suff+k+var+'ratio'].GetBinError(xx)/self.histos[suff+k+'full'+var].GetBinContent(xx) )
                            self.histos[suff+k+var+'ratioBand'].SetBinError(xx,0)
                        self.histos[suff+k+var+'ratioBand'].Divide(self.histos[suff+k+'full'+var])
                    else : 
                        self.histos[suff+k+var+'ratio'].Divide(self.histos[suff+k+'full'+var])
                        self.histos[suff+k+var+'ratioBand'] =  self.histos[suff+k+var+'ratio'].Clone(self.histos[suff+k+var+'ratio'].GetName()+'Band')

                                
        #-----------------------------------------------------------------------------------------------------------------------------------------#
        #-----------------------------------------------------------------------------------------------------------------------------------------#
        #-----------------------------------------------------------------------------------------------------------------------------------------#
        #-----------------------------------------------------------------------------------------------------------------------------------------#
        #----------------------------------- CANVAS PREPARATION --------------------------------------------------------------------------------- #  
        #-----------------------------------------------------------------------------------------------------------------------------------------#
        #-----------------------------------------------------------------------------------------------------------------------------------------#
        
        cmslab = "#bf{CMS} #scale[0.7]{#it{Simulation Work in progress}}"
        lumilab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
        cmsLatex = ROOT.TLatex()
        qtLab = "q_{T}^{W}<60 GeV"
        yLab = "|Y_{W}|<2.4"
            
        #--------------- Ai 1D canvas -------------------#
        for c in self.coeffDict:
            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                self.canvas[suff+'FitAC'+'qt'+str(i)+c] = ROOT.TCanvas(suff+"_c_projQt{}_{}".format(i,c),suff+"_c_projQt{}_{}".format(i,c),1200,1050)
                self.leg[suff+'FitAC'+'qt'+str(i)+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
                self.canvas[suff+'FitAC'+'qt'+str(i)+c].cd()
                self.canvas['ph'+suff+'FitAC'+'qt'+str(i)+c] = ROOT.TPad('ph'+suff+"_c_projQt{}_{}".format(i,c), 'ph'+suff+"_c_projQt{}_{}".format(i,c),0,0.3,1,1)
                self.canvas['pr'+suff+'FitAC'+'qt'+str(i)+c] = ROOT.TPad('pr'+suff+"_c_projQt{}_{}".format(i,c),  'pr'+suff+"_c_projQt{}_{}".format(i,c), 0,0,1,0.3)
                self.canvas['ph'+suff+'FitAC'+'qt'+str(i)+c].SetBottomMargin(0.02)
                self.canvas['ph'+suff+'FitAC'+'qt'+str(i)+c].Draw()
                self.canvas['pr'+suff+'FitAC'+'qt'+str(i)+c].SetTopMargin(0.01)
                self.canvas['pr'+suff+'FitAC'+'qt'+str(i)+c].SetBottomMargin(0.4)
                self.canvas['pr'+suff+'FitAC'+'qt'+str(i)+c].Draw()
                
                self.canvas['ph'+suff+'FitAC'+'qt'+str(i)+c].cd()
                self.canvas['ph'+suff+'FitAC'+'qt'+str(i)+c].SetGridx()
                self.canvas['ph'+suff+'FitAC'+'qt'+str(i)+c].SetGridy()
                self.histos[suff+'FitAC'+'qt'+str(i)+c] = self.histos[suff+'FitAC'+c].ProjectionY("projQt{}_{}".format(i,c),i,i)
                self.histos[suff+'FitBand'+'qt'+str(i)+c] = self.histos[suff+'FitBand'+c].ProjectionY("projQt{}_{}_err".format(i,c),i,i)
                # self.histos[suff+'FitBandPDF'+'qt'+str(i)+c] = self.histos[suff+'FitBandPDF'+c].ProjectionY("projQt{}_{}_errPDF".format(i,c),i,i)
                # self.histos[suff+'FitBandScale'+'qt'+str(i)+c] = self.histos[suff+'FitBandScale'+c].ProjectionY("projQt{}_{}_errScale".format(i,c),i,i)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].SetTitle(self.coeffDict[c][2]+", transverse momentum distribution, "+str(self.yArr[i-1])+"<|Y_{W}|<"+str(self.yArr[i])+", "+self.signDict[self.sign])
                self.histos[suff+'FitAC'+'qt'+str(i)+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
                self.histos[suff+'FitAC'+'qt'+str(i)+c].SetLineWidth(2)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].SetLineColor(ROOT.kBlack)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].SetStats(0)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].SetMarkerStyle(20)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].SetMarkerSize(0.5)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].Draw()
                self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetTitleSize(0.06)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].SetLabelSize(0.05,'y')
                if not 'unpol' in c and self.anaKind['name']=='angCoeff':
                    # self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetRangeUser(-4,4)
                    self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetTitle(self.coeffDict[c][2])
                elif self.anaKind['name']=='helXsec':
                     self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetTitle('d#sigma_'+c+'/dq_{T}^{W}d|Y_{W}| [fb/GeV]')
                elif self.anaKind['name']=='normXsec':
                    self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetTitle('d#sigma_'+c+'/dq_{T}^{W}d|Y_{W}|/#sigma_{tot} [1/GeV]')
                elif self.anaKind['name']=='mu':
                    self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetTitle('#mu_'+c)
                else :
                    self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetTitle('d#sigma^{U+L}/dq_{T}^{W}d|Y_{W}| [fb/GeV]')
                    # maxvalMain = self.histos[suff+'FitAC'+'qt'+str(i)+c].GetMaximum()
                    # self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetRangeUser(0,maxvalMain*1.5)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].SetFillColor(ROOT.kOrange)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].SetFillStyle(0)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].SetLineColor(ROOT.kOrange)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].SetLineWidth(2)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetTitleOffset(0.8)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].GetXaxis().SetTitleOffset(3)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].GetXaxis().SetLabelOffset(3)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].DrawCopy("E2 same")
                self.histos[suff+'FitBand'+'qt'+str(i)+c].SetFillStyle(3001)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].Draw("E2 same")
                self.histos[suff+'FitAC'+'qt'+str(i)+c].DrawCopy("E1 same") #to have foreground
                
                if aposteriori!='' :
                    self.histos[suff+'apo'+'qt'+str(i)+c] = self.histos[suff+'apo'+c].ProjectionY('apo'+"projQt{}_{}".format(i,c),i,i)
                    self.histos[suff+'apo'+'qt'+str(i)+c].SetLineColor(ROOT.kAzure+1)
                    self.histos[suff+'apo'+'qt'+str(i)+c].SetMarkerColor(ROOT.kBlue-4)
                    self.histos[suff+'apo'+'qt'+str(i)+c].SetLineWidth(5)
                    self.histos[suff+'apo'+'qt'+str(i)+c].SetMarkerStyle(20)
                    self.histos[suff+'apo'+'qt'+str(i)+c].Draw("EX0 same")
                    self.histos[suff+'apo'+'qt'+str(i)+c].SetMarkerSize(1)
                
                ROOT.gPad.Update()
                if 'unpol' in c and self.anaKind['angNames']: minvalMain = 0
                else : minvalMain = self.canvas['ph'+suff+'FitAC'+'qt'+str(i)+c].GetFrame().GetY1()
                self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetRangeUser(minvalMain,1.3*self.canvas['ph'+suff+'FitAC'+'qt'+str(i)+c].GetFrame().GetY2())
                self.leg[suff+'FitAC'+'qt'+str(i)+c].Draw("same")
                
                self.canvas['pr'+suff+'FitAC'+'qt'+str(i)+c].cd()
                self.canvas['pr'+suff+'FitAC'+'qt'+str(i)+c].SetGridx()
                self.canvas['pr'+suff+'FitAC'+'qt'+str(i)+c].SetGridy()
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c] = self.histos[suff+'FitRatioAC'+c].ProjectionY("projQt{}_{}_RatioAC".format(i,c),i,i)
                self.histos[suff+'FitRatio'+'qt'+str(i)+c] = self.histos[suff+'FitRatio'+c].ProjectionY("projQt{}_{}_Ratio".format(i,c),i,i)
                for d in self.dirList :
                    self.histos[suff+'FitRatioPDF'+'qt'+str(i)+d+c] = self.histos[suff+'FitRatioPDF'+d+c].ProjectionY("projQt{}_{}_RatioPDF".format(i,c)+d,i,i)
                    self.histos[suff+'FitRatioScale'+'qt'+str(i)+d+c] = self.histos[suff+'FitRatioScale'+d+c].ProjectionY("projQt{}_{}_RatioScale".format(i,c)+d,i,i)
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
                if ('unpol' in c or (not self.anaKind['differential'])): 
                    self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].GetYaxis().SetTitle('Pred./Theory')
                else :
                    self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].GetYaxis().SetTitle('Pred.-Theory')
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].SetMarkerStyle(20)
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].SetMarkerSize(0.5)
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].SetLineWidth(2)
                self.histos[suff+'FitRatio'+'qt'+str(i)+c].SetLineWidth(2)
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].SetTitle("")
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].GetYaxis().SetTitleOffset(0.32)
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].SetTitleSize(0.13,'y')
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].SetLabelSize(0.10,'y')
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].GetXaxis().SetTitleOffset(1.2)
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].SetTitleSize(0.12,'x')
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].SetLabelSize(0.12,'x')
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].SetStats(0)
                for d in self.dirList :
                    self.histos[suff+'FitRatioPDF'+'qt'+str(i)+d+c].SetLineWidth(2)
                    self.histos[suff+'FitRatioScale'+'qt'+str(i)+d+c].SetLineWidth(2)
                self.histos[suff+'FitRatio'+'qt'+str(i)+c].SetLineColor(ROOT.kOrange)
                for d in self.dirList :
                    self.histos[suff+'FitRatioPDF'+'qt'+str(i)+d+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                    self.histos[suff+'FitRatioScale'+'qt'+str(i)+d+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
                self.histos[suff+'FitRatio'+'qt'+str(i)+c].SetFillColor(ROOT.kOrange)
                self.histos[suff+'FitRatio'+'qt'+str(i)+c].SetFillStyle(3001)
                
                maxvalRatio=0
                for xx in range(1,self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].GetNbinsX()+1) :
                    if self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].GetBinError(xx)>maxvalRatio :
                        maxvalRatio= self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].GetBinError(xx)
                    # self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].SetBinError(xx,0)
                if self.anaKind['angNames'] :
                    PadYCut =  self.RatioPadYcutDict[c][0]
                else : PadYCut = self.RatioPadYcut 
                if maxvalRatio>PadYCut :
                    maxvalRatio=PadYCut
                # self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].GetYaxis().SetRangeUser(-maxvalRatio*1.1+1,maxvalRatio*1.1+1)
                
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].Draw("EX0")
                self.histos[suff+'FitRatio'+'qt'+str(i)+c].DrawCopy("same E2")
                for d in self.dirList :
                    self.histos[suff+'FitRatioPDF'+'qt'+str(i)+d+c].Draw("hist same")
                    self.histos[suff+'FitRatioScale'+'qt'+str(i)+d+c].Draw("hist same")
                self.histos[suff+'FitRatioAC'+'qt'+str(i)+c].Draw("EX0 same") #to have foreground
                
                if aposteriori!='' :
                    self.histos[suff+'apoRatio'+'qt'+str(i)+c] = self.histos[suff+'apo'+'Ratio'+c].ProjectionY("projQt{}_{}_apoRatio".format(i,c),i,i)
                    self.histos[suff+'apoRatio'+'qt'+str(i)+c].SetLineColor(ROOT.kAzure+1)
                    self.histos[suff+'apoRatio'+'qt'+str(i)+c].SetMarkerColor(ROOT.kBlue-4)
                    self.histos[suff+'apoRatio'+'qt'+str(i)+c].SetLineWidth(4)
                    self.histos[suff+'apoRatio'+'qt'+str(i)+c].SetMarkerStyle(20)
                    self.histos[suff+'apoRatio'+'qt'+str(i)+c].Draw("EX0 same")
                    self.histos[suff+'apoRatio'+'qt'+str(i)+c].SetMarkerSize(1)
                
                self.leg[suff+'FitAC'+'qt'+str(i)+c].AddEntry(self.histos[suff+'FitAC'+'qt'+str(i)+c], self.signDict[self.sign]+' Fit')
                self.leg[suff+'FitAC'+'qt'+str(i)+c].AddEntry(self.histos[suff+'FitBand'+'qt'+str(i)+c], "MC Syst. Unc.")
                self.leg[suff+'FitAC'+'qt'+str(i)+c].AddEntry(self.histos[suff+'FitRatioPDF'+'qt'+str(i)+'up'+c], self.groupedSystColors['LHEPdfWeightVars'][1])
                self.leg[suff+'FitAC'+'qt'+str(i)+c].AddEntry(self.histos[suff+'FitRatioScale'+'qt'+str(i)+'up'+c], self.groupedSystColors['LHEScaleWeightVars'][1])
                self.leg[suff+'FitAC'+'qt'+str(i)+c].SetNColumns(2)
                if aposteriori!='' :
                    self.leg[suff+'FitAC'+'qt'+str(i)+c].AddEntry(self.histos[suff+'apo'+'qt'+str(i)+c], "Post-fit regularized")
                
                #Additional block for scale-pdf splitted
                # self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                # self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].SetFillStyle(0)
                # self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                # self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].SetLineWidth(2)
                # self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].DrawCopy("E2 same")
                # self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].SetFillStyle(3004)
                # self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].Draw("E2 same")
                # self.histos[suff+'FitBandScale'+'qt'+str(i)+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])
                # self.histos[suff+'FitBandScale'+'qt'+str(i)+c].SetFillStyle(0)
                # self.histos[suff+'FitBandScale'+'qt'+str(i)+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
                # self.histos[suff+'FitBandScale'+'qt'+str(i)+c].SetLineWidth(2)
                # self.histos[suff+'FitBandScale'+'qt'+str(i)+c].DrawCopy("E2 same")
                # self.histos[suff+'FitBandScale'+'qt'+str(i)+c].SetFillStyle(3005)
                # self.histos[suff+'FitBandScale'+'qt'+str(i)+c].Draw("E2 same")
                
                
            for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                self.canvas[suff+'FitAC'+'y'+str(j)+c] = ROOT.TCanvas(suff+"_c_projY{}_{}".format(j,c),suff+"_c_projY{}_{}".format(j,c),1200,1050)
                self.leg[suff+'FitAC'+'y'+str(j)+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
                self.canvas[suff+'FitAC'+'y'+str(j)+c].cd()
                self.canvas['ph'+suff+'FitAC'+'y'+str(j)+c] = ROOT.TPad('ph'+suff+"_c_projY{}_{}".format(i,c), 'ph'+suff+"_c_projY{}_{}".format(i,c),0,0.3,1,1)
                self.canvas['pr'+suff+'FitAC'+'y'+str(j)+c] = ROOT.TPad('pr'+suff+"_c_projY{}_{}".format(i,c),  'pr'+suff+"_c_projY{}_{}".format(i,c), 0,0,1,0.3)
                self.canvas['ph'+suff+'FitAC'+'y'+str(j)+c].SetBottomMargin(0.02)
                self.canvas['ph'+suff+'FitAC'+'y'+str(j)+c].Draw()
                self.canvas['pr'+suff+'FitAC'+'y'+str(j)+c].SetTopMargin(0.01)
                self.canvas['pr'+suff+'FitAC'+'y'+str(j)+c].SetBottomMargin(0.4)
                self.canvas['pr'+suff+'FitAC'+'y'+str(j)+c].Draw()
                
                self.canvas['ph'+suff+'FitAC'+'y'+str(j)+c].cd()
                self.canvas['ph'+suff+'FitAC'+'y'+str(j)+c].SetGridx()
                self.canvas['ph'+suff+'FitAC'+'y'+str(j)+c].SetGridy()
                self.histos[suff+'FitAC'+'y'+str(j)+c] = self.histos[suff+'FitAC'+c].ProjectionX("projY{}_{}".format(j,c),j,j)
                self.histos[suff+'FitBand'+'y'+str(j)+c] = self.histos[suff+'FitBand'+c].ProjectionX("projY{}_{}_err".format(j,c),j,j)
                # self.histos[suff+'FitBandPDF'+'y'+str(j)+c] = self.histos[suff+'FitBandPDF'+c].ProjectionX("projY{}_{}_errPDF".format(j,c),j,j)
                # self.histos[suff+'FitBandScale'+'y'+str(j)+c] = self.histos[suff+'FitBandScale'+c].ProjectionX("projY{}_{}_errScale".format(j,c),j,j)
                self.histos[suff+'FitAC'+'y'+str(j)+c].SetTitle(self.coeffDict[c][2]+", rapidity distribution, "+str(self.qtArr[j-1])+"<q_{T}^{W}<"+str(self.qtArr[j])+", "+self.signDict[self.sign])
                self.histos[suff+'FitAC'+'y'+str(j)+c].GetXaxis().SetTitle('|Y_{W}|')
                self.histos[suff+'FitAC'+'y'+str(j)+c].SetLineWidth(2)
                self.histos[suff+'FitAC'+'y'+str(j)+c].SetLineColor(ROOT.kBlack)
                self.histos[suff+'FitAC'+'y'+str(j)+c].SetStats(0)
                self.histos[suff+'FitAC'+'y'+str(j)+c].SetMarkerStyle(20)
                self.histos[suff+'FitAC'+'y'+str(j)+c].SetMarkerSize(0.5)
                self.histos[suff+'FitAC'+'y'+str(j)+c].Draw()
                self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetTitleSize(0.06)
                self.histos[suff+'FitAC'+'y'+str(j)+c].SetLabelSize(0.05,'y')
                if not 'unpol' in c and self.anaKind['name']=='angCoeff':
                    # self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetRangeUser(-4,4)
                    self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetTitle(self.coeffDict[c][2])
                elif self.anaKind['name']=='helXsec':
                     self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetTitle('d#sigma_'+c+'/dq_{T}^{W}d|Y_{W}| [fb/GeV]')
                elif self.anaKind['name']=='normXsec':
                    self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetTitle('d#sigma_'+c+'/dq_{T}^{W}d|Y_{W}|/#sigma_{tot} [1/GeV]')
                elif self.anaKind['name']=='mu':
                    self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetTitle('#mu_'+c)
                else :
                    self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetTitle('d#sigma^{U+L}/dq_{T}^{W}d|Y_{W}| [fb/GeV]')
                    # maxvalMain = self.histos[suff+'FitAC'+'y'+str(j)+c].GetMaximum()
                    # self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetRangeUser(0,maxvalMain*1.5)
                self.histos[suff+'FitBand'+'y'+str(j)+c].SetFillColor(ROOT.kOrange)#kMagenta-7)
                self.histos[suff+'FitBand'+'y'+str(j)+c].SetFillStyle(0)
                self.histos[suff+'FitBand'+'y'+str(j)+c].SetLineColor(ROOT.kOrange)#kMagenta-7)
                self.histos[suff+'FitBand'+'y'+str(j)+c].SetLineWidth(2)
                self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetTitleOffset(0.8)
                self.histos[suff+'FitAC'+'y'+str(j)+c].GetXaxis().SetTitleOffset(3)
                self.histos[suff+'FitAC'+'y'+str(j)+c].GetXaxis().SetLabelOffset(3)
                self.histos[suff+'FitBand'+'y'+str(j)+c].DrawCopy("E2 same")
                self.histos[suff+'FitBand'+'y'+str(j)+c].SetFillStyle(3001)
                self.histos[suff+'FitBand'+'y'+str(j)+c].Draw("E2 same")
                self.histos[suff+'FitAC'+'y'+str(j)+c].DrawCopy("E1 same") #to have foreground   
                
                if aposteriori!='' :
                    self.histos[suff+'apo'+'y'+str(j)+c] = self.histos[suff+'apo'+c].ProjectionX('apo'+"projY{}_{}".format(j,c),j,j)
                    self.histos[suff+'apo'+'y'+str(j)+c].SetLineColor(ROOT.kAzure+1)
                    self.histos[suff+'apo'+'y'+str(j)+c].SetMarkerColor(ROOT.kBlue-4)
                    self.histos[suff+'apo'+'y'+str(j)+c].SetLineWidth(5)
                    self.histos[suff+'apo'+'y'+str(j)+c].SetMarkerStyle(20)
                    self.histos[suff+'apo'+'y'+str(j)+c].Draw("EX0 same")
                    self.histos[suff+'apo'+'y'+str(j)+c].SetMarkerSize(1)
                    
                ROOT.gPad.Update()
                if 'unpol' in c and self.anaKind['angNames']: minvalMain = 0
                else : minvalMain = self.canvas['ph'+suff+'FitAC'+'y'+str(j)+c].GetFrame().GetY1()
                self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetRangeUser(minvalMain,1.3*self.canvas['ph'+suff+'FitAC'+'y'+str(j)+c].GetFrame().GetY2())
                self.leg[suff+'FitAC'+'y'+str(j)+c].Draw("same")
                
                self.canvas['pr'+suff+'FitAC'+'y'+str(j)+c].cd()
                self.canvas['pr'+suff+'FitAC'+'y'+str(j)+c].SetGridx()
                self.canvas['pr'+suff+'FitAC'+'y'+str(j)+c].SetGridy()
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c] = self.histos[suff+'FitRatioAC'+c].ProjectionX("projY{}_{}_RatioAC".format(j,c),j,j)
                self.histos[suff+'FitRatio'+'y'+str(j)+c] = self.histos[suff+'FitRatio'+c].ProjectionX("projY{}_{}_Ratio".format(j,c),j,j)
                for d in self.dirList :
                    self.histos[suff+'FitRatioPDF'+'y'+str(j)+d+c] = self.histos[suff+'FitRatioPDF'+d+c].ProjectionX("projY{}_{}_RatioPDF".format(j,c)+d,j,j)
                    self.histos[suff+'FitRatioScale'+'y'+str(j)+d+c] = self.histos[suff+'FitRatioScale'+d+c].ProjectionX("projY{}_{}_RatioScale".format(j,c)+d,j,j)
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].GetXaxis().SetTitle('|Y_{W}|')
                if ('unpol' in c or (not self.anaKind['differential'])):
                    self.histos[suff+'FitRatioAC'+'y'+str(j)+c].GetYaxis().SetTitle('Pred./Theory')
                else :
                    self.histos[suff+'FitRatioAC'+'y'+str(j)+c].GetYaxis().SetTitle('Pred.-Theory')
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].SetMarkerStyle(20)
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].SetMarkerSize(0.5)
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].SetLineWidth(2)
                self.histos[suff+'FitRatio'+'y'+str(j)+c].SetLineWidth(2)
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].SetTitle("")
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].GetYaxis().SetTitleOffset(0.32)
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].SetTitleSize(0.13,'y')
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].SetLabelSize(0.10,'y')
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].GetXaxis().SetTitleOffset(1.2)
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].SetTitleSize(0.12,'x')
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].SetLabelSize(0.12,'x')
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].SetStats(0)
                for d in self.dirList :
                    self.histos[suff+'FitRatioPDF'+'y'+str(j)+d+c].SetLineWidth(2)
                    self.histos[suff+'FitRatioScale'+'y'+str(j)+d+c].SetLineWidth(2)
                self.histos[suff+'FitRatio'+'y'+str(j)+c].SetLineColor(ROOT.kOrange)
                for d in self.dirList :
                    self.histos[suff+'FitRatioPDF'+'y'+str(j)+d+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                    self.histos[suff+'FitRatioScale'+'y'+str(j)+d+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
                self.histos[suff+'FitRatio'+'y'+str(j)+c].SetFillColor(ROOT.kOrange)
                self.histos[suff+'FitRatio'+'y'+str(j)+c].SetFillStyle(3001)
                
                maxvalRatio=0
                for xx in range(1,self.histos[suff+'FitRatioAC'+'y'+str(j)+c].GetNbinsX()+1) :
                    if self.histos[suff+'FitRatioAC'+'y'+str(j)+c].GetBinError(xx)>maxvalRatio :
                        maxvalRatio= self.histos[suff+'FitRatioAC'+'y'+str(j)+c].GetBinError(xx)
                    # self.histos[suff+'FitRatioAC'+'y'+str(j)+c].SetBinError(xx,0)
                if self.anaKind['angNames'] :
                    PadYCut =  self.RatioPadYcutDict[c][0]
                else : PadYCut = self.RatioPadYcut 
                if maxvalRatio>PadYCut :
                    maxvalRatio=PadYCut
                # self.histos[suff+'FitRatioAC'+'y'+str(j)+c].GetYaxis().SetRangeUser(-maxvalRatio*1.1+1,maxvalRatio*1.1+1)
                
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].Draw("EX0")
                self.histos[suff+'FitRatio'+'y'+str(j)+c].DrawCopy("same E2")
                for d in self.dirList :
                    self.histos[suff+'FitRatioPDF'+'y'+str(j)+d+c].Draw("hist same")
                    self.histos[suff+'FitRatioScale'+'y'+str(j)+d+c].Draw("hist same")
                self.histos[suff+'FitRatioAC'+'y'+str(j)+c].Draw("EX0 same") #to have foreground
                
                if aposteriori!='' :
                    self.histos[suff+'apoRatio'+'y'+str(j)+c] = self.histos[suff+'apo'+'Ratio'+c].ProjectionX("projY{}_{}_apoRatio".format(j,c),j,j)
                    self.histos[suff+'apoRatio'+'y'+str(j)+c].SetLineColor(ROOT.kAzure+1)
                    self.histos[suff+'apoRatio'+'y'+str(j)+c].SetMarkerColor(ROOT.kBlue-4)
                    self.histos[suff+'apoRatio'+'y'+str(j)+c].SetLineWidth(4)
                    self.histos[suff+'apoRatio'+'y'+str(j)+c].SetMarkerStyle(20)
                    self.histos[suff+'apoRatio'+'y'+str(j)+c].Draw("EX0 same")
                    self.histos[suff+'apoRatio'+'y'+str(j)+c].SetMarkerSize(1)
                
                self.leg[suff+'FitAC'+'y'+str(j)+c].AddEntry(self.histos[suff+'FitAC'+'y'+str(j)+c], self.signDict[self.sign]+' Fit')
                self.leg[suff+'FitAC'+'y'+str(j)+c].AddEntry(self.histos[suff+'FitBand'+'y'+str(j)+c], "MC Syst. Unc.")
                self.leg[suff+'FitAC'+'y'+str(j)+c].AddEntry(self.histos[suff+'FitRatioPDF'+'y'+str(j)+'up'+c], self.groupedSystColors['LHEPdfWeightVars'][1])
                self.leg[suff+'FitAC'+'y'+str(j)+c].AddEntry(self.histos[suff+'FitRatioScale'+'y'+str(j)+'up'+c], self.groupedSystColors['LHEScaleWeightVars'][1])
                self.leg[suff+'FitAC'+'y'+str(j)+c].SetNColumns(2)
                if aposteriori!='' :
                    self.leg[suff+'FitAC'+'y'+str(j)+c].AddEntry(self.histos[suff+'apo'+'y'+str(j)+c], "Post-fit regularized")
                
                
                #Additional block for scale-pdf splitted
                # self.histos[suff+'FitBandPDF'+'y'+str(j)+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
                # self.histos[suff+'FitBandPDF'+'y'+str(j)+c].SetFillStyle(0)
                # self.histos[suff+'FitBandPDF'+'y'+str(j)+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
                # self.histos[suff+'FitBandPDF'+'y'+str(j)+c].SetLineWidth(2)
                # self.histos[suff+'FitBandPDF'+'y'+str(j)+c].DrawCopy("E2 same")
                # self.histos[suff+'FitBandPDF'+'y'+str(j)+c].SetFillStyle(3001)
                # self.histos[suff+'FitBandPDF'+'y'+str(j)+c].Draw("E2 same")
                # self.histos[suff+'FitBandScale'+'y'+str(j)+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
                # self.histos[suff+'FitBandScale'+'y'+str(j)+c].SetFillStyle(0)
                # self.histos[suff+'FitBandScale'+'y'+str(j)+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
                # self.histos[suff+'FitBandScale'+'y'+str(j)+c].SetLineWidth(2)
                # self.histos[suff+'FitBandScale'+'y'+str(j)+c].DrawCopy("E2 same")
                # self.histos[suff+'FitBandScale'+'y'+str(j)+c].SetFillStyle(3001)
                # self.histos[suff+'FitBandScale'+'y'+str(j)+c].Draw("E2 same")
                
                     
        
          

            
            #-------------- unrolled: qt large, y small ----------------- #
            self.histos[suff+'FitAC'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_UNRqty'+c,suff+'_coeff_UNRqty'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            self.histos[suff+'FitBand'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_UNRqty_syst'+c,suff+'_coeff_UNRqty_syst'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            # self.histos[suff+'FitBandScale'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_UNRqty_systPDF'+c,suff+'_coeff_UNRqty_systPDF'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            # self.histos[suff+'FitBandPDF'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_UNRqty_systScale'+c,suff+'_coeff_UNRqty_systScale'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            self.histos[suff+'FitRatio'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_UNRqty_Ratio'+c,suff+'_coeff_UNRqty_Ratio'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            self.histos[suff+'FitRatioAC'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_UNRqty_RatioAC'+c,suff+'_coeff_UNRqty_RatioAC'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'UNRqty'+d+c] = ROOT.TH1F(suff+'_coeff_UNRqty_RatioPDF'+d+c,suff+'_coeff_UNRqty_RatioPDF'+d+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
                self.histos[suff+'FitRatioScale'+'UNRqty'+d+c] = ROOT.TH1F(suff+'_coeff_UNRqty_RatioScale'+d+c,suff+'_coeff_UNRqty_RatioScale'+d+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            if aposteriori!='' :
                self.histos[suff+'apo'+'UNRqty'+c] = ROOT.TH1F(suff+'_apo_UNRqty'+c,suff+'_apo_UNRqty'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
                self.histos[suff+'apoRatio'+'UNRqty'+c] = ROOT.TH1F(suff+'_apoRatio_UNRqty'+c,suff+'_apoRatio_UNRqty'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            if toy !='' :
                self.histos[suff+'FitAC'+'UNRqty'+c+'toy'] = ROOT.TH1F(suff+'FitAC'+'UNRqty'+c+'toy',suff+'FitAC'+'UNRqty'+c+'toy',len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
                
            
            for q in self.qtArr[:-1] :
                for y in self.yArr[:-1] :
                    indexUNRqty = self.qtArr.index(float(q))*(len(self.yArr)-1)+self.yArr.index(float(y))
                    self.histos[suff+'FitAC'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitAC'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitAC'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitAC'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBand'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitBand'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBand'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitBand'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))   
                    # self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitBandPDF'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    # self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitBandPDF'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    # self.histos[suff+'FitBandScale'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitBandScale'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    # self.histos[suff+'FitBandScale'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitBandScale'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitRatio'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitRatio'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitRatio'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitRatio'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitRatioAC'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitRatioAC'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))

                    for d in self.dirList :
                        self.histos[suff+'FitRatioPDF'+'UNRqty'+d+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitRatioPDF'+d+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'FitRatioPDF'+'UNRqty'+d+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitRatioPDF'+d+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'FitRatioScale'+'UNRqty'+d+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitRatioScale'+d+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'FitRatioScale'+'UNRqty'+d+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitRatioScale'+d+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    
                    if aposteriori!='' :
                        self.histos[suff+'apo'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'apo'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'apo'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'apo'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'apoRatio'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'apoRatio'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'apoRatio'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'apoRatio'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    
                    if toy !='' :
                        self.histos[suff+'FitAC'+'UNRqty'+c+'toy'].SetBinContent(indexUNRqty+1,self.histos[suff+'FitAC'+c+'toy'].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'FitAC'+'UNRqty'+c+'toy'].SetBinError(indexUNRqty+1,self.histos[suff+'FitAC'+c+'toy'].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))

                    if self.yArr.index(y)==0 :
                        # print("unrbin=",indexUNRqty+1, "binq=",q)
                        self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetBinLabel(indexUNRqty+1,"q_{T}^{W}#in[%.0f,%.0f]" % (q, self.qtArr[self.qtArr.index(q)+1]))
                        self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().ChangeLabel(indexUNRqty+1,340,0.03)
                        
                        self.histos[suff+'FitRatioAC'+'UNRqty'+c].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+'FitRatioAC'+'UNRqty'+c].GetXaxis().SetBinLabel(indexUNRqty+1,"q_{T}^{W}#in[%.0f,%.0f]" % (q, self.qtArr[self.qtArr.index(q)+1]))
                        self.histos[suff+'FitRatioAC'+'UNRqty'+c].GetXaxis().ChangeLabel(indexUNRqty+1,340,0.08)
    
            # self.histos[suff+'FitAC'+'UNRqty'+c].SetTitle(self.coeffDict[c][2]+", unrolled q_{T}(Y), "+self.signDict[self.sign])
            self.histos[suff+'FitAC'+'UNRqty'+c].SetTitle("")
            # self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetTitle('fine binning: |Y_{W}| 0#rightarrow 2.4')
            # self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetTitleOffset(1.45)
            self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetTitleSize(0.06)
            if not 'unpol' in c and self.anaKind['name']=='angCoeff':
                # self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetRangeUser(-4,4)
                self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetTitle(self.coeffDict[c][2])
            elif self.anaKind['name']=='helXsec':
                self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetTitle('d#sigma_'+c+'/dq_{T}^{W}d|Y_{W}| [fb/GeV]')
            elif self.anaKind['name']=='normXsec':
                self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetTitle('d#sigma_'+c+'/dq_{T}^{W}d|Y_{W}|/#sigma_{tot} [1/GeV]')
            elif self.anaKind['name']=='mu':
                self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetTitle('#mu_'+c)
            else :
                self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetTitle('d#sigma^{U+L}/dq_{T}^{W}d|Y_{W}| [fb/GeV]')
                # maxvalMain = self.histos[suff+'FitAC'+'UNRqty'+c].GetMaximum()
                # self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetRangeUser(0,maxvalMain*1.5)
            self.histos[suff+'FitAC'+'UNRqty'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'UNRqty'+c].SetStats(0)
            self.histos[suff+'FitBand'+'UNRqty'+c].SetFillColor(ROOT.kOrange)#kMagenta-7)
            self.histos[suff+'FitBand'+'UNRqty'+c].SetFillStyle(0)
            self.histos[suff+'FitBand'+'UNRqty'+c].SetLineColor(ROOT.kOrange)#kMagenta-7)
            self.histos[suff+'FitBand'+'UNRqty'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'UNRqty'+c].SetLineColor(ROOT.kBlack)
            self.histos[suff+'FitAC'+'UNRqty'+c].SetMarkerStyle(20)
            self.histos[suff+'FitAC'+'UNRqty'+c].SetMarkerSize(0.5)
            self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetTitleOffset(0.8)
            self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetTitleOffset(3)
            self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetLabelOffset(3)
            self.histos[suff+'FitAC'+'UNRqty'+c].SetLabelSize(0.05,'y')
            
            if ('unpol' in c or (not self.anaKind['differential'])) : 
                self.histos[suff+'FitRatioAC'+'UNRqty'+c].GetYaxis().SetTitle('Pred./Theory')
            else :
                self.histos[suff+'FitRatioAC'+'UNRqty'+c].GetYaxis().SetTitle('Pred.-Theory')
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetMarkerStyle(20)
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetMarkerSize(0.5)
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetLineWidth(2)
            self.histos[suff+'FitRatio'+'UNRqty'+c].SetLineWidth(2)
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetTitle("")
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].GetYaxis().SetTitleOffset(0.32)
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetTitleSize(0.13,'y')
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetLabelSize(0.10,'y')
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].GetXaxis().SetTitleOffset(1.65)
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetTitleSize(0.12,'x')
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetLabelSize(0.15,'x')
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetLabelOffset(0.03,'x')
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].LabelsOption("d")
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetStats(0)
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].GetXaxis().SetTitle('fine binning: |Y_{W}| 0#rightarrow 2.4')
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'UNRqty'+d+c].SetLineWidth(2)
                self.histos[suff+'FitRatioScale'+'UNRqty'+d+c].SetLineWidth(2)
            self.histos[suff+'FitRatio'+'UNRqty'+c].SetLineColor(ROOT.kOrange)
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'UNRqty'+d+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                self.histos[suff+'FitRatioScale'+'UNRqty'+d+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            self.histos[suff+'FitRatio'+'UNRqty'+c].SetFillColor(ROOT.kOrange)
            self.histos[suff+'FitRatio'+'UNRqty'+c].SetFillStyle(3001)
            
            maxvalRatio=0
            for xx in range(1,len(self.unrolledQtY)) :
                if not 'unpol' in c and not 'UL' in c: ACstring = '' #for unpol integrated in qt the uncertainity of MC is higher than prediction--> used the former as limit
                else : ACstring ='AC'
                if self.histos[suff+'FitRatio'+ACstring+'UNRqty'+c].GetBinError(xx)>maxvalRatio :
                    maxvalRatio= self.histos[suff+'FitRatio'+ACstring+'UNRqty'+c].GetBinError(xx)
                # if not 'unpol' in c and not 'UL' in c : self.histos[suff+'FitRatioAC'+'UNRqty'+c].SetBinError(xx,0)
            if self.anaKind['angNames'] :
                    PadYCut =  self.RatioPadYcutDict[c][0]
            else : PadYCut = self.RatioPadYcut 
            if maxvalRatio>PadYCut :
                    maxvalRatio=PadYCut
            # self.histos[suff+'FitRatioAC'+'UNRqty'+c].GetYaxis().SetRangeUser(-maxvalRatio*1.1+1,maxvalRatio*1.1+1)

            # self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
            # self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetFillStyle(0)
            # self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
            # self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetLineWidth(2)    
            # self.histos[suff+'FitBandScale'+'UNRqty'+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
            # self.histos[suff+'FitBandScale'+'UNRqty'+c].SetFillStyle(0)
            # self.histos[suff+'FitBandScale'+'UNRqty'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
            # self.histos[suff+'FitBandScale'+'UNRqty'+c].SetLineWidth(2)
            
            self.canvas[suff+'FitAC'+'UNRqty'+c] = ROOT.TCanvas(suff+'_c_UNRqty'+c,suff+'_c_UNRqty'+c,1200,1050)
            self.leg[suff+'FitAC'+'UNRqty'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            self.canvas[suff+'FitAC'+'UNRqty'+c].cd()
            self.canvas['ph'+suff+'FitAC'+'UNRqty'+c] = ROOT.TPad('ph'+suff+'_c_UNRqty'+c, 'ph'+suff+'_c_UNRqty'+c,0,0.3,1,1)
            self.canvas['pr'+suff+'FitAC'+'UNRqty'+c] = ROOT.TPad('pr'+suff+'_c_UNRqty'+c,  'pr'+suff+'_c_UNRqty'+c, 0,0,1,0.3)
            self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].SetBottomMargin(0.02)
            self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].Draw()
            self.canvas['pr'+suff+'FitAC'+'UNRqty'+c].SetTopMargin(0.01)
            self.canvas['pr'+suff+'FitAC'+'UNRqty'+c].SetBottomMargin(0.4)
            self.canvas['pr'+suff+'FitAC'+'UNRqty'+c].Draw()
            
            self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].cd()
            self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].SetGridx()
            self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].SetGridy()
            self.histos[suff+'FitAC'+'UNRqty'+c].Draw()
            self.histos[suff+'FitBand'+'UNRqty'+c].DrawCopy('E2 same')
            self.histos[suff+'FitBand'+'UNRqty'+c].SetFillStyle(3001)
            self.histos[suff+'FitBand'+'UNRqty'+c].Draw("E2 same")
            self.histos[suff+'FitAC'+'UNRqty'+c].DrawCopy("E1 same") #to have foreground 
            if aposteriori!='' :
                self.histos[suff+'apo'+'UNRqty'+c].SetLineColor(ROOT.kAzure+1)
                self.histos[suff+'apo'+'UNRqty'+c].SetMarkerColor(ROOT.kBlue-4)
                self.histos[suff+'apo'+'UNRqty'+c].SetLineWidth(5)
                self.histos[suff+'apo'+'UNRqty'+c].SetMarkerStyle(20)
                self.histos[suff+'apo'+'UNRqty'+c].Draw("EX0 same")
            
            if toy!='' :
                self.histos[suff+'FitAC'+'UNRqty'+c+'toy'].SetLineColor(ROOT.kViolet+1)
                self.histos[suff+'FitAC'+'UNRqty'+c+'toy'].SetFillStyle(0)
                self.histos[suff+'FitAC'+'UNRqty'+c+'toy'].SetLineWidth(3)
                self.histos[suff+'FitAC'+'UNRqty'+c+'toy'].SetMarkerStyle(2)
                self.histos[suff+'FitAC'+'UNRqty'+c+'toy'].SetMarkerColor(ROOT.kViolet)
                self.histos[suff+'FitAC'+'UNRqty'+c+'toy'].Draw("E2 same")
            
            ROOT.gPad.Update()
            if 'unpol' in c and self.anaKind['angNames']: minvalMain = 0
            else : minvalMain = self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetFrame().GetY1()
            self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetRangeUser(minvalMain,1.3*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetFrame().GetY2())
            self.leg[suff+'FitAC'+'UNRqty'+c].Draw("same")                   

            self.canvas['pr'+suff+'FitAC'+'UNRqty'+c].cd()
            self.canvas['pr'+suff+'FitAC'+'UNRqty'+c].SetGridx()
            self.canvas['pr'+suff+'FitAC'+'UNRqty'+c].SetGridy()
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].Draw("EX0")
            self.histos[suff+'FitRatio'+'UNRqty'+c].DrawCopy("same E2")
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'UNRqty'+d+c].Draw("hist same")
                self.histos[suff+'FitRatioScale'+'UNRqty'+d+c].Draw("hist same")
            self.histos[suff+'FitRatioAC'+'UNRqty'+c].Draw("EX0 same") #to have foreground
            
            if aposteriori!='' :
                    self.histos[suff+'apoRatio'+'UNRqty'+c].SetLineColor(ROOT.kAzure+1)
                    self.histos[suff+'apoRatio'+'UNRqty'+c].SetMarkerColor(ROOT.kBlue-4)
                    self.histos[suff+'apoRatio'+'UNRqty'+c].SetLineWidth(4)
                    self.histos[suff+'apoRatio'+'UNRqty'+c].SetMarkerStyle(20)
                    self.histos[suff+'apoRatio'+'UNRqty'+c].SetMarkerSize(1)
                    self.histos[suff+'apoRatio'+'UNRqty'+c].Draw("EX0 same")
            
            self.leg[suff+'FitAC'+'UNRqty'+c].AddEntry(self.histos[suff+'FitAC'+'UNRqty'+c], self.signDict[self.sign]+' Fit')
            self.leg[suff+'FitAC'+'UNRqty'+c].AddEntry(self.histos[suff+'FitBand'+'UNRqty'+c], "MC Syst. Unc.")
            self.leg[suff+'FitAC'+'UNRqty'+c].AddEntry(self.histos[suff+'FitRatioPDF'+'UNRqty'+'up'+c], self.groupedSystColors['LHEPdfWeightVars'][1])
            self.leg[suff+'FitAC'+'UNRqty'+c].AddEntry(self.histos[suff+'FitRatioScale'+'UNRqty'+'up'+c], self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitAC'+'UNRqty'+c].SetNColumns(2)
            if aposteriori!='' :
                    self.leg[suff+'FitAC'+'UNRqty'+c].AddEntry(self.histos[suff+'apo'+'UNRqty'+c], "Post-fit regularized")
            if toy!='' :
                    self.leg[suff+'FitAC'+'UNRqty'+c].AddEntry(self.histos[suff+'FitAC'+'UNRqty'+c+'toy'], "MC toys  #mu#pm#sigma")
            
            self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].cd()
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetRightMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetLeftMargin()+0.07,1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetTopMargin(),cmslab)

                    
            # self.histos[suff+'FitBandPDF'+'UNRqty'+c].DrawCopy('E2 same')
            # self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetFillStyle(3004)
            # self.histos[suff+'FitBandPDF'+'UNRqty'+c].Draw("E2 same")
            # self.histos[suff+'FitBandScale'+'UNRqty'+c].DrawCopy('E2 same')
            # self.histos[suff+'FitBandScale'+'UNRqty'+c].SetFillStyle(3005)
            # self.histos[suff+'FitBandScale'+'UNRqty'+c].Draw("E2 same")
            



            #------------- unrolled: y large, qt small -----------------                    
            self.histos[suff+'FitAC'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_UNRyqt'+c,suff+'_coeff_UNRyqt'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            self.histos[suff+'FitBand'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_UNRyqt_syst'+c,suff+'_coeff_UNRyqt_syst'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            # self.histos[suff+'FitBandPDF'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_UNRyqt_systPDF'+c,suff+'_coeff_UNRyqt_systPDF'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            # self.histos[suff+'FitBandScale'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_UNRyqt_systScale'+c,suff+'_coeff_UNRyqt_systScale'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            self.histos[suff+'FitRatio'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_UNRyqt_Ratio'+c,suff+'_coeff_UNRyqt_Ratio'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_UNRyqt_RatioAC'+c,suff+'_coeff_UNRyqt_RatioAC'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'UNRyqt'+d+c] = ROOT.TH1F(suff+'_coeff_UNRyqt_RatioPDF'+d+c,suff+'_coeff_UNRyqt_RatioPDF'+d+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
                self.histos[suff+'FitRatioScale'+'UNRyqt'+d+c] = ROOT.TH1F(suff+'_coeff_UNRyqt_RatioScale'+d+c,suff+'_coeff_UNRyqt_RatioScale'+d+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            if aposteriori!='' :
                self.histos[suff+'apo'+'UNRyqt'+c] = ROOT.TH1F(suff+'_apo_UNRyqt'+c,suff+'_apo_UNRyqt'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
                self.histos[suff+'apoRatio'+'UNRyqt'+c] = ROOT.TH1F(suff+'_apoRatio_UNRyqt'+c,suff+'_apoRatio_UNRyqt'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            if toy !='' :
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toy'] = ROOT.TH1F(suff+'FitAC'+'UNRyqt'+c+'toy',suff+'FitAC'+'UNRyqt'+c+'toy',len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            
            for y in self.yArr[:-1] :
                for q in self.qtArr[:-1] :
                    indexUNRyqt = self.yArr.index(float(y))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                    self.histos[suff+'FitAC'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitAC'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitAC'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitAC'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBand'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitBand'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBand'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitBand'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    # self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitBandPDF'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    # self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitBandPDF'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    # self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitBandScale'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    # self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitBandScale'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))                    
                    self.histos[suff+'FitRatio'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitRatio'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitRatio'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitRatio'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitRatioAC'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitRatioAC'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))

                    for d in self.dirList :
                        self.histos[suff+'FitRatioPDF'+'UNRyqt'+d+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitRatioPDF'+d+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'FitRatioPDF'+'UNRyqt'+d+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitRatioPDF'+d+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'FitRatioScale'+'UNRyqt'+d+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitRatioScale'+d+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'FitRatioScale'+'UNRyqt'+d+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitRatioScale'+d+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    
                    if aposteriori!='' :
                        self.histos[suff+'apo'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'apo'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'apo'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'apo'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'apoRatio'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'apoRatio'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'apoRatio'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'apoRatio'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        
                    if toy !='' :
                        self.histos[suff+'FitAC'+'UNRyqt'+c+'toy'].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitAC'+c+'toy'].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'FitAC'+'UNRyqt'+c+'toy'].SetBinError(indexUNRyqt+1,self.histos[suff+'FitAC'+c+'toy'].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))

                    if self.qtArr.index(q)==0 :
                        self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetBinLabel(indexUNRyqt+1,"|Y_{W}|#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
                        self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().ChangeLabel(indexUNRyqt+1,340,0.03)
                        
                        self.histos[suff+'FitRatioAC'+'UNRyqt'+c].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+'FitRatioAC'+'UNRyqt'+c].GetXaxis().SetBinLabel(indexUNRyqt+1,"|Y_{W}|#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
                        self.histos[suff+'FitRatioAC'+'UNRyqt'+c].GetXaxis().ChangeLabel(indexUNRyqt+1,340,0.08)
                        
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetTitle(self.coeffDict[c][2]+", unrolled Y(q_{T}), "+self.signDict[self.sign])
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetTitle("")
            # self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 60 GeV')
            # self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetTitleOffset(1.45) 
            self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetTitleSize(0.06)
            if not 'unpol' in c and self.anaKind['name']=='angCoeff':
                # self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetRangeUser(-4,4)
                self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetTitle(self.coeffDict[c][2])
            elif self.anaKind['name']=='helXsec':
                self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetTitle('d#sigma_'+c+'/dq_{T}^{W}d|Y_{W}| [fb/GeV]')
            elif self.anaKind['name']=='normXsec':
                self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetTitle('d#sigma_'+c+'/dq_{T}^{W}d|Y_{W}|/#sigma_{tot} [1/GeV]')
            elif self.anaKind['name']=='mu':
                self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetTitle('#mu_'+c)
            else :
                self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetTitle('d#sigma^{U+L}/dq_{T}^{W}d|Y_{W}| [fb/GeV]')
                # maxvalMain = self.histos[suff+'FitAC'+'UNRyqt'+c].GetMaximum()
                # self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetRangeUser(0,maxvalMain*1.5)
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetStats(0)
            self.histos[suff+'FitBand'+'UNRyqt'+c].SetFillColor(ROOT.kOrange)
            self.histos[suff+'FitBand'+'UNRyqt'+c].SetFillStyle(0)
            self.histos[suff+'FitBand'+'UNRyqt'+c].SetLineColor(ROOT.kOrange)
            self.histos[suff+'FitBand'+'UNRyqt'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetLineColor(ROOT.kBlack)
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetMarkerStyle(20)
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetMarkerSize(0.5)
            self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetTitleOffset(0.8)
            self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetTitleOffset(3)
            self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetLabelOffset(3)
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetLabelSize(0.05,'y')
            
            if ('unpol' in c or (not self.anaKind['differential'])) : 
                self.histos[suff+'FitRatioAC'+'UNRyqt'+c].GetYaxis().SetTitle('Pred./Theory')
            else :
                self.histos[suff+'FitRatioAC'+'UNRyqt'+c].GetYaxis().SetTitle('Pred.-Theory')
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetMarkerStyle(20)
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetMarkerSize(0.5)
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetLineWidth(2)
            self.histos[suff+'FitRatio'+'UNRyqt'+c].SetLineWidth(2)
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetTitle("")
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].GetYaxis().SetTitleOffset(0.32)
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetTitleSize(0.13,'y')
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetLabelSize(0.10,'y')
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].GetXaxis().SetTitleOffset(1.65)
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetTitleSize(0.12,'x')
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetLabelSize(0.15,'x')
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetLabelOffset(0.03,'x')
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].LabelsOption("d")
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetStats(0)
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 60 GeV')
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'UNRyqt'+d+c].SetLineWidth(2)
                self.histos[suff+'FitRatioScale'+'UNRyqt'+d+c].SetLineWidth(2)
            self.histos[suff+'FitRatio'+'UNRyqt'+c].SetLineColor(ROOT.kOrange)
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'UNRyqt'+d+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                self.histos[suff+'FitRatioScale'+'UNRyqt'+d+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            self.histos[suff+'FitRatio'+'UNRyqt'+c].SetFillColor(ROOT.kOrange)
            self.histos[suff+'FitRatio'+'UNRyqt'+c].SetFillStyle(3001)
            
            maxvalRatio=0
            for xx in range(1,len(self.unrolledYQt)) :
                if not 'unpol' in c and not 'UL' in c: ACstring = '' #for unpol integrated in qt the uncertainity of MC is higher than prediction--> used the former as limit
                else : ACstring ='AC'
                if self.histos[suff+'FitRatio'+ACstring+'UNRyqt'+c].GetBinError(xx)>maxvalRatio :
                    maxvalRatio= self.histos[suff+'FitRatio'+ACstring+'UNRyqt'+c].GetBinError(xx)
                # if not 'unpol' in c and not 'UL' in c :self.histos[suff+'FitRatioAC'+'UNRyqt'+c].SetBinError(xx,0)
            if self.anaKind['angNames'] :
                    PadYCut =  self.RatioPadYcutDict[c][0]
            else : PadYCut = self.RatioPadYcut 
            if maxvalRatio>PadYCut :
                    maxvalRatio=PadYCut
            if not self.asimov and c=='A4' : self.histos[suff+'FitRatioAC'+'UNRyqt'+c].GetYaxis().SetRangeUser(-1.5,1.5)
            
            # self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            # self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetFillStyle(0)
            # self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            # self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetLineWidth(2)
            # self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            # self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetFillStyle(0)
            # self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            # self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetLineWidth(2)
            
            self.canvas[suff+'FitAC'+'UNRyqt'+c] = ROOT.TCanvas(suff+'_c_UNRyqt'+c,suff+'_c_UNRyqt'+c,1200,1050)
            self.leg[suff+'FitAC'+'UNRyqt'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            self.canvas[suff+'FitAC'+'UNRyqt'+c].cd()
            self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c] = ROOT.TPad('ph'+suff+'_c_UNRyqt'+c, 'ph'+suff+'_c_UNRyqt'+c,0,0.3,1,1)
            self.canvas['pr'+suff+'FitAC'+'UNRyqt'+c] = ROOT.TPad('pr'+suff+'_c_UNRyqt'+c,  'pr'+suff+'_c_UNRyqt'+c, 0,0,1,0.3)
            self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].SetBottomMargin(0.02)
            self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].Draw()
            self.canvas['pr'+suff+'FitAC'+'UNRyqt'+c].SetTopMargin(0.01)
            self.canvas['pr'+suff+'FitAC'+'UNRyqt'+c].SetBottomMargin(0.4)
            self.canvas['pr'+suff+'FitAC'+'UNRyqt'+c].Draw()
            
            self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].cd()
            self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].SetGridx()
            self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].SetGridy()
            self.histos[suff+'FitAC'+'UNRyqt'+c].Draw()
            self.histos[suff+'FitBand'+'UNRyqt'+c].DrawCopy('E2 same')
            self.histos[suff+'FitBand'+'UNRyqt'+c].SetFillStyle(3001)
            self.histos[suff+'FitBand'+'UNRyqt'+c].Draw("E2 same")
            self.histos[suff+'FitAC'+'UNRyqt'+c].DrawCopy("E1 same") #to have foreground 
            
            if aposteriori!='' :
                self.histos[suff+'apo'+'UNRyqt'+c].SetLineColor(ROOT.kAzure+1)
                self.histos[suff+'apo'+'UNRyqt'+c].SetMarkerColor(ROOT.kBlue-4)
                self.histos[suff+'apo'+'UNRyqt'+c].SetLineWidth(5)
                self.histos[suff+'apo'+'UNRyqt'+c].SetMarkerStyle(20)
                self.histos[suff+'apo'+'UNRyqt'+c].Draw("EX0 same")
            
            if toy!='' :
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toy'].SetLineColor(ROOT.kViolet+1)
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toy'].SetFillStyle(0)
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toy'].SetLineWidth(3)
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toy'].SetMarkerStyle(2)
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toy'].SetMarkerColor(ROOT.kViolet)
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toy'].Draw("E2 same")
            
            ROOT.gPad.Update()
            if 'unpol' in c and self.anaKind['angNames']: minvalMain = 0
            else : minvalMain = self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].GetFrame().GetY1()
            if ((not self.asimov) and c=='A4') : 
                self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetRangeUser(-2.3,1.5)                 
            else :
                self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetRangeUser(minvalMain,1.3*self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].GetFrame().GetY2())
            self.leg[suff+'FitAC'+'UNRyqt'+c].Draw("same")
            self.leg[suff+'FitAC'+'UNRyqt'+c].SetNColumns(2)  

            self.canvas['pr'+suff+'FitAC'+'UNRyqt'+c].cd()
            self.canvas['pr'+suff+'FitAC'+'UNRyqt'+c].SetGridx()
            self.canvas['pr'+suff+'FitAC'+'UNRyqt'+c].SetGridy()
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].Draw("EX0")
            self.histos[suff+'FitRatio'+'UNRyqt'+c].DrawCopy("same E2")
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'UNRyqt'+d+c].Draw("hist same")
                self.histos[suff+'FitRatioScale'+'UNRyqt'+d+c].Draw("hist same")
            self.histos[suff+'FitRatioAC'+'UNRyqt'+c].Draw("EX0 same") #to have foreground
            
            if aposteriori!='' :
                self.histos[suff+'apoRatio'+'UNRyqt'+c].SetLineColor(ROOT.kAzure+1)
                self.histos[suff+'apoRatio'+'UNRyqt'+c].SetMarkerColor(ROOT.kBlue-4)
                self.histos[suff+'apoRatio'+'UNRyqt'+c].SetLineWidth(4)
                self.histos[suff+'apoRatio'+'UNRyqt'+c].SetMarkerStyle(20)
                self.histos[suff+'apoRatio'+'UNRyqt'+c].Draw("EX0 same")
                self.histos[suff+'apoRatio'+'UNRyqt'+c].SetMarkerSize(1)
            
            self.leg[suff+'FitAC'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitAC'+'UNRyqt'+c], self.signDict[self.sign]+' Fit')
            self.leg[suff+'FitAC'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitBand'+'UNRyqt'+c], "MC Syst. Unc.")
            self.leg[suff+'FitAC'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitRatioPDF'+'UNRyqt'+'up'+c], self.groupedSystColors['LHEPdfWeightVars'][1])
            self.leg[suff+'FitAC'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitRatioScale'+'UNRyqt'+'up'+c], self.groupedSystColors['LHEScaleWeightVars'][1])
            if aposteriori!='' :
                self.leg[suff+'FitAC'+'UNRyqt'+c].AddEntry(self.histos[suff+'apo'+'UNRyqt'+c], "Post-fit regularized")
            if toy!='' :
                self.leg[suff+'FitAC'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitAC'+'UNRyqt'+c+'toy'], "MC toys  #mu#pm#sigma")
                
            self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].cd()
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].GetRightMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].GetLeftMargin()+0.07,1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRyqt'+c].GetTopMargin(),cmslab)

            
            # self.histos[suff+'FitBandPDF'+'UNRyqt'+c].DrawCopy('E2 same')
            # self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetFillStyle(3004)
            # self.histos[suff+'FitBandPDF'+'UNRyqt'+c].Draw("E2 same")
            # self.histos[suff+'FitBandScale'+'UNRyqt'+c].DrawCopy('E2 same')
            # self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetFillStyle(3005)
            # self.histos[suff+'FitBandScale'+'UNRyqt'+c].Draw("E2 same")
            
            
            
            #--------------------- Ai QT only canvas -------------------------#                
            self.canvas[suff+'FitAC'+'qt'+c] = ROOT.TCanvas(suff+"_c_QT_{}".format(c),suff+"_c_QT_{}".format(c),1200,1050)
            self.leg[suff+'FitAC'+'qt'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            self.canvas[suff+'FitAC'+'qt'+c].cd()
            self.canvas['ph'+suff+'FitAC'+'qt'+c] = ROOT.TPad('ph'+suff+"_c_QT_{}".format(c), 'ph'+suff+"_c_QT_{}".format(c),0,0.3,1,1)
            self.canvas['pr'+suff+'FitAC'+'qt'+c] = ROOT.TPad('pr'+suff+"_c_QT_{}".format(c),  'pr'+suff+"_c_QT_{}".format(c), 0,0,1,0.3)
            self.canvas['ph'+suff+'FitAC'+'qt'+c].SetBottomMargin(0.02)
            self.canvas['ph'+suff+'FitAC'+'qt'+c].Draw()
            self.canvas['pr'+suff+'FitAC'+'qt'+c].SetTopMargin(0.01)
            self.canvas['pr'+suff+'FitAC'+'qt'+c].SetBottomMargin(0.4)
            self.canvas['pr'+suff+'FitAC'+'qt'+c].Draw()
            
            self.canvas['ph'+suff+'FitAC'+'qt'+c].cd()
            self.canvas['ph'+suff+'FitAC'+'qt'+c].SetGridx()
            self.canvas['ph'+suff+'FitAC'+'qt'+c].SetGridy()
            # self.histos[suff+'FitAC'+'qt'+c].SetTitle(self.coeffDict[c][2]+", transverse momentum distribution, Y integrated"+", "+self.signDict[self.sign])
            self.histos[suff+'FitAC'+'qt'+c].SetTitle("")
            self.histos[suff+'FitAC'+'qt'+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
            self.histos[suff+'FitAC'+'qt'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'qt'+c].SetLineColor(ROOT.kBlack)
            self.histos[suff+'FitAC'+'qt'+c].SetStats(0)
            self.histos[suff+'FitAC'+'qt'+c].SetMarkerStyle(20)
            self.histos[suff+'FitAC'+'qt'+c].SetMarkerSize(0.5)
            self.histos[suff+'FitAC'+'qt'+c].Draw()
            self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetTitleSize(0.06)
            self.histos[suff+'FitAC'+'qt'+c].SetLabelSize(0.05,'y')
            if not 'unpol' in c and self.anaKind['name']=='angCoeff':
                # self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetRangeUser(-4,4)
                self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetTitle(self.coeffDict[c][2]+', '+yLab)
            elif self.anaKind['name']=='helXsec':
                    self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetTitle('d#sigma_'+c+'/dq_{T}^{W} [fb/GeV]'+', '+yLab)
            elif self.anaKind['name']=='normXsec':
                self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetTitle('d#sigma_'+c+'/dq_{T}^{W}/#sigma_{tot} [1/GeV]'+', '+yLab)
            elif self.anaKind['name']=='mu':
                self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetTitle('#mu_'+c+', '+yLab)
            else :
                 self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetTitle('d#sigma^{U+L}/dq_{T} [fb/GeV]'+', '+yLab)
                #  maxvalMain = self.histos[suff+'FitAC'+'qt'+c].GetMaximum()
                #  self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetRangeUser(0,maxvalMain*1.5)
            
            self.histos[suff+'FitBand'+'qt'+c].SetFillColor(ROOT.kOrange)
            self.histos[suff+'FitBand'+'qt'+c].SetFillStyle(0)
            self.histos[suff+'FitBand'+'qt'+c].SetLineColor(ROOT.kOrange)
            self.histos[suff+'FitBand'+'qt'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetTitleOffset(0.8)
            self.histos[suff+'FitAC'+'qt'+c].GetXaxis().SetTitleOffset(3)
            self.histos[suff+'FitAC'+'qt'+c].GetXaxis().SetLabelOffset(3)
            self.histos[suff+'FitBand'+'qt'+c].DrawCopy("E2 same")
            self.histos[suff+'FitBand'+'qt'+c].SetFillStyle(3001)
            self.histos[suff+'FitBand'+'qt'+c].Draw("E2 same")
            self.histos[suff+'FitAC'+'qt'+c].DrawCopy("E1 same") #to have foreground
            
            if aposteriori!='' :
                self.histos[suff+'apo'+'qt'+c].SetLineColor(ROOT.kAzure+1)
                self.histos[suff+'apo'+'qt'+c].SetMarkerColor(ROOT.kBlue-4)
                self.histos[suff+'apo'+'qt'+c].SetMarkerStyle(20)
                if not 'unpol' in c : 
                    self.histos[suff+'apo'+'qt'+c].SetLineWidth(1)
                    self.histos[suff+'apo'+'qt'+c].DrawCopy("Lhist same")
                    self.histos[suff+'apo'+'qt'+c].SetFillColor(ROOT.kAzure+1)
                    self.histos[suff+'apo'+'qt'+c].SetFillStyle(3002)
                    self.histos[suff+'apo'+'qt'+c].Draw("E3 same")
                    self.histos[suff+'apo'+'qt'+c].SetMarkerSize(1)
                else :
                    self.histos[suff+'apo'+'qt'+c].SetLineWidth(5)
                    self.histos[suff+'apo'+'qt'+c].Draw("EX0 same")
            
            if toy!='' :
                self.histos[suff+'FitAC'+'qt'+c+'toy'].SetLineColor(ROOT.kViolet+1)
                self.histos[suff+'FitAC'+'qt'+c+'toy'].SetFillStyle(0)
                self.histos[suff+'FitAC'+'qt'+c+'toy'].SetLineWidth(3)
                self.histos[suff+'FitAC'+'qt'+c+'toy'].SetMarkerStyle(2)
                self.histos[suff+'FitAC'+'qt'+c+'toy'].SetMarkerColor(ROOT.kViolet)
                self.histos[suff+'FitAC'+'qt'+c+'toy'].Draw("E2 same")
            
            ROOT.gPad.Update()
            if 'unpol' in c and self.anaKind['angNames']: minvalMain = 0
            else : minvalMain = self.canvas['ph'+suff+'FitAC'+'qt'+c].GetFrame().GetY1()
            self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetRangeUser(minvalMain,1.3*self.canvas['ph'+suff+'FitAC'+'qt'+c].GetFrame().GetY2())
            self.leg[suff+'FitAC'+'qt'+c].Draw("same")

            self.canvas['pr'+suff+'FitAC'+'qt'+c].cd()
            self.canvas['pr'+suff+'FitAC'+'qt'+c].SetGridx()
            self.canvas['pr'+suff+'FitAC'+'qt'+c].SetGridy()
            self.histos[suff+'FitRatioAC'+'qt'+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
            if ('unpol' in c or (not self.anaKind['differential'])): 
                self.histos[suff+'FitRatioAC'+'qt'+c].GetYaxis().SetTitle('Pred./Theory')
            else :
                self.histos[suff+'FitRatioAC'+'qt'+c].GetYaxis().SetTitle('Pred.-Theory')
            self.histos[suff+'FitRatioAC'+'qt'+c].SetMarkerStyle(20)
            self.histos[suff+'FitRatioAC'+'qt'+c].SetMarkerSize(0.5)
            self.histos[suff+'FitRatioAC'+'qt'+c].SetLineWidth(2)
            self.histos[suff+'FitRatio'+'qt'+c].SetLineWidth(2)
            self.histos[suff+'FitRatioAC'+'qt'+c].SetTitle("")
            self.histos[suff+'FitRatioAC'+'qt'+c].GetYaxis().SetTitleOffset(0.32)
            self.histos[suff+'FitRatioAC'+'qt'+c].SetTitleSize(0.13,'y')
            self.histos[suff+'FitRatioAC'+'qt'+c].SetLabelSize(0.10,'y')
            self.histos[suff+'FitRatioAC'+'qt'+c].GetXaxis().SetTitleOffset(1.2)
            self.histos[suff+'FitRatioAC'+'qt'+c].SetTitleSize(0.12,'x')
            self.histos[suff+'FitRatioAC'+'qt'+c].SetLabelSize(0.12,'x')
            self.histos[suff+'FitRatioAC'+'qt'+c].SetStats(0)
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'qt'+d+c].SetLineWidth(2)
                self.histos[suff+'FitRatioScale'+'qt'+d+c].SetLineWidth(2)
            self.histos[suff+'FitRatio'+'qt'+c].SetLineColor(ROOT.kOrange)
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'qt'+d+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                self.histos[suff+'FitRatioScale'+'qt'+d+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            self.histos[suff+'FitRatio'+'qt'+c].SetFillColor(ROOT.kOrange)
            self.histos[suff+'FitRatio'+'qt'+c].SetFillStyle(3001)
            
            maxvalRatio=0
            for xx in range(1,self.histos[suff+'FitRatioAC'+'qt'+c].GetNbinsX()+1) :
                if self.histos[suff+'FitRatioAC'+'qt'+c].GetBinError(xx)>maxvalRatio :
                    maxvalRatio= self.histos[suff+'FitRatioAC'+'qt'+c].GetBinError(xx)
                # self.histos[suff+'FitRatioAC'+'qt'+c].SetBinError(xx,0)
            if self.anaKind['angNames'] :
                    PadYCut =  self.RatioPadYcutDict[c][2]
            else : PadYCut = self.RatioPadYcut 
            if maxvalRatio>PadYCut :
                    maxvalRatio=PadYCut
            # self.histos[suff+'FitRatioAC'+'qt'+c].GetYaxis().SetRangeUser(-maxvalRatio*1.1+1,maxvalRatio*1.1+1)
            
            self.histos[suff+'FitRatioAC'+'qt'+c].Draw("EX0")
            self.histos[suff+'FitRatio'+'qt'+c].DrawCopy("same E2")
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'qt'+d+c].Draw("hist same")
                self.histos[suff+'FitRatioScale'+'qt'+d+c].Draw("hist same")
            self.histos[suff+'FitRatioAC'+'qt'+c].Draw("E1 same") #to have foreground
            
            if aposteriori!='' :
                self.histos[suff+'apoRatio'+'qt'+c].SetLineColor(ROOT.kAzure+1)
                self.histos[suff+'apoRatio'+'qt'+c].SetMarkerColor(ROOT.kBlue-4)
                self.histos[suff+'apoRatio'+'qt'+c].SetLineWidth(4)
                self.histos[suff+'apoRatio'+'qt'+c].SetMarkerStyle(20)
                self.histos[suff+'apoRatio'+'qt'+c].Draw("EX0 same")
                self.histos[suff+'apoRatio'+'qt'+c].SetMarkerSize(1)
            
            self.leg[suff+'FitAC'+'qt'+c].AddEntry(self.histos[suff+'FitAC'+'qt'+c], self.signDict[self.sign]+' Fit')
            self.leg[suff+'FitAC'+'qt'+c].AddEntry(self.histos[suff+'FitBand'+'qt'+c], "MC Syst. Unc.")
            self.leg[suff+'FitAC'+'qt'+c].AddEntry(self.histos[suff+'FitRatioPDF'+'qt'+'up'+c], self.groupedSystColors['LHEPdfWeightVars'][1])
            self.leg[suff+'FitAC'+'qt'+c].AddEntry(self.histos[suff+'FitRatioScale'+'qt'+'up'+c], self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitAC'+'qt'+c].SetNColumns(2) 
            if toy!='' :  
                self.leg[suff+'FitAC'+'qt'+c].AddEntry(self.histos[suff+'FitAC'+'qt'+c+'toy'], "MC toys  #mu#pm#sigma")
            if aposteriori!='' :
                    self.leg[suff+'FitAC'+'qt'+c].AddEntry(self.histos[suff+'apo'+'qt'+c], "Post-fit regularized")
                    
            self.canvas['ph'+suff+'FitAC'+'qt'+c].cd()
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas['ph'+suff+'FitAC'+'qt'+c].GetRightMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'qt'+c].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas['ph'+suff+'FitAC'+'qt'+c].GetLeftMargin()+0.07,1-0.8*self.canvas['ph'+suff+'FitAC'+'qt'+c].GetTopMargin(),cmslab)
            # cmsLatex.DrawLatex(self.canvas['ph'+suff+'FitAC'+'qt'+c].GetLeftMargin()+0.07,1-self.canvas['ph'+suff+'FitAC'+'qt'+c].GetTopMargin()-0.02,yLab)

        
            # self.histos[suff+'FitBandPDF'+'qt'+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            # self.histos[suff+'FitBandPDF'+'qt'+c].SetFillStyle(0)
            # self.histos[suff+'FitBandPDF'+'qt'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            # self.histos[suff+'FitBandPDF'+'qt'+c].SetLineWidth(2)
            # self.histos[suff+'FitBandPDF'+'qt'+c].DrawCopy("E2 same")
            # self.histos[suff+'FitBandPDF'+'qt'+c].SetFillStyle(3004)
            # self.histos[suff+'FitBandPDF'+'qt'+c].Draw("E2 same")
            # self.histos[suff+'FitBandScale'+'qt'+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            # self.histos[suff+'FitBandScale'+'qt'+c].SetFillStyle(0)
            # self.histos[suff+'FitBandScale'+'qt'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            # self.histos[suff+'FitBandScale'+'qt'+c].SetLineWidth(2)
            # self.histos[suff+'FitBandScale'+'qt'+c].DrawCopy("E2 same")
            # self.histos[suff+'FitBandScale'+'qt'+c].SetFillStyle(3005)
            # self.histos[suff+'FitBandScale'+'qt'+c].Draw("E2 same")
            
            
            #--------------------- Ai Y only canvas -------------------------#          
            self.canvas[suff+'FitAC'+'y'+c] = ROOT.TCanvas(suff+"_c_Y_{}".format(c),suff+"_c_Y_{}".format(c),1200,1050)
            self.leg[suff+'FitAC'+'y'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            self.canvas[suff+'FitAC'+'y'+c].cd()
            self.canvas['ph'+suff+'FitAC'+'y'+c] = ROOT.TPad('ph'+suff+"_c_Y_{}".format(c), 'ph'+suff+"_c_Y_{}".format(c),0,0.3,1,1)
            self.canvas['pr'+suff+'FitAC'+'y'+c] = ROOT.TPad('pr'+suff+"_c_Y_{}".format(c),  'pr'+suff+"_c_Y_{}".format(c), 0,0,1,0.3)
            self.canvas['ph'+suff+'FitAC'+'y'+c].SetBottomMargin(0.02)
            self.canvas['ph'+suff+'FitAC'+'y'+c].Draw()
            self.canvas['pr'+suff+'FitAC'+'y'+c].SetTopMargin(0.01)
            self.canvas['pr'+suff+'FitAC'+'y'+c].SetBottomMargin(0.4)
            self.canvas['pr'+suff+'FitAC'+'y'+c].Draw()
            
            self.canvas['ph'+suff+'FitAC'+'y'+c].cd()
            self.canvas['ph'+suff+'FitAC'+'y'+c].SetGridx()
            self.canvas['ph'+suff+'FitAC'+'y'+c].SetGridy()
            # self.histos[suff+'FitAC'+'y'+c].SetTitle(self.coeffDict[c][2]+", rapidity distribution,  q_{T} integrated"+", "+self.signDict[self.sign])
            self.histos[suff+'FitAC'+'y'+c].SetTitle("")
            self.histos[suff+'FitAC'+'y'+c].GetXaxis().SetTitle('|Y_{W}|')
            self.histos[suff+'FitAC'+'y'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'y'+c].SetLineColor(ROOT.kBlack)
            self.histos[suff+'FitAC'+'y'+c].SetStats(0)
            self.histos[suff+'FitAC'+'y'+c].SetMarkerStyle(20)
            self.histos[suff+'FitAC'+'y'+c].SetMarkerSize(0.5)
            self.histos[suff+'FitAC'+'y'+c].Draw()
            self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetTitleSize(0.06)
            self.histos[suff+'FitAC'+'y'+c].SetLabelSize(0.05,'y')
            if not 'unpol' in c and self.anaKind['name']=='angCoeff':
                # self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetRangeUser(-4,4)
                self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetTitle(self.coeffDict[c][2]+', '+qtLab)
            elif self.anaKind['name']=='helXsec':
                    self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetTitle('d#sigma_'+c+'/d|Y_{W}| [fb]'+', '+qtLab)
            elif self.anaKind['name']=='normXsec':
                self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetTitle('d#sigma_'+c+'/d|Y_{W}|/#sigma_{tot}'+', '+qtLab)
            elif self.anaKind['name']=='mu':
                self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetTitle('#mu_'+c+', '+qtLab)
            else :
                self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetTitle('d#sigma^{U+L}/d|Y_{W}| [fb]'+', '+qtLab)
                # maxvalMain = self.histos[suff+'FitAC'+'y'+c].GetMaximum()
                # self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetRangeUser(0,maxvalMain*1.5)
            self.histos[suff+'FitBand'+'y'+c].SetFillColor(ROOT.kOrange)#kMagenta-7)
            self.histos[suff+'FitBand'+'y'+c].SetFillStyle(0)
            self.histos[suff+'FitBand'+'y'+c].SetLineColor(ROOT.kOrange)#kMagenta-7)
            self.histos[suff+'FitBand'+'y'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetTitleOffset(0.8)
            self.histos[suff+'FitAC'+'y'+c].GetXaxis().SetTitleOffset(3)
            self.histos[suff+'FitAC'+'y'+c].GetXaxis().SetLabelOffset(3)
            self.histos[suff+'FitBand'+'y'+c].DrawCopy("E2 same")
            self.histos[suff+'FitBand'+'y'+c].SetFillStyle(3001)
            self.histos[suff+'FitBand'+'y'+c].Draw("E2 same")
            self.histos[suff+'FitAC'+'y'+c].DrawCopy("E1 same") #to have foreground
            
            if aposteriori!='' :
                self.histos[suff+'apo'+'y'+c].SetLineColor(ROOT.kAzure+1)
                self.histos[suff+'apo'+'y'+c].SetMarkerColor(ROOT.kBlue-4)
                self.histos[suff+'apo'+'y'+c].SetMarkerStyle(20)
                if not 'unpol' in c : 
                    self.histos[suff+'apo'+'y'+c].SetLineWidth(1)
                    self.histos[suff+'apo'+'y'+c].DrawCopy("Lhist same")
                    self.histos[suff+'apo'+'y'+c].SetFillColor(ROOT.kAzure+1)
                    self.histos[suff+'apo'+'y'+c].SetFillStyle(3002)
                    self.histos[suff+'apo'+'y'+c].Draw("E3 same")
                    self.histos[suff+'apo'+'y'+c].SetMarkerSize(1)
                else :
                    self.histos[suff+'apo'+'y'+c].SetLineWidth(5)
                    self.histos[suff+'apo'+'y'+c].Draw("EX0 same")

                        
            if toy!='' :
                self.histos[suff+'FitAC'+'y'+c+'toy'].SetLineColor(ROOT.kViolet+1)
                self.histos[suff+'FitAC'+'y'+c+'toy'].SetFillStyle(0)
                self.histos[suff+'FitAC'+'y'+c+'toy'].SetLineWidth(3)
                self.histos[suff+'FitAC'+'y'+c+'toy'].SetMarkerStyle(2)
                self.histos[suff+'FitAC'+'y'+c+'toy'].SetMarkerColor(ROOT.kViolet)
                self.histos[suff+'FitAC'+'y'+c+'toy'].Draw("E2 same")
            
            ROOT.gPad.Update()
            if 'unpol' in c and self.anaKind['angNames']: minvalMain = 0
            else : minvalMain = self.canvas['ph'+suff+'FitAC'+'y'+c].GetFrame().GetY1()
            self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetRangeUser(minvalMain,1.3*self.canvas['ph'+suff+'FitAC'+'y'+c].GetFrame().GetY2())
            self.leg[suff+'FitAC'+'y'+c].Draw("same")
            
            self.canvas['pr'+suff+'FitAC'+'y'+c].cd()
            self.canvas['pr'+suff+'FitAC'+'y'+c].SetGridx()
            self.canvas['pr'+suff+'FitAC'+'y'+c].SetGridy()
            self.histos[suff+'FitRatioAC'+'y'+c].GetXaxis().SetTitle('|Y_{W}|')
            if ('unpol' in c or (not self.anaKind['differential'])): 
                self.histos[suff+'FitRatioAC'+'y'+c].GetYaxis().SetTitle('Pred./Theory')
            else :
                self.histos[suff+'FitRatioAC'+'y'+c].GetYaxis().SetTitle('Pred.-Theory')
            self.histos[suff+'FitRatioAC'+'y'+c].SetMarkerStyle(20)
            self.histos[suff+'FitRatioAC'+'y'+c].SetMarkerSize(0.5)
            self.histos[suff+'FitRatioAC'+'y'+c].SetLineWidth(2)
            self.histos[suff+'FitRatio'+'y'+c].SetLineWidth(2)
            self.histos[suff+'FitRatioAC'+'y'+c].SetTitle("")
            self.histos[suff+'FitRatioAC'+'y'+c].GetYaxis().SetTitleOffset(0.32)
            self.histos[suff+'FitRatioAC'+'y'+c].SetTitleSize(0.13,'y')
            self.histos[suff+'FitRatioAC'+'y'+c].SetLabelSize(0.10,'y')
            self.histos[suff+'FitRatioAC'+'y'+c].GetXaxis().SetTitleOffset(1.2)
            self.histos[suff+'FitRatioAC'+'y'+c].SetTitleSize(0.12,'x')
            self.histos[suff+'FitRatioAC'+'y'+c].SetLabelSize(0.12,'x')
            self.histos[suff+'FitRatioAC'+'y'+c].SetStats(0)
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'y'+d+c].SetLineWidth(2)
                self.histos[suff+'FitRatioScale'+'y'+d+c].SetLineWidth(2)
            self.histos[suff+'FitRatio'+'y'+c].SetLineColor(ROOT.kOrange)
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'y'+d+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                self.histos[suff+'FitRatioScale'+'y'+d+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            self.histos[suff+'FitRatio'+'y'+c].SetFillColor(ROOT.kOrange)
            self.histos[suff+'FitRatio'+'y'+c].SetFillStyle(3001)
            
            maxvalRatio=0
            for xx in range(1,self.histos[suff+'FitRatioAC'+'y'+c].GetNbinsX()+1) :
                if 'unpol' in c : ACstring = '' #for unpol integrated in qt the uncertainity of MC is higher than prediction--> used the former as limit
                else : ACstring ='AC'
                if self.histos[suff+'FitRatio'+ACstring+'y'+c].GetBinError(xx)>maxvalRatio :
                    maxvalRatio= self.histos[suff+'FitRatio'+ACstring+'y'+c].GetBinError(xx)
                # self.histos[suff+'FitRatioAC'+'y'+c].SetBinError(xx,0)
            if self.anaKind['angNames'] :
                    PadYCut =  self.RatioPadYcutDict[c][1]
            else : PadYCut = self.RatioPadYcut 
            if maxvalRatio>PadYCut :
                    maxvalRatio=PadYCut
            # if 'unpol' in c : #not needed if used ratio for everything 
            if 'unpol' in c :    
                self.histos[suff+'FitRatioAC'+'y'+c].GetYaxis().SetRangeUser(-maxvalRatio*1.1+1,maxvalRatio*1.1+1)
            
            self.histos[suff+'FitRatioAC'+'y'+c].Draw("EX0")
            self.histos[suff+'FitRatio'+'y'+c].DrawCopy("same E2")
            for d in self.dirList :
                self.histos[suff+'FitRatioPDF'+'y'+d+c].Draw("hist same")
                self.histos[suff+'FitRatioScale'+'y'+d+c].Draw("hist same")
            self.histos[suff+'FitRatioAC'+'y'+c].Draw("E1 same") #to have foreground
            
            if aposteriori!='' :
                self.histos[suff+'apoRatio'+'y'+c].SetLineColor(ROOT.kAzure+1)
                self.histos[suff+'apoRatio'+'y'+c].SetMarkerColor(ROOT.kBlue-4)
                self.histos[suff+'apoRatio'+'y'+c].SetLineWidth(4)
                self.histos[suff+'apoRatio'+'y'+c].SetMarkerStyle(20)
                self.histos[suff+'apoRatio'+'y'+c].Draw("EX0 same")
                self.histos[suff+'apoRatio'+'y'+c].SetMarkerSize(1)
            
            self.leg[suff+'FitAC'+'y'+c].AddEntry(self.histos[suff+'FitAC'+'y'+c], self.signDict[self.sign]+' Fit')
            self.leg[suff+'FitAC'+'y'+c].AddEntry(self.histos[suff+'FitBand'+'y'+c], "MC Syst. Unc.")
            self.leg[suff+'FitAC'+'y'+c].AddEntry(self.histos[suff+'FitRatioPDF'+'y'+'up'+c], self.groupedSystColors['LHEPdfWeightVars'][1])
            self.leg[suff+'FitAC'+'y'+c].AddEntry(self.histos[suff+'FitRatioScale'+'y'+'up'+c], self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitAC'+'y'+c].SetNColumns(2)    
            if toy!='' :  
                self.leg[suff+'FitAC'+'y'+c].AddEntry(self.histos[suff+'FitAC'+'y'+c+'toy'], "MC toys  #mu#pm#sigma")
            if aposteriori!='' :
                    self.leg[suff+'FitAC'+'y'+c].AddEntry(self.histos[suff+'apo'+'y'+c], "Post-fit regularized")
            
            self.canvas['ph'+suff+'FitAC'+'y'+c].cd()
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas['ph'+suff+'FitAC'+'y'+c].GetRightMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'y'+c].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas['ph'+suff+'FitAC'+'y'+c].GetLeftMargin()+0.07,1-0.8*self.canvas['ph'+suff+'FitAC'+'y'+c].GetTopMargin(),cmslab)
            # cmsLatex.DrawLatex(self.canvas['ph'+suff+'FitAC'+'qt'+c].GetLeftMargin()+0.07,1-self.canvas['ph'+suff+'FitAC'+'qt'+c].GetTopMargin()-0.02,qtLab)

            
            # self.histos[suff+'FitBandPDF'+'y'+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
            # self.histos[suff+'FitBandPDF'+'y'+c].SetFillStyle(0)
            # self.histos[suff+'FitBandPDF'+'y'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
            # self.histos[suff+'FitBandPDF'+'y'+c].SetLineWidth(2)
            # self.histos[suff+'FitBandPDF'+'y'+c].DrawCopy("E2 same")
            # self.histos[suff+'FitBandPDF'+'y'+c].SetFillStyle(3004)
            # self.histos[suff+'FitBandPDF'+'y'+c].Draw("E2 same")
            # self.histos[suff+'FitBandScale'+'y'+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
            # self.histos[suff+'FitBandScale'+'y'+c].SetFillStyle(0)
            # self.histos[suff+'FitBandScale'+'y'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
            # self.histos[suff+'FitBandScale'+'y'+c].SetLineWidth(2)
            # self.histos[suff+'FitBandScale'+'y'+c].DrawCopy("E2 same")
            # self.histos[suff+'FitBandScale'+'y'+c].SetFillStyle(3005)
            # self.histos[suff+'FitBandScale'+'y'+c].Draw("E2 same")
            
            
        
        
             #--------------- Ai 1D relative syst. band canvas -------------------#
            for i in range(1, self.histos[suff+'FitErr'+c].GetNbinsX()+1): #loop over rapidity bins
                self.canvas[suff+'FitErr'+'qt'+str(i)+c] = ROOT.TCanvas(suff+"_c_projQt{}_{}_Err".format(i,c),suff+"_c_ErrprojQt{}_{}_Err".format(i,c),1200,900)
                self.canvas[suff+'FitErr'+'qt'+str(i)+c].SetGridx()
                self.canvas[suff+'FitErr'+'qt'+str(i)+c].SetGridy()
                self.leg[suff+'FitErr'+'qt'+str(i)+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
                self.histos[suff+'FitErr'+'qt'+str(i)+c] = self.histos[suff+'FitErr'+c].ProjectionY("projQt{}_{}_Err".format(i,c),i,i)
                self.histos[suff+'FitErrPDF'+'qt'+str(i)+c] = self.histos[suff+'FitErrPDF'+c].ProjectionY("projQt{}_{}_ErrPDF".format(i,c),i,i)
                self.histos[suff+'FitErrScale'+'qt'+str(i)+c] = self.histos[suff+'FitErrScale'+c].ProjectionY("projQt{}_{}_ErrScale".format(i,c),i,i)
                self.histos[suff+'FitErr'+'qt'+str(i)+c].SetTitle(self.coeffDict[c][2]+" relative uncertainity, transverse momentum distribution, "+str(self.yArr[i-1])+"<|Y_{W}|<"+str(self.yArr[i])+", "+self.signDict[self.sign])
                self.canvas[suff+'FitErr'+'qt'+str(i)+c].cd()
                self.histos[suff+'FitErr'+'qt'+str(i)+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
                self.histos[suff+'FitErr'+'qt'+str(i)+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
                self.histos[suff+'FitErr'+'qt'+str(i)+c].GetYaxis().SetTitle('Relative uncertainty on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                self.histos[suff+'FitErr'+'qt'+str(i)+c].GetYaxis().SetRangeUser(0.001,7000)
                self.canvas[suff+'FitErr'+'qt'+str(i)+c].SetLogy()
                self.histos[suff+'FitErr'+'qt'+str(i)+c].SetLineWidth(3)
                self.histos[suff+'FitErr'+'qt'+str(i)+c].SetStats(0)
                self.histos[suff+'FitErr'+'qt'+str(i)+c].Draw('hist')
                self.histos[suff+'FitErr'+'qt'+str(i)+c].SetLineColor(self.groupedSystColors['Nominal'][0])
                self.histos[suff+'FitErrPDF'+'qt'+str(i)+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                self.histos[suff+'FitErrScale'+'qt'+str(i)+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
                self.histos[suff+'FitErrPDF'+'qt'+str(i)+c].SetLineWidth(3)
                self.histos[suff+'FitErrScale'+'qt'+str(i)+c].SetLineWidth(3)
                self.histos[suff+'FitErrPDF'+'qt'+str(i)+c].Draw("hist same")
                self.histos[suff+'FitErrScale'+'qt'+str(i)+c].Draw("hist same")
                self.leg[suff+'FitErr'+'qt'+str(i)+c].AddEntry(self.histos[suff+'FitErr'+'qt'+str(i)+c],self.groupedSystColors['Nominal'][1])
                self.leg[suff+'FitErr'+'qt'+str(i)+c].AddEntry(self.histos[suff+'FitErrPDF'+'qt'+str(i)+c],self.groupedSystColors['LHEPdfWeightVars'][1])
                self.leg[suff+'FitErr'+'qt'+str(i)+c].AddEntry(self.histos[suff+'FitErrScale'+'qt'+str(i)+c],self.groupedSystColors['LHEScaleWeightVars'][1])
                self.leg[suff+'FitErr'+'qt'+str(i)+c].Draw("Same")
            for j in range(1, self.histos[suff+'FitErr'+c].GetNbinsY()+1): #loop over pt bins
                self.canvas[suff+'FitErr'+'y'+str(j)+c] = ROOT.TCanvas(suff+"_c_projY{}_{}_Err".format(j,c),suff+"_c_projY{}_{}_Err".format(j,c),1200,900)
                self.canvas[suff+'FitErr'+'y'+str(j)+c].SetGridx()
                self.canvas[suff+'FitErr'+'y'+str(j)+c].SetGridy()
                self.leg[suff+'FitErr'+'y'+str(j)+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
                self.histos[suff+'FitErr'+'y'+str(j)+c] = self.histos[suff+'FitErr'+c].ProjectionX("projY{}_{}".format(j,c),j,j)
                self.histos[suff+'FitErrPDF'+'y'+str(j)+c] = self.histos[suff+'FitErrPDF'+c].ProjectionX("projY{}_{}_ErrPDF".format(j,c),j,j)
                self.histos[suff+'FitErrScale'+'y'+str(j)+c] = self.histos[suff+'FitErrScale'+c].ProjectionX("projY{}_{}_ErrScale".format(j,c),j,j)
                self.histos[suff+'FitErr'+'y'+str(j)+c].SetTitle(self.coeffDict[c][2]+"relative uncertainity, rapidity distribution, "+str(self.qtArr[j-1])+"<q_{T}^{W}<"+str(self.qtArr[j])+", "+self.signDict[self.sign])
                self.canvas[suff+'FitErr'+'y'+str(j)+c].cd()
                self.histos[suff+'FitErr'+'y'+str(j)+c].GetXaxis().SetTitle('|Y_{W}|')
                # self.histos[suff+'FitErr'+'y'+str(j)+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
                self.histos[suff+'FitErr'+'y'+str(j)+c].GetYaxis().SetTitle('Relative uncertainty on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                self.histos[suff+'FitErr''y'+str(j)+c].GetYaxis().SetRangeUser(0.001,7000)
                self.canvas[suff+'FitErr''y'+str(j)+c].SetLogy()
                self.histos[suff+'FitErr'+'y'+str(j)+c].SetLineWidth(3)
                self.histos[suff+'FitErr'+'y'+str(j)+c].SetStats(0)
                self.histos[suff+'FitErr'+'y'+str(j)+c].Draw('hist')
                self.histos[suff+'FitErr'+'y'+str(j)+c].SetLineColor(self.groupedSystColors['Nominal'][0])
                self.histos[suff+'FitErrPDF'+'y'+str(j)+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                self.histos[suff+'FitErrScale'+'y'+str(j)+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
                self.histos[suff+'FitErrPDF'+'y'+str(j)+c].SetLineWidth(3)
                self.histos[suff+'FitErrScale'+'y'+str(j)+c].SetLineWidth(3)
                self.histos[suff+'FitErrPDF'+'y'+str(j)+c].Draw("hist same")
                self.histos[suff+'FitErrScale'+'y'+str(j)+c].Draw("hist same")
                self.leg[suff+'FitErr'+'y'+str(j)+c].AddEntry(self.histos[suff+'FitErr'+'y'+str(j)+c],self.groupedSystColors['Nominal'][1])
                self.leg[suff+'FitErr'+'y'+str(j)+c].AddEntry(self.histos[suff+'FitErrPDF'+'y'+str(j)+c],self.groupedSystColors['LHEPdfWeightVars'][1])
                self.leg[suff+'FitErr'+'y'+str(j)+c].AddEntry(self.histos[suff+'FitErrScale'+'y'+str(j)+c],self.groupedSystColors['LHEScaleWeightVars'][1])
                self.leg[suff+'FitErr'+'y'+str(j)+c].Draw("Same")
                
                
            #-------------- unrolled: qt large, y small - relative syst. band canvas ----------------- #
            self.histos[suff+'FitErr'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_Err_UNRqty'+c,suff+'_coeff_Err_UNRqty'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            self.histos[suff+'FitErrPDF'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_ErrPDF_UNRqty'+c,suff+'_coeff_ErrPDF_UNRqty'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            self.histos[suff+'FitErrScale'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_ErrScale_UNRqty'+c,suff+'_coeff_ErrScale_UNRqty'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            for q in self.qtArr[:-1] :
                for y in self.yArr[:-1] :
                    indexUNRqty = self.qtArr.index(float(q))*(len(self.yArr)-1)+self.yArr.index(float(y))
                    self.histos[suff+'FitErr'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitErr'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitErrPDF'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitErrPDF'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitErrScale'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitErrScale'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    if self.yArr.index(y)==0 :
                        self.histos[suff+'FitErr'+'UNRqty'+c].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+'FitErr'+'UNRqty'+c].GetXaxis().SetBinLabel(indexUNRqty+1,"q_{T}^{W}#in[%.0f,%.0f]" % (q, self.qtArr[self.qtArr.index(q)+1]))
                        self.histos[suff+'FitErr'+'UNRqty'+c].GetXaxis().ChangeLabel(indexUNRqty+1,340,0.03)    
                        self.histos[suff+'FitErr'+'UNRqty'+c].GetXaxis().LabelsOption("d")
            # self.histos[suff+'FitErr'+'UNRqty'+c].SetTitle(self.coeffDict[c][2]+", unrolled q_{T}(Y), "+self.signDict[self.sign])
            self.histos[suff+'FitErr'+'UNRqty'+c].SetTitle('')
            self.histos[suff+'FitErr'+'UNRqty'+c].GetXaxis().SetTitle('fine binning: |Y_{W}| 0#rightarrow 2.4')
            self.histos[suff+'FitErr'+'UNRqty'+c].GetXaxis().SetTitleOffset(1.45)
            self.histos[suff+'FitErr'+'UNRqty'+c].GetYaxis().SetTitle('Relative uncertainty on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
            self.histos[suff+'FitErr'+'UNRqty'+c].SetLineWidth(3)
            self.histos[suff+'FitErrPDF'+'UNRqty'+c].SetLineWidth(3)
            self.histos[suff+'FitErrScale'+'UNRqty'+c].SetLineWidth(3)
            self.histos[suff+'FitErr'+'UNRqty'+c].SetStats(0)
            self.histos[suff+'FitErr'+'UNRqty'+c].SetLineColor(self.groupedSystColors['Nominal'][0])
            self.histos[suff+'FitErrPDF'+'UNRqty'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            self.histos[suff+'FitErrScale'+'UNRqty'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])

            self.canvas[suff+'FitErr'+'UNRqty'+c] = ROOT.TCanvas(suff+'_c_Err_UNRqty'+c,suff+'_c_Err_UNRqty'+c,1200,900)
            self.leg[suff+'FitErr'+'UNRqty'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            self.canvas[suff+'FitErr'+'UNRqty'+c].cd()
            self.canvas[suff+'FitErr'+'UNRqty'+c].SetGridx()
            self.canvas[suff+'FitErr'+'UNRqty'+c].SetGridy()
            self.histos[suff+'FitErr'+'UNRqty'+c].GetYaxis().SetRangeUser(0.001,7000)
            self.canvas[suff+'FitErr'+'UNRqty'+c].SetLogy()
            self.histos[suff+'FitErr'+'UNRqty'+c].Draw()
            if not impact : self.histos[suff+'FitErrPDF'+'UNRqty'+c].Draw("same")
            if not impact : self.histos[suff+'FitErrScale'+'UNRqty'+c].Draw("same")
                
            self.leg[suff+'FitErr'+'UNRqty'+c].AddEntry(self.histos[suff+'FitErr'+'UNRqty'+c],self.groupedSystColors['Nominal'][1])
            if not impact : self.leg[suff+'FitErr'+'UNRqty'+c].AddEntry(self.histos[suff+'FitErrPDF'+'UNRqty'+c],self.groupedSystColors['LHEPdfWeightVars'][1])
            if not impact : self.leg[suff+'FitErr'+'UNRqty'+c].AddEntry(self.histos[suff+'FitErrScale'+'UNRqty'+c],self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitErr'+'UNRqty'+c].SetLineWidth(0)
            self.leg[suff+'FitErr'+'UNRqty'+c].SetFillStyle(0)
            # self.leg[suff+'FitErr'+'UNRqty'+c].SetNColumns(2)
            self.leg[suff+'FitErr'+'UNRqty'+c].Draw("Same")

            if 'unpol' in c :
                if impact : self.histos[suff+'impact'+'UNRqty'+c+'stat'].Draw("samehist")
                if impact : self.leg[suff+'FitErr'+'UNRqty'+c].AddEntry(self.histos[suff+'impact'+'UNRqty'+c+'stat'],'Stat. only')
                self.histos[suff+'impact'+'UNRqty'+c+'poiss'].Draw("same")
                self.histos[suff+'impact'+'UNRqty'+c+'poiss'].SetLineWidth(3)
                self.histos[suff+'impact'+'UNRqty'+c+'poiss'].SetLineColor(self.groupedSystColors['poiss'][0])
                self.leg[suff+'FitErr'+'UNRqty'+c].AddEntry(self.histos[suff+'impact'+'UNRqty'+c+'poiss'],self.groupedSystColors['poiss'][1])
                self.histos[suff+'FitErr'+'UNRqty'+c].GetYaxis().SetRangeUser(0.0001,1)

            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetRightMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetLeftMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetTopMargin(),cmslab)    



            #------------- unrolled: y large, qt small - relative syst. band canvas -----------------                    
            self.histos[suff+'FitErr'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_Err_UNRyqt'+c,suff+'_coeff_Err_UNRyqt'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            self.histos[suff+'FitErrPDF'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_ErrPDF_UNRyqt'+c,suff+'_coeff_ErrPDF_UNRyqt'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            self.histos[suff+'FitErrScale'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_ErrScale_UNRyqt'+c,suff+'_coeff_ErrScale_UNRyqt'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            for y in self.yArr[:-1] :
                for q in self.qtArr[:-1] :
                    indexUNRyqt = self.yArr.index(float(y))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                    self.histos[suff+'FitErr'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitErr'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitErrPDF'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitErrPDF'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitErrScale'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitErrScale'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    if self.qtArr.index(q)==0 :
                        self.histos[suff+'FitErr'+'UNRyqt'+c].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+'FitErr'+'UNRyqt'+c].GetXaxis().SetBinLabel(indexUNRyqt+1,"|Y_{W}|#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
                        self.histos[suff+'FitErr'+'UNRyqt'+c].GetXaxis().ChangeLabel(indexUNRyqt+1,340,0.03)
                        self.histos[suff+'FitErr'+'UNRyqt'+c].GetXaxis().LabelsOption("d")
            # self.histos[suff+'FitErr'+'UNRyqt'+c].SetTitle(self.coeffDict[c][2]+", unrolled Y(q_{T}), "+self.signDict[self.sign])
            self.histos[suff+'FitErr'+'UNRyqt'+c].SetTitle('')
            self.histos[suff+'FitErr'+'UNRyqt'+c].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 60 GeV')
            self.histos[suff+'FitErr'+'UNRyqt'+c].GetXaxis().SetTitleOffset(1.45)
            # self.histos[suff+'FitErr'+'UNRyqt'+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
            self.histos[suff+'FitErr'+'UNRyqt'+c].GetYaxis().SetTitle('Relative uncertainty on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
            self.histos[suff+'FitErr'+'UNRyqt'+c].SetLineWidth(3)
            self.histos[suff+'FitErrPDF'+'UNRyqt'+c].SetLineWidth(3)
            self.histos[suff+'FitErrScale'+'UNRyqt'+c].SetLineWidth(3)
            self.histos[suff+'FitErr'+'UNRyqt'+c].SetStats(0)
            self.histos[suff+'FitErr'+'UNRyqt'+c].SetLineColor(self.groupedSystColors['Nominal'][0])
            self.histos[suff+'FitErrPDF'+'UNRyqt'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            self.histos[suff+'FitErrScale'+'UNRyqt'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            
            self.canvas[suff+'FitErr'+'UNRyqt'+c] = ROOT.TCanvas(suff+'_c_Err_UNRyqt'+c,suff+'_c_Err_UNRyqt'+c,1200,900)
            self.leg[suff+'FitErr'+'UNRyqt'+c] = ROOT.TLegend(0.3,0.75,0.7,0.9)
            self.canvas[suff+'FitErr'+'UNRyqt'+c].cd()
            # self.canvas[suff+'FitErr'+'UNRyqt'+c].SetGridx()
            self.canvas[suff+'FitErr'+'UNRyqt'+c].SetGridy()
            self.histos[suff+'FitErr'+'UNRyqt'+c].GetYaxis().SetRangeUser(0.001,7000)
            self.canvas[suff+'FitErr'+'UNRyqt'+c].SetLogy()
            self.histos[suff+'FitErr'+'UNRyqt'+c].Draw()
            if not impact : self.histos[suff+'FitErrPDF'+'UNRyqt'+c].Draw("same")
            if not impact : self.histos[suff+'FitErrScale'+'UNRyqt'+c].Draw("same")
            self.leg[suff+'FitErr'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitErr'+'UNRyqt'+c],self.groupedSystColors['Nominal'][1])
            if not impact : self.leg[suff+'FitErr'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitErrPDF'+'UNRyqt'+c],self.groupedSystColors['LHEPdfWeightVars'][1])
            if not impact : self.leg[suff+'FitErr'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitErrScale'+'UNRyqt'+c],self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitErr'+'UNRyqt'+c].SetLineWidth(0)
            self.leg[suff+'FitErr'+'UNRyqt'+c].SetFillStyle(0)
            # self.leg[suff+'FitErr'+'UNRyqt'+c].SetNColumns(2)
            self.leg[suff+'FitErr'+'UNRyqt'+c].Draw("Same")
            
            if 'unpol' in c :
                if impact : self.histos[suff+'impact'+'UNRyqt'+c+'stat'].Draw("samehist")
                if impact : self.leg[suff+'FitErr'+'UNRyqt'+c].AddEntry(self.histos[suff+'impact'+'UNRyqt'+c+'stat'],'Stat. only')
                self.histos[suff+'impact'+'UNRyqt'+c+'poiss'].Draw("same")
                self.histos[suff+'impact'+'UNRyqt'+c+'poiss'].SetLineWidth(3)
                self.histos[suff+'impact'+'UNRyqt'+c+'poiss'].SetLineColor(self.groupedSystColors['poiss'][0])
                self.leg[suff+'FitErr'+'UNRyqt'+c].AddEntry(self.histos[suff+'impact'+'UNRyqt'+c+'poiss'],self.groupedSystColors['poiss'][1])
                self.histos[suff+'FitErr'+'UNRyqt'+c].GetYaxis().SetRangeUser(0.0001,1)
                bottomY4lineErr = 0.0001
                topY4lineErr = 1
                vErrLines = []
                for yqt in self.unrolledYQt :
                    if yqt%(self.qtArr[-1])==0 :
                        if yqt==0 :continue
                        if yqt==len(self.unrolledYQt)-1 : continue 
                        vErrLines.append(ROOT.TLine(yqt-2,bottomY4lineErr,yqt-2,topY4lineErr))
                        vErrLines[-1].SetLineStyle(ROOT.kDotted)
                        vErrLines[-1].DrawLine(yqt-2,bottomY4lineErr,yqt-2,topY4lineErr)

                
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetRightMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetLeftMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetTopMargin(),cmslab)
            
            
            #--------------- Ai QT only relative syst. band canvas -------------------#
            self.canvas[suff+'FitErr'+'qt'+c] = ROOT.TCanvas(suff+"_c_QT_{}_Err".format(c),suff+"_c_QT_{}_Err".format(c),1200,900)
            self.canvas[suff+'FitErr'+'qt'+c].SetGridx()
            self.canvas[suff+'FitErr'+'qt'+c].SetGridy()
            self.leg[suff+'FitErr'+'qt'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            # self.histos[suff+'FitErr'+'qt'+c].SetTitle(self.coeffDict[c][2]+" relative uncertainity, transverse momentum distribution, Y integrated"+", "+self.signDict[self.sign])
            self.histos[suff+'FitErr'+'qt'+c].SetTitle('')
            self.canvas[suff+'FitErr'+'qt'+c].cd()
            self.histos[suff+'FitErr'+'qt'+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
            # self.histos[suff+'FitErr'+'qt'+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
            self.histos[suff+'FitErr'+'qt'+c].GetYaxis().SetTitle('Relative uncertainty on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
            self.histos[suff+'FitErr'+'qt'+c].GetYaxis().SetRangeUser(0.001,7000)
            self.canvas[suff+'FitErr'+'qt'+c].SetLogy()
            self.histos[suff+'FitErr'+'qt'+c].SetLineWidth(3)
            self.histos[suff+'FitErr'+'qt'+c].SetStats(0)
            self.histos[suff+'FitErr'+'qt'+c].Draw('hist')
            self.histos[suff+'FitErr'+'qt'+c].SetLineColor(self.groupedSystColors['Nominal'][0])
            self.histos[suff+'FitErrPDF'+'qt'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            self.histos[suff+'FitErrScale'+'qt'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            self.histos[suff+'FitErrPDF'+'qt'+c].SetLineWidth(3)
            self.histos[suff+'FitErrScale'+'qt'+c].SetLineWidth(3)
            if not impact : self.histos[suff+'FitErrPDF'+'qt'+c].Draw("hist same")
            if not impact : self.histos[suff+'FitErrScale'+'qt'+c].Draw("hist same")
            self.leg[suff+'FitErr'+'qt'+c].AddEntry(self.histos[suff+'FitErr'+'qt'+c],self.groupedSystColors['Nominal'][1])
            if not impact : self.leg[suff+'FitErr'+'qt'+c].AddEntry(self.histos[suff+'FitErrPDF'+'qt'+c],self.groupedSystColors['LHEPdfWeightVars'][1])
            if not impact : self.leg[suff+'FitErr'+'qt'+c].AddEntry(self.histos[suff+'FitErrScale'+'qt'+c],self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitErr'+'qt'+c].SetLineWidth(0)
            self.leg[suff+'FitErr'+'qt'+c].SetFillStyle(0)
            # self.leg[suff+'FitErr'+'qt'+c].SetNColumns(2)
            self.leg[suff+'FitErr'+'qt'+c].Draw("Same")
            
            if 'unpol' in c :
                if impact : self.histos[suff+'impact'+'qt'+c+'stat'].Draw("samehist")
                if impact : self.leg[suff+'FitErr'+'qt'+c].AddEntry(self.histos[suff+'impact'+'qt'+c+'stat'],'Stat. only')
                self.histos[suff+'impact'+'qt'+c+'poiss'].Draw("same")
                self.histos[suff+'impact'+'qt'+c+'poiss'].SetLineWidth(3)
                self.histos[suff+'impact'+'qt'+c+'poiss'].SetLineColor(self.groupedSystColors['poiss'][0])
                self.leg[suff+'FitErr'+'qt'+c].AddEntry(self.histos[suff+'impact'+'qt'+c+'poiss'],self.groupedSystColors['poiss'][1])
                self.histos[suff+'FitErr'+'qt'+c].GetYaxis().SetRangeUser(0.0001,1)
                
                
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetRightMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetLeftMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetTopMargin(),cmslab)
            
            
            #--------------- Ai Y only relative syst. band canvas -------------------#
            self.canvas[suff+'FitErr'+'y'+c] = ROOT.TCanvas(suff+"_c_Y_{}_Err".format(c),suff+"_c_Y_{}_Err".format(c),1200,900)
            self.canvas[suff+'FitErr'+'y'+c].SetGridx()
            self.canvas[suff+'FitErr'+'y'+c].SetGridy()
            self.leg[suff+'FitErr'+'y'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            # self.histos[suff+'FitErr'+'y'+c].SetTitle(self.coeffDict[c][2]+"relative uncertainity, rapidity distribution,  q_{T} integrated"+", "+self.signDict[self.sign])
            self.histos[suff+'FitErr'+'y'+c].SetTitle('')
            self.canvas[suff+'FitErr'+'y'+c].cd()
            self.histos[suff+'FitErr'+'y'+c].GetXaxis().SetTitle('|Y_{W}|')
            self.histos[suff+'FitErr'+'y'+c].GetYaxis().SetTitle('Relative uncertainty on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
            # self.histos[suff+'FitErr'+'y'+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
            self.histos[suff+'FitErr'+'y'+c].GetYaxis().SetRangeUser(0.001,7000)
            self.canvas[suff+'FitErr'+'y'+c].SetLogy()
            self.histos[suff+'FitErr'+'y'+c].SetLineWidth(3)
            self.histos[suff+'FitErr'+'y'+c].SetStats(0)
            self.histos[suff+'FitErr'+'y'+c].Draw('hist')
            self.histos[suff+'FitErr'+'y'+c].SetLineColor(self.groupedSystColors['Nominal'][0])
            self.histos[suff+'FitErrPDF'+'y'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            self.histos[suff+'FitErrScale'+'y'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            self.histos[suff+'FitErrPDF'+'y'+c].SetLineWidth(3)
            self.histos[suff+'FitErrScale'+'y'+c].SetLineWidth(3)
            if not impact : self.histos[suff+'FitErrPDF'+'y'+c].Draw("hist same")
            if not impact : self.histos[suff+'FitErrScale'+'y'+c].Draw("hist same")
            self.leg[suff+'FitErr'+'y'+c].AddEntry(self.histos[suff+'FitErr'+'y'+c],self.groupedSystColors['Nominal'][1])
            if not impact : self.leg[suff+'FitErr'+'y'+c].AddEntry(self.histos[suff+'FitErrPDF'+'y'+c],self.groupedSystColors['LHEPdfWeightVars'][1])
            if not impact : self.leg[suff+'FitErr'+'y'+c].AddEntry(self.histos[suff+'FitErrScale'+'y'+c],self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitErr'+'y'+c].SetLineWidth(0)
            self.leg[suff+'FitErr'+'y'+c].SetFillStyle(0)
            # self.leg[suff+'FitErr'+'y'+c].SetNColumns(2)
            self.leg[suff+'FitErr'+'y'+c].Draw("Same")
            
            if 'unpol' in c :
                if impact : self.histos[suff+'impact'+'y'+c+'stat'].Draw("samehist")
                if impact : self.leg[suff+'FitErr'+'y'+c].AddEntry(self.histos[suff+'impact'+'y'+c+'stat'],'Stat. only')
                self.histos[suff+'impact'+'y'+c+'poiss'].Draw("same")
                self.histos[suff+'impact'+'y'+c+'poiss'].SetLineWidth(3)
                self.histos[suff+'impact'+'y'+c+'poiss'].SetLineColor(self.groupedSystColors['poiss'][0])
                self.leg[suff+'FitErr'+'y'+c].AddEntry(self.histos[suff+'impact'+'y'+c+'poiss'],self.groupedSystColors['poiss'][1])
                self.histos[suff+'FitErr'+'y'+c].GetYaxis().SetRangeUser(0.0001,1)
            
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetRightMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetLeftMargin(),1-0.8*self.canvas['ph'+suff+'FitAC'+'UNRqty'+c].GetTopMargin(),cmslab)
        
        
        #----------------prefit and postfit canvas -----------------------#
        if postfit : 
            for k in ['prefit','postfit'] :
                if k == 'prefit' : kName = 'pre-fit'
                if k == 'postfit' : kName = 'post-fit'
                for var,varInfo in self.postFitVars.items() :
                    if var=='' :
                        self.canvas[suff+k+var] = ROOT.TCanvas(suff+'_c_'+k+'_'+var,suff+'_c_'+k+'_'+var,2400,800)   
                        self.canvas['ph'+suff+k+var] = ROOT.TPad('ph'+suff+k+var, 'ph'+suff+k+var,0,0.3,1,1)
                        self.canvas['pr'+suff+k+var] = ROOT.TPad('pr'+suff+k+var, 'pr'+suff+k+var,0,0,1,0.3)  
                        self.canvas['pr'+suff+k].SetBottomMargin(0.5)
                        self.canvas['pr'+suff+k].SetRightMargin(0.02)
                        self.canvas['pr'+suff+k].SetLeftMargin(0.05)
                        self.canvas['ph'+suff+k].SetRightMargin(0.02)
                        self.canvas['ph'+suff+k].SetLeftMargin(0.05)
                   
                    else :
                        self.canvas[suff+k+var] = ROOT.TCanvas(suff+'_c_'+k+'_'+var,suff+'_c_'+k+'_'+var,800,700)
                        self.canvas['ph'+suff+k+var] = ROOT.TPad('ph'+suff+k+var, 'ph'+suff+k+var,0,0.2,1,1)
                        self.canvas['pr'+suff+k+var] = ROOT.TPad('pr'+suff+k+var, 'pr'+suff+k+var,0,0,1,0.2)
                        self.canvas['pr'+suff+k+var].SetBottomMargin(0.32)

                    self.canvas[suff+k+var].cd()
                    # self.canvas['ph'+suff+k+var] = ROOT.TPad('ph'+suff+k+var, 'ph'+suff+k+var,0,0.2,1,1)
                    # self.canvas['pr'+suff+k+var] = ROOT.TPad('pr'+suff+k+var, 'pr'+suff+k+var,0,0,1,0.2)
                    self.canvas['ph'+suff+k+var].SetBottomMargin(0.02)
                    self.canvas['ph'+suff+k+var].Draw()
                    self.canvas['pr'+suff+k+var].SetTopMargin(0.01)

                    self.canvas['pr'+suff+k+var].Draw()
                    
                    self.canvas['ph'+suff+k+var].cd()
                    self.canvas['ph'+suff+k+var].SetGridx()
                    self.canvas['ph'+suff+k+var].SetGridy()
                    # self.histos[suff+k+var+'stack'].SetTitle(varInfo[0]+", "+kName+", "+self.signDict[self.sign])
                    self.histos[suff+k+var+'stack'].SetTitle()
                    self.histos[suff+k+var+'stack'].Draw("HIST")
                    # self.histos[suff+k+var+'stack'].Draw("nostack")
                    self.histos[suff+k+'obs'+var].Draw("PE1SAME")
                    self.leg[suff+k+var].Draw("SAME")
                    
                    self.canvas['pr'+suff+k+var].cd()
                    self.canvas['pr'+suff+k+var].SetGridx()
                    self.canvas['pr'+suff+k+var].SetGridy()
                    self.histos[suff+k+var+'ratioBand'].Draw('E2')
                    if k=='prefit' and data!='' :
                         self.histos[suff+k+var+'ratio'].Draw("PE1same")
                    else :
                        self.histos[suff+k+var+'ratio'].Draw("histSAME")
                    self.histos[suff+k+var+'ratioBand'].SetStats(0)
                    self.histos[suff+k+var+'ratio'].SetStats(0)
                    
                    #aesthetic features
                    self.leg[suff+k+var].AddEntry(self.histos[suff+k+var+'ratioBand'], "Tot. Syst. Unc.")
                    self.leg[suff+k+var].SetFillStyle(0)
                    self.leg[suff+k+var].SetBorderSize(0)
                    self.leg[suff+k+var].SetNColumns(2)
                    
                    self.histos[suff+k+var+'stack'].SetMaximum(1.7*max(self.histos[suff+k+'obs'+var].GetMaximum(),self.histos[suff+k+var+'stack'].GetMaximum())) 
                    self.histos[suff+k+var+'stack'].GetYaxis().SetTitle(varInfo[1]+str(round(self.histos[suff+k+var+'stack'].GetXaxis().GetBinWidth(1),1))+varInfo[3])
                    self.histos[suff+k+var+'stack'].GetYaxis().SetTitleOffset(1.3)
                    self.histos[suff+k+var+'stack'].GetXaxis().SetTitle(varInfo[2])
                    self.histos[suff+k+var+'stack'].GetXaxis().SetTitleOffset(3)
                    self.histos[suff+k+var+'stack'].GetXaxis().SetLabelOffset(3)
                    
                    self.histos[suff+k+'obs'+var].SetMarkerStyle(20)
                    if var=='' : self.histos[suff+k+'obs'+var].SetMarkerSize(0.5)
                    self.histos[suff+k+'obs'+var].SetMarkerColor(self.processDict['obs'][0])
                    
                    self.histos[suff+k+var+'ratio'].SetMarkerStyle(1)#20
                    self.histos[suff+k+var+'ratio'].SetLineColor(ROOT.kBlack)
                    if var=='' and k=='prefit' and data!='' : 
                        self.histos[suff+k+var+'ratio'].SetLineWidth(1)
                    else :
                        self.histos[suff+k+var+'ratio'].SetLineWidth(2)
                    
                    self.histos[suff+k+var+'ratioBand'].SetLineWidth(0)
                    self.histos[suff+k+var+'ratioBand'].SetFillColor(ROOT.kOrange)#(ROOT.kCyan-4)
                    # self.histos[suff+k+var+'ratioBand'].SetFillStyle(1)
                    # self.histos[suff+k+var+'ratioBand'].SetFillStyle(3001)
                    self.histos[suff+k+var+'ratioBand'].SetMarkerStyle(1)
                    self.histos[suff+k+var+'ratioBand'].SetTitle("")
                    self.histos[suff+k+var+'ratioBand'].GetYaxis().SetTitle('Data/Pred.')
                    self.histos[suff+k+var+'ratioBand'].GetYaxis().SetTitleOffset(0.25)
                    self.histos[suff+k+var+'ratioBand'].GetYaxis().SetNdivisions(506)
                    self.histos[suff+k+var+'ratioBand'].SetTitleSize(0.15,'y')
                    self.histos[suff+k+var+'ratioBand'].SetLabelSize(0.12,'y')
                    self.histos[suff+k+var+'ratioBand'].GetXaxis().SetTitle(varInfo[2])
                    self.histos[suff+k+var+'ratioBand'].GetXaxis().SetTitleOffset(0.8)
                    self.histos[suff+k+var+'ratioBand'].SetTitleSize(0.18,'x')
                    self.histos[suff+k+var+'ratioBand'].SetLabelSize(0.15,'x')
                    # elif var.startswith('Mu1_pt'):
                    #     self.histos[suff+k+var+'ratioBand'].GetYaxis().SetRangeUser(0.88,1.18)
                    # else :
                    self.histos[suff+k+var+'ratioBand'].GetYaxis().SetRangeUser(0.98,1.019)

                    
                    if var =='' :
                        self.histos[suff+k+var+'ratioBand'].GetYaxis().SetRangeUser(0.85,1.149)
                        if k=='postfit' :
                            self.histos[suff+k+var+'ratioBand'].GetYaxis().SetRangeUser(0.98,1.019)
                        for ieta in range(1,self.histos[suff+k+p+'2D'].GetNbinsX()+1) :
                            for ipt in range(1,self.histos[suff+k+p+'2D'].GetNbinsY()+1) :
                                if ipt==30 and (ieta-1)%2==0 : #in the middle of the pt-bin for each even eta bin
                                    indexUNRetpt = (ieta-1)*(len(self.ptArr)-1)+(ipt-1)
                                    self.histos[suff+k+var+'ratioBand'].GetXaxis().SetBinLabel(indexUNRetpt+1,"#eta^{#mu}#in[%.1f,%.1f]" % (self.etaArr[ieta-1],self.etaArr[ieta]))
                                    self.histos[suff+k+var+'stack'].GetXaxis().SetBinLabel(indexUNRetpt+1,"#eta^{#mu}#in[%.1f,%.1f]" % (self.etaArr[ieta-1],self.etaArr[ieta]))
                        self.histos[suff+k+var+'ratioBand'].GetXaxis().SetTitle('fine binning: p_{T}^{#mu} 25#rightarrow 55 GeV')
                        # self.histos[suff+k+var+'ratioBand'].GetXaxis().SetTitleOffset(1.45)
                        self.histos[suff+k+var+'ratioBand'].LabelsOption("v")
                        self.histos[suff+k+var+'ratioBand'].GetXaxis().SetTitleOffset(2.8)
                        self.histos[suff+k+var+'ratioBand'].SetTitleSize(0.1,'x')
                        self.histos[suff+k+var+'ratioBand'].SetLabelSize(0.12,'x')
                        self.histos[suff+k+var+'ratioBand'].SetLabelOffset(0.005,'x')
                        self.histos[suff+k+var+'stack'].GetXaxis().SetTickLength(0)
                        self.histos[suff+k+var+'stack'].GetYaxis().SetTitleOffset(0.45)
                        self.histos[suff+k+var+'stack'].GetYaxis().SetTitleSize(0.05)
                        self.histos[suff+k+var+'stack'].GetYaxis().SetLabelSize(0.05)
                        self.histos[suff+k+var+'ratioBand'].GetYaxis().SetTitleOffset(0.45)
                        self.histos[suff+k+var+'ratioBand'].SetTitleSize(0.12,'y')
                        self.histos[suff+k+var+'ratioBand'].SetTitleOffset(0.17,'y')
                        self.histos[suff+k+var+'ratioBand'].SetLabelSize(0.1,'y')
                        self.histos[suff+k+var+'ratioBand'].GetXaxis().SetTickLength(0)
                        
                        self.canvas['ph'+suff+k+var].cd()
                        cmsLatex.SetNDC()
                        cmsLatex.SetTextFont(42)
                        cmsLatex.SetTextColor(ROOT.kBlack)
                        cmsLatex.SetTextAlign(31) 
                        cmsLatex.DrawLatex(1-self.canvas['ph'+suff+k+var].GetRightMargin(),1-0.8*self.canvas['ph'+suff+k+var].GetTopMargin(),lumilab)
                        cmsLatex.SetTextAlign(11) 
                        if 'prefit' in k and data!='': 
                            cmslabmod = cmslab.replace('Simulation','')
                        else :
                            cmslabmod = cmslab
                        cmsLatex.DrawLatex(self.canvas['ph'+suff+k+var].GetLeftMargin()+0.03,1-0.8*self.canvas['ph'+suff+k+var].GetTopMargin(),cmslabmod)


                        
                        
                        # for etapt in range(1,self.histos[suff+k+var+'ratioBand'].GetNbinsX()+1) :
                        #     print("out", (etapt-1)%len(self.ptArr))
                        #     if (etapt-1)%len(self.ptArr)==0 and etapt%2==0:
                        #         print(etapt)
                        #         print((etapt-1)%len(self.ptArr))
                        #         indPt = math.ceil(etapt/len(self.ptArr))
                        #         self.histos[suff+k+var+'ratioBand'].GetXaxis().SetBinLabel(etapt,"p_{T}^{#mu}#in[%.1f,%.1f]" % (self.ptArr[indPt],self.ptArr[indPt+1]))

            
        #---------------------------- Canvas correlation matrices ------------------------
        for mtx in ['corr','cov'] :
            if mtx == 'corr' : htitle = 'Correlation' 
            if mtx == 'cov' : htitle = 'Covariance'
            
            for c in self.coeffDict:
                self.canvas[suff+mtx+'Mat'+c] = ROOT.TCanvas(suff+'_c_'+mtx+'Mat_'+c,suff+'_'+mtx+'Mat_'+c,1200,900)
                self.canvas[suff+mtx+'Mat'+c].cd()
                # self.canvas[suff+mtx+'Mat'+c].SetGridx()
                # self.canvas[suff+mtx+'Mat'+c].SetGridy()
                self.histos[suff+mtx+'Mat'+c].Draw("colz")
                self.histos[suff+mtx+'Mat'+c].SetStats(0)
                # self.histos[suff+mtx+'Mat'+c].SetTitle('POI ' +htitle+' Matrix, '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                self.histos[suff+mtx+'Mat'+c].SetTitle('')
                self.histos[suff+mtx+'Mat'+c].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 60 GeV')
                self.histos[suff+mtx+'Mat'+c].GetXaxis().SetTitleOffset(1.48)
                self.histos[suff+mtx+'Mat'+c].GetXaxis().SetTickLength(0)
                self.histos[suff+mtx+'Mat'+c].GetYaxis().SetTickLength(0)
                self.histos[suff+mtx+'Mat'+c].GetYaxis().SetLabelOffset(0.001)
                
                self.histos[suff+mtx+'Mat'+c].GetZaxis().SetTitle(htitle)
                self.canvas[suff+mtx+'Mat'+c].SetRightMargin(0.15)
                self.canvas[suff+mtx+'Mat'+c].Update()
                palette = self.histos[suff+mtx+'Mat'+c].GetListOfFunctions().FindObject("palette")
                palette.SetX1NDC(0.875)
                
                if mtx == 'corr' : 
                    self.histos[suff+mtx+'Mat'+c].GetZaxis().SetCanExtend(1)
                    self.histos[suff+mtx+'Mat'+c].GetZaxis().SetRangeUser(-1,1)
                    
                for y in self.yArr[:-1] :
                    for q in self.qtArr[:-1] :
                        if self.qtArr.index(q)==0 :
                            indexUNR= self.yArr.index(float(y))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                            self.histos[suff+mtx+'Mat'+c].GetXaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+c].GetXaxis().SetBinLabel(indexUNR+1,"|Y_{W}|#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
                            self.histos[suff+mtx+'Mat'+c].GetXaxis().ChangeLabel(indexUNR+1,340,0.025)
                            self.histos[suff+mtx+'Mat'+c].GetXaxis().LabelsOption("d")
                            self.histos[suff+mtx+'Mat'+c].GetYaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+c].GetYaxis().SetBinLabel(indexUNR+1,"|Y_{W}|#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
                            self.histos[suff+mtx+'Mat'+c].GetYaxis().ChangeLabel(indexUNR+1,20,0.025)
                            self.histos[suff+mtx+'Mat'+c].GetYaxis().LabelsOption("u")

                cmsLatex.SetNDC()
                cmsLatex.SetTextFont(42)
                cmsLatex.SetTextColor(ROOT.kBlack)
                cmsLatex.SetTextAlign(31) 
                cmsLatex.DrawLatex(1-self.canvas[suff+mtx+'Mat'+c].GetRightMargin(),1-0.8*self.canvas[suff+mtx+'Mat'+c].GetTopMargin(),lumilab)
                cmsLatex.SetTextAlign(11) 
                cmsLatex.DrawLatex(self.canvas[suff+mtx+'Mat'+c].GetLeftMargin(),1-0.8*self.canvas[suff+mtx+'Mat'+c].GetTopMargin(),cmslab)
                # cmsLatex.SetTextColor(ROOT.kWhite)
                
                self.leg[suff+mtx+'Mat'+c] = ROOT.TPaveText(self.canvas[suff+mtx+'Mat'+c].GetLeftMargin()+0.02, 1-self.canvas[suff+mtx+'Mat'+c].GetTopMargin()-0.02,self.canvas[suff+mtx+'Mat'+c].GetLeftMargin()+0.2,1-self.canvas[suff+mtx+'Mat'+c].GetTopMargin()-0.1,"NDC")
                self.leg[suff+mtx+'Mat'+c].AddText(self.coeffDict[c][2]+", "+self.signDict[self.sign])
                self.leg[suff+mtx+'Mat'+c].SetFillColorAlpha(ROOT.kWhite,0.6)
                self.leg[suff+mtx+'Mat'+c].Draw("same")
                # cmsLatex.DrawLatex(self.canvas[suff+mtx+'Mat'+c].GetLeftMargin()+0.02,1-self.canvas[suff+mtx+'Mat'+c].GetTopMargin()-0.05,self.coeffDict[c][2]+", "+self.signDict[self.sign])
    
                
            for i in range(1, self.histos[suff+'FitAC'+self.coeffList[0]].GetNbinsX()+1): #loop over y bins
                self.canvas[suff+mtx+'Mat'+'y'+str(i)] = ROOT.TCanvas(suff+'_c_'+mtx+'Mat_'+'y'+str(i),suff+'_'+mtx+'Mat_'+'y'+str(i),1200,900)
                self.canvas[suff+mtx+'Mat'+'y'+str(i)].cd()
                # self.canvas[suff+mtx+'Mat'+'y'+str(i)].SetGridx()
                # self.canvas[suff+mtx+'Mat'+'y'+str(i)].SetGridy()
                self.histos[suff+mtx+'Mat'+'y'+str(i)].Draw("colz")
                self.histos[suff+mtx+'Mat'+'y'+str(i)].SetStats(0)
                self.histos[suff+mtx+'Mat'+'y'+str(i)].SetTitle('POI ' +htitle+' Matrix, '+str(self.yArr[i-1])+"<|Y_{W}|<"+str(self.yArr[i])+", "+self.signDict[self.sign])
                self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 60 GeV')
                self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().SetTitleOffset(1.45)
                self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().SetTickLength(0)
                self.histos[suff+mtx+'Mat'+'y'+str(i)].GetYaxis().SetTickLength(0)
                if mtx == 'corr' : 
                    self.histos[suff+mtx+'Mat'+'y'+str(i)].GetZaxis().SetCanExtend(1)
                    self.histos[suff+mtx+'Mat'+'y'+str(i)].GetZaxis().SetRangeUser(-1,1)
                for c in self.coeffArr[:-1] :
                    for q in self.qtArr[:-1] :
                        if self.qtArr.index(q)==0 :
                            indexUNR= self.coeffArr.index(float(c))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                            # if 'unpol' not in self.coeffList[self.coeffArr.index(c)] :    
                            #     coeffName = self.coeffList[self.coeffArr.index(c)]  
                            # else :
                            #     coeffName = 'unpol'
                                
                            coeffName = self.coeffDict[self.coeffList[self.coeffArr.index(c)]][2]
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().SetBinLabel(indexUNR+1,coeffName)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().ChangeLabel(indexUNR+1,340,0.04)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().LabelsOption("h")
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetYaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetYaxis().SetBinLabel(indexUNR+1,coeffName)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetYaxis().ChangeLabel(indexUNR+1,340,0.04)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetYaxis().LabelsOption("h")
            
            
            for j in range(1, self.histos[suff+'FitAC'+self.coeffList[0]].GetNbinsY()+1): #loop over qt bins
                self.canvas[suff+mtx+'Mat'+'qt'+str(j)] = ROOT.TCanvas(suff+'_c_'+mtx+'Mat_'+'qt'+str(j),suff+'_'+mtx+'Mat_'+'qt'+str(j),1200,900)
                self.canvas[suff+mtx+'Mat'+'qt'+str(j)].cd()
                self.canvas[suff+mtx+'Mat'+'qt'+str(j)].SetGridx()
                self.canvas[suff+mtx+'Mat'+'qt'+str(j)].SetGridy()
                self.histos[suff+mtx+'Mat'+'qt'+str(j)].Draw("colz")
                self.histos[suff+mtx+'Mat'+'qt'+str(j)].SetStats(0)
                self.histos[suff+mtx+'Mat'+'qt'+str(j)].SetTitle('POI ' +htitle+' Matrix, '+str(self.qtArr[j-1])+"<q_{T}^{W}<"+str(self.qtArr[j])+", "+self.signDict[self.sign])
                self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().SetTitle('fine binning: |Y_{W}| 0#rightarrow 2.4')
                self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().SetTitleOffset(1.45)
                if mtx == 'corr' : 
                    self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetZaxis().SetCanExtend(1)
                    self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetZaxis().SetRangeUser(-1,1)
                for c in self.coeffArr[:-1] :
                    for y in self.yArr[:-1] :
                        if self.yArr.index(y)==0 :
                            indexUNR= self.coeffArr.index(float(c))*(len(self.yArr)-1)+self.yArr.index(float(y))
                            # if 'unpol' not in self.coeffList[self.coeffArr.index(c)] :    
                            #     coeffName = self.coeffList[self.coeffArr.index(c)]  
                            # else :
                            #     coeffName = 'unpol'
                            coeffName = self.coeffDict[self.coeffList[self.coeffArr.index(c)]][2]
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().SetBinLabel(indexUNR+1,coeffName)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().ChangeLabel(indexUNR+1,340,0.04)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().LabelsOption("h")
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetYaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetYaxis().SetBinLabel(indexUNR+1,coeffName)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetYaxis().ChangeLabel(indexUNR+1,340,0.04)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetYaxis().LabelsOption("h")
            
            
            #---------------------- qt integrated-----------------------#
            self.canvas[suff+mtx+'Mat'+'Integrated'+'y'] = ROOT.TCanvas(suff+'_c_'+mtx+'Mat_'+'Integrated_'+'y',suff+'_'+mtx+'Mat_'+'Integrated_'+'y',1200,900)
            self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].cd()
            self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].SetGridx()
            self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].SetGridy()
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].Draw("colz")
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].SetStats(0)
            # self.histos[suff+mtx+'Mat'+'Integrated'+'y'].SetTitle('POI ' +htitle+' Matrix, q_{T} integrated'+", "+self.signDict[self.sign])
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].SetTitle('')
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().SetTitle('fine binning: |Y_{W}| 0#rightarrow 2.4')
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().SetTitleOffset(1.45)
                
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetZaxis().SetTitle(htitle)
            self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].SetRightMargin(0.15)
            self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].Update()
            try : 
                palette = self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetListOfFunctions().FindObject("palette")
                palette.SetX1NDC(0.875)
            except :
                print("no palette")
            
            if mtx == 'corr' : 
                self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetZaxis().SetCanExtend(1)
                self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetZaxis().SetRangeUser(-1,1)
            for c in self.coeffArr[:-1] :
                for y in self.yArr[:-1] :
                    if self.yArr.index(y)==0 :
                        indexUNR= self.coeffArr.index(float(c))*(len(self.yArr)-1)+self.yArr.index(float(y))
                        # if 'unpol' not in self.coeffList[self.coeffArr.index(c)] :    
                        #     coeffName = self.coeffList[self.coeffArr.index(c)]  
                        # else :
                        #     coeffName = 'unpol'
                        coeffName = self.coeffDict[self.coeffList[self.coeffArr.index(c)]][2]
                        halfrange = int(len(self.yArr)/2)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().SetBinLabel(indexUNR+1,coeffName)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().ChangeLabel(indexUNR+1,340,0.04)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().LabelsOption("h")
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetYaxis().SetNdivisions(-1)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetYaxis().SetBinLabel(indexUNR+1,coeffName)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetYaxis().LabelsOption("h")
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetYaxis().ChangeLabel(indexUNR+1,340,0.04)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetYaxis().SetLabelSize(0.06)
            
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].GetRightMargin(),1-0.8*self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].GetLeftMargin(),1-0.8*self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].GetTopMargin(),cmslab)
            # cmsLatex.SetTextColor(ROOT.kWhite)
            # cmsLatex.DrawLatex(self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].GetLeftMargin()+0.02,1-self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].GetTopMargin()-0.05,"q_{T}^{W} integrated, "+self.signDict[self.sign])
            
            self.leg[suff+mtx+'Mat'+'Integrated'+'y'] = ROOT.TPaveText(self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].GetLeftMargin()+0.02, 1-self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].GetTopMargin()-0.02,self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].GetLeftMargin()+0.3,1-self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].GetTopMargin()-0.1,"NDC")
            self.leg[suff+mtx+'Mat'+'Integrated'+'y'].AddText("q_{T}^{W} integrated, "+self.signDict[self.sign])
            self.leg[suff+mtx+'Mat'+'Integrated'+'y'].SetFillColorAlpha(ROOT.kWhite,0.6)
            self.leg[suff+mtx+'Mat'+'Integrated'+'y'].Draw("same")
    
            #------------------- y integrated-----------------------#
            self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'] = ROOT.TCanvas(suff+'_c_'+mtx+'Mat_'+'Integrated_'+'qt',suff+'_'+mtx+'Mat_'+'Integrated_'+'qt',1200,900)
            self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].cd()
            # self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].SetGridx()
            # self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].SetGridy()
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].Draw("colz")
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].SetStats(0)
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].SetTitle('POI ' +htitle+' Matrix, Y integrated'+", "+self.signDict[self.sign])
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].SetTitle('')
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 60 GeV')
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().SetTitleOffset(1.45)
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().SetTickLength(0)
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetYaxis().SetTickLength(0)
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetZaxis().SetTitle(htitle)
            
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetZaxis().SetTitle(htitle)
            self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].SetRightMargin(0.15)
            self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].Update()
            try :
                palette = self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetListOfFunctions().FindObject("palette")
                palette.SetX1NDC(0.875)
            except :
                print("no palette")
            
            if mtx == 'corr' : 
                self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetZaxis().SetCanExtend(1)
                self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetZaxis().SetRangeUser(-1,1)
            for c in self.coeffArr[:-1] :
                for q in self.qtArr[:-1] :
                    if self.qtArr.index(q)==0 :
                        indexUNR= self.coeffArr.index(float(c))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                        # if 'unpol' not in self.coeffList[self.coeffArr.index(c)] :    
                        #     coeffName = self.coeffList[self.coeffArr.index(c)]  
                        # else :
                        #     coeffName = 'unpol'
                        coeffName = self.coeffDict[self.coeffList[self.coeffArr.index(c)]][2]
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().SetBinLabel(indexUNR+1,coeffName)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().ChangeLabel(indexUNR+1,340,0.04)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().LabelsOption("h")
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetYaxis().SetBinLabel(indexUNR+1,coeffName)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetYaxis().SetNdivisions(-1)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetYaxis().LabelsOption("h")
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetYaxis().ChangeLabel(indexUNR+1,340,0.04)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetYaxis().SetLabelSize(0.06)
            
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].GetRightMargin(),1-0.8*self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].GetLeftMargin(),1-0.8*self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].GetTopMargin(),cmslab)
            # cmsLatex.SetTextColor(ROOT.kWhite)
            # cmsLatex.DrawLatex(self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].GetLeftMargin()+0.02,1-self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].GetTopMargin()-0.05,"|Y_{W}| integrated, "+self.signDict[self.sign])
            
            self.leg[suff+mtx+'Mat'+'Integrated'+'qt'] = ROOT.TPaveText(self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].GetLeftMargin()+0.02, 1-self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].GetTopMargin()-0.02,self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].GetLeftMargin()+0.3,1-self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].GetTopMargin()-0.1,"NDC")
            self.leg[suff+mtx+'Mat'+'Integrated'+'qt'].AddText("|Y_{W}| integrated, "+self.signDict[self.sign])
            self.leg[suff+mtx+'Mat'+'Integrated'+'qt'].SetFillColorAlpha(ROOT.kWhite,0.6)
            self.leg[suff+mtx+'Mat'+'Integrated'+'qt'].Draw("same")
        
        #---------------------------- Canvas impact plots ------------------------
        if impact :
            for c in self.coeffDict:

                self.histos[suff+'impact'+'UNRqty'+c]=ROOT.THStack(suff+'_impact_UNRqty'+c,suff+'_impact_UNRqty'+c)
                hempty = ROOT.TH1D(suff+'_impact_UNRqty'+c+'empty', suff+'_impact_UNRqty'+c+'empty', len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
                self.histos[suff+'impact'+'UNRqty'+c].SetHistogram(hempty)        
                self.leg[suff+'impact'+'UNRqty'+c] = ROOT.TLegend(0.3,0.15,0.85,0.35)
                if cleanNuisance :
                    # if 'unpol' in c :
                    self.leg[suff+'impact'+'UNRqty'+c] = ROOT.TLegend(0.25,0.72,0.6,0.88)
                    # else :
                    #     self.leg[suff+'impact'+'UNRqty'+c] = ROOT.TLegend(0.50,0.72,0.85,0.88)
                    self.leg[suff+'impact'+'UNRqty'+c].SetNColumns(2)
                    self.leg[suff+'impact'+'UNRqty'+c].SetFillStyle(0)
                    self.leg[suff+'impact'+'UNRqty'+c].SetLineWidth(0)
                else :
                    self.leg[suff+'impact'+'UNRqty'+c].SetNColumns(5)

                for nui in self.nuisanceDict: 
                    self.histos[suff+'impact'+'UNRqty'+c+nui].SetLineColor(self.nuisanceDict[nui][0])
                    self.histos[suff+'impact'+'UNRqty'+c+nui].SetMarkerColor(self.nuisanceDict[nui][0])
                    self.histos[suff+'impact'+'UNRqty'+c+nui].SetMarkerStyle(self.nuisanceDict[nui][2])
                    self.histos[suff+'impact'+'UNRqty'+c+nui].SetMarkerSize(1.)
                    self.histos[suff+'impact'+'UNRqty'+c+nui].SetLineWidth(2)
                    self.leg[suff+'impact'+'UNRqty'+c].AddEntry( self.histos[suff+'impact'+'UNRqty'+c+nui], self.nuisanceDict[nui][1])
                    self.histos[suff+'impact'+'UNRqty'+c].Add( self.histos[suff+'impact'+'UNRqty'+c+nui], 'lp')
                
                # self.histos[suff+'impact'+'UNRqty'+c].SetTitle('Impacts on '+self.coeffDict[c][2]+" unrolled q_{T}(Y)"+", "+self.signDict[self.sign])
                self.histos[suff+'impact'+'UNRqty'+c].SetTitle('')
                self.histos[suff+'impact'+'UNRqty'+c].GetXaxis().SetTitle('fine binning: |Y_{W}| 0#rightarrow 2.4')
                self.histos[suff+'impact'+'UNRqty'+c].GetXaxis().SetTitleOffset(1.47)
                if ('unpol' in c or (not self.anaKind['differential'])):  
                    self.histos[suff+'impact'+'UNRqty'+c].GetYaxis().SetTitle('Relative Impact on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                else :
                    self.histos[suff+'impact'+'UNRqty'+c].GetYaxis().SetTitle('Absolute Impact on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                hempty.SetStats(0)
                for q in self.qtArr[:-1] :
                    for y in self.yArr[:-1] :
                        indexUNRqty = self.qtArr.index(float(q))*(len(self.yArr)-1)+self.yArr.index(float(y))
                        if self.yArr.index(y)==0 :
                            self.histos[suff+'impact'+'UNRqty'+c].GetXaxis().SetNdivisions(-1)
                            self.histos[suff+'impact'+'UNRqty'+c].GetXaxis().SetBinLabel(indexUNRqty+1,"q_{T}^{W}#in[%.0f,%.0f]" % (q, self.qtArr[self.qtArr.index(q)+1]))
                            self.histos[suff+'impact'+'UNRqty'+c].GetXaxis().ChangeLabel(indexUNRqty+1,340,0.025)    
                            self.histos[suff+'impact'+'UNRqty'+c].GetXaxis().LabelsOption("d")
                self.canvas[suff+'impact'+'UNRqty'+c]=ROOT.TCanvas(suff+'_c_impact_UNRqty'+c,suff+'_c_impact_UNRqty'+c,1200,900)
                self.canvas[suff+'impact'+'UNRqty'+c].cd()
                self.canvas[suff+'impact'+'UNRqty'+c].SetLogy()
                self.canvas[suff+'impact'+'UNRqty'+c].SetGridy()
                self.canvas[suff+'impact'+'UNRqty'+c].SetGridx()
                self.histos[suff+'impact'+'UNRqty'+c].Draw('nostack')
                self.leg[suff+'impact'+'UNRqty'+c].Draw("same")
                if 'unpol' in c :  
                    self.histos[suff+'impact'+'UNRqty'+c].SetMinimum(0.002)
                    self.histos[suff+'impact'+'UNRqty'+c].SetMaximum(0.4)
                elif 'A4' in c:
                    self.histos[suff+'impact'+'UNRqty'+c].SetMinimum(0.0005)
                    self.histos[suff+'impact'+'UNRqty'+c].SetMaximum(3)
                else :
                    self.histos[suff+'impact'+'UNRqty'+c].SetMinimum(0.0005)
                    self.histos[suff+'impact'+'UNRqty'+c].SetMaximum(9)
                
                cmsLatex.SetNDC()
                cmsLatex.SetTextFont(42)
                cmsLatex.SetTextColor(ROOT.kBlack)
                cmsLatex.SetTextAlign(31) 
                cmsLatex.DrawLatex(1-self.canvas[suff+'impact'+'UNRqty'+c].GetRightMargin(),1-0.8*self.canvas[suff+'impact'+'UNRqty'+c].GetTopMargin(),lumilab)
                cmsLatex.SetTextAlign(11) 
                cmsLatex.DrawLatex(self.canvas[suff+'impact'+'UNRqty'+c].GetLeftMargin(),1-0.8*self.canvas[suff+'impact'+'UNRqty'+c].GetTopMargin(),cmslab)
                
                
                self.histos[suff+'impact'+'UNRyqt'+c]=ROOT.THStack(suff+'_impact_UNRyqt'+c,suff+'_impact_UNRyqt'+c)
                hempty = ROOT.TH1D(suff+'_impact_UNRyqt'+c+'empty', suff+'_impact_UNRyqt'+c+'empty', len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
                self.histos[suff+'impact'+'UNRyqt'+c].SetHistogram(hempty)        
                self.leg[suff+'impact'+'UNRyqt'+c] = ROOT.TLegend(0.3,0.15,0.85,0.35)
                if cleanNuisance :
                    #  if 'unpol' in c :
                     self.leg[suff+'impact'+'UNRyqt'+c] = ROOT.TLegend(0.25,0.72,0.6,0.88)
                    #  else : 
                        # self.leg[suff+'impact'+'UNRyqt'+c] = ROOT.TLegend(0.50,0.72,0.85,0.88)
                     self.leg[suff+'impact'+'UNRyqt'+c].SetNColumns(2)
                     self.leg[suff+'impact'+'UNRyqt'+c].SetFillStyle(0)
                     self.leg[suff+'impact'+'UNRyqt'+c].SetLineWidth(0)
                else :
                    self.leg[suff+'impact'+'UNRyqt'+c].SetNColumns(5)

                for nui in self.nuisanceDict: 
                    self.histos[suff+'impact'+'UNRyqt'+c+nui].SetLineColor(self.nuisanceDict[nui][0])
                    self.histos[suff+'impact'+'UNRyqt'+c+nui].SetMarkerColor(self.nuisanceDict[nui][0])
                    self.histos[suff+'impact'+'UNRyqt'+c+nui].SetMarkerStyle(self.nuisanceDict[nui][2])
                    self.histos[suff+'impact'+'UNRyqt'+c+nui].SetMarkerSize(1.)
                    self.histos[suff+'impact'+'UNRyqt'+c+nui].SetLineWidth(2)
                    self.leg[suff+'impact'+'UNRyqt'+c].AddEntry( self.histos[suff+'impact'+'UNRyqt'+c+nui], self.nuisanceDict[nui][1])
                    self.histos[suff+'impact'+'UNRyqt'+c].Add( self.histos[suff+'impact'+'UNRyqt'+c+nui], 'lp')
                
                # self.histos[suff+'impact'+'UNRyqt'+c].SetTitle('Impacts on '+self.coeffDict[c][2]+" unrolled Y(q_{T})"+", "+self.signDict[self.sign])
                self.histos[suff+'impact'+'UNRyqt'+c].SetTitle('')
                self.histos[suff+'impact'+'UNRyqt'+c].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 60 GeV')
                self.histos[suff+'impact'+'UNRyqt'+c].GetXaxis().SetTitleOffset(1.47)
                if ('unpol' in c or (not self.anaKind['differential'])):
                    self.histos[suff+'impact'+'UNRyqt'+c].GetYaxis().SetTitle('Relative Impact on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                else :
                    self.histos[suff+'impact'+'UNRyqt'+c].GetYaxis().SetTitle('Absolute Impact on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                self.histos[suff+'impact'+'UNRyqt'+c].GetXaxis().SetTickLength(0)
                hempty.SetStats(0)
                for y in self.yArr[:-1] :
                    for q in self.qtArr[:-1] :
                        indexUNRyqt = self.yArr.index(float(y))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                        if self.qtArr.index(q)==0 :
                            self.histos[suff+'impact'+'UNRyqt'+c].GetXaxis().SetNdivisions(-1)
                            self.histos[suff+'impact'+'UNRyqt'+c].GetXaxis().SetBinLabel(indexUNRyqt+1,"|Y_{W}|#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
                            self.histos[suff+'impact'+'UNRyqt'+c].GetXaxis().ChangeLabel(indexUNRyqt+1,340,0.025)    
                            self.histos[suff+'impact'+'UNRyqt'+c].GetXaxis().LabelsOption("d")
                self.canvas[suff+'impact'+'UNRyqt'+c]=ROOT.TCanvas(suff+'_c_impact_UNRyqt'+c,suff+'_c_impact_UNRyqt'+c,1200,900)
                self.canvas[suff+'impact'+'UNRyqt'+c].cd()
                self.canvas[suff+'impact'+'UNRyqt'+c].SetLogy()
                self.canvas[suff+'impact'+'UNRyqt'+c].SetGridy()
                # self.canvas[suff+'impact'+'UNRyqt'+c].SetGridx()
                self.histos[suff+'impact'+'UNRyqt'+c].Draw('nostack')
                if 'unpol' in c :  
                    self.histos[suff+'impact'+'UNRyqt'+c].SetMinimum(0.002)
                    self.histos[suff+'impact'+'UNRyqt'+c].SetMaximum(0.4)
                    bottomY4line = 0.002
                    topY4line = 0.4
                elif 'A4' in c :
                    self.histos[suff+'impact'+'UNRyqt'+c].SetMinimum(0.0005)
                    self.histos[suff+'impact'+'UNRyqt'+c].SetMaximum(3)
                    bottomY4line = 0.0005
                    topY4line = 3
                else :
                    self.histos[suff+'impact'+'UNRyqt'+c].SetMinimum(0.0005)
                    self.histos[suff+'impact'+'UNRyqt'+c].SetMaximum(9)
                    bottomY4line = 0.0005
                    topY4line = 9
                
                vImpLines = []
                # ROOT.gPad.Update()
                # bottomY4line = self.histos[suff+'impact'+'UNRyqt'+c].GetMinimum()
                # topY4line = self.histos[suff+'impact'+'UNRyqt'+c].GetMaximum()
                # self.canvas[suff+'impact'+'UNRyqt'+c].Range(0,0.001,1,200)
                
                # print(c, bottomY4line, topY4line)
                # if 'A4' in c or 'unpol' in c : 
                for yqt in self.unrolledYQt :
                        if yqt%(self.qtArr[-1])==0 :
                            if yqt==0 :continue
                            if yqt==len(self.unrolledYQt)-1 : continue 
                            # print(yqt-2,bottomY4line,yqt-2,topY4line)
                            vImpLines.append(ROOT.TLine(yqt-2,bottomY4line,yqt-2,topY4line))
                            vImpLines[-1].SetLineStyle(ROOT.kDotted)
                            # self.canvas[suff+'impact'+'UNRyqt'+c].cd()
                            vImpLines[-1].DrawLine(yqt-2,bottomY4line,yqt-2,topY4line)
                            # ll = ROOT.TLine(60,-2,60,2)
                            # ll.DrawLine(60,-2,60,2)
                            # ll.Draw("same")
                            # print(vImpLines)
                
                self.leg[suff+'impact'+'UNRyqt'+c].Draw("same")
                
                cmsLatex.SetNDC()
                cmsLatex.SetTextFont(42)
                cmsLatex.SetTextColor(ROOT.kBlack)
                cmsLatex.SetTextAlign(31) 
                cmsLatex.DrawLatex(1-self.canvas[suff+'impact'+'UNRyqt'+c].GetRightMargin(),1-0.8*self.canvas[suff+'impact'+'UNRyqt'+c].GetTopMargin(),lumilab)
                cmsLatex.SetTextAlign(11) 
                cmsLatex.DrawLatex(self.canvas[suff+'impact'+'UNRyqt'+c].GetLeftMargin(),1-0.8*self.canvas[suff+'impact'+'UNRyqt'+c].GetTopMargin(),cmslab)
                
                
                
                
                
                self.histos[suff+'impact'+'y'+c]=ROOT.THStack(suff+'_impact_y'+c,suff+'_impact_y'+c)
                hempty = ROOT.TH1D(suff+'_impact_y'+c+'empty', suff+'_impact_y'+c+'empty', len(self.yArr)-1, array('f',self.yArr))
                self.histos[suff+'impact'+'y'+c].SetHistogram(hempty)        
                self.leg[suff+'impact'+'y'+c] = ROOT.TLegend(0.3,0.15,0.85,0.35)
                if cleanNuisance :
                    # if 'unpol' in c :
                    self.leg[suff+'impact'+'y'+c] = ROOT.TLegend(0.25,0.72,0.6,0.88)
                    # else :
                        # self.leg[suff+'impact'+'y'+c] = ROOT.TLegend(0.50,0.72,0.85,0.88)
                    self.leg[suff+'impact'+'y'+c].SetNColumns(2)
                    self.leg[suff+'impact'+'y'+c].SetFillStyle(0)
                    self.leg[suff+'impact'+'y'+c].SetLineWidth(0)
                else :
                    self.leg[suff+'impact'+'y'+c].SetNColumns(5)
                for nui in self.nuisanceDict: 
                    self.histos[suff+'impact'+'y'+c+nui].SetLineColor(self.nuisanceDict[nui][0])
                    self.histos[suff+'impact'+'y'+c+nui].SetMarkerColor(self.nuisanceDict[nui][0])
                    self.histos[suff+'impact'+'y'+c+nui].SetMarkerStyle(self.nuisanceDict[nui][2])
                    self.histos[suff+'impact'+'y'+c+nui].SetMarkerSize(1.)
                    self.histos[suff+'impact'+'y'+c+nui].SetLineWidth(2)
                    self.leg[suff+'impact'+'y'+c].AddEntry( self.histos[suff+'impact'+'y'+c+nui], self.nuisanceDict[nui][1])
                    self.histos[suff+'impact'+'y'+c].Add( self.histos[suff+'impact'+'y'+c+nui], 'lp')
                # self.histos[suff+'impact'+'y'+c].SetTitle('Impacts on '+self.coeffDict[c][2]+" rapidity distribution, q_{T}^{W} integrated"+", "+self.signDict[self.sign])
                self.histos[suff+'impact'+'y'+c].SetTitle('')
                self.histos[suff+'impact'+'y'+c].GetXaxis().SetTitle('|Y_{W}|')
                # self.histos[suff+'impact'+'y'+c].GetXaxis().SetTitleOffset(1.45)
                if ('unpol' in c or (not self.anaKind['differential'])): 
                    self.histos[suff+'impact'+'y'+c].GetYaxis().SetTitle('Relative Impact on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                else :
                    self.histos[suff+'impact'+'y'+c].GetYaxis().SetTitle('Absolute Impact on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                hempty.SetStats(0)
                self.canvas[suff+'impact'+'y'+c]=ROOT.TCanvas(suff+'_c_impact_y'+c,suff+'_c_impact_y'+c,1200,900)
                self.canvas[suff+'impact'+'y'+c].cd()
                self.canvas[suff+'impact'+'y'+c].SetLogy()
                self.canvas[suff+'impact'+'y'+c].SetGridy()
                self.canvas[suff+'impact'+'y'+c].SetGridx()
                self.histos[suff+'impact'+'y'+c].Draw('nostack')
                self.leg[suff+'impact'+'y'+c].Draw("same")
                if 'unpol' in c :  
                    self.histos[suff+'impact'+'y'+c].SetMinimum(0.002)
                    self.histos[suff+'impact'+'y'+c].SetMaximum(0.4)
                elif 'A4' in c : 
                    self.histos[suff+'impact'+'y'+c].SetMinimum(0.0005)
                    self.histos[suff+'impact'+'y'+c].SetMaximum(3)
                else :
                   self.histos[suff+'impact'+'y'+c].SetMinimum(0.0005)
                   self.histos[suff+'impact'+'y'+c].SetMaximum(9)  
                                  
                cmsLatex.SetNDC()
                cmsLatex.SetTextFont(42)
                cmsLatex.SetTextColor(ROOT.kBlack)
                cmsLatex.SetTextAlign(31) 
                cmsLatex.DrawLatex(1-self.canvas[suff+'impact'+'y'+c].GetRightMargin(),1-0.8*self.canvas[suff+'impact'+'y'+c].GetTopMargin(),lumilab)
                cmsLatex.SetTextAlign(11) 
                cmsLatex.DrawLatex(self.canvas[suff+'impact'+'y'+c].GetLeftMargin(),1-0.8*self.canvas[suff+'impact'+'y'+c].GetTopMargin(),cmslab)
                
                
                self.histos[suff+'impact'+'qt'+c]=ROOT.THStack(suff+'_impact_qt'+c,suff+'_impact_qt'+c)
                hempty = ROOT.TH1D(suff+'_impact_qt'+c+'empty', suff+'_impact_qt'+c+'empty', len(self.qtArr)-1, array('f',self.qtArr))
                self.histos[suff+'impact'+'qt'+c].SetHistogram(hempty)           
                self.leg[suff+'impact'+'qt'+c] = ROOT.TLegend(0.3,0.15,0.85,0.35)
                if cleanNuisance :
                    # if 'unpol' in c :
                    self.leg[suff+'impact'+'qt'+c] = ROOT.TLegend(0.25,0.72,0.6,0.88)
                    # else :
                        # self.leg[suff+'impact'+'qt'+c] = ROOT.TLegend(0.50,0.72,0.85,0.88)
                    self.leg[suff+'impact'+'qt'+c].SetNColumns(2)
                    self.leg[suff+'impact'+'qt'+c].SetFillStyle(0)
                    self.leg[suff+'impact'+'qt'+c].SetLineWidth(0)
                else :
                    self.leg[suff+'impact'+'qt'+c].SetNColumns(5)
                for nui in self.nuisanceDict: 
                    self.histos[suff+'impact'+'qt'+c+nui].SetLineColor(self.nuisanceDict[nui][0])
                    self.histos[suff+'impact'+'qt'+c+nui].SetMarkerColor(self.nuisanceDict[nui][0])
                    self.histos[suff+'impact'+'qt'+c+nui].SetMarkerStyle(self.nuisanceDict[nui][2])
                    self.histos[suff+'impact'+'qt'+c+nui].SetMarkerSize(1.)
                    self.histos[suff+'impact'+'qt'+c+nui].SetLineWidth(2)
                    self.leg[suff+'impact'+'qt'+c].AddEntry( self.histos[suff+'impact'+'qt'+c+nui], self.nuisanceDict[nui][1])
                    self.histos[suff+'impact'+'qt'+c].Add( self.histos[suff+'impact'+'qt'+c+nui], 'lp')
                # self.histos[suff+'impact'+'qt'+c].SetTitle('Impacts on '+self.coeffDict[c][2]+" trasverse momentum distirbution, |Y_{W}| integrated"+", "+self.signDict[self.sign])
                self.histos[suff+'impact'+'qt'+c].SetTitle('')
                self.histos[suff+'impact'+'qt'+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
                # self.histos[suff+'impact'+'qt'+c].GetXaxis().SetTitleOffset(1.45)
                if ('unpol' in c or (not self.anaKind['differential'])):  
                    self.histos[suff+'impact'+'qt'+c].GetYaxis().SetTitle('Relative Impact on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                else :
                    self.histos[suff+'impact'+'qt'+c].GetYaxis().SetTitle('Absolute Impact on '+self.coeffDict[c][2]+", "+self.signDict[self.sign])
                hempty.SetStats(0)
                self.canvas[suff+'impact'+'qt'+c]=ROOT.TCanvas(suff+'_c_impact_qt'+c,suff+'_c_impact_qt'+c,1200,900)
                self.canvas[suff+'impact'+'qt'+c].cd()
                self.canvas[suff+'impact'+'qt'+c].SetLogy()
                self.canvas[suff+'impact'+'qt'+c].SetGridy()
                self.canvas[suff+'impact'+'qt'+c].SetGridx()
                self.histos[suff+'impact'+'qt'+c].Draw('nostack')
                self.leg[suff+'impact'+'qt'+c].Draw("same")
                if 'unpol' in c :  
                    self.histos[suff+'impact'+'qt'+c].SetMinimum(0.002)
                    self.histos[suff+'impact'+'qt'+c].SetMaximum(0.4)
                elif 'A4' in c :
                    self.histos[suff+'impact'+'qt'+c].SetMinimum(0.0005)
                    self.histos[suff+'impact'+'qt'+c].SetMaximum(3)
                else :
                   self.histos[suff+'impact'+'qt'+c].SetMinimum(0.0005)
                   self.histos[suff+'impact'+'qt'+c].SetMaximum(9) 
                
                cmsLatex.SetNDC()
                cmsLatex.SetTextFont(42)
                cmsLatex.SetTextColor(ROOT.kBlack)
                cmsLatex.SetTextAlign(31) 
                cmsLatex.DrawLatex(1-self.canvas[suff+'impact'+'qt'+c].GetRightMargin(),1-0.8*self.canvas[suff+'impact'+'qt'+c].GetTopMargin(),lumilab)
                cmsLatex.SetTextAlign(11) 
                cmsLatex.DrawLatex(self.canvas[suff+'impact'+'qt'+c].GetLeftMargin(),1-0.8*self.canvas[suff+'impact'+'qt'+c].GetTopMargin(),cmslab)
            
            
            self.histos[suff+'impact'+'mass']=ROOT.THStack(suff+'_impact_mass',suff+'_impact_mass')
            hempty = ROOT.TH1D(suff+'_impact_mass'+'empty', suff+'_impact_mass'+'empty', 1, 0,1)
            self.histos[suff+'impact'+'mass'].SetHistogram(hempty)        
            # self.leg[suff+'impact'+'mass'] = ROOT.TLegend(0.6,0.1,0.98,0.9)
            self.leg[suff+'impact'+'mass'] = ROOT.TLegend(0.59,0.2,0.98,0.8)
            self.leg[suff+'impact'+'mass'].SetLineStyle(2)
            # self.leg[suff+'impact'+'mass'].SetNColumns(1)
            for nui in self.nuisanceDict: 
                if nui=='mass' : continue 
                self.histos[suff+'impact'+'mass'+nui].SetLineColor(self.nuisanceDict[nui][0])
                self.histos[suff+'impact'+'mass'+nui].SetMarkerColor(self.nuisanceDict[nui][0])
                self.histos[suff+'impact'+'mass'+nui].SetMarkerStyle(self.nuisanceDict[nui][2])
                # self.histos[suff+'impact'+'mass'+nui].SetMarkerSize(1.)
                self.histos[suff+'impact'+'mass'+nui].SetLineWidth(4)
                # self.leg[suff+'impact'+'mass'].AddEntry( self.histos[suff+'impact'+'mass'+nui], self.nuisanceDict[nui][1]+'', self.histos[suff+'impact'+'mass'+nui].GetBinContent(1) )
                self.histos[suff+'impact'+'mass'].Add( self.histos[suff+'impact'+'mass'+nui], 'lp')
            
            #order the legend
            sortNuiMass = [] 
            for nui in self.nuisanceDict: 
                if nui=='mass' : continue 
                sortNuiMass.append([self.histos[suff+'impact'+'mass'+nui].GetBinContent(1), self.histos[suff+'impact'+'mass'+nui], self.nuisanceDict[nui][1]])
            sortNuiMass.sort(reverse=True)
            for nui in sortNuiMass: 
                if nui=='mass' : continue 
                val = '%.1f' % nui[0]
                self.leg[suff+'impact'+'mass'].AddEntry( nui[1], val+' MeV, '+nui[2])
                
            
            self.histos[suff+'impact'+'mass'].GetXaxis().SetBinLabel(1,'')
            self.histos[suff+'impact'+'mass'].GetXaxis().SetLabelSize(0.08)
            self.histos[suff+'impact'+'mass'].GetYaxis().SetMoreLogLabels(1)
            # self.histos[suff+'impact'+'mass'].SetTitle('mass uncertainty'+", "+self.signDict[self.sign])
            self.histos[suff+'impact'+'mass'].SetTitle('')
            self.histos[suff+'impact'+'mass'].GetXaxis().SetTitle('')
            # self.histos[suff+'impact'+'mass'].GetXaxis().SetTitleOffset(1.45)
            self.histos[suff+'impact'+'mass'].SetMinimum(0.025)
            self.histos[suff+'impact'+'mass'].SetMaximum(14)
            if self.sign=='plus' :
                self.histos[suff+'impact'+'mass'].GetYaxis().SetTitle('Impact on m_{W^{+}} [MeV]')
            else :
                self.histos[suff+'impact'+'mass'].GetYaxis().SetTitle('Impact on m_{W^{-}} [MeV]')
            hempty.SetStats(0)
            self.canvas[suff+'impact'+'mass']=ROOT.TCanvas(suff+'_c_impact_mass',suff+'_c_impact_mass',800,900)
            self.canvas[suff+'impact'+'mass'].cd()
            self.canvas[suff+'impact'+'mass'].SetLogy()
            self.canvas[suff+'impact'+'mass'].SetGridy()
            self.canvas[suff+'impact'+'mass'].SetGridx()
            self.canvas[suff+'impact'+'mass'].SetLeftMargin(0.13)
            self.canvas[suff+'impact'+'mass'].SetRightMargin(0.45)
            self.canvas[suff+'impact'+'mass'].SetBottomMargin(0.02)
            self.histos[suff+'impact'+'mass'].Draw('nostack E1')
            self.leg[suff+'impact'+'mass'].Draw("same")
            
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(0.98,1-0.8*self.canvas[suff+'impact'+'mass'].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas[suff+'impact'+'mass'].GetLeftMargin(),1-0.8*self.canvas[suff+'impact'+'mass'].GetTopMargin(),cmslab)
            
            
            # if cleanNuisance : #clone of the impacts for unpol only with the impacts + poissonian error aposteriori fit error
            poissColor = ROOT.kAzure+10
            apoColor = ROOT.kBlue-4
            c = 'unpolarizedxsec'
            
            for kkind in ['UNRqty','UNRyqt','y','qt'] :
            
                self.canvas[suff+'impact'+kkind+c+'poiss'] = self.canvas[suff+'impact'+kkind+c].Clone(self.canvas[suff+'impact'+kkind+c].GetName()+'_poiss')
                self.canvas[suff+'impact'+kkind+c+'poiss'].SetName(self.canvas[suff+'impact'+kkind+c+'poiss'].GetName()+'_poiss')
                self.canvas[suff+'impact'+kkind+c+'poiss'].cd()
                self.histos[suff+'impact'+kkind+c+'poiss'].SetLineColor(poissColor)
                self.histos[suff+'impact'+kkind+c+'poiss'].SetLineWidth(4)
                self.histos[suff+'impact'+kkind+c+'poiss'].Draw("Lhist same")
                self.leg[suff+'impact'+kkind+c+'poiss'] = self.canvas[suff+'impact'+kkind+c+'poiss'].GetPrimitive('TPave')
                self.leg[suff+'impact'+kkind+c+'poiss'].AddEntry(self.histos[suff+'impact'+kkind+c+'poiss'], 'Poiss. Unc.')
                if aposteriori!= '' :
                    self.histos[suff+'impact'+kkind+c+'apo'] = self.histos[suff+'impact'+kkind+c+'poiss'].Clone(self.histos[suff+'impact'+kkind+c+'poiss'].GetName().replace('poiss','apo'))
                    for xx in range(1,self.histos[suff+'impact'+kkind+c+'apo'].GetNbinsX()+1) :
                        self.histos[suff+'impact'+kkind+c+'apo'].SetBinContent(xx,self.histos[suff+'apoRatio'+kkind+c].GetBinError(xx))
                    self.histos[suff+'impact'+kkind+c+'apo'].SetLineColor(apoColor)
                    self.histos[suff+'impact'+kkind+c+'apo'].SetLineWidth(4)
                    self.histos[suff+'impact'+kkind+c+'apo'].Draw("Lhist same")
                    self.leg[suff+'impact'+kkind+c+'poiss'].AddEntry(self.histos[suff+'impact'+kkind+c+'apo'], 'Post-fit reg.')

            
            
        
        #------------------------canvas nuisance plots-------------------#  
        for nuiDict_key, nuiDict_val in self.NuiConstrDict.items() :
            self.histos[suff+'NuiConstr'+nuiDict_key].SetTitle("Nuisance parameters post-fit uncertainty "+ nuiDict_key+", "+self.signDict[self.sign])
            # self.histos[suff+'NuiConstr'+nuiDict_key].GetYaxis().SetTitle("#sigma(#theta-#theta^{0})")
            self.histos[suff+'NuiConstr'+nuiDict_key].GetYaxis().SetTitle("#sigma(#theta) post-fit, "+self.signDict[self.sign])
            self.histos[suff+'NuiConstr'+nuiDict_key].SetLineWidth(3)
            self.histos[suff+'NuiConstr'+nuiDict_key].SetStats(0)  
            if 'other' in nuiDict_key or 'all' in nuiDict_key :
                self.histos[suff+'NuiConstr'+nuiDict_key].SetLabelSize(0.05,'x')
            else :
                self.histos[suff+'NuiConstr'+nuiDict_key].SetLabelSize(0.015,'x')
            self.canvas[suff+'NuiConstr'+nuiDict_key]=ROOT.TCanvas(suff+'_c_NuiConstr_'+nuiDict_key,suff+'_c_NuiConstr_'+nuiDict_key,2400,800)
            self.canvas[suff+'NuiConstr'+nuiDict_key].cd()
            # self.canvas[suff+'NuiConstr'+nuiDict_key].SetGridx()
            self.canvas[suff+'NuiConstr'+nuiDict_key].SetGridy()
            self.histos[suff+'NuiConstr'+nuiDict_key].SetMinimum(-0.001)
            self.histos[suff+'NuiConstr'+nuiDict_key].GetYaxis().SetNdivisions(522)
            self.histos[suff+'NuiConstr'+nuiDict_key].GetYaxis().SetTickLength(0)
            self.histos[suff+'NuiConstr'+nuiDict_key].GetXaxis().SetTickLength(0)
            if 'all' in nuiDict_key :
                self.canvas[suff+'NuiConstr'+nuiDict_key].SetBottomMargin(0.22)
                self.canvas[suff+'NuiConstr'+nuiDict_key].SetLeftMargin(0.05)
                self.canvas[suff+'NuiConstr'+nuiDict_key].SetRightMargin(0.01)
                # self.histos[suff+'NuiConstr'+nuiDict_key].SetTitle("Nuisance parameters post-fit uncertainty, "+self.signDict[self.sign])
                self.histos[suff+'NuiConstr'+nuiDict_key].SetTitle("")
                self.histos[suff+'NuiConstr'+nuiDict_key].GetYaxis().SetTitleOffset(0.5)
                self.histos[suff+'NuiConstr'+nuiDict_key].GetYaxis().SetTitleSize(0.05)
            self.histos[suff+'NuiConstr'+nuiDict_key].Draw("EX01")
            for xx in range(1, self.histos[suff+'NuiConstr'+nuiDict_key].GetNbinsX()+1) :
                oldlab = self.histos[suff+'NuiConstr'+nuiDict_key].GetXaxis().GetBinLabel(xx)
                try : 
                    newlab = self.NuiConstrLabels[oldlab]
                except :
                    newlab = oldlab
                self.histos[suff+'NuiConstr'+nuiDict_key].GetXaxis().SetBinLabel(xx,newlab)
            
            cmsLatex.SetNDC()
            cmsLatex.SetTextFont(42)
            cmsLatex.SetTextColor(ROOT.kBlack)
            cmsLatex.SetTextAlign(31) 
            cmsLatex.DrawLatex(1-self.canvas[suff+'NuiConstr'+nuiDict_key].GetRightMargin(),1-0.8*self.canvas[suff+'NuiConstr'+nuiDict_key].GetTopMargin(),lumilab)
            cmsLatex.SetTextAlign(11) 
            cmsLatex.DrawLatex(self.canvas[suff+'NuiConstr'+nuiDict_key].GetLeftMargin(),1-0.8*self.canvas[suff+'NuiConstr'+nuiDict_key].GetTopMargin(),cmslab)

        
        if toy !='' :
            
            if self.toyPullFit :
                #------------fit the toy pulls --------#
                for c in self.coeffDict:       
                    self.histos[suff+'FitAC'+c+'toyPull'+'Fit'] = ROOT.TH2D(suff+'FitAC{c}_toyPull_Fit'.format(c=c), suff+'FitAC{c}_toyPull_Fit'.format(c=c), len(self.yArr)-1, array('f',self.yArr), len(self.qtArr)-1, array('f',self.qtArr))
                    self.histos[suff+'FitACqt'+c+'toyPull'+'Fit'] = ROOT.TH1D(suff+'FitACqt{c}_toyPull_Fit'.format(c=c), suff+'FitACqt{c}_toyPull_Fit'.format(c=c), len(self.qtArr)-1, array('f',self.qtArr))
                    self.histos[suff+'FitACy'+c+'toyPull'+'Fit'] = ROOT.TH1D(suff+'FitACy{c}_toyPull_Fit'.format(c=c), suff+'FitACy{c}_toyPull_Fit'.format(c=c), len(self.yArr)-1, array('f',self.yArr))
    
                    for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                        for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                            self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)].Fit('gaus','Q')
                            func = self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)].GetFunction('gaus')
                            # print(func.GetChisquare(), func.GetNDF())
                            self.histos[suff+'FitAC'+c+'toyPull'+'Fit'].SetBinContent(i,j,func.GetParameter(1))
                            self.histos[suff+'FitAC'+c+'toyPull'+'Fit'].SetBinError(i,j,func.GetParameter(2))
                    
                    for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                        self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)].Fit('gaus','Q')
                        func = self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)].GetFunction('gaus')
                        # print(func.GetChisquare())
                        self.histos[suff+'FitACqt'+c+'toyPull'+'Fit'].SetBinContent(j,func.GetParameter(1))
                        self.histos[suff+'FitACqt'+c+'toyPull'+'Fit'].SetBinError(j,func.GetParameter(2))
                    
                    for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                        self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)].Fit('gaus','Q')
                        func = self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)].GetFunction('gaus')
                        # print(func.GetChisquare())
                        self.histos[suff+'FitACy'+c+'toyPull'+'Fit'].SetBinContent(i,func.GetParameter(1))
                        self.histos[suff+'FitACy'+c+'toyPull'+'Fit'].SetBinError(i,func.GetParameter(2))
                    
                    self.histos[suff+'FitAC'+c+'toyPull'] = self.histos[suff+'FitAC'+c+'toyPull'+'Fit'] 
                    self.histos[suff+'FitACqt'+c+'toyPull'] = self.histos[suff+'FitACqt'+c+'toyPull'+'Fit'] 
                    self.histos[suff+'FitACy'+c+'toyPull'] = self.histos[suff+'FitACy'+c+'toyPull'+'Fit']            
            
            
            #---------------------------- Canvas pulls - unrolled: qt large, y small canvas (only unrolled produced for pulls) ------------------------------------
            for c in self.coeffDict:
                self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'] = ROOT.TH1F(suff+'_coeff_toyPull_UNRqty'+c,suff+'_coeff_toyPull_UNRqty'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
                for q in self.qtArr[:-1] :
                    for y in self.yArr[:-1] :
                        indexUNRqty = self.qtArr.index(float(q))*(len(self.yArr)-1)+self.yArr.index(float(y))
                        self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].SetBinContent(indexUNRqty+1,self.histos[suff+'FitAC'+c+'toyPull'].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].SetBinError(indexUNRqty+1,self.histos[suff+'FitAC'+c+'toyPull'].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        if self.yArr.index(y)==0 :
                            self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].GetXaxis().SetNdivisions(-1)
                            self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].GetXaxis().SetBinLabel(indexUNRqty+1,"q_{T}^{W}#in[%.0f,%.0f]" % (q, self.qtArr[self.qtArr.index(q)+1]))
                            self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].GetXaxis().ChangeLabel(indexUNRqty+1,340,0.02)    
                            self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].GetXaxis().LabelsOption("d")  
                self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].SetTitle("Toy pulls mean and RMS  "+c+", unrolled q_{T}(Y)"+", "+self.signDict[self.sign])
                self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].GetXaxis().SetTitle('fine binning: |Y_{W}| 0#rightarrow 2.4')
                self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].GetXaxis().SetTitleOffset(1.45)
                # self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].GetYaxis().SetTitle("(toy_{i}-gen)/#sigma_{toy_{i}}")
                self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].GetYaxis().SetTitle("#mu_{pull} #pm #sigma_{pull}")
                self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].SetLineWidth(3)
                self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].SetStats(0)
                self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].SetLineColor(self.groupedSystColors['Nominal'][0])

                self.canvas[suff+'FitAC'+'UNRqty'+c+'toyPull'] = ROOT.TCanvas(suff+'_c_toyPull_UNRqty'+c,suff+'_c_toyPull_UNRqty'+c,1200,900)
                self.canvas[suff+'FitAC'+'UNRqty'+c+'toyPull'].cd()
                self.canvas[suff+'FitAC'+'UNRqty'+c+'toyPull'].SetGridx()
                self.canvas[suff+'FitAC'+'UNRqty'+c+'toyPull'].SetGridy()
                self.histos[suff+'FitAC'+'UNRqty'+c+'toyPull'].DrawCopy()
                
                #---------------------------- Canvas pulls - unrolled: y large, qt small canvas  (only unrolled produced for pulls) ------------------------------------
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'] = ROOT.TH1F(suff+'_coeff_toyPull_UNRyqt'+c,suff+'_coeff_toyPull_UNRyqt'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
                for y in self.yArr[:-1] :
                    for q in self.qtArr[:-1] :
                        indexUNRyqt = self.yArr.index(float(y))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                        self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitAC'+c+'toyPull'].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].SetBinError(indexUNRyqt+1,self.histos[suff+'FitAC'+c+'toyPull'].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                        if self.qtArr.index(q)==0 :
                            self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].GetXaxis().SetNdivisions(-1)
                            self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].GetXaxis().SetBinLabel(indexUNRyqt+1,"|Y_{W}|#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
                            self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].GetXaxis().ChangeLabel(indexUNRyqt+1,340,0.02)
                            self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].GetXaxis().LabelsOption("d")
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].SetTitle("Toy pulls mean and RMS  "+c+", unrolled Y(q_{T})"+", "+self.signDict[self.sign])
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 60 GeV')
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].GetXaxis().SetTitleOffset(1.45)
                # self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].GetYaxis().SetTitle("(toy_{i}-gen)/#sigma_{toy_{i}}")
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].GetYaxis().SetTitle("#mu_{pull} #pm #sigma_{pull}")
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].SetLineWidth(3)
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].SetStats(0)
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].SetLineColor(self.groupedSystColors['Nominal'][0])

                self.canvas[suff+'FitAC'+'UNRyqt'+c+'toyPull'] = ROOT.TCanvas(suff+'_c_toyPull_UNRyqt'+c,suff+'_c_toyPull_UNRyqt'+c,1200,900)
                self.canvas[suff+'FitAC'+'UNRyqt'+c+'toyPull'].cd()
                self.canvas[suff+'FitAC'+'UNRyqt'+c+'toyPull'].SetGridx()
                self.canvas[suff+'FitAC'+'UNRyqt'+c+'toyPull'].SetGridy()
                self.histos[suff+'FitAC'+'UNRyqt'+c+'toyPull'].DrawCopy()
                
            
                #--------------- Canvas pulls QT only canvas -------------------#
                self.canvas[suff+'FitAC'+'qt'+c+'toyPull'] = ROOT.TCanvas(suff+"_c_QT_{}_toyPull".format(c),suff+"_c_QT_{}_toyPull".format(c),1200,900)
                self.canvas[suff+'FitAC'+'qt'+c+'toyPull'].SetGridx()
                self.canvas[suff+'FitAC'+'qt'+c+'toyPull'].SetGridy()
                self.histos[suff+'FitAC'+'qt'+c+'toyPull'].SetTitle("Toy pulls mean and RMS  "+c+", Y integrated"+", "+self.signDict[self.sign])
                self.canvas[suff+'FitAC'+'qt'+c+'toyPull'].cd()
                self.histos[suff+'FitAC'+'qt'+c+'toyPull'].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
                # self.histos[suff+'FitAC'+'qt'+c+'toyPull'].GetYaxis().SetTitle("(toy_{i}-gen)/#sigma_{toy_{i}}")
                self.histos[suff+'FitAC'+'qt'+c+'toyPull'].GetYaxis().SetTitle("#mu_{pull} #pm #sigma_{pull}")
                self.histos[suff+'FitAC'+'qt'+c+'toyPull'].SetLineWidth(3)
                self.histos[suff+'FitAC'+'qt'+c+'toyPull'].SetStats(0)
                self.histos[suff+'FitAC'+'qt'+c+'toyPull'].DrawCopy()
                self.histos[suff+'FitAC'+'qt'+c+'toyPull'].SetLineColor(self.groupedSystColors['Nominal'][0])
                    
                    
                #--------------- Canvas pulls Y only canvas -------------------#
                self.canvas[suff+'FitAC'+'y'+c+'toyPull'] = ROOT.TCanvas(suff+"_c_Y_{}_toyPull".format(c),suff+"_c_Y_{}_toyPull".format(c),1200,900)
                self.canvas[suff+'FitAC'+'y'+c+'toyPull'].SetGridx()
                self.canvas[suff+'FitAC'+'y'+c+'toyPull'].SetGridy()
                self.histos[suff+'FitAC'+'y'+c+'toyPull'].SetTitle("Toy pulls mean and RMS  "+c+",  q_{T} integrated"+", "+self.signDict[self.sign])
                self.canvas[suff+'FitAC'+'y'+c+'toyPull'].cd()
                self.histos[suff+'FitAC'+'y'+c+'toyPull'].GetXaxis().SetTitle('|Y_{W}|')
                self.histos[suff+'FitAC'+'y'+c+'toyPull'].GetYaxis().SetTitle("(toy_{i}-gen)/#sigma_{toy_{i}}")
                self.histos[suff+'FitAC'+'y'+c+'toyPull'].GetYaxis().SetTitle("#mu_{pull} #pm #sigma_{pull}")
                self.histos[suff+'FitAC'+'y'+c+'toyPull'].SetLineWidth(3)
                self.histos[suff+'FitAC'+'y'+c+'toyPull'].SetStats(0)
                self.histos[suff+'FitAC'+'y'+c+'toyPull'].DrawCopy()
                self.histos[suff+'FitAC'+'y'+c+'toyPull'].SetLineColor(self.groupedSystColors['Nominal'][0])
            
            
            #build up/down band with extra histos
            for kindpull in ['UNRqty','UNRyqt','qt','y'] :
                for c in self.coeffList : 
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'+'Up'] = self.histos[suff+'FitAC'+kindpull+c+'toyPull'].Clone(self.histos[suff+'FitAC'+kindpull+c+'toyPull'].GetName()+'Up')
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'+'Down'] = self.histos[suff+'FitAC'+kindpull+c+'toyPull'].Clone(self.histos[suff+'FitAC'+kindpull+c+'toyPull'].GetName()+'Down')
                    for xp in range(1, self.histos[suff+'FitAC'+kindpull+c+'toyPull'+'Up'].GetNbinsX()+1) :
                        self.histos[suff+'FitAC'+kindpull+c+'toyPull'+'Up'].SetBinContent(xp, self.histos[suff+'FitAC'+kindpull+c+'toyPull'].GetBinContent(xp)+self.histos[suff+'FitAC'+kindpull+c+'toyPull'].GetBinError(xp))
                        self.histos[suff+'FitAC'+kindpull+c+'toyPull'+'Down'].SetBinContent(xp, self.histos[suff+'FitAC'+kindpull+c+'toyPull'].GetBinContent(xp)-self.histos[suff+'FitAC'+kindpull+c+'toyPull'].GetBinError(xp))
                        self.histos[suff+'FitAC'+kindpull+c+'toyPull'+'Up'].SetBinError(xp, 0)
                        self.histos[suff+'FitAC'+kindpull+c+'toyPull'+'Down'].SetBinError(xp, 0)
            
            # #pulls with all coeff canvas 
            for kindpull in ['UNRqty','UNRyqt','qt','y'] :
                self.canvas[suff+'FitAC'+kindpull+'toyPull'] = ROOT.TCanvas(suff+'_c_toyPull_'+kindpull,suff+'_c_toyPull_'+kindpull,1200,900)
                self.canvas[suff+'FitAC'+kindpull+'toyPull'].cd()
                self.canvas[suff+'FitAC'+kindpull+'toyPull'].SetGridx()
                self.canvas[suff+'FitAC'+kindpull+'toyPull'].SetGridy()
                self.leg[suff+'FitAC'+kindpull+'toyPull'] = ROOT.TLegend(0.15, 0.72, 0.4, 0.9) 
                self.leg[suff+'FitAC'+kindpull+'toyPull'].SetFillStyle(0)
                self.leg[suff+'FitAC'+kindpull+'toyPull'].SetBorderSize(0)
                self.leg[suff+'FitAC'+kindpull+'toyPull'].SetNColumns(2)
                sameFlag = ''
                for c in self.coeffList : 
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetMaximum(2)
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetMinimum(-2)
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetLineColor(self.coeffDict[c][3])
                    # self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetFillColorAlpha(self.coeffDict[c][3],0.1)
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetFillColor(self.coeffDict[c][3])
                    # self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetFillStyle(3002+self.coeffList.index(c))
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetFillStyle(3003)
                    # self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetFillStyle(1001)
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetTitle(self.histos[suff+'FitAC'+kindpull+c+'toyPull'].GetTitle().replace(" "+c,''))
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetTitle('')
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].DrawCopy(sameFlag+'E2')
                    sameFlag = 'same'
                    self.leg[suff+'FitAC'+kindpull+'toyPull'].AddEntry(self.histos[suff+'FitAC'+kindpull+c+'toyPull'], self.coeffDict[c][2])
                for c in self.coeffList : 
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetMarkerColor(self.coeffDict[c][3])
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetMarkerStyle(20)
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetMarkerSize(0.7)
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].SetLineWidth(0)
                    self.histos[suff+'FitAC'+kindpull+c+'toyPull'].Draw("same")
                    
                    for ud in ['Up', 'Down'] :
                        self.histos[suff+'FitAC'+kindpull+c+'toyPull'+ud].SetLineColor(self.coeffDict[c][3])
                        self.histos[suff+'FitAC'+kindpull+c+'toyPull'+ud].SetLineWidth(2)
                        self.histos[suff+'FitAC'+kindpull+c+'toyPull'+ud].Draw("same hist")
                    
    
                self.leg[suff+'FitAC'+kindpull+'toyPull'].Draw("same")
                # self.canvas[suff+'FitAC'+kindpull+'toyPull'].SaveAs('TEST'+suff+'FitAC'+kindpull+'toyPull'+'.png')
                
                cmsLatex.SetNDC()
                cmsLatex.SetTextFont(42)
                cmsLatex.SetTextColor(ROOT.kBlack)
                cmsLatex.SetTextAlign(31) 
                cmsLatex.DrawLatex(1-self.canvas[suff+'FitAC'+kindpull+'toyPull'].GetRightMargin(),1-0.8*self.canvas[suff+'FitAC'+kindpull+'toyPull'].GetTopMargin(),lumilab)
                cmsLatex.SetTextAlign(11) 
                cmsLatex.DrawLatex(self.canvas[suff+'FitAC'+kindpull+'toyPull'].GetLeftMargin(),1-0.8*self.canvas[suff+'FitAC'+kindpull+'toyPull'].GetTopMargin(),cmslab)
                cmsLatex.SetTextAlign(31) 
                if kindpull == 'qt' :
                    extralab = '|Y_{W}| integrated, '+self.signDict[self.sign]
                elif kindpull == 'y' :
                    extralab = 'q_{T}^{W} integrated, '+self.signDict[self.sign]
                else :
                    extralab = self.signDict[self.sign]
                cmsLatex.DrawLatex(1-self.canvas[suff+'FitAC'+kindpull+'toyPull'].GetRightMargin()-0.02,1-self.canvas[suff+'FitAC'+kindpull+'toyPull'].GetTopMargin()-0.05,extralab)

            
            
            
            
            
            
            
        
        
    
    
    def GenRecoComparison(self,suffGen, suffReco) :
        
        for c in self.coeffDict:
            
            # ------------------ Ai plots -------------------------------
            
            for i in range(1, self.histos[suffGen+'FitAC'+c].GetNbinsX()+1):
                self.canvas['comp'+'FitAC'+'qt'+str(i)+c] = self.canvas[suffGen+'FitAC'+'qt'+str(i)+c].Clone(self.canvas[suffGen+'FitAC'+'qt'+str(i)+c].GetName().replace(suffGen,'comp'))
                self.canvas['comp'+'FitAC'+'qt'+str(i)+c].SetTitle(self.canvas[suffGen+'FitAC'+'qt'+str(i)+c].GetTitle().replace(suffGen,'comp'))
                self.canvas['comp'+'FitAC'+'qt'+str(i)+c].SetName(self.canvas[suffGen+'FitAC'+'qt'+str(i)+c].GetName().replace(suffGen,'comp'))
                self.canvas['comp'+'FitAC'+'qt'+str(i)+c].cd()
                self.histos[suffReco+'FitAC'+'qt'+str(i)+c+'4comp'] = self.histos[suffReco+'FitAC'+'qt'+str(i)+c].Clone(suffReco+'FitAC'+'qt'+str(i)+c+'4comp')
                self.histos[suffReco+'FitAC'+'qt'+str(i)+c+'4comp'].SetLineColor(ROOT.kGreen-3)
                self.histos[suffReco+'FitAC'+'qt'+str(i)+c+'4comp'].Draw("same")
                self.histos[suffGen+'FitAC'+'qt'+str(i)+c].DrawCopy("same") #to have foreground
                
                self.leg['comp'+'FitAC'+'qt'+str(i)+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
                self.leg['comp'+'FitAC'+'qt'+str(i)+c].AddEntry(self.histos[suffReco+'FitAC'+'qt'+str(i)+c+'4comp'],'Reco Asimov')
                self.leg['comp'+'FitAC'+'qt'+str(i)+c].AddEntry(self.histos[suffGen+'FitAC'+'qt'+str(i)+c],'Gen Asimov')
                self.leg['comp'+'FitAC'+'qt'+str(i)+c].AddEntry(self.histos[suffGen+'FitBand'+'qt'+str(i)+c],'Gen Syst')
                self.leg['comp'+'FitAC'+'qt'+str(i)+c].Draw("same")
                
            for j in range(1, self.histos[suffGen+'FitErr'+c].GetNbinsY()+1):
                self.canvas['comp'+'FitAC'+'y'+str(j)+c] = self.canvas[suffGen+'FitAC'+'y'+str(j)+c].Clone(self.canvas[suffGen+'FitAC'+'y'+str(j)+c].GetName().replace(suffGen,'comp'))
                self.canvas['comp'+'FitAC'+'y'+str(j)+c].SetTitle(self.canvas[suffGen+'FitAC'+'y'+str(j)+c].GetTitle().replace(suffGen,'comp'))
                self.canvas['comp'+'FitAC'+'y'+str(j)+c].SetName(self.canvas[suffGen+'FitAC'+'y'+str(j)+c].GetName().replace(suffGen,'comp'))
                self.canvas['comp'+'FitAC'+'y'+str(j)+c].cd()
                self.histos[suffReco+'FitAC'+'y'+str(j)+c+'4comp'] = self.histos[suffReco+'FitAC'+'y'+str(j)+c].Clone(suffReco+'FitAC'+'y'+str(j)+c+'4comp')
                self.histos[suffReco+'FitAC'+'y'+str(j)+c+'4comp'].SetLineColor(ROOT.kGreen-3)
                self.histos[suffReco+'FitAC'+'y'+str(j)+c+'4comp'].Draw("same")
                self.histos[suffGen+'FitAC'+'y'+str(j)+c].DrawCopy("same") #to have foreground
                
                self.leg['comp'+'FitAC'+'y'+str(j)+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
                self.leg['comp'+'FitAC'+'y'+str(j)+c].AddEntry(self.histos[suffReco+'FitAC'+'y'+str(j)+c+'4comp'],'Reco Asimov')
                self.leg['comp'+'FitAC'+'y'+str(j)+c].AddEntry(self.histos[suffGen+'FitAC'+'y'+str(j)+c],'Gen Asimov')
                self.leg['comp'+'FitAC'+'y'+str(j)+c].AddEntry(self.histos[suffGen+'FitBand'+'y'+str(j)+c],'Gen Syst')
                self.leg['comp'+'FitAC'+'y'+str(j)+c].Draw("same")
            
            self.canvas['comp'+'FitAC'+'UNRqty'+c] = self.canvas[suffGen+'FitAC'+'UNRqty'+c].Clone(self.canvas[suffGen+'FitAC'+'UNRqty'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'UNRqty'+c].SetTitle(self.canvas[suffGen+'FitAC'+'UNRqty'+c].GetTitle().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'UNRqty'+c].SetName(self.canvas[suffGen+'FitAC'+'UNRqty'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'UNRqty'+c].cd()
            self.histos[suffReco+'FitAC'+'UNRqty'+c+'4comp'] = self.histos[suffReco+'FitAC'+'UNRqty'+c].Clone(suffReco+'FitAC'+'y'+str(j)+c+'4comp')
            self.histos[suffReco+'FitAC'+'UNRqty'+c+'4comp'].SetLineColor(ROOT.kGreen-3)
            self.histos[suffReco+'FitAC'+'UNRqty'+c+'4comp'].Draw("same")
            self.histos[suffGen+'FitAC'+'UNRqty'+c].DrawCopy("same") #to have foreground
            self.leg['comp'+'FitAC'+'UNRqty'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            self.leg['comp'+'FitAC'+'UNRqty'+c].AddEntry(self.histos[suffReco+'FitAC'+'UNRqty'+c+'4comp'],'Reco Asimov')
            self.leg['comp'+'FitAC'+'UNRqty'+c].AddEntry(self.histos[suffGen+'FitAC'+'UNRqty'+c],'Gen Asimov')
            self.leg['comp'+'FitAC'+'UNRqty'+c].AddEntry(self.histos[suffGen+'FitBand'+'UNRqty'+c],'Gen Syst')
            self.leg['comp'+'FitAC'+'UNRqty'+c].Draw("same")
            
            self.canvas['comp'+'FitAC'+'UNRyqt'+c] = self.canvas[suffGen+'FitAC'+'UNRyqt'+c].Clone(self.canvas[suffGen+'FitAC'+'UNRyqt'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'UNRyqt'+c].SetTitle(self.canvas[suffGen+'FitAC'+'UNRyqt'+c].GetTitle().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'UNRyqt'+c].SetName(self.canvas[suffGen+'FitAC'+'UNRyqt'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'UNRyqt'+c].cd()
            self.histos[suffReco+'FitAC'+'UNRyqt'+c+'4comp'] = self.histos[suffReco+'FitAC'+'UNRyqt'+c].Clone(suffReco+'FitAC'+'y'+str(j)+c+'4comp')
            self.histos[suffReco+'FitAC'+'UNRyqt'+c+'4comp'].SetLineColor(ROOT.kGreen-3)
            self.histos[suffReco+'FitAC'+'UNRyqt'+c+'4comp'].Draw("same")
            self.histos[suffGen+'FitAC'+'UNRyqt'+c].DrawCopy("same") #to have foreground    
            self.leg['comp'+'FitAC'+'UNRyqt'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            self.leg['comp'+'FitAC'+'UNRyqt'+c].AddEntry(self.histos[suffReco+'FitAC'+'UNRyqt'+c+'4comp'],'Reco Asimov')
            self.leg['comp'+'FitAC'+'UNRyqt'+c].AddEntry(self.histos[suffGen+'FitAC'+'UNRyqt'+c],'Gen Asimov')
            self.leg['comp'+'FitAC'+'UNRyqt'+c].AddEntry(self.histos[suffGen+'FitBand'+'UNRyqt'+c],'Gen Syst')
            self.leg['comp'+'FitAC'+'UNRyqt'+c].Draw("same")
            
            #integrated
            self.canvas['comp'+'FitAC'+'qt'+c] = self.canvas[suffGen+'FitAC'+'qt'+c].Clone(self.canvas[suffGen+'FitAC'+'qt'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'qt'+c].SetTitle(self.canvas[suffGen+'FitAC'+'qt'+c].GetTitle().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'qt'+c].SetName(self.canvas[suffGen+'FitAC'+'qt'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'qt'+c].cd()
            self.histos[suffReco+'FitAC'+'qt'+c+'4comp'] = self.histos[suffReco+'FitAC'+'qt'+c].Clone(suffReco+'FitAC'+'qt'+c+'4comp')
            self.histos[suffReco+'FitAC'+'qt'+c+'4comp'].SetLineColor(ROOT.kGreen-3)
            self.histos[suffReco+'FitAC'+'qt'+c+'4comp'].Draw("same")
            self.histos[suffGen+'FitAC'+'qt'+c].DrawCopy("same") #to have foreground
            self.leg['comp'+'FitAC'+'qt'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            self.leg['comp'+'FitAC'+'qt'+c].AddEntry(self.histos[suffReco+'FitAC'+'qt'+c+'4comp'],'Reco Asimov')
            self.leg['comp'+'FitAC'+'qt'+c].AddEntry(self.histos[suffGen+'FitAC'+'qt'+c],'Gen Asimov')
            self.leg['comp'+'FitAC'+'qt'+c].AddEntry(self.histos[suffGen+'FitBand'+'qt'+c],'Gen Syst')
            self.leg['comp'+'FitAC'+'qt'+c].Draw("same")
            
            self.canvas['comp'+'FitAC'+'y'+c] = self.canvas[suffGen+'FitAC'+'y'+c].Clone(self.canvas[suffGen+'FitAC'+'y'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'y'+c].SetTitle(self.canvas[suffGen+'FitAC'+'y'+c].GetTitle().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'y'+c].SetName(self.canvas[suffGen+'FitAC'+'y'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitAC'+'y'+c].cd()
            self.histos[suffReco+'FitAC'+'y'+c+'4comp'] = self.histos[suffReco+'FitAC'+'y'+c].Clone(suffReco+'FitAC'+'y'+c+'4comp')
            self.histos[suffReco+'FitAC'+'y'+c+'4comp'].SetLineColor(ROOT.kGreen-3)
            self.histos[suffReco+'FitAC'+'y'+c+'4comp'].Draw("same")
            self.histos[suffGen+'FitAC'+'y'+c].DrawCopy("same") #to have foreground
            self.leg['comp'+'FitAC'+'y'+c] = ROOT.TLegend(0.5,0.75,0.9,0.9)
            self.leg['comp'+'FitAC'+'y'+c].AddEntry(self.histos[suffReco+'FitAC'+'y'+c+'4comp'],'Reco Asimov')
            self.leg['comp'+'FitAC'+'y'+c].AddEntry(self.histos[suffGen+'FitAC'+'y'+c],'Gen Asimov')
            self.leg['comp'+'FitAC'+'y'+c].AddEntry(self.histos[suffGen+'FitBand'+'y'+c],'Gen Syst')
            self.leg['comp'+'FitAC'+'y'+c].Draw("same")
            
            
            #-------------------- Ai errors plots -------------------------------
            
            for i in range(1, self.histos[suffGen+'FitErr'+c].GetNbinsX()+1):
                self.canvas['comp'+'FitErr'+'qt'+str(i)+c] = self.canvas[suffGen+'FitErr'+'qt'+str(i)+c].Clone(self.canvas[suffGen+'FitErr'+'qt'+str(i)+c].GetName().replace(suffGen,'comp'))
                self.canvas['comp'+'FitErr'+'qt'+str(i)+c].SetTitle(self.canvas[suffGen+'FitErr'+'qt'+str(i)+c].GetTitle().replace(suffGen,'comp'))
                self.canvas['comp'+'FitErr'+'qt'+str(i)+c].SetName(self.canvas[suffGen+'FitErr'+'qt'+str(i)+c].GetName().replace(suffGen,'comp'))
                self.canvas['comp'+'FitErr'+'qt'+str(i)+c].cd()
                self.histos[suffReco+'FitErr'+'qt'+str(i)+c+'4comp'] = self.histos[suffReco+'FitErr'+'qt'+str(i)+c].Clone(suffReco+'FitErr'+'qt'+str(i)+c+'4comp')
                self.histos[suffReco+'FitErr'+'qt'+str(i)+c+'4comp'].SetLineColor(ROOT.kGreen-3)
                self.histos[suffReco+'FitErr'+'qt'+str(i)+c+'4comp'].Draw("same")
                self.histos[suffGen+'FitErr'+'qt'+str(i)+c].DrawCopy("same") #to have foreground
                
                self.leg['comp'+'FitErr'+'qt'+str(i)+c] = self.leg[suffGen+'FitErr'+'qt'+str(i)+c].Clone()
                self.leg['comp'+'FitErr'+'qt'+str(i)+c].AddEntry(self.histos[suffReco+'FitErr'+'qt'+str(i)+c+'4comp'],'Reco Fit Err.')
                self.leg['comp'+'FitErr'+'qt'+str(i)+c].Draw("same")
                
            for j in range(1, self.histos[suffGen+'FitErr'+c].GetNbinsY()+1):
                self.canvas['comp'+'FitErr'+'y'+str(j)+c] = self.canvas[suffGen+'FitErr'+'y'+str(j)+c].Clone(self.canvas[suffGen+'FitErr'+'y'+str(j)+c].GetName().replace(suffGen,'comp'))
                self.canvas['comp'+'FitErr'+'y'+str(j)+c].SetTitle(self.canvas[suffGen+'FitErr'+'y'+str(j)+c].GetTitle().replace(suffGen,'comp'))
                self.canvas['comp'+'FitErr'+'y'+str(j)+c].SetName(self.canvas[suffGen+'FitErr'+'y'+str(j)+c].GetName().replace(suffGen,'comp'))
                self.canvas['comp'+'FitErr'+'y'+str(j)+c].cd()
                self.histos[suffReco+'FitErr'+'y'+str(j)+c+'4comp'] = self.histos[suffReco+'FitErr'+'y'+str(j)+c].Clone(suffReco+'FitErr'+'y'+str(j)+c+'4comp')
                self.histos[suffReco+'FitErr'+'y'+str(j)+c+'4comp'].SetLineColor(ROOT.kGreen-3)
                self.histos[suffReco+'FitErr'+'y'+str(j)+c+'4comp'].Draw("same")
                self.histos[suffGen+'FitErr'+'y'+str(j)+c].DrawCopy("same") #to have foreground
                
                self.leg['comp'+'FitErr'+'y'+str(j)+c] = self.leg[suffGen+'FitErr'+'y'+str(j)+c].Clone()
                self.leg['comp'+'FitErr'+'y'+str(j)+c].AddEntry(self.histos[suffReco+'FitErr'+'y'+str(j)+c+'4comp'],'Reco Fit Err.')
                self.leg['comp'+'FitErr'+'y'+str(j)+c].Draw("same")
            
            self.canvas['comp'+'FitErr'+'UNRqty'+c] = self.canvas[suffGen+'FitErr'+'UNRqty'+c].Clone(self.canvas[suffGen+'FitErr'+'UNRqty'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'UNRqty'+c].SetTitle(self.canvas[suffGen+'FitErr'+'UNRqty'+c].GetTitle().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'UNRqty'+c].SetName(self.canvas[suffGen+'FitErr'+'UNRqty'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'UNRqty'+c].cd()
            self.histos[suffReco+'FitErr'+'UNRqty'+c+'4comp'] = self.histos[suffReco+'FitErr'+'UNRqty'+c].Clone(suffReco+'FitErr'+'y'+str(j)+c+'4comp')
            self.histos[suffReco+'FitErr'+'UNRqty'+c+'4comp'].SetLineColor(ROOT.kGreen-3)
            self.histos[suffReco+'FitErr'+'UNRqty'+c+'4comp'].Draw("same")
            self.histos[suffGen+'FitErr'+'UNRqty'+c].DrawCopy("same") #to have foreground
            self.leg['comp'+'FitErr'+'UNRqty'+c] = self.leg[suffGen+'FitErr'+'UNRqty'+c].Clone()
            self.leg['comp'+'FitErr'+'UNRqty'+c].AddEntry(self.histos[suffReco+'FitErr'+'UNRqty'+c+'4comp'],'Reco Fit Err.')
            self.leg['comp'+'FitErr'+'UNRqty'+c].Draw("same")
            
            self.canvas['comp'+'FitErr'+'UNRyqt'+c] = self.canvas[suffGen+'FitErr'+'UNRyqt'+c].Clone(self.canvas[suffGen+'FitErr'+'UNRyqt'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'UNRyqt'+c].SetTitle(self.canvas[suffGen+'FitErr'+'UNRyqt'+c].GetTitle().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'UNRyqt'+c].SetName(self.canvas[suffGen+'FitErr'+'UNRyqt'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'UNRyqt'+c].cd()
            self.histos[suffReco+'FitErr'+'UNRyqt'+c+'4comp'] = self.histos[suffReco+'FitErr'+'UNRyqt'+c].Clone(suffReco+'FitErr'+'y'+str(j)+c+'4comp')
            self.histos[suffReco+'FitErr'+'UNRyqt'+c+'4comp'].SetLineColor(ROOT.kGreen-3)
            self.histos[suffReco+'FitErr'+'UNRyqt'+c+'4comp'].Draw("same")
            self.histos[suffGen+'FitErr'+'UNRyqt'+c].DrawCopy("same") #to have foreground    
            self.leg['comp'+'FitErr'+'UNRyqt'+c] = self.leg[suffGen+'FitErr'+'UNRyqt'+c].Clone()
            self.leg['comp'+'FitErr'+'UNRyqt'+c].AddEntry(self.histos[suffReco+'FitErr'+'UNRyqt'+c+'4comp'],'Reco Fit Err.')
            self.leg['comp'+'FitErr'+'UNRyqt'+c].Draw("same")
            
            #integrated
            self.canvas['comp'+'FitErr'+'qt'+c] = self.canvas[suffGen+'FitErr'+'qt'+c].Clone(self.canvas[suffGen+'FitErr'+'qt'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'qt'+c].SetTitle(self.canvas[suffGen+'FitErr'+'qt'+c].GetTitle().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'qt'+c].SetName(self.canvas[suffGen+'FitErr'+'qt'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'qt'+c].cd()
            self.histos[suffReco+'FitErr'+'qt'+c+'4comp'] = self.histos[suffReco+'FitErr'+'qt'+c].Clone(suffReco+'FitErr'+'qt'+c+'4comp')
            self.histos[suffReco+'FitErr'+'qt'+c+'4comp'].SetLineColor(ROOT.kGreen-3)
            self.histos[suffReco+'FitErr'+'qt'+c+'4comp'].Draw("same")
            self.histos[suffGen+'FitErr'+'qt'+c].DrawCopy("same") #to have foreground
            
            self.leg['comp'+'FitErr'+'qt'+c] = self.leg[suffGen+'FitErr'+'qt'+c].Clone()
            self.leg['comp'+'FitErr'+'qt'+c].AddEntry(self.histos[suffReco+'FitErr'+'qt'+c+'4comp'],'Reco Fit Err.')
            self.leg['comp'+'FitErr'+'qt'+c].Draw("same")
            
            self.canvas['comp'+'FitErr'+'y'+c] = self.canvas[suffGen+'FitErr'+'y'+c].Clone(self.canvas[suffGen+'FitErr'+'y'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'y'+c].SetTitle(self.canvas[suffGen+'FitErr'+'y'+c].GetTitle().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'y'+c].SetName(self.canvas[suffGen+'FitErr'+'y'+c].GetName().replace(suffGen,'comp'))
            self.canvas['comp'+'FitErr'+'y'+c].cd()
            self.histos[suffReco+'FitErr'+'y'+c+'4comp'] = self.histos[suffReco+'FitErr'+'y'+c].Clone(suffReco+'FitErr'+'y'+c+'4comp')
            self.histos[suffReco+'FitErr'+'y'+c+'4comp'].SetLineColor(ROOT.kGreen-3)
            self.histos[suffReco+'FitErr'+'y'+c+'4comp'].Draw("same")
            self.histos[suffGen+'FitErr'+'y'+c].DrawCopy("same") #to have foreground
            
            self.leg['comp'+'FitErr'+'y'+c] = self.leg[suffGen+'FitErr'+'y'+c].Clone()
            self.leg['comp'+'FitErr'+'y'+c].AddEntry(self.histos[suffReco+'FitErr'+'y'+c+'4comp'],'Reco Fit Err.')
            self.leg['comp'+'FitErr'+'y'+c].Draw("same")
            

    def makeRootOutput(self,outFileName,SAVE, suffList, comparison,toy,impact, postfit,saver, aposteriori) :
        
        outFile = ROOT.TFile(outFileName+"_"+self.anaKind['name']+".root", "recreate")
        outFile.cd()
        dirFinalDict = {}

        for suff in suffList : 
            dirFinalDict['coeff2D'+suff] = outFile.mkdir('coeff2D_'+suff)
            dirFinalDict['coeff2D'+suff].cd()
            # for hel in self.hels:
            #     self.histos[suff+'FitHel'+hel].Write()
            #     self.histos[suff+'FitHel'+hel+'norm'].Write()
            for c in self.coeffDict:
                self.histos[suff+'FitAC'+c].Write()
                self.histos[suff+'FitACqt'+c].Write()
                self.histos[suff+'FitACy'+c].Write()
                self.histos[suff+'FitBand'+c].Write()
                        
            if self.anaKind['angNames'] : 
                for cat in self.category :
                    self.histos[suff+'FitAC'+cat+'unpolarizedxsec'+'4apo'].Write()
                
                for c in self.coeffDict:
                    for kat in ['y','qt', 'UNRyqt'] :
                        self.histos[suff+'FitRatioAC'+kat+c].Write()
                        if aposteriori :
                            self.histos[suff+'apoRatio'+kat+c].Write()
                    
                    
            dirFinalDict['coeff'+suff] = outFile.mkdir('coeff_'+suff)
            dirFinalDict['coeff'+suff].cd()
            for c in self.coeffDict:
                self.canvas[suff+'FitAC'+'qt'+c].Write()
                self.canvas[suff+'FitErr'+'qt'+c].Write()
                self.canvas[suff+'FitAC'+'y'+c].Write()
                self.canvas[suff+'FitErr'+'y'+c].Write()
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1):
                    self.canvas[suff+'FitAC'+'qt'+str(i)+c].Write()
                    self.canvas[suff+'FitErr'+'qt'+str(i)+c].Write()
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1):
                    self.canvas[suff+'FitAC'+'y'+str(j)+c].Write()
                    self.canvas[suff+'FitErr'+'y'+str(j)+c].Write()
        
            dirFinalDict['coeff_unrolled'+suff] = outFile.mkdir('coeff_unrolled_'+suff)
            dirFinalDict['coeff_unrolled'+suff].cd()
            for c in self.coeffDict:
                self.canvas[suff+'FitAC'+'UNRqty'+c].Write()
                self.canvas[suff+'FitErr'+'UNRqty'+c].Write()
                self.canvas[suff+'FitAC'+'UNRyqt'+c].Write()
                self.canvas[suff+'FitErr'+'UNRyqt'+c].Write()
            
            dirFinalDict['matrices'+suff] = outFile.mkdir('matrices_'+suff)
            dirFinalDict['matrices'+suff].cd()
            for mtx in ['corr','cov'] :
                try : 
                    self.histos[suff+mtx+'Mat'].Write()
                    self.histos[suff+mtx+'MatIntegrated'].Write()
                except :
                    print("missing total matrix in writing")
                for c in self.coeffDict:
                    self.canvas[suff+mtx+'Mat'+c].Write()
                for i in range(1, self.histos[suff+'FitAC'+self.coeffList[0]].GetNbinsX()+1):
                    self.canvas[suff+mtx+'Mat'+'y'+str(i)].Write()
                for j in range(1, self.histos[suff+'FitAC'+self.coeffList[0]].GetNbinsY()+1):
                    self.canvas[suff+mtx+'Mat'+'qt'+str(j)].Write()
                self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].Write()
                self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].Write()
            
            dirFinalDict['nuisance'+suff] = outFile.mkdir('nuisance_'+suff)
            dirFinalDict['nuisance'+suff].cd()
            for nuiDict_key, nuiDict_val in self.NuiConstrDict.items() :
                self.canvas[suff+'NuiConstr'+nuiDict_key].Write()
                self.histos[suff+'NuiConstr'+nuiDict_key].Write()
            self.histos[suff+'mass'].Write()
            
            if postfit : 
                dirFinalDict['prefit_postfit'+suff] = outFile.mkdir('prefit_postfit_'+suff)
                dirFinalDict['prefit_postfit'+suff].cd()
                for k in ['prefit','postfit'] :
                    for var,varInfo in self.postFitVars.items() :
                        self.canvas[suff+k+var].Write()       
                
            if impact :
                dirFinalDict['impact'+suff] = outFile.mkdir('impact_'+suff)
                dirFinalDict['impact'+suff].cd()
                for c in self.coeffDict :
                    self.canvas[suff+'impact'+'UNRqty'+c].Write()
                    self.canvas[suff+'impact'+'UNRyqt'+c].Write()
                    self.canvas[suff+'impact'+'y'+c].Write()
                    self.canvas[suff+'impact'+'qt'+c].Write()
                self.canvas[suff+'impact'+'mass'].Write()
                
                if self.anaKind['angNames'] :
                    cp = 'unpolarizedxsec'+'poiss'
                    self.canvas[suff+'impact'+'UNRqty'+cp].Write()
                    self.canvas[suff+'impact'+'UNRyqt'+cp].Write()
                    self.canvas[suff+'impact'+'y'+cp].Write()
                    self.canvas[suff+'impact'+'qt'+cp].Write()
                
                
                
            
            if toy!='' :
                dirFinalDict['toys'+suff] = outFile.mkdir('toys'+suff)
                dirFinalDict['toys'+suff].cd()
                self.histos[suff+'mass'+'toy'].Write()
                self.histos[suff+'mass'+'toyPull'].Write()
                for c in self.coeffDict:
                    for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1):
                        for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1):
                            self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)].Write()
                            self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)].Write()
                    for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1):
                        self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)].Write()
                        self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)].Write()
                    for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1):
                        self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)].Write()
                        self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)].Write()
                
                dirFinalDict['toys_summary_'+suff] = outFile.mkdir('toys_summary_'+suff)
                dirFinalDict['toys_summary_'+suff].cd()
                self.histos[suff+'mass'+'toy'+'mean'].Write()
                self.histos[suff+'mass'+'toyPull'+'mean'].Write()
                for c in self.coeffDict:
                    self.canvas[suff+'FitAC'+'UNRqty'+c+'toyPull'].Write()
                    self.canvas[suff+'FitAC'+'UNRyqt'+c+'toyPull'].Write()
                    self.canvas[suff+'FitAC'+'y'+c+'toyPull'].Write()
                    self.canvas[suff+'FitAC'+'qt'+c+'toyPull'].Write()
                for kindpull in ['UNRqty','UNRyqt','qt','y'] :
                    self.canvas[suff+'FitAC'+kindpull+'toyPull'].Write()
            
          
        if comparison :
            
            dirFinalDict['comparison_coeff'] = outFile.mkdir('comparison_coeff')
            dirFinalDict['comparison_coeff'].cd()
            for c in self.coeffDict:
                self.canvas['comp'+'FitAC'+'qt'+c].Write()
                self.canvas['comp'+'FitErr'+'qt'+c].Write()
                self.canvas['comp'+'FitAC'+'y'+c].Write()
                self.canvas['comp'+'FitErr'+'y'+c].Write()
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1):
                    self.canvas['comp'+'FitAC'+'qt'+str(i)+c].Write()
                    self.canvas['comp'+'FitErr'+'qt'+str(i)+c].Write()
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1):
                     self.canvas['comp'+'FitAC'+'y'+str(j)+c].Write()
                     self.canvas['comp'+'FitErr'+'y'+str(j)+c].Write()
                     
            dirFinalDict['comparison_coeff_unrolled'] = outFile.mkdir('comparison_coeff_unrolled')
            dirFinalDict['comparison_coeff_unrolled'].cd()
            for c in self.coeffDict:
                self.canvas['comp'+'FitAC'+'UNRqty'+c].Write()                
                self.canvas['comp'+'FitErr'+'UNRqty'+c].Write()
                self.canvas['comp'+'FitAC'+'UNRyqt'+c].Write()
                self.canvas['comp'+'FitErr'+'UNRyqt'+c].Write()
        
        
        
        if toy!= '' :
            outFileToy = ROOT.TFile(outFileName+"_TOYonly_"+self.anaKind['name']+".root", "recreate")
            outFileToy.cd()
            self.histos[suff+'mass'+'toy'].Write()
            self.histos[suff+'mass'+'toyPull'].Write()
            for c in self.coeffDict:
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1):
                    for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1):
                        self.histos[suff+'FitAC'+c+'toy'+'y'+str(i)+'_qt'+str(j)].Write()
                        self.histos[suff+'FitAC'+c+'toyPull'+'y'+str(i)+'_qt'+str(j)].Write()
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1):
                    self.histos[suff+'FitACy'+c+'toy'+'y'+str(i)].Write()
                    self.histos[suff+'FitACy'+c+'toyPull'+'y'+str(i)].Write()
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1):
                    self.histos[suff+'FitACqt'+c+'toy'+'qt'+str(j)].Write()
                    self.histos[suff+'FitACqt'+c+'toyPull'+'qt'+str(j)].Write()
            
        #################### saver###########################
        if saver :
            foldername = "histo_" + outFile.GetName().replace('.root','')
            if not os.path.exists(foldername): os.system("mkdir -p "+foldername)
            
            if 'plus' in foldername : sign = 'plus'
            elif 'minus' in foldername : sign = 'minus'
            else : sign = ''
            
            if self.asimov :
                asimovString = 'asimov'
            else :
                asimovString = 'toy'
            if aposteriori :
                asimovString=asimovString+'_aposteriori'


            nameDict = { #which canvas must be saved, folder : [pe coeff?, canvasname1,...canvasnameN]
                'AngCoeff' : [1, 'FitACUNRyqt', 'FitACqt','FitACy'],
                'matrices' : [1, 'corrMat'],
                'matricesInt': [0, 'corrMatIntegratedqt', 'corrMatIntegratedy'],
                'nuisance' : [0, 'NuiConstrall'],
                # 'impact' : [0, 'impactUNRyqtA4', 'impactyA4','impactqtA4', 'impactUNRyqtunpolarizedxsec', 'impactyunpolarizedxsec','impactqtunpolarizedxsec', 'impactmass'],
                'impact' : [1, 'impactUNRyqt', 'impacty','impactqt'],
                'impactmass' : [0, 'impactmass'],
                'prefit' : [0, 'prefit', 'postfit'],
                'toy' : [0,'FitAC'+'UNRyqt'+'toyPull','FitAC'+'qt'+'toyPull','FitAC'+'y'+'toyPull'],
                'coeffErr' : [0,'FitErr'+'UNRyqt'+'unpolarizedxsec', 'FitErr'+'y'+'unpolarizedxsec','FitErr'+'qt'+'unpolarizedxsec'],
            }
                        
            if not postfit : del nameDict['prefit']
            if not impact : del nameDict['impact']
            if not impact : del nameDict['impactmass']
            if toy=='' : 
                del nameDict['toy']
            
            for suff in suffList :
                for ind, nameList in nameDict.items() :
                    if nameList[0] :
                        for c in self.coeffDict:
                            for name in nameList[1:] :
                                self.canvas[suff+name+c].SaveAs(foldername+'/'+self.canvas[suff+name+c].GetName().replace(suff+'_c','Fit_'+asimovString+'_'+s)+'.pdf')
                                self.canvas[suff+name+c].SaveAs(foldername+'/'+self.canvas[suff+name+c].GetName().replace(suff+'_c','Fit_'+asimovString+'_'+s)+'.png')
                    else :
                        for name in nameList[1:] :
                            self.canvas[suff+name].SaveAs(foldername+'/'+self.canvas[suff+name].GetName().replace(suff+'_c','Fit_'+asimovString+'_'+s)+'.pdf')
                            self.canvas[suff+name].SaveAs(foldername+'/'+self.canvas[suff+name].GetName().replace(suff+'_c','Fit_'+asimovString+'_'+s)+'.png')
            
        
        
        outFile.Close()

        
                
               
    def getCoeffDict(self) :
        return self.coeffDict
    
    def getYArr(self) :
        return self.yArr
    
    def getQtArr(self) :
        return self.qtArr
    
# def saver(rootInput, suffList, comparison,coeffDict,yArr,qtArr,impact, postfit, anaKind) : #not used, integrated in the writer
#     # if anaKind['angNames']!='angCoeff' : 
#         # print("saved plot not implemented for helicity x section")
#         # return 
    
#     if not os.path.exists(rootInput): os.system("mkdir -p histo_" + rootInput)
#     FitFile = ROOT.TFile.Open(rootInput+'.root')
    
#     nameDict = { #which canvas must be saved, folder : [pe coeff?, canvasname1,...canvasnameN]
#         'coeff_unrolled_' : [1, '_c_UNRyqt'],
#         'coeff_' : [1, '_c_QT_','_c_Y_'],
#         'matrices_' : [1, '_c_corrMat_'],
#         'matrices_': [0, '_c_corrMat_Integrated_qt', '_c_corrMat_Integrated_y'],
#         'nuisance_' : [0, '_c_NuiConstr_all'],
#         'impact_' : [0, '_c_impact_UNRqtyA4', '_c_impact_yA4','_c_impact_qtA4', '_c_impact_UNRqtyunpolarizedxsec', '_c_impact_yunpolarizedxsec','_c_impact_qtunpolarizedxsec'],
#         'prefit_postfit_' : [0, '_prefit_', '_postfit_']
#     }
    
#     for suff in suffList :
#         for fold, nameList in nameDict.items() :
#             if nameList[0] :
#                 for name in nameList[1:] :
#                     for c in coeffDict:
#                         can = FitFile.Get(fold+suff+'/'+suff+name+c)
#                         print(can, fold+suff+'/'+suff+name+c, FitFile)
#                         can.SaveAs('histo_'+rootInput+'/'+can.GetName().replace(suff,'Fit')+'.pdf')
#                         can.SaveAs('histo_'+rootInput+'/'+can.GetName().replace(suff,'Fit')+'.png')
#             else :
#                 for name in nameList[1:] :
#                     can = FitFile.Get(fold+suff+'/'+suff+name)
#                     print(can, fold+suff+'/'+suff+name)
#                     can.SaveAs('histo_'+rootInput+'/'+can.GetName().replace(suff,'Fit')+'.pdf')
#                     can.SaveAs('histo_'+rootInput+'/'+can.GetName().replace(suff,'Fit')+'.png')
            
        
        
    
    
#     if (0) : #old saver
#         unrList = ['UNRqty','Err_UNRqty','UNRyqt','Err_UNRyqt']
#         intList = ['QT_','Y_']

#         if comparison :
#             for c in coeffDict:
#                 for name in unrList :
#                     can = FitFile.Get('comparison_coeff_unrolled/comp_c_'+name+c)
#                     can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
#                     can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
        
#         for suff in suffList : 
            
#             if not comparison :
#                 for c in coeffDict:
#                     for name in unrList :
#                         can = FitFile.Get('coeff_unrolled_'+suff+'/'+suff+'_c_'+name+c)
#                         can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
#                         can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
#                     for name in intList :
#                         for ee in ['','_Err'] :
#                             can = FitFile.Get('coeff_'+suff+'/'+suff+'_c_'+name+c+ee)
#                             can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
#                             can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
                
            
#             for mtx in ['corr','cov'] :
#                 for c in coeffDict:
#                     can = FitFile.Get('matrices_'+suff+'/'+suff+'_'+mtx+'Mat_'+c)
#                     can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
#                     can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
#                 for i in range(1, len(yArr)):
#                     can = FitFile.Get('matrices_'+suff+'/'+suff+'_'+mtx+'Mat_y'+str(i))
#                     can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
#                     can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
#                 for j in range(1, len(qtArr)):
#                     can = FitFile.Get('matrices_'+suff+'/'+suff+'_'+mtx+'Mat_qt'+str(j))
#                     can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
#                     can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
#                 for nn in ['qt','y'] :
#                     can = FitFile.Get('matrices_'+suff+'/'+suff+'_'+mtx+'Mat_Integrated_'+nn)
#                     can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
#                     can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
            
#             if impact :
#                 for c in coeffDict:
#                     for name in ['UNR','y','qt'] :
#                         can = FitFile.Get('impact_'+suff+'/'+suff+'impact_'+name+c)
#                         can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
#                         can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
#                 can = FitFile.Get('impact_'+suff+'/'+suff+'impact_mass')
#                 can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
#                 can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')        
            
                
                
                
###############################################################################################################################
#
#  usage: python plotter_fitResult.py --output fitPlots_Wplus --fitFile fit_Wplus_reco.root --suff reco --input ../../../analysisOnGen/genInput_Wplus.root --impact 1
#
#  all the parameters are described below, but some note:
#  --uncorrelate: it should be set to True accordingly to theory prescription (31 Scale variation in total)
#  --comparison: if True expected two files called --fitFile and --FiteFile_reco
#
################################################################################################################################
       
parser = argparse.ArgumentParser("")

parser.add_argument('-o','--output', type=str, default='fitResult_Wplus',help="name of the output file")
parser.add_argument('-f','--fitFile', type=str, default='fit_Wplus.root',help="name of the fit result root file. Always plus charge must be given. If comparison active the name must be: NAME.root, NAME_reco.root")
parser.add_argument('-i','--input', type=str, default='../analysisOnGen/genInput_Wplus.root',help="name of the input root file")
parser.add_argument('-u','--uncorrelate', type=int, default=True,help="if true uncorrelate num and den of Angular Coeff in MC scale variation")
parser.add_argument('-c','--comparison', type=int, default=False,help="comparison between reco and gen fit")
parser.add_argument('-s','--save', type=int, default=False,help="save .png and .pdf canvas")
parser.add_argument('-l','--suffList', type=str, default='',nargs='*', help="list of suff to be processed in the form: gen,reco")
parser.add_argument('-a','--aposteriori', type=str, default='',help="name of the aposteriori fit file, if empty not plotted")
parser.add_argument('-t','--toy', type=str, default='',help="name of the toys fit file, if empty not plotted")
parser.add_argument('-m','--impact', type=int, default=False,help="use the fitFile to produce also the impact plots (the required histos must be filled in the fit)")
parser.add_argument('-p','--postfit', type=int, default=False,help="use the fitFile to produce also the postfit plots (the required histos must be filled in the fit)")
parser.add_argument('-cl','--cleanNuisance', type=int, default=False,help="plot only the most relevant nuisances in the impacts")
parser.add_argument('-mass','--mass', type=str, default='',help="plot the mass impact using external stat using this file, if empty not done")
parser.add_argument('-data','--data', type=str, default='../Wplus_reco.root',help="data file used in the prefit, used MC if not provided")
parser.add_argument('-asimov','--asimov', type=int, default=True,help="if asimov fit (for plotting ranges and naming)")
parser.add_argument('-acc','--acceptanceFile', type=str, default='/scratch/bertacch/wmass/wproperties-analysis/Fit/OUTPUT_25Apr_noCF_qtMax60_fixLowAcc/accMap.root',help="file with acceptance info")





args = parser.parse_args()
OUTPUT = args.output
FITFILE = args.fitFile
INPUT = args.input
UNCORR = args.uncorrelate
COMP = args.comparison
SAVE= args.save
SUFFL =args.suffList
APO = args.aposteriori
TOY = args.toy
IMPACT = args.impact
POSTFIT = args.postfit
CLEANNUI = args.cleanNuisance
MASS = args.mass
DATA = args.data 
ASIMOV = args.asimov
ACCEPTANCE = args.acceptanceFile

plusOnly=False
if plusOnly :
    signList = ['plus']
else :
    signList = ['plus', 'minus']
# signList = ['minus']    
kindDict = {
    'angCoeff' :{
        'name' : 'angCoeff',
        'apoFlag' : True,
        'angNames' : True, #A0,A1...
        'differential' : True,
        'diffString' : 'helpois',
        'intString' : 'helmetapois'
    },
    'helXsec' : {
        'name' : 'helXsec',
        'apoFlag' : False,
        'angNames' : False, #L,I,T...
        'differential' : False ,
        'diffString' : 'pmaskedexp' ,
        'intString' : 'sumpois'
    },
    'mu' : {
        'name' : 'mu',
        'apoFlag' : False,
        'angNames' : False,
        'differential' : True ,
        'diffString' : 'mu',
        'intString' : 'sumpois' #placeholder
    },
    # 'normXsec' : { #this should not be used, since the normalization is not a quantity with a physical meaning (sum of the xsection in ALL the processes, NOT the UNPOL)
    #     'name' : 'normXsec',
    #     'apoFlag' : False,
    #     'angNames' : False,
    #     'differential' : False,
    #     'diffString' :  'pmaskedexpnorm',
    #     'intString' : 'sumpoisnorm'
    # }
}

for s in signList : 
    FITFILE_s = FITFILE.replace('plus',s)
    INPUT_s = INPUT.replace('plus',s)
    OUTPUT_s = OUTPUT.replace('plus',s)
    APO_s = APO.replace('plus',s)
    MASS_s = MASS.replace('plus',s)
    DATA_s = DATA.replace('plus',s)
    
    # for helXsec in [False,True] :
    for kk, kVal in kindDict.items() :
        print("kind=",kk)
        # if helXsec : continue 
        if not kVal['apoFlag'] : APO_s = '' #we do not regularize helXsec, but Ai only

        p=plotter(anaKind=kVal, sign=s, asimov=ASIMOV)
        p.AngCoeffPlots(inputFile=INPUT_s, fitFile=FITFILE_s, uncorrelate=UNCORR,suff=SUFFL[0],aposteriori=APO_s,toy=TOY,impact=IMPACT, postfit=POSTFIT, cleanNuisance=CLEANNUI, massImp=MASS_s,data=DATA_s,acceptance=ACCEPTANCE)
        if COMP :
            recoFit = FITFILE_s.replace('.root', '_'+str(SUFFL[1])+'.root')
            p.AngCoeffPlots(inputFile=INPUT_s, fitFile=recoFit, uncorrelate=UNCORR, suff=SUFFL[1])
            p.GenRecoComparison(suffGen='gen', suffReco='reco')
        p.makeRootOutput(outFileName=OUTPUT_s, SAVE=SAVE,suffList=SUFFL,comparison=COMP,toy=TOY, impact = IMPACT, postfit=POSTFIT,saver=SAVE,aposteriori=APO_s)
        # if SAVE : #old version, now insite makeRootOutput
        #     saverInput = OUTPUT_s+'_'+kk
        #     saver(rootInput=saverInput,suffList=SUFFL,comparison=COMP,coeffDict=p.getCoeffDict(),yArr=p.getYArr(),qtArr=p.getQtArr(), impact=IMPACT, postfit=POSTFIT, anaKind=kVal)