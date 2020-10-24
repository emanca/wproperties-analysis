import ROOT
from array import array
import math
from termcolor import colored
import copy
import argparse
import os


ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

class plotter :
    
    def __init__(self):
            
        self.yArr = [0, 0.4, 0.8, 1.2 ,1.6, 2.0, 2.4]
        self.qtArr = [0., 4., 8., 12., 16., 20., 24., 28., 32.]
        self.coeffArr = [0,1,2,3,4,5,6]
        self.nuisArr = ['mass']
        
        self.systDict = {
            "_LHEScaleWeight" : ["_LHEScaleWeight_muR0p5_muF0p5", "_LHEScaleWeight_muR0p5_muF1p0","_LHEScaleWeight_muR1p0_muF0p5","_LHEScaleWeight_muR1p0_muF2p0","_LHEScaleWeight_muR2p0_muF1p0", "_LHEScaleWeight_muR2p0_muF2p0"],
            "_LHEPdfWeight" : ["_LHEPdfWeightHess" + str(i)  for i in range(1, 61)],
        }
        
        self.vetoScaleList = []
        for scNum in self.systDict['_LHEScaleWeight'] :
            for scDen in self.systDict['_LHEScaleWeight'] :
                if ('0p5' in scNum and '2p0' in scDen) or ('2p0' in scNum and '0p5' in scDen): 
                    self.vetoScaleList.append([scNum,scDen])        
        
        self.hels = ['L', 'I', 'T', 'A', 'P', 'UL']
        self.coeffDict = {
            'A0' : 1.,
            'A1' : 5.,
            'A2' : 20.,
            'A3' : 4.,
            'A4' : 4.,
            'unpolarizedxsec' : 1
        }
        self.coeffList = []
        for ind,val in self.coeffDict.iteritems() :
            self.coeffList.append(ind)
        
        self.histos = {}
        self.canvas = {}
        self.leg = {}
        
        self.groupedSystColors = {
            # "WHSFVars"  : [ROOT.kGreen+1, 'Scale Factors'],
            "LHEScaleWeightVars" : [ROOT.kViolet-2, 'MC Scale'],
            # "LHEScaleWeightVars" : [ROOT.kGreen-2, 'MC Scale'],
            # "ptScaleVars" : [ROOT.kBlue-4, 'pT Scale'],
            # "jmeVars" : [ROOT.kAzure+10, 'MET'],
            "LHEPdfWeightVars" : [ROOT.kRed+1, 'PDF'],
            # "LHEPdfWeightVars" : [ROOT.kBlue+1, 'PDF'],
            "Nominal" : [1, 'Fit Unc.']
        }
        
        
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
        
        # print "unrolledQtY", self.unrolledQtY
        # print "unrolledYQt", self.unrolledYQt
        # print "unrolledAY", self.unrolledAY
        # print "unrolledAQt", self.unrolledAQt
        
        
        
    def getHistos(self,inFile, FitFile, uncorrelate,suff) :
        
        resFit = FitFile.fitresults

        #MC (Pre Fit) angular coefficients and variations
        self.histos[suff+'MC'+'mapTot'] =  inFile.Get('angularCoefficients/mapTot')
        self.histos[suff+'MCy'+'mapTot'] =  inFile.Get('angularCoefficients/Y')
        self.histos[suff+'MCqt'+'mapTot'] =  inFile.Get('angularCoefficients/Pt')
        for coeff,div in self.coeffDict.iteritems() :
            self.histos[suff+'MC'+coeff] =  inFile.Get('angularCoefficients/harmonics'+coeff+'_nom_nom')
            self.histos[suff+'MC'+'y'+coeff] =  inFile.Get('angularCoefficients/harmonicsY'+coeff+'_nom_nom')
            self.histos[suff+'MC'+'qt'+coeff] =  inFile.Get('angularCoefficients/harmonicsPt'+coeff+'_nom_nom')

        for sKind, sList in self.systDict.iteritems():

            sListMod = copy.deepcopy(sList)
            if sKind=='_LHEScaleWeight' and UNCORR :
                sListMod.append("_nom") #add nominal variation
            
            for sName in sListMod :
                self.histos[suff+'MC'+sName+'mapTot'] =  inFile.Get('angularCoefficients'+sKind+'/mapTot'+sName)
                self.histos[suff+'MCy'+sName+'mapTot'] =  inFile.Get('angularCoefficients'+sKind+'/Y'+sName)
                self.histos[suff+'MCqt'+sName+'mapTot'] =  inFile.Get('angularCoefficients'+sKind+'/Pt'+sName)
                for sNameDen in sListMod :
                    if sNameDen!=sName and not (sKind=='_LHEScaleWeight' and UNCORR) : #PDF or correlated Scale
                        continue 
                    if sName=='_nom' and sNameDen=='_nom' : 
                        continue
                    for coeff,div in self.coeffDict.iteritems() :
                        if "unpolarizedxsec" in coeff: continue
                        if UNCORR :
                            if sKind=='_LHEScaleWeight':
                                self.histos[suff+'MC'+sName+sNameDen+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonics'+coeff+sName+sNameDen)
                                self.histos[suff+'MCy'+sName+sNameDen+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonicsY'+coeff+sName+sNameDen)
                                self.histos[suff+'MCqt'+sName+sNameDen+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonicsPt'+coeff+sName+sNameDen)
                            else :
                                self.histos[suff+'MC'+sName+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonics'+coeff+sName+sNameDen)
                                self.histos[suff+'MCy'+sName+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonicsY'+coeff+sName+sNameDen)
                                self.histos[suff+'MCqt'+sName+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonicsPt'+coeff+sName+sNameDen)
                        else :
                            self.histos[suff+'MC'+sName+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonics'+coeff+sName) 
                            self.histos[suff+'MCy'+sName+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonicsY'+coeff+sName) 
                            self.histos[suff+'MCqt'+sName+coeff] =  inFile.Get('angularCoefficients'+sKind+'/harmonicsPt'+coeff+sName) 
        
         
        #fit - helicity xsecs histos 
        for hel in self.hels:
            self.histos[suff+'FitHel'+hel] = ROOT.TH2D(suff+'FitHel{c}'.format(c=hel), suff+'FitHel{c}'.format(c=hel), len(self.yArr)-1, array('f',self.yArr), len(self.qtArr)-1, array('f',self.qtArr))
            self.histos[suff+'FitHel'+hel+'norm'] = ROOT.TH2D(suff+'FitHel{c}norm'.format(c=hel), suff+'FitHel{c}norm'.format(c=hel), len(self.yArr)-1, array('f',self.yArr), len(self.qtArr)-1, array('f',self.qtArr))
            self.histos[suff+'FitHel'+hel].GetXaxis().SetTitle('W rapidity')
            self.histos[suff+'FitHel'+hel].GetYaxis().SetTitle('W q_{T}')
            self.histos[suff+'FitHel'+hel+'norm'].GetXaxis().SetTitle('W rapidity')
            self.histos[suff+'FitHel'+hel+'norm'].GetYaxis().SetTitle('W q_{T}')
            for ev in resFit: #dummy because there's one event only
                for i in range(1, self.histos[suff+'FitHel'+hel].GetNbinsX()+1): #loop over rapidity bins
                    for j in range(1, self.histos[suff+'FitHel'+hel].GetNbinsY()+1): #loop over pt bins
                        try:
                            coeff = eval('ev.helXsecs{}_y_{}_qt_{}_pmaskedexp'.format(hel, i, j))
                            coeff_err = eval('ev.helXsecs{}_y_{}_qt_{}_pmaskedexp_err'.format(hel, i, j))
                            self.histos[suff+'FitHel'+hel].SetBinContent(i,j,coeff)
                            self.histos[suff+'FitHel'+hel].SetBinError(i,j,coeff_err)

                            coeffnorm = eval('ev.helXsecs{}_y_{}_qt_{}_pmaskedexpnorm'.format(hel, i, j))
                            coeffnorm_err = eval('ev.helXsecs{}_y_{}_qt_{}_pmaskedexpnorm_err'.format(hel, i, j))
                            self.histos[suff+'FitHel'+hel+'norm'].SetBinContent(i,j,coeffnorm)
                            self.histos[suff+'FitHel'+hel+'norm'].SetBinError(i,j,coeffnorm_err)
                            
                        except AttributeError:
                            pass

        #fit - angular coefficients with band
        for c in self.coeffDict:
            self.histos[suff+'FitAC'+c] = ROOT.TH2D(suff+'FitAC{c}'.format(c=c), suff+'FitAC{c}'.format(c=c), len(self.yArr)-1, array('f',self.yArr), len(self.qtArr)-1, array('f',self.qtArr))
            self.histos[suff+'FitACqt'+c] = ROOT.TH1D(suff+'FitACqt{c}'.format(c=c), suff+'FitACqt{c}'.format(c=c), len(self.qtArr)-1, array('f',self.qtArr))
            self.histos[suff+'FitACy'+c] = ROOT.TH1D(suff+'FitACy{c}'.format(c=c), suff+'FitACy{c}'.format(c=c), len(self.yArr)-1, array('f',self.yArr))
            for ev in resFit: #dummy because there's one event only
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                    for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                        try:
                            coeff = eval('ev.y_{i}_qt_{j}_{c}'.format(c=c, j=j, i=i))
                            coeff_err = eval('ev.y_{i}_qt_{j}_{c}_err'.format(c=c, j=j, i=i))
                            if 'unpol' in c:
                                coeff = coeff/(3./16./math.pi)
                                coeff_err = coeff_err/(3./16./math.pi)
                            self.histos[suff+'FitAC'+c].SetBinContent(i,j,coeff)
                            self.histos[suff+'FitAC'+c].SetBinError(i,j,coeff_err)
                            
                        except AttributeError: 
                            pass
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                    try:
                        coeff = eval('ev.qt_{j}_helmeta_{c}'.format(c=c, j=j))
                        coeff_err = eval('ev.qt_{j}_helmeta_{c}_err'.format(c=c, j=j))
                        if 'unpol' in c:
                                coeff = coeff/(3./16./math.pi)
                                coeff_err = coeff_err/(3./16./math.pi)
                        self.histos[suff+'FitACqt'+c].SetBinContent(j,coeff)
                        self.histos[suff+'FitACqt'+c].SetBinError(j,coeff_err)


                    except AttributeError: 
                        pass
                for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                    try:
                        coeff = eval('ev.y_{i}_helmeta_{c}'.format(c=c, i=i))
                        coeff_err = eval('ev.y_{i}_helmeta_{c}_err'.format(c=c, i=i))
                        if 'unpol' in c:
                                coeff = coeff/(3./16./math.pi)
                                coeff_err = coeff_err/(3./16./math.pi)
                        self.histos[suff+'FitACy'+c].SetBinContent(i,coeff)
                        self.histos[suff+'FitACy'+c].SetBinError(i,coeff_err)

                    except AttributeError: 
                        pass
            
        #covariance and correlation matrices 
        self.histos[suff+'corrMat'] = FitFile.Get('correlation_matrix_channelhelpois')
        self.histos[suff+'covMat'] = FitFile.Get('covariance_matrix_channelhelpois')  
        self.histos[suff+'corrMat'+'Integrated'] = FitFile.Get('correlation_matrix_channelhelmetapois')
        self.histos[suff+'covMat'+'Integrated'] = FitFile.Get('covariance_matrix_channelhelmetapois')  
              
        
    def AngCoeffPlots(self,inputFile, fitFile, uncorrelate,suff) :
        
        FitFile = ROOT.TFile.Open(fitFile)
        inFile = ROOT.TFile.Open(inputFile)
        self.getHistos(inFile=inFile, FitFile=FitFile, uncorrelate=uncorrelate,suff=suff)
        # print "WARNING: syst.band done with the fit central value (ok if asimov only)"
        
        # ------------- build the bands for the angular coefficient ---------------------------# 
        for c in self.coeffDict:
            # self.histos[suff+'FitBand'+c] = self.histos[suff+'FitAC'+c].Clone('FitBand'+c) #from fit (good in case of asimov)
            # self.histos[suff+'FitBandPDF'+c] = self.histos[suff+'FitAC'+c].Clone('FitBandPDF'+c) #from fit (good in case of asimov) NB: NB: not used in the plots now
            # self.histos[suff+'FitBandScale'+c] = self.histos[suff+'FitAC'+c].Clone('FitBandScale'+c) #from fit (good in case of asimov) NB: not used in the plots now
            
            # self.histos[suff+'FitBandy'+c] = self.histos[suff+'FitACy'+c].Clone('FitBandy'+c)
            # self.histos[suff+'FitBandPDFy'+c] = self.histos[suff+'FitACy'+c].Clone('FitBandPDFy'+c)
            # self.histos[suff+'FitBandScaley'+c] = self.histos[suff+'FitACy'+c].Clone('FitBandScaley'+c)
            
            # self.histos[suff+'FitBandqt'+c] = self.histos[suff+'FitACqt'+c].Clone('FitBandqt'+c)
            # self.histos[suff+'FitBandPDFqt'+c] = self.histos[suff+'FitACqt'+c].Clone('FitBandPDFqt'+c)
            # self.histos[suff+'FitBandScaleqt'+c] = self.histos[suff+'FitACqt'+c].Clone('FitBandScaleqt'+c)
            if not "unpol" in c:    
                self.histos[suff+'FitBand'+c] = self.histos[suff+'MC'+c].Clone('FitBand'+c)
                self.histos[suff+'FitBandPDF'+c] = self.histos[suff+'MC'+c].Clone('FitBandPDF'+c) 
                self.histos[suff+'FitBandScale'+c] = self.histos[suff+'MC'+c].Clone('FitBandScale'+c) 
                
                self.histos[suff+'FitBandy'+c] = self.histos[suff+'MC'+'y'+c].Clone('FitBandy'+c)
                self.histos[suff+'FitBandPDFy'+c] = self.histos[suff+'MC'+'y'+c].Clone('FitBandPDFy'+c)
                self.histos[suff+'FitBandScaley'+c] = self.histos[suff+'MC'+'y'+c].Clone('FitBandScaley'+c)
                
                self.histos[suff+'FitBandqt'+c] = self.histos[suff+'MC'+'qt'+c].Clone('FitBandqt'+c)
                self.histos[suff+'FitBandPDFqt'+c] = self.histos[suff+'MC'+'qt'+c].Clone('FitBandPDFqt'+c)
                self.histos[suff+'FitBandScaleqt'+c] = self.histos[suff+'MC'+'qt'+c].Clone('FitBandScaleqt'+c)
            else :
                self.histos[suff+'FitBand'+c] = self.histos[suff+'MC'+'mapTot'].Clone('FitBand'+c)
                self.histos[suff+'FitBandPDF'+c] = self.histos[suff+'MC'+'mapTot'].Clone('FitBandPDF'+c) 
                self.histos[suff+'FitBandScale'+c] = self.histos[suff+'MC'+'mapTot'].Clone('FitBandScale'+c) 
                
                self.histos[suff+'FitBandy'+c] = self.histos[suff+'MC'+'y'+'mapTot'].Clone('FitBandy'+c)
                self.histos[suff+'FitBandPDFy'+c] = self.histos[suff+'MC'+'y'+'mapTot'].Clone('FitBandPDFy'+c)
                self.histos[suff+'FitBandScaley'+c] = self.histos[suff+'MC'+'y'+'mapTot'].Clone('FitBandScaley'+c)
                
                self.histos[suff+'FitBandqt'+c] = self.histos[suff+'MC'+'qt'+'mapTot'].Clone('FitBandqt'+c)
                self.histos[suff+'FitBandPDFqt'+c] = self.histos[suff+'MC'+'qt'+'mapTot'].Clone('FitBandPDFqt'+c)
                self.histos[suff+'FitBandScaleqt'+c] = self.histos[suff+'MC'+'qt'+'mapTot'].Clone('FitBandScaleqt'+c)
                self.histos[suff+'FitBand'+c].Scale(1/35.9)
                self.histos[suff+'FitBandPDF'+c].Scale(1/35.9)
                self.histos[suff+'FitBandScale'+c].Scale(1/35.9)
                self.histos[suff+'FitBandy'+c].Scale(1/35.9)
                self.histos[suff+'FitBandPDFy'+c].Scale(1/35.9)
                self.histos[suff+'FitBandScaley'+c].Scale(1/35.9)
                self.histos[suff+'FitBandqt'+c].Scale(1/35.9)
                self.histos[suff+'FitBandPDFqt'+c].Scale(1/35.9)
                self.histos[suff+'FitBandScaleqt'+c].Scale(1/35.9)


            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                
                    #debug unclousure
                    if abs(self.histos[suff+'FitBand'+c].GetBinContent(i,j)-self.histos[suff+'FitAC'+c].GetBinContent(i,j))/self.histos[suff+'FitBand'+c].GetBinContent(i,j)>0.0000001 :
                            print "not clousure of", c, i, j , ",   (fitted-mc)/fitted=", (self.histos[suff+'FitBand'+c].GetBinContent(i,j)-self.histos[suff+'FitAC'+c].GetBinContent(i,j))/self.histos[suff+'FitAC'+c].GetBinContent(i,j)
                    if i==1 :
                        if abs(self.histos[suff+'FitBandqt'+c].GetBinContent(j)-self.histos[suff+'FitACqt'+c].GetBinContent(j))/self.histos[suff+'FitBandqt'+c].GetBinContent(j)>0.0000001 :
                            print "not clousure of", c, j , "(qT),   (fitted-mc)/fitted=", (self.histos[suff+'FitBandqt'+c].GetBinContent(j)-self.histos[suff+'FitACqt'+c].GetBinContent(j))/self.histos[suff+'FitBandqt'+c].GetBinContent(j)
                    if j==1 : 
                        if abs(self.histos[suff+'FitBandy'+c].GetBinContent(i)-self.histos[suff+'FitACy'+c].GetBinContent(i))/self.histos[suff+'FitBandqt'+c].GetBinContent(i)>0.0000001 :
                            print "not clousure of", c, i , "(y),   (fitted-mc)/fitted=", (self.histos[suff+'FitBandy'+c].GetBinContent(i)-self.histos[suff+'FitACy'+c].GetBinContent(i))/self.histos[suff+'FitBandy'+c].GetBinContent(i)
                    
                    errPDF = 0.
                    if i==1 : errPDFqt = 0.
                    if j==1 : errPDFy = 0.
                    MCVal = self.histos[suff+'FitBand'+c].GetBinContent(i,j) #like MC, but already lumi scaled
                    # if 'unpol' in c:
                    #     MCVal=MCVal/35.9
                    if i==1 : MCValqt = self.histos[suff+'FitBandqt'+c].GetBinContent(j)
                    if j==1 : MCValy = self.histos[suff+'FitBandy'+c].GetBinContent(i)
                    
                    for sName in self.systDict['_LHEPdfWeight']:
                        if 'unpol' in c:
                            systVal=self.histos[suff+'MC'+sName+'mapTot'].GetBinContent(i,j)/35.9
                            if i==1 : systValqt=self.histos[suff+'MCqt'+sName+'mapTot'].GetBinContent(j)/35.9
                            if j==1 : systValy=self.histos[suff+'MCy'+sName+'mapTot'].GetBinContent(i)/35.9
                        else:
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
                            if sNameDen!=sName and not UNCORR : continue
                            if sNameDen!=sName and 'unpol' in c : continue
                            if sName=='_nom' and sNameDen=='_nom' : continue
                            if ([sName,sNameDen] in self.vetoScaleList) : continue  #extremal cases
                            if 'unpol' in c:
                                systVal=self.histos[suff+'MC'+sName+'mapTot'].GetBinContent(i,j)/35.9
                                if i==1 : systValqt =self.histos[suff+'MCqt'+sName+'mapTot'].GetBinContent(j)/35.9
                                if j==1 : systValy =self.histos[suff+'MCy'+sName+'mapTot'].GetBinContent(i)/35.9
                            else:
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
            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): 
                for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1):
                    if not "unpol" in c:   
                        valCentral = self.histos[suff+'MC'+c].GetBinContent(i,j)
                    else :
                        valCentral = self.histos[suff+'MC'+'mapTot'].GetBinContent(i,j)
                        valCentral= valCentral/35.9
                    self.histos[suff+'FitRatioAC'+c].SetBinContent(i,j, self.histos[suff+'FitRatioAC'+c].GetBinContent(i,j)/valCentral)
                    self.histos[suff+'FitRatio'+c].SetBinContent(i,j, self.histos[suff+'FitRatio'+c].GetBinContent(i,j)/valCentral)
                    self.histos[suff+'FitRatioPDF'+c].SetBinContent(i,j,self.histos[suff+'FitRatioPDF'+c].GetBinContent(i,j)/valCentral)
                    self.histos[suff+'FitRatioScale'+c].SetBinContent(i,j,self.histos[suff+'FitRatioScale'+c].GetBinContent(i,j)/valCentral)
                    self.histos[suff+'FitRatioAC'+c].SetBinError(i,j,self.histos[suff+'FitRatioAC'+c].GetBinError(i,j)/valCentral)
                    self.histos[suff+'FitRatio'+c].SetBinError(i,j,self.histos[suff+'FitRatio'+c].GetBinError(i,j)/valCentral)
                    self.histos[suff+'FitRatioPDF'+c].SetBinError(i,j,self.histos[suff+'FitRatioPDF'+c].GetBinError(i,j)/valCentral)
                    self.histos[suff+'FitRatioScale'+c].SetBinError(i,j,self.histos[suff+'FitRatioScale'+c].GetBinError(i,j)/valCentral) 
                
            
            #--------------- build the ratio plots  Y trends -----------------------#
            self.histos[suff+'FitRatioACy'+c] = self.histos[suff+'FitAC'+'y'+c].Clone('FitRatioACy'+c)
            self.histos[suff+'FitRatioy'+c] = self.histos[suff+'FitBand'+'y'+c].Clone('FitRatioy'+c)
            self.histos[suff+'FitRatioPDFy'+c] = self.histos[suff+'FitBandPDF'+'y'+c].Clone('FitRatioPDFy'+c)
            self.histos[suff+'FitRatioScaley'+c] = self.histos[suff+'FitBandScale'+'y'+c].Clone('FitRatioScaley'+c)
            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): 
                if not "unpol" in c:   
                    valCentral = self.histos[suff+'MC'+'y'+c].GetBinContent(i)
                else :
                    valCentral = self.histos[suff+'MC'+'y'+'mapTot'].GetBinContent(i)
                    valCentral= valCentral/35.9
                self.histos[suff+'FitRatioAC'+'y'+c].SetBinContent(i, self.histos[suff+'FitRatioAC'+'y'+c].GetBinContent(i)/valCentral)
                self.histos[suff+'FitRatio'+'y'+c].SetBinContent(i, self.histos[suff+'FitRatio'+'y'+c].GetBinContent(i)/valCentral)
                self.histos[suff+'FitRatioPDF'+'y'+c].SetBinContent(i,self.histos[suff+'FitRatioPDF'+'y'+c].GetBinContent(i)/valCentral)
                self.histos[suff+'FitRatioScale'+'y'+c].SetBinContent(i,self.histos[suff+'FitRatioScale'+'y'+c].GetBinContent(i)/valCentral)
                self.histos[suff+'FitRatioAC'+'y'+c].SetBinError(i,self.histos[suff+'FitRatioAC'+'y'+c].GetBinError(i)/valCentral)
                self.histos[suff+'FitRatio'+'y'+c].SetBinError(i,self.histos[suff+'FitRatio'+'y'+c].GetBinError(i)/valCentral)
                self.histos[suff+'FitRatioPDF'+'y'+c].SetBinError(i,self.histos[suff+'FitRatioPDF'+'y'+c].GetBinError(i)/valCentral)
                self.histos[suff+'FitRatioScale'+'y'+c].SetBinError(i,self.histos[suff+'FitRatioScale'+'y'+c].GetBinError(i,j)/valCentral) 
            
            
            #--------------- build the ratio plots  qt trends -----------------------#
            self.histos[suff+'FitRatioACqt'+c] = self.histos[suff+'FitAC'+'qt'+c].Clone('FitRatioACqt'+c)
            self.histos[suff+'FitRatioqt'+c] = self.histos[suff+'FitBand'+'qt'+c].Clone('FitRatioqt'+c)
            self.histos[suff+'FitRatioPDFqt'+c] = self.histos[suff+'FitBandPDF'+'qt'+c].Clone('FitRatioPDFqt'+c)
            self.histos[suff+'FitRatioScaleqt'+c] = self.histos[suff+'FitBandScale'+'qt'+c].Clone('FitBandScaleqt'+c)
            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): 
                if not "unpol" in c:   
                    valCentral = self.histos[suff+'MC'+'qt'+c].GetBinContent(i)
                else :
                    valCentral = self.histos[suff+'MC'+'qt'+'mapTot'].GetBinContent(i)
                    valCentral= valCentral/35.9
                self.histos[suff+'FitRatioAC'+'qt'+c].SetBinContent(i, self.histos[suff+'FitRatioAC'+'qt'+c].GetBinContent(i)/valCentral)
                self.histos[suff+'FitRatio'+'qt'+c].SetBinContent(i, self.histos[suff+'FitRatio'+'qt'+c].GetBinContent(i)/valCentral)
                self.histos[suff+'FitRatioPDF'+'qt'+c].SetBinContent(i,self.histos[suff+'FitRatioPDF'+'qt'+c].GetBinContent(i)/valCentral)
                self.histos[suff+'FitRatioScale'+'qt'+c].SetBinContent(i,self.histos[suff+'FitRatioScale'+'qt'+c].GetBinContent(i)/valCentral)
                self.histos[suff+'FitRatioAC'+'qt'+c].SetBinError(i,self.histos[suff+'FitRatioAC'+'qt'+c].GetBinError(i)/valCentral)
                self.histos[suff+'FitRatio'+'qt'+c].SetBinError(i,self.histos[suff+'FitRatio'+'qt'+c].GetBinError(i)/valCentral)
                self.histos[suff+'FitRatioPDF'+'qt'+c].SetBinError(i,self.histos[suff+'FitRatioPDF'+'qt'+c].GetBinError(i)/valCentral)
                self.histos[suff+'FitRatioScale'+'qt'+c].SetBinError(i,self.histos[suff+'FitRatioScale'+'qt'+c].GetBinError(i,j)/valCentral) 
        
                
                    
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
                            nameX = 'y_'+str(self.yArr.index(y1)+1)+'_qt_'+str(self.qtArr.index(q1)+1)+'_'+c
                            nameY = 'y_'+str(self.yArr.index(y2)+1)+'_qt_'+str(self.qtArr.index(q2)+1)+'_'+c
                            corrVal = self.histos[suff+'corrMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                            covVal = self.histos[suff+'covMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                            self.histos[suff+'corrMat'+c].SetBinContent(indexUNR_X,indexUNR_Y,corrVal)
                            self.histos[suff+'covMat'+c].SetBinContent(indexUNR_X,indexUNR_Y,covVal)
                            
        for i in range(1, self.histos[suff+'FitACA0'].GetNbinsX()+1): #loop over y bins
            self.histos[suff+'corrMat'+'y'+str(i)] = ROOT.TH2D(suff+'corrMat_y{c}'.format(c=i), suff+'corrMat_y{c}'.format(c=i), len(self.unrolledAQt)-1, array('f',self.unrolledAQt), len(self.unrolledAQt)-1, array('f',self.unrolledAQt))
            self.histos[suff+'covMat'+'y'+str(i)] = ROOT.TH2D(suff+'covMat_y{c}'.format(c=i), suff+'covMat_y{c}'.format(c=i), len(self.unrolledAQt)-1, array('f',self.unrolledAQt), len(self.unrolledAQt)-1, array('f',self.unrolledAQt))
            for q1 in self.qtArr[:-1]: 
                for c1 in self.coeffArr[:-1] :
                    for q2 in self.qtArr[:-1]: 
                        for c2 in self.coeffArr[:-1] : 
                            indexUNR_X = self.coeffArr.index(float(c1))*(len(self.qtArr)-1)+self.qtArr.index(float(q1))
                            indexUNR_Y = self.coeffArr.index(float(c2))*(len(self.qtArr)-1)+self.qtArr.index(float(q2))
                            nameX = 'y_'+str(i)+'_qt_'+str(self.qtArr.index(q1)+1)+'_'+self.coeffList[self.coeffArr.index(c1)]
                            nameY = 'y_'+str(i)+'_qt_'+str(self.qtArr.index(q2)+1)+'_'+self.coeffList[self.coeffArr.index(c2)]
                            corrVal = self.histos[suff+'corrMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                            covVal = self.histos[suff+'covMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                            self.histos[suff+'corrMat'+'y'+str(i)].SetBinContent(indexUNR_X,indexUNR_Y,corrVal)
                            self.histos[suff+'covMat'+'y'+str(i)].SetBinContent(indexUNR_X,indexUNR_Y,covVal)
        
        for j in range(1, self.histos[suff+'FitACA0'].GetNbinsY()+1): #loop over qt bins
            self.histos[suff+'corrMat'+'qt'+str(j)] = ROOT.TH2D(suff+'corrMat_qt{c}'.format(c=j), suff+'corrMat_qt{c}'.format(c=j), len(self.unrolledAY)-1, array('f',self.unrolledAY), len(self.unrolledAY)-1, array('f',self.unrolledAY))
            self.histos[suff+'covMat'+'qt'+str(j)] = ROOT.TH2D(suff+'covMat_qt{c}'.format(c=j), suff+'covMat_qt{c}'.format(c=j), len(self.unrolledAY)-1, array('f',self.unrolledAY), len(self.unrolledAY)-1, array('f',self.unrolledAY))
            for y1 in self.yArr[:-1]: 
                for c1 in self.coeffArr[:-1] :
                    for y2 in self.yArr[:-1]: 
                        for c2 in self.coeffArr[:-1] : 
                            indexUNR_X = self.coeffArr.index(float(c1))*(len(self.yArr)-1)+self.yArr.index(float(y1))
                            indexUNR_Y = self.coeffArr.index(float(c2))*(len(self.yArr)-1)+self.yArr.index(float(y2))
                            nameX = 'y_'+str(self.yArr.index(y1)+1)+'_qt_'+str(j)+'_'+self.coeffList[self.coeffArr.index(c1)]
                            nameY = 'y_'+str(self.yArr.index(y2)+1)+'_qt_'+str(j)+'_'+self.coeffList[self.coeffArr.index(c2)]
                            corrVal = self.histos[suff+'corrMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                            covVal = self.histos[suff+'covMat'].GetBinContent(self.histos[suff+'corrMat'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'].GetYaxis().FindBin(nameY))
                            self.histos[suff+'corrMat'+'qt'+str(j)].SetBinContent(indexUNR_X,indexUNR_Y,corrVal)
                            self.histos[suff+'covMat'+'qt'+str(j)].SetBinContent(indexUNR_X,indexUNR_Y,covVal)        
        
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
                        nameX = 'y_'+str(self.yArr.index(y1)+1)+'_helmeta_'+self.coeffList[self.coeffArr.index(c1)]
                        nameY = 'y_'+str(self.yArr.index(y2)+1)+'_helmeta_'+self.coeffList[self.coeffArr.index(c2)]
                        corrVal = self.histos[suff+'corrMat'+'Integrated'].GetBinContent(self.histos[suff+'corrMat'+'Integrated'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'+'Integrated'].GetYaxis().FindBin(nameY))
                        covVal = self.histos[suff+'covMat'+'Integrated'].GetBinContent(self.histos[suff+'corrMat'+'Integrated'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'+'Integrated'].GetYaxis().FindBin(nameY))
                        self.histos[suff+'corrMat'+'Integrated'+'y'].SetBinContent(indexUNR_X,indexUNR_Y,corrVal)
                        self.histos[suff+'covMat'+'Integrated'+'y'].SetBinContent(indexUNR_X,indexUNR_Y,covVal)
        
        
        self.histos[suff+'corrMat'+'Integrated'+'qt'] = ROOT.TH2D(suff+'corrMat_Integrated_qt', suff+'corrMat_Integrated_qt', len(self.unrolledAQt)-1, array('f',self.unrolledAQt), len(self.unrolledAQt)-1, array('f',self.unrolledAQt))
        self.histos[suff+'covMat'+'Integrated'+'qt'] = ROOT.TH2D(suff+'covMat_Integrated_qt', suff+'covMat_Integrated_qt', len(self.unrolledAQt)-1, array('f',self.unrolledAQt), len(self.unrolledAQt)-1, array('f',self.unrolledAQt))
        for q1 in self.qtArr[:-1]: 
            for c1 in self.coeffArr[:-1] :
                for q2 in self.qtArr[:-1]: 
                    for c2 in self.coeffArr[:-1] : 
                        indexUNR_X = self.coeffArr.index(float(c1))*(len(self.qtArr)-1)+self.qtArr.index(float(q1))
                        indexUNR_Y = self.coeffArr.index(float(c2))*(len(self.qtArr)-1)+self.qtArr.index(float(q2))
                        nameX = 'qt_'+str(self.qtArr.index(q1)+1)+'_helmeta_'+self.coeffList[self.coeffArr.index(c1)]
                        nameY = 'qt_'+str(self.qtArr.index(q2)+1)+'_helmeta_'+self.coeffList[self.coeffArr.index(c2)]
                        corrVal = self.histos[suff+'corrMat'+'Integrated'].GetBinContent(self.histos[suff+'corrMat'+'Integrated'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'+'Integrated'].GetYaxis().FindBin(nameY))
                        covVal = self.histos[suff+'covMat'+'Integrated'].GetBinContent(self.histos[suff+'corrMat'+'Integrated'].GetXaxis().FindBin(nameX),self.histos[suff+'corrMat'+'Integrated'].GetYaxis().FindBin(nameY))
                        self.histos[suff+'corrMat'+'Integrated'+'qt'].SetBinContent(indexUNR_X,indexUNR_Y,corrVal)
                        self.histos[suff+'covMat'+'Integrated'+'qt'].SetBinContent(indexUNR_X,indexUNR_Y,covVal)
        
          
        #----------------------------------- CANVAS PREPARATION ---------------------------------------------------------------------------------          
            
        #--------------- Ai 1D canvas -------------------#
        for c in self.coeffDict:
            for i in range(1, self.histos[suff+'FitAC'+c].GetNbinsX()+1): #loop over rapidity bins
                self.canvas[suff+'FitAC'+'qt'+str(i)+c] = ROOT.TCanvas(suff+"_c_projQt{}_{}".format(i,c),suff+"_c_projQt{}_{}".format(i,c))
                self.canvas[suff+'FitAC'+'qt'+str(i)+c].SetGridx()
                self.canvas[suff+'FitAC'+'qt'+str(i)+c].SetGridy()
                self.histos[suff+'FitAC'+'qt'+str(i)+c] = self.histos[suff+'FitAC'+c].ProjectionY("projQt{}_{}".format(i,c),i,i)
                self.histos[suff+'FitBand'+'qt'+str(i)+c] = self.histos[suff+'FitBand'+c].ProjectionY("projQt{}_{}_err".format(i,c),i,i)
                self.histos[suff+'FitBandPDF'+'qt'+str(i)+c] = self.histos[suff+'FitBandPDF'+c].ProjectionY("projQt{}_{}_errPDF".format(i,c),i,i)
                self.histos[suff+'FitBandScale'+'qt'+str(i)+c] = self.histos[suff+'FitBandScale'+c].ProjectionY("projQt{}_{}_errScale".format(i,c),i,i)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].SetTitle(suff+' '+c+" vs W transverse momentum, "+str(self.yArr[i-1])+"<Y_{W}<"+str(self.yArr[i]))
                self.canvas[suff+'FitAC'+'qt'+str(i)+c].cd()
                self.histos[suff+'FitAC'+'qt'+str(i)+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
                self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetTitle(c)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].SetLineWidth(2)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].SetStats(0)
                self.histos[suff+'FitAC'+'qt'+str(i)+c].Draw()
                if not 'unpol' in c :
                    self.histos[suff+'FitAC'+'qt'+str(i)+c].GetYaxis().SetRangeUser(-4,4)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].SetFillColor(ROOT.kOrange-3)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].SetFillStyle(0)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].SetLineColor(ROOT.kOrange-3)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].SetLineWidth(2)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].DrawCopy("E2 same")
                self.histos[suff+'FitBand'+'qt'+str(i)+c].SetFillStyle(3001)
                self.histos[suff+'FitBand'+'qt'+str(i)+c].Draw("E2 same")
                
                #Additional block for scale-pdf splitted
                self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].SetFillStyle(0)
                self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
                self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].SetLineWidth(2)
                self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].DrawCopy("E2 same")
                self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].SetFillStyle(3004)
                self.histos[suff+'FitBandPDF'+'qt'+str(i)+c].Draw("E2 same")
                self.histos[suff+'FitBandScale'+'qt'+str(i)+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])
                self.histos[suff+'FitBandScale'+'qt'+str(i)+c].SetFillStyle(0)
                self.histos[suff+'FitBandScale'+'qt'+str(i)+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
                self.histos[suff+'FitBandScale'+'qt'+str(i)+c].SetLineWidth(2)
                self.histos[suff+'FitBandScale'+'qt'+str(i)+c].DrawCopy("E2 same")
                self.histos[suff+'FitBandScale'+'qt'+str(i)+c].SetFillStyle(3005)
                self.histos[suff+'FitBandScale'+'qt'+str(i)+c].Draw("E2 same")
                
                self.histos[suff+'FitAC'+'qt'+str(i)+c].DrawCopy("same") #to have foreground
            for j in range(1, self.histos[suff+'FitAC'+c].GetNbinsY()+1): #loop over pt bins
                self.canvas[suff+'FitAC'+'y'+str(j)+c] = ROOT.TCanvas(suff+"_c_projY{}_{}".format(j,c),suff+"_c_projY{}_{}".format(j,c))
                self.canvas[suff+'FitAC'+'y'+str(j)+c].SetGridx()
                self.canvas[suff+'FitAC'+'y'+str(j)+c].SetGridy()
                self.histos[suff+'FitAC'+'y'+str(j)+c] = self.histos[suff+'FitAC'+c].ProjectionX("projY{}_{}".format(j,c),j,j)
                self.histos[suff+'FitBand'+'y'+str(j)+c] = self.histos[suff+'FitBand'+c].ProjectionX("projY{}_{}_err".format(j,c),j,j)
                self.histos[suff+'FitBandPDF'+'y'+str(j)+c] = self.histos[suff+'FitBandPDF'+c].ProjectionX("projY{}_{}_errPDF".format(j,c),j,j)
                self.histos[suff+'FitBandScale'+'y'+str(j)+c] = self.histos[suff+'FitBandScale'+c].ProjectionX("projY{}_{}_errScale".format(j,c),j,j)
                self.histos[suff+'FitAC'+'y'+str(j)+c].SetTitle(suff+' '+c+" vs W rapidity, "+str(self.qtArr[j-1])+"<q_{T}^{W}<"+str(self.qtArr[j]))
                self.canvas[suff+'FitAC'+'y'+str(j)+c].cd()
                self.histos[suff+'FitAC'+'y'+str(j)+c].GetXaxis().SetTitle('Y_{W}')
                self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetTitle(c)
                self.histos[suff+'FitAC'+'y'+str(j)+c].SetLineWidth(2)
                self.histos[suff+'FitAC'+'y'+str(j)+c].SetStats(0)
                self.histos[suff+'FitAC'+'y'+str(j)+c].Draw()
                if not 'unpol' in c :
                    self.histos[suff+'FitAC'+'y'+str(j)+c].GetYaxis().SetRangeUser(-4,4)
                self.histos[suff+'FitBand'+'y'+str(j)+c].SetFillColor(ROOT.kOrange-3)#kMagenta-7)
                self.histos[suff+'FitBand'+'y'+str(j)+c].SetFillStyle(0)
                self.histos[suff+'FitBand'+'y'+str(j)+c].SetLineColor(ROOT.kOrange-3)#kMagenta-7)
                self.histos[suff+'FitBand'+'y'+str(j)+c].SetLineWidth(2)
                self.histos[suff+'FitBand'+'y'+str(j)+c].DrawCopy("E2 same")
                self.histos[suff+'FitBand'+'y'+str(j)+c].SetFillStyle(3001)
                self.histos[suff+'FitBand'+'y'+str(j)+c].Draw("E2 same")
                
                #Additional block for scale-pdf splitted
                self.histos[suff+'FitBandPDF'+'y'+str(j)+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
                self.histos[suff+'FitBandPDF'+'y'+str(j)+c].SetFillStyle(0)
                self.histos[suff+'FitBandPDF'+'y'+str(j)+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
                self.histos[suff+'FitBandPDF'+'y'+str(j)+c].SetLineWidth(2)
                self.histos[suff+'FitBandPDF'+'y'+str(j)+c].DrawCopy("E2 same")
                self.histos[suff+'FitBandPDF'+'y'+str(j)+c].SetFillStyle(3001)
                self.histos[suff+'FitBandPDF'+'y'+str(j)+c].Draw("E2 same")
                self.histos[suff+'FitBandScale'+'y'+str(j)+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
                self.histos[suff+'FitBandScale'+'y'+str(j)+c].SetFillStyle(0)
                self.histos[suff+'FitBandScale'+'y'+str(j)+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
                self.histos[suff+'FitBandScale'+'y'+str(j)+c].SetLineWidth(2)
                self.histos[suff+'FitBandScale'+'y'+str(j)+c].DrawCopy("E2 same")
                self.histos[suff+'FitBandScale'+'y'+str(j)+c].SetFillStyle(3001)
                self.histos[suff+'FitBandScale'+'y'+str(j)+c].Draw("E2 same")
                
                self.histos[suff+'FitAC'+'y'+str(j)+c].DrawCopy("same") #to have foreground        
        
          

            
            #-------------- unrolled: qt large, y small ----------------- #
            self.histos[suff+'FitAC'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_UNRqty'+c,suff+'_coeff_UNRqty'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            self.histos[suff+'FitBand'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_UNRqty_syst'+c,suff+'_coeff_UNRqty_syst'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            self.histos[suff+'FitBandScale'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_UNRqty_systPDF'+c,suff+'_coeff_UNRqty_systPDF'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            self.histos[suff+'FitBandPDF'+'UNRqty'+c] = ROOT.TH1F(suff+'_coeff_UNRqty_systScale'+c,suff+'_coeff_UNRqty_systScale'+c,len(self.unrolledQtY)-1, array('f',self.unrolledQtY))
            for q in self.qtArr[:-1] :
                for y in self.yArr[:-1] :
                    indexUNRqty = self.qtArr.index(float(q))*(len(self.yArr)-1)+self.yArr.index(float(y))
                    self.histos[suff+'FitAC'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitAC'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitAC'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitAC'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBand'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitBand'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBand'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitBand'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    
                    self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitBandPDF'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitBandPDF'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBandScale'+'UNRqty'+c].SetBinContent(indexUNRqty+1,self.histos[suff+'FitBandScale'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBandScale'+'UNRqty'+c].SetBinError(indexUNRqty+1,self.histos[suff+'FitBandScale'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    
                    # if self.yArr.index(y)==len(self.yArr)-2 :
                        # self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetBinLabel(indexUNRqty+1,"q_{T}="+str(q)+", Y="+str(y))
                    if self.yArr.index(y)==0 :
                        self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetBinLabel(indexUNRqty+1,"q_{T}^{W}#in[%.0f,%.0f]" % (q, self.qtArr[self.qtArr.index(q)+1]))
                        self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().ChangeLabel(indexUNRqty+1,340,0.03)
                    # else :
                    #     self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetBinLabel(indexUNRqty+1,str(y))        
            self.histos[suff+'FitAC'+'UNRqty'+c].SetTitle(suff+' '+c+", unrolled q_{T}(Y) ")
            # self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetTitle('q_{T}^{W} large, Y_{W} small')
            self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetTitle('fine binning: Y_{W} 0#rightarrow 2.4')
            self.histos[suff+'FitAC'+'UNRqty'+c].GetXaxis().SetTitleOffset(1.45)
            self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetTitle(c)
            if not 'unpol' in c :
                self.histos[suff+'FitAC'+'UNRqty'+c].GetYaxis().SetRangeUser(-4,4)
            self.histos[suff+'FitAC'+'UNRqty'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'UNRqty'+c].SetStats(0)
            self.histos[suff+'FitBand'+'UNRqty'+c].SetFillColor(ROOT.kOrange-3)#kMagenta-7)
            self.histos[suff+'FitBand'+'UNRqty'+c].SetFillStyle(0)
            self.histos[suff+'FitBand'+'UNRqty'+c].SetLineColor(ROOT.kOrange-3)#kMagenta-7)
            self.histos[suff+'FitBand'+'UNRqty'+c].SetLineWidth(2)
            
            self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
            self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetFillStyle(0)
            self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
            self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetLineWidth(2)    
            self.histos[suff+'FitBandScale'+'UNRqty'+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
            self.histos[suff+'FitBandScale'+'UNRqty'+c].SetFillStyle(0)
            self.histos[suff+'FitBandScale'+'UNRqty'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
            self.histos[suff+'FitBandScale'+'UNRqty'+c].SetLineWidth(2)
            
            self.canvas[suff+'FitAC'+'UNRqty'+c] = ROOT.TCanvas(suff+'_c_UNRqty'+c,suff+'_c_UNRqty'+c,800,600)
            self.canvas[suff+'FitAC'+'UNRqty'+c].cd()
            self.canvas[suff+'FitAC'+'UNRqty'+c].SetGridx()
            self.canvas[suff+'FitAC'+'UNRqty'+c].SetGridy()
            self.histos[suff+'FitAC'+'UNRqty'+c].Draw()
            self.histos[suff+'FitBand'+'UNRqty'+c].DrawCopy('E2 same')
            self.histos[suff+'FitBand'+'UNRqty'+c].SetFillStyle(3001)
            self.histos[suff+'FitBand'+'UNRqty'+c].Draw("E2 same")
            
            self.histos[suff+'FitBandPDF'+'UNRqty'+c].DrawCopy('E2 same')
            self.histos[suff+'FitBandPDF'+'UNRqty'+c].SetFillStyle(3004)
            self.histos[suff+'FitBandPDF'+'UNRqty'+c].Draw("E2 same")
            self.histos[suff+'FitBandScale'+'UNRqty'+c].DrawCopy('E2 same')
            self.histos[suff+'FitBandScale'+'UNRqty'+c].SetFillStyle(3005)
            self.histos[suff+'FitBandScale'+'UNRqty'+c].Draw("E2 same")
            
            self.histos[suff+'FitAC'+'UNRqty'+c].DrawCopy("same") #to have foreground                    



            #------------- unrolled: y large, qt small -----------------                    
            self.histos[suff+'FitAC'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_UNRyqt'+c,suff+'_coeff_UNRyqt'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            self.histos[suff+'FitBand'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_UNRyqt_syst'+c,suff+'_coeff_UNRyqt_syst'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            self.histos[suff+'FitBandPDF'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_UNRyqt_systPDF'+c,suff+'_coeff_UNRyqt_systPDF'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            self.histos[suff+'FitBandScale'+'UNRyqt'+c] = ROOT.TH1F(suff+'_coeff_UNRyqt_systScale'+c,suff+'_coeff_UNRyqt_systScale'+c,len(self.unrolledYQt)-1, array('f',self.unrolledYQt))
            for y in self.yArr[:-1] :
                for q in self.qtArr[:-1] :
                    indexUNRyqt = self.yArr.index(float(y))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                    self.histos[suff+'FitAC'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitAC'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitAC'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitAC'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBand'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitBand'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBand'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitBand'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    
                    self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitBandPDF'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitBandPDF'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetBinContent(indexUNRyqt+1,self.histos[suff+'FitBandScale'+c].GetBinContent(self.yArr.index(y)+1,self.qtArr.index(q)+1))
                    self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetBinError(indexUNRyqt+1,self.histos[suff+'FitBandScale'+c].GetBinError(self.yArr.index(y)+1,self.qtArr.index(q)+1))                    
                    
                    if self.qtArr.index(q)==0 :
                        # self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetBinLabel(indexUNRyqt+1,"Y_{W}=%.1f" % y)
                        self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetBinLabel(indexUNRyqt+1,"Y_{W}#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
                        # self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().ChangeLabel(indexUNRyqt+1,340,0.03)
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetTitle(suff+' '+c+", unrolled Y(q_{T}) ")
            # self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetTitle('Y_{W} large, q_{T}^{W}small')
            self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 32 GeV')
            self.histos[suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetTitleOffset(1.45)
            self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetTitle(c)
            if not 'unpol' in c :
                self.histos[suff+'FitAC'+'UNRyqt'+c].GetYaxis().SetRangeUser(-4,4)
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'UNRyqt'+c].SetStats(0)
            self.histos[suff+'FitBand'+'UNRyqt'+c].SetFillColor(ROOT.kOrange-3)
            self.histos[suff+'FitBand'+'UNRyqt'+c].SetFillStyle(0)
            self.histos[suff+'FitBand'+'UNRyqt'+c].SetLineColor(ROOT.kOrange-3)
            self.histos[suff+'FitBand'+'UNRyqt'+c].SetLineWidth(2)
            
            self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetFillStyle(0)
            self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetLineWidth(2)
            self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetFillStyle(0)
            self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetLineWidth(2)
            
            self.canvas[suff+'FitAC'+'UNRyqt'+c] = ROOT.TCanvas(suff+'_c_UNRyqt'+c,suff+'_c_UNRyqt'+c,800,600)
            self.canvas[suff+'FitAC'+'UNRyqt'+c].cd()
            self.canvas[suff+'FitAC'+'UNRyqt'+c].SetGridx()
            self.canvas[suff+'FitAC'+'UNRyqt'+c].SetGridy()
            self.histos[suff+'FitAC'+'UNRyqt'+c].Draw()
            self.histos[suff+'FitBand'+'UNRyqt'+c].DrawCopy('E2 same')
            self.histos[suff+'FitBand'+'UNRyqt'+c].SetFillStyle(3001)
            self.histos[suff+'FitBand'+'UNRyqt'+c].Draw("E2 same")
            
            self.histos[suff+'FitBandPDF'+'UNRyqt'+c].DrawCopy('E2 same')
            self.histos[suff+'FitBandPDF'+'UNRyqt'+c].SetFillStyle(3004)
            self.histos[suff+'FitBandPDF'+'UNRyqt'+c].Draw("E2 same")
            self.histos[suff+'FitBandScale'+'UNRyqt'+c].DrawCopy('E2 same')
            self.histos[suff+'FitBandScale'+'UNRyqt'+c].SetFillStyle(3005)
            self.histos[suff+'FitBandScale'+'UNRyqt'+c].Draw("E2 same")
            
            self.histos[suff+'FitAC'+'UNRyqt'+c].DrawCopy("same") #to have foreground
            
            
            #--------------------- Ai QT only canvas -------------------------#                
            self.canvas[suff+'FitAC'+'qt'+c] = ROOT.TCanvas(suff+"_c_QT_{}".format(c),suff+"_c_QT_{}".format(c))
            self.canvas[suff+'FitAC'+'qt'+c].SetGridx()
            self.canvas[suff+'FitAC'+'qt'+c].SetGridy()
            self.histos[suff+'FitAC'+'qt'+c].SetTitle(suff+' '+c+" vs W transverse momentum, Y integrated")
            self.canvas[suff+'FitAC'+'qt'+c].cd()
            self.histos[suff+'FitAC'+'qt'+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
            self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetTitle(c)
            self.histos[suff+'FitAC'+'qt'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'qt'+c].SetStats(0)
            self.histos[suff+'FitAC'+'qt'+c].Draw()
            if not 'unpol' in c :
                self.histos[suff+'FitAC'+'qt'+c].GetYaxis().SetRangeUser(-4,4)
            self.histos[suff+'FitBand'+'qt'+c].SetFillColor(ROOT.kOrange-3)
            self.histos[suff+'FitBand'+'qt'+c].SetFillStyle(0)
            self.histos[suff+'FitBand'+'qt'+c].SetLineColor(ROOT.kOrange-3)
            self.histos[suff+'FitBand'+'qt'+c].SetLineWidth(2)
            self.histos[suff+'FitBand'+'qt'+c].DrawCopy("E2 same")
            self.histos[suff+'FitBand'+'qt'+c].SetFillStyle(3001)
            self.histos[suff+'FitBand'+'qt'+c].Draw("E2 same")
            
            self.histos[suff+'FitBandPDF'+'qt'+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            self.histos[suff+'FitBandPDF'+'qt'+c].SetFillStyle(0)
            self.histos[suff+'FitBandPDF'+'qt'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            self.histos[suff+'FitBandPDF'+'qt'+c].SetLineWidth(2)
            self.histos[suff+'FitBandPDF'+'qt'+c].DrawCopy("E2 same")
            self.histos[suff+'FitBandPDF'+'qt'+c].SetFillStyle(3004)
            self.histos[suff+'FitBandPDF'+'qt'+c].Draw("E2 same")
            self.histos[suff+'FitBandScale'+'qt'+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            self.histos[suff+'FitBandScale'+'qt'+c].SetFillStyle(0)
            self.histos[suff+'FitBandScale'+'qt'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            self.histos[suff+'FitBandScale'+'qt'+c].SetLineWidth(2)
            self.histos[suff+'FitBandScale'+'qt'+c].DrawCopy("E2 same")
            self.histos[suff+'FitBandScale'+'qt'+c].SetFillStyle(3005)
            self.histos[suff+'FitBandScale'+'qt'+c].Draw("E2 same")
            
            self.histos[suff+'FitAC'+'qt'+c].DrawCopy("same") #to have foreground
            
            #--------------------- Ai Y only canvas -------------------------#          
            self.canvas[suff+'FitAC'+'y'+c] = ROOT.TCanvas(suff+"_c_Y_{}".format(c),suff+"_c_Y_{}".format(c))
            self.canvas[suff+'FitAC'+'y'+c].SetGridx()
            self.canvas[suff+'FitAC'+'y'+c].SetGridy()
            self.histos[suff+'FitAC'+'y'+c].SetTitle(suff+' '+c+" vs W rapidity,  q_{T} integrated")
            self.canvas[suff+'FitAC'+'y'+c].cd()
            self.histos[suff+'FitAC'+'y'+c].GetXaxis().SetTitle('Y_{W}')
            self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetTitle(c)
            self.histos[suff+'FitAC'+'y'+c].SetLineWidth(2)
            self.histos[suff+'FitAC'+'y'+c].SetStats(0)
            self.histos[suff+'FitAC'+'y'+c].Draw()
            if not 'unpol' in c :
                self.histos[suff+'FitAC'+'y'+c].GetYaxis().SetRangeUser(-4,4)
            self.histos[suff+'FitBand'+'y'+c].SetFillColor(ROOT.kOrange-3)#kMagenta-7)
            self.histos[suff+'FitBand'+'y'+c].SetFillStyle(0)
            self.histos[suff+'FitBand'+'y'+c].SetLineColor(ROOT.kOrange-3)#kMagenta-7)
            self.histos[suff+'FitBand'+'y'+c].SetLineWidth(2)
            self.histos[suff+'FitBand'+'y'+c].DrawCopy("E2 same")
            self.histos[suff+'FitBand'+'y'+c].SetFillStyle(3001)
            self.histos[suff+'FitBand'+'y'+c].Draw("E2 same")
            
            self.histos[suff+'FitBandPDF'+'y'+c].SetFillColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
            self.histos[suff+'FitBandPDF'+'y'+c].SetFillStyle(0)
            self.histos[suff+'FitBandPDF'+'y'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])#kMagenta-7)
            self.histos[suff+'FitBandPDF'+'y'+c].SetLineWidth(2)
            self.histos[suff+'FitBandPDF'+'y'+c].DrawCopy("E2 same")
            self.histos[suff+'FitBandPDF'+'y'+c].SetFillStyle(3004)
            self.histos[suff+'FitBandPDF'+'y'+c].Draw("E2 same")
            self.histos[suff+'FitBandScale'+'y'+c].SetFillColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
            self.histos[suff+'FitBandScale'+'y'+c].SetFillStyle(0)
            self.histos[suff+'FitBandScale'+'y'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])#kMagenta-7)
            self.histos[suff+'FitBandScale'+'y'+c].SetLineWidth(2)
            self.histos[suff+'FitBandScale'+'y'+c].DrawCopy("E2 same")
            self.histos[suff+'FitBandScale'+'y'+c].SetFillStyle(3005)
            self.histos[suff+'FitBandScale'+'y'+c].Draw("E2 same")
            
            self.histos[suff+'FitAC'+'y'+c].DrawCopy("same") #to have foreground
        
        
             #--------------- Ai 1D relative syst. band canvas -------------------#
            for i in range(1, self.histos[suff+'FitErr'+c].GetNbinsX()+1): #loop over rapidity bins
                self.canvas[suff+'FitErr'+'qt'+str(i)+c] = ROOT.TCanvas(suff+"_c_projQt{}_{}_Err".format(i,c),suff+"_c_ErrprojQt{}_{}_Err".format(i,c))
                self.canvas[suff+'FitErr'+'qt'+str(i)+c].SetGridx()
                self.canvas[suff+'FitErr'+'qt'+str(i)+c].SetGridy()
                self.leg[suff+'FitErr'+'qt'+str(i)+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
                self.histos[suff+'FitErr'+'qt'+str(i)+c] = self.histos[suff+'FitErr'+c].ProjectionY("projQt{}_{}_Err".format(i,c),i,i)
                self.histos[suff+'FitErrPDF'+'qt'+str(i)+c] = self.histos[suff+'FitErrPDF'+c].ProjectionY("projQt{}_{}_ErrPDF".format(i,c),i,i)
                self.histos[suff+'FitErrScale'+'qt'+str(i)+c] = self.histos[suff+'FitErrScale'+c].ProjectionY("projQt{}_{}_ErrScale".format(i,c),i,i)
                self.histos[suff+'FitErr'+'qt'+str(i)+c].SetTitle(suff+' '+c+" relative uncertainity vs W transverse momentum, "+str(self.yArr[i-1])+"<Y_{W}<"+str(self.yArr[i]))
                self.canvas[suff+'FitErr'+'qt'+str(i)+c].cd()
                self.histos[suff+'FitErr'+'qt'+str(i)+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
                self.histos[suff+'FitErr'+'qt'+str(i)+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
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
                self.canvas[suff+'FitErr'+'y'+str(j)+c] = ROOT.TCanvas(suff+"_c_projY{}_{}_Err".format(j,c),suff+"_c_projY{}_{}_Err".format(j,c))
                self.canvas[suff+'FitErr'+'y'+str(j)+c].SetGridx()
                self.canvas[suff+'FitErr'+'y'+str(j)+c].SetGridy()
                self.leg[suff+'FitErr'+'y'+str(j)+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
                self.histos[suff+'FitErr'+'y'+str(j)+c] = self.histos[suff+'FitErr'+c].ProjectionX("projY{}_{}".format(j,c),j,j)
                self.histos[suff+'FitErrPDF'+'y'+str(j)+c] = self.histos[suff+'FitErrPDF'+c].ProjectionY("projY{}_{}_ErrPDF".format(j,c),j,j)
                self.histos[suff+'FitErrScale'+'y'+str(j)+c] = self.histos[suff+'FitErrScale'+c].ProjectionY("projY{}_{}_ErrScale".format(j,c),j,j)
                self.histos[suff+'FitErr'+'y'+str(j)+c].SetTitle(suff+' '+c+"relative uncertainity vs W rapidity, "+str(self.qtArr[j-1])+"<q_{T}^{W}<"+str(self.qtArr[j]))
                self.canvas[suff+'FitErr'+'y'+str(j)+c].cd()
                self.histos[suff+'FitErr'+'y'+str(j)+c].GetXaxis().SetTitle('Y_{W}')
                self.histos[suff+'FitErr'+'y'+str(j)+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
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
            self.histos[suff+'FitErr'+'UNRqty'+c].SetTitle(suff+' '+c+", unrolled q_{T}(Y) ")
            self.histos[suff+'FitErr'+'UNRqty'+c].GetXaxis().SetTitle('fine binning: Y_{W} 0#rightarrow 2.4')
            self.histos[suff+'FitErr'+'UNRqty'+c].GetXaxis().SetTitleOffset(1.45)
            self.histos[suff+'FitErr'+'UNRqty'+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
            self.histos[suff+'FitErr'+'UNRqty'+c].SetLineWidth(3)
            self.histos[suff+'FitErrPDF'+'UNRqty'+c].SetLineWidth(3)
            self.histos[suff+'FitErrScale'+'UNRqty'+c].SetLineWidth(3)
            self.histos[suff+'FitErr'+'UNRqty'+c].SetStats(0)
            self.histos[suff+'FitErr'+'UNRqty'+c].SetLineColor(self.groupedSystColors['Nominal'][0])
            self.histos[suff+'FitErrPDF'+'UNRqty'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            self.histos[suff+'FitErrScale'+'UNRqty'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])

            self.canvas[suff+'FitErr'+'UNRqty'+c] = ROOT.TCanvas(suff+'_c_Err_UNRqty'+c,suff+'_c_Err_UNRqty'+c,800,600)
            self.leg[suff+'FitErr'+'UNRqty'+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
            self.canvas[suff+'FitErr'+'UNRqty'+c].cd()
            self.canvas[suff+'FitErr'+'UNRqty'+c].SetGridx()
            self.canvas[suff+'FitErr'+'UNRqty'+c].SetGridy()
            self.histos[suff+'FitErr'+'UNRqty'+c].GetYaxis().SetRangeUser(0.001,7000)
            self.canvas[suff+'FitErr'+'UNRqty'+c].SetLogy()
            self.histos[suff+'FitErr'+'UNRqty'+c].Draw()
            self.histos[suff+'FitErrPDF'+'UNRqty'+c].Draw("same")
            self.histos[suff+'FitErrScale'+'UNRqty'+c].Draw("same")
            self.leg[suff+'FitErr'+'UNRqty'+c].AddEntry(self.histos[suff+'FitErr'+'UNRqty'+c],self.groupedSystColors['Nominal'][1])
            self.leg[suff+'FitErr'+'UNRqty'+c].AddEntry(self.histos[suff+'FitErrPDF'+'UNRqty'+c],self.groupedSystColors['LHEPdfWeightVars'][1])
            self.leg[suff+'FitErr'+'UNRqty'+c].AddEntry(self.histos[suff+'FitErrScale'+'UNRqty'+c],self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitErr'+'UNRqty'+c].Draw("Same")



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
                        self.histos[suff+'FitErr'+'UNRyqt'+c].GetXaxis().SetBinLabel(indexUNRyqt+1,"Y_{W}#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
            self.histos[suff+'FitErr'+'UNRyqt'+c].SetTitle(suff+' '+c+", unrolled Y(q_{T}) ")
            self.histos[suff+'FitErr'+'UNRyqt'+c].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 32 GeV')
            self.histos[suff+'FitErr'+'UNRyqt'+c].GetXaxis().SetTitleOffset(1.45)
            self.histos[suff+'FitErr'+'UNRyqt'+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
            self.histos[suff+'FitErr'+'UNRyqt'+c].SetLineWidth(3)
            self.histos[suff+'FitErrPDF'+'UNRyqt'+c].SetLineWidth(3)
            self.histos[suff+'FitErrScale'+'UNRyqt'+c].SetLineWidth(3)
            self.histos[suff+'FitErr'+'UNRyqt'+c].SetStats(0)
            self.histos[suff+'FitErr'+'UNRyqt'+c].SetLineColor(self.groupedSystColors['Nominal'][0])
            self.histos[suff+'FitErrPDF'+'UNRyqt'+c].SetLineColor(self.groupedSystColors['LHEPdfWeightVars'][0])
            self.histos[suff+'FitErrScale'+'UNRyqt'+c].SetLineColor(self.groupedSystColors['LHEScaleWeightVars'][0])
            
            self.canvas[suff+'FitErr'+'UNRyqt'+c] = ROOT.TCanvas(suff+'_c_Err_UNRyqt'+c,suff+'_c_Err_UNRyqt'+c,800,600)
            self.leg[suff+'FitErr'+'UNRyqt'+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
            self.canvas[suff+'FitErr'+'UNRyqt'+c].cd()
            self.canvas[suff+'FitErr'+'UNRyqt'+c].SetGridx()
            self.canvas[suff+'FitErr'+'UNRyqt'+c].SetGridy()
            self.histos[suff+'FitErr'+'UNRyqt'+c].GetYaxis().SetRangeUser(0.001,7000)
            self.canvas[suff+'FitErr'+'UNRyqt'+c].SetLogy()
            self.histos[suff+'FitErr'+'UNRyqt'+c].Draw()
            self.histos[suff+'FitErrPDF'+'UNRyqt'+c].Draw("same")
            self.histos[suff+'FitErrScale'+'UNRyqt'+c].Draw("same")
            self.leg[suff+'FitErr'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitErr'+'UNRyqt'+c],self.groupedSystColors['Nominal'][1])
            self.leg[suff+'FitErr'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitErrPDF'+'UNRyqt'+c],self.groupedSystColors['LHEPdfWeightVars'][1])
            self.leg[suff+'FitErr'+'UNRyqt'+c].AddEntry(self.histos[suff+'FitErrScale'+'UNRyqt'+c],self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitErr'+'UNRyqt'+c].Draw("Same")
            
            
            #--------------- Ai QT only relative syst. band canvas -------------------#
            self.canvas[suff+'FitErr'+'qt'+c] = ROOT.TCanvas(suff+"_c_QT_{}_Err".format(c),suff+"_c_QT_{}_Err".format(c))
            self.canvas[suff+'FitErr'+'qt'+c].SetGridx()
            self.canvas[suff+'FitErr'+'qt'+c].SetGridy()
            self.leg[suff+'FitErr'+'qt'+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
            self.histos[suff+'FitErr'+'qt'+c].SetTitle(suff+' '+c+" relative uncertainity vs W transverse momentum, Y integrated")
            self.canvas[suff+'FitErr'+'qt'+c].cd()
            self.histos[suff+'FitErr'+'qt'+c].GetXaxis().SetTitle('q_{T}^{W} [GeV]')
            self.histos[suff+'FitErr'+'qt'+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
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
            self.histos[suff+'FitErrPDF'+'qt'+c].Draw("hist same")
            self.histos[suff+'FitErrScale'+'qt'+c].Draw("hist same")
            self.leg[suff+'FitErr'+'qt'+c].AddEntry(self.histos[suff+'FitErr'+'qt'+c],self.groupedSystColors['Nominal'][1])
            self.leg[suff+'FitErr'+'qt'+c].AddEntry(self.histos[suff+'FitErrPDF'+'qt'+c],self.groupedSystColors['LHEPdfWeightVars'][1])
            self.leg[suff+'FitErr'+'qt'+c].AddEntry(self.histos[suff+'FitErrScale'+'qt'+c],self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitErr'+'qt'+c].Draw("Same")
            
            #--------------- Ai Y only relative syst. band canvas -------------------#
            self.canvas[suff+'FitErr'+'y'+c] = ROOT.TCanvas(suff+"_c_Y_{}_Err".format(c),suff+"_c_Y_{}_Err".format(c))
            self.canvas[suff+'FitErr'+'y'+c].SetGridx()
            self.canvas[suff+'FitErr'+'y'+c].SetGridy()
            self.leg[suff+'FitErr'+'y'+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
            self.histos[suff+'FitErr'+'y'+c].SetTitle(suff+' '+c+"relative uncertainity vs W rapidity,  q_{T} integrated")
            self.canvas[suff+'FitErr'+'y'+c].cd()
            self.histos[suff+'FitErr'+'y'+c].GetXaxis().SetTitle('Y_{W}')
            self.histos[suff+'FitErr'+'y'+c].GetYaxis().SetTitle('#Delta'+c+'/'+c)
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
            self.histos[suff+'FitErrPDF'+'y'+c].Draw("hist same")
            self.histos[suff+'FitErrScale'+'y'+c].Draw("hist same")
            self.leg[suff+'FitErr'+'y'+c].AddEntry(self.histos[suff+'FitErr'+'y'+c],self.groupedSystColors['Nominal'][1])
            self.leg[suff+'FitErr'+'y'+c].AddEntry(self.histos[suff+'FitErrPDF'+'y'+c],self.groupedSystColors['LHEPdfWeightVars'][1])
            self.leg[suff+'FitErr'+'y'+c].AddEntry(self.histos[suff+'FitErrScale'+'y'+c],self.groupedSystColors['LHEScaleWeightVars'][1])
            self.leg[suff+'FitErr'+'y'+c].Draw("Same")
        
        #---------------------------- Canvas correlation matrices ------------------------
        for mtx in ['corr','cov'] :
            if mtx == 'corr' : htitle = 'Correlation' 
            if mtx == 'cov' : htitle = 'Covariance'
            
            for c in self.coeffDict:
                self.canvas[suff+mtx+'Mat'+c] = ROOT.TCanvas(suff+'_'+mtx+'Mat_'+c,suff+'_'+mtx+'Mat_'+c,800,600)
                self.canvas[suff+mtx+'Mat'+c].cd()
                self.canvas[suff+mtx+'Mat'+c].SetGridx()
                self.canvas[suff+mtx+'Mat'+c].SetGridy()
                self.histos[suff+mtx+'Mat'+c].Draw("colz")
                self.histos[suff+mtx+'Mat'+c].SetStats(0)
                self.histos[suff+mtx+'Mat'+c].SetTitle(suff+ ' '+ htitle+' Matrix, '+c)
                self.histos[suff+mtx+'Mat'+c].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 32 GeV')
                self.histos[suff+mtx+'Mat'+c].GetXaxis().SetTitleOffset(1.45)
                if mtx == 'corr' : 
                    self.histos[suff+mtx+'Mat'+c].GetZaxis().SetCanExtend(1)
                    self.histos[suff+mtx+'Mat'+c].GetZaxis().SetRangeUser(-1,1)
                    
                for y in self.yArr[:-1] :
                    for q in self.qtArr[:-1] :
                        if self.qtArr.index(q)==0 :
                            indexUNR= self.yArr.index(float(y))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                            self.histos[suff+mtx+'Mat'+c].GetXaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+c].GetXaxis().SetBinLabel(indexUNR+1,"Y_{W}#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
                            self.histos[suff+mtx+'Mat'+c].GetXaxis().ChangeLabel(indexUNR+1,340,0.03)
                            self.histos[suff+mtx+'Mat'+c].GetYaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+c].GetYaxis().SetBinLabel(indexUNR+1,"Y_{W}#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))
                            self.histos[suff+mtx+'Mat'+c].GetYaxis().ChangeLabel(indexUNR+1,340,0.03)
                
            for i in range(1, self.histos[suff+'FitACA0'].GetNbinsX()+1): #loop over y bins
                self.canvas[suff+mtx+'Mat'+'y'+str(i)] = ROOT.TCanvas(suff+'_'+mtx+'Mat_'+'y'+str(i),suff+'_'+mtx+'Mat_'+'y'+str(i),800,600)
                self.canvas[suff+mtx+'Mat'+'y'+str(i)].cd()
                self.canvas[suff+mtx+'Mat'+'y'+str(i)].SetGridx()
                self.canvas[suff+mtx+'Mat'+'y'+str(i)].SetGridy()
                self.histos[suff+mtx+'Mat'+'y'+str(i)].Draw("colz")
                self.histos[suff+mtx+'Mat'+'y'+str(i)].SetStats(0)
                self.histos[suff+mtx+'Mat'+'y'+str(i)].SetTitle(suff+ ' '+ htitle+' Matrix, '+str(self.yArr[i-1])+"<Y_{W}<"+str(self.yArr[i]))
                self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 32 GeV')
                self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().SetTitleOffset(1.45)
                if mtx == 'corr' : 
                    self.histos[suff+mtx+'Mat'+'y'+str(i)].GetZaxis().SetCanExtend(1)
                    self.histos[suff+mtx+'Mat'+'y'+str(i)].GetZaxis().SetRangeUser(-1,1)
                for c in self.coeffArr[:-1] :
                    for q in self.qtArr[:-1] :
                        if self.qtArr.index(q)==0 :
                            indexUNR= self.coeffArr.index(float(c))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                            if 'unpol' not in self.coeffList[self.coeffArr.index(c)] :    
                                coeffName = self.coeffList[self.coeffArr.index(c)]  
                            else :
                                coeffName = 'unpol'
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().SetBinLabel(indexUNR+1,coeffName)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetXaxis().ChangeLabel(indexUNR+1,340,0.03)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetYaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetYaxis().SetBinLabel(indexUNR+1,coeffName)
                            self.histos[suff+mtx+'Mat'+'y'+str(i)].GetYaxis().ChangeLabel(indexUNR+1,340,0.03)
            
            
            for j in range(1, self.histos[suff+'FitACA0'].GetNbinsY()+1): #loop over qt bins
                self.canvas[suff+mtx+'Mat'+'qt'+str(j)] = ROOT.TCanvas(suff+'_'+mtx+'Mat_'+'qt'+str(j),suff+'_'+mtx+'Mat_'+'qt'+str(j),800,600)
                self.canvas[suff+mtx+'Mat'+'qt'+str(j)].cd()
                self.canvas[suff+mtx+'Mat'+'qt'+str(j)].SetGridx()
                self.canvas[suff+mtx+'Mat'+'qt'+str(j)].SetGridy()
                self.histos[suff+mtx+'Mat'+'qt'+str(j)].Draw("colz")
                self.histos[suff+mtx+'Mat'+'qt'+str(j)].SetStats(0)
                self.histos[suff+mtx+'Mat'+'qt'+str(j)].SetTitle(suff+ ' '+ htitle+' Matrix, '+str(self.qtArr[j-1])+"<q_{T}^{W}<"+str(self.qtArr[j]))
                self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().SetTitle('fine binning: Y_{W} 0#rightarrow 2.4')
                self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().SetTitleOffset(1.45)
                if mtx == 'corr' : 
                    self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetZaxis().SetCanExtend(1)
                    self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetZaxis().SetRangeUser(-1,1)
                for c in self.coeffArr[:-1] :
                    for y in self.yArr[:-1] :
                        if self.yArr.index(y)==0 :
                            indexUNR= self.coeffArr.index(float(c))*(len(self.yArr)-1)+self.yArr.index(float(y))
                            if 'unpol' not in self.coeffList[self.coeffArr.index(c)] :    
                                coeffName = self.coeffList[self.coeffArr.index(c)]  
                            else :
                                coeffName = 'unpol'
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().SetBinLabel(indexUNR+1,coeffName)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetXaxis().ChangeLabel(indexUNR+1,340,0.03)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetYaxis().SetNdivisions(-1)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetYaxis().SetBinLabel(indexUNR+1,coeffName)
                            self.histos[suff+mtx+'Mat'+'qt'+str(j)].GetYaxis().ChangeLabel(indexUNR+1,340,0.03)
            
            
            #---------------------------- Canvas integrated correlation matrices ------------------------ 
            self.canvas[suff+mtx+'Mat'+'Integrated'+'y'] = ROOT.TCanvas(suff+'_'+mtx+'Mat_'+'Integrated_'+'y',suff+'_'+mtx+'Mat_'+'Integrated_'+'y',800,600)
            self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].cd()
            self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].SetGridx()
            self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].SetGridy()
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].Draw("colz")
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].SetStats(0)
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].SetTitle(suff+ ' '+ htitle+' Matrix, q_{T} integrated')
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().SetTitle('fine binning: Y_{W} 0#rightarrow 2.4')
            self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().SetTitleOffset(1.45)
            if mtx == 'corr' : 
                self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetZaxis().SetCanExtend(1)
                self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetZaxis().SetRangeUser(-1,1)
            for c in self.coeffArr[:-1] :
                for y in self.yArr[:-1] :
                    if self.yArr.index(y)==0 :
                        indexUNR= self.coeffArr.index(float(c))*(len(self.yArr)-1)+self.yArr.index(float(y))
                        if 'unpol' not in self.coeffList[self.coeffArr.index(c)] :    
                            coeffName = self.coeffList[self.coeffArr.index(c)]  
                        else :
                            coeffName = 'unpol'
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().SetBinLabel(indexUNR+1,coeffName)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetXaxis().ChangeLabel(indexUNR+1,340,0.03)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetYaxis().SetNdivisions(-1)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetYaxis().SetBinLabel(indexUNR+1,coeffName)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'y'].GetYaxis().ChangeLabel(indexUNR+1,340,0.03)
            
                    
            self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'] = ROOT.TCanvas(suff+'_'+mtx+'Mat_'+'Integrated_'+'qt',suff+'_'+mtx+'Mat_'+'Integrated_'+'qt',800,600)
            self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].cd()
            self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].SetGridx()
            self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].SetGridy()
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].Draw("colz")
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].SetStats(0)
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].SetTitle(suff+ ' '+ htitle+' Matrix, Y integrated')
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 32 GeV')
            self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().SetTitleOffset(1.45)
            if mtx == 'corr' : 
                self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetZaxis().SetCanExtend(1)
                self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetZaxis().SetRangeUser(-1,1)
            for c in self.coeffArr[:-1] :
                for q in self.qtArr[:-1] :
                    if self.qtArr.index(q)==0 :
                        indexUNR= self.coeffArr.index(float(c))*(len(self.qtArr)-1)+self.qtArr.index(float(q))
                        if 'unpol' not in self.coeffList[self.coeffArr.index(c)] :    
                            coeffName = self.coeffList[self.coeffArr.index(c)]  
                        else :
                            coeffName = 'unpol'
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().SetNdivisions(-1)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().SetBinLabel(indexUNR+1,coeffName)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetXaxis().ChangeLabel(indexUNR+1,340,0.03)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetYaxis().SetBinLabel(indexUNR+1,coeffName)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetYaxis().SetNdivisions(-1)
                        self.histos[suff+mtx+'Mat'+'Integrated'+'qt'].GetYaxis().ChangeLabel(indexUNR+1,340,0.03)
        

            
        
        
    
    
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
                
                self.leg['comp'+'FitAC'+'qt'+str(i)+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
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
                
                self.leg['comp'+'FitAC'+'y'+str(j)+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
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
            self.leg['comp'+'FitAC'+'UNRqty'+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
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
            self.leg['comp'+'FitAC'+'UNRyqt'+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
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
            self.leg['comp'+'FitAC'+'qt'+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
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
            self.leg['comp'+'FitAC'+'y'+c] = ROOT.TLegend(0.6,0.75,0.9,0.9)
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
            

    def makeRootOutput(self,outFileName,SAVE, suffList, comparison) :

        outFile = ROOT.TFile(outFileName+".root", "recreate")
        outFile.cd()
        dirFinalDict = {}

        for suff in suffList : 
            dirFinalDict['coeff2D'+suff] = outFile.mkdir('coeff2D_'+suff)
            dirFinalDict['coeff2D'+suff].cd()
            for hel in self.hels:
                self.histos[suff+'FitHel'+hel].Write()
                self.histos[suff+'FitHel'+hel+'norm'].Write()
            for c in self.coeffDict:
                self.histos[suff+'FitAC'+c].Write()
                self.histos[suff+'FitACqt'+c].Write()
                self.histos[suff+'FitACy'+c].Write()
                self.histos[suff+'FitBand'+c].Write()
                    
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
                self.histos[suff+mtx+'Mat'].Write()
                for c in self.coeffDict:
                    self.canvas[suff+mtx+'Mat'+c].Write()
                for i in range(1, self.histos[suff+'FitACA0'].GetNbinsX()+1):
                    self.canvas[suff+mtx+'Mat'+'y'+str(i)].Write()
                for j in range(1, self.histos[suff+'FitACA0'].GetNbinsY()+1):
                    self.canvas[suff+mtx+'Mat'+'qt'+str(j)].Write()
                self.canvas[suff+mtx+'Mat'+'Integrated'+'qt'].Write()
                self.canvas[suff+mtx+'Mat'+'Integrated'+'y'].Write()
                
                  
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
        
        outFile.Close()
                
               
    def getCoeffDict(self) :
        return self.coeffDict
    
    def getYArr(self) :
        return self.yArr
    
    def getQtArr(self) :
        return self.qtArr
    
def saver(rootInput, suffList, comparison,coeffDict,yArr,qtArr) :
    
    if not os.path.exists(rootInput): os.system("mkdir -p " + rootInput)
    FitFile = ROOT.TFile.Open(rootInput+'.root')
        
    unrList = ['UNRqty','Err_UNRqty','UNRyqt','Err_UNRyqt']
    intList = ['QT','Y']

    if comparison :
        for c in coeffDict:
            for name in unrList :
                can = FitFile.Get('comparison_coeff_unrolled/comp_c_'+name+c)
                can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
                can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
    
    for suff in suffList : 
        
        if not comparison :
             for c in coeffDict:
                for name in unrList :
                    can = FitFile.Get('coeff_unrolled'+suff+'/'+suff+'_c_'+name+c)
                    can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
                    can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
                for name in intList :
                    for ee in ['','_Err'] :
                        can = FitFile.Get('coeff_'+suff+'/'+suff+'_c_'+name+c+ee)
                        can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
                        can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
            
        
        for mtx in ['corr','cov'] :
            for c in coeffDict:
                can = FitFile.Get('matrices_'+suff+'/'+suff+'_'+mtx+'Mat_'+c)
                can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
                can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
            for i in range(1, len(yArr)):
                can = FitFile.Get('matrices_'+suff+'/'+suff+'_'+mtx+'Mat_y'+str(i))
                can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
                can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
            for j in range(1, len(qtArr)):
                can = FitFile.Get('matrices_'+suff+'/'+suff+'_'+mtx+'Mat_qt'+str(j))
                can.SaveAs(rootInput+'/'+can.GetTitle()+'.pdf')
                can.SaveAs(rootInput+'/'+can.GetTitle()+'.png')
                
                
###############################################################################################################################
#
#  usage: python plotter_fitResult.py -o OUTNAME -f /scratchnvme/emanca/wproperties-analysis/Fit/fit_Wplus.root -u 1 -c 1 -s 1
#
#  all the parameters are described below, but some note:
#  --uncorrelate: it should be set to True accordingly to theory prescription (31 Scale variation in total)
#  --comparison: if True expected two files called --fitFile and --FiteFile_reco
#
################################################################################################################################
       
parser = argparse.ArgumentParser("")

parser.add_argument('-o','--output', type=str, default='fitResult',help="name of the output file")
parser.add_argument('-f','--fitFile', type=str, default='fit_Wplus.root',help="name of the fit result root file. If comparison active the name must be: NAME.root, NAME_reco.root")
parser.add_argument('-i','--input', type=str, default='../analysisOnGen/genInput_Wplus.root',help="name of the input root file")
parser.add_argument('-u','--uncorrelate', type=int, default=True,help="if true uncorrelate num and den of Angular Coeff in MC scale variation")
parser.add_argument('-c','--comparison', type=int, default=True,help="comparison between reco and gen fit")
parser.add_argument('-s','--save', type=int, default=False,help="save .png and .pdf canvas")
parser.add_argument('-l','--suffList', type=str, default='',nargs='*', help="list of suff to be processed in the form: gen,reco")

args = parser.parse_args()
OUTPUT = args.output
FITFILE = args.fitFile
INPUT = args.input
UNCORR = args.uncorrelate
COMP = args.comparison
SAVE= args.save
SUFFL =args.suffList


p=plotter()
p.AngCoeffPlots(inputFile=INPUT, fitFile=FITFILE, uncorrelate=UNCORR,suff=SUFFL[0])
if COMP :
    recoFit = FITFILE.replace('.root', '_'+str(SUFFL[1])+'.root')
    p.AngCoeffPlots(inputFile=INPUT, fitFile=recoFit, uncorrelate=UNCORR, suff=SUFFL[1])
    p.GenRecoComparison(suffGen='gen', suffReco='reco')
p.makeRootOutput(outFileName=OUTPUT, SAVE=SAVE,suffList=SUFFL,comparison=COMP)
if SAVE :
    saver(rootInput=OUTPUT,suffList=SUFFL,comparison=COMP,coeffDict=p.getCoeffDict(),yArr=p.getYArr(),qtArr=p.getQtArr())