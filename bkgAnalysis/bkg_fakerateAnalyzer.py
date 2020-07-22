import ROOT
from array import array
import math
import sys
import copy as copy
from scipy.special import erf
import numpy as np
import ctypes

import bkg_utils


class bkg_analyzer:
    def __init__(self, ptBinning, etaBinning, outdir='./bkg', inputDir='./data', systKind='Nominal',systName='',correlatedFit=False,statAna=False, nameSuff = '')  :

        self.outdir = outdir
        self.inputDir = inputDir
        self.ptBinning = ptBinning
        self.etaBinning = etaBinning
        self.nameSuff = nameSuff
        self.systKind = systKind
        self.systName = systName
        self.correlatedFit = correlatedFit
        self.statAna = statAna
        
        if self.correlatedFit :
            self.corrFitSuff = '_CF'  
        else :
            self.corrFitSuff = ''
        if self.statAna :
            self.corrFitSuff = self.corrFitSuff+'statAna'
        self.signList = ['Plus','Minus']
        self.sampleList =  ['WToMuNu','Data'] #WToMuNu= All EWK samples
        self.maxPt_linearFit = 55
        print "WARNING: pt fit range 25-55 GeV"
        self.rootFiles = []
        for f in self.sampleList : 
             self.rootFiles.append(ROOT.TFile.Open(self.inputDir+'/'+f+'.root'))    

        self.ptBinningS = ['{:.2g}'.format(x) for x in self.ptBinning[:-1]]
        self.etaBinningS = ['{:.2g}'.format(x) for x in self.etaBinning[:-1]]
        
        print "WARNING: hardcoded LHE dict"
        self.LHEdict = {
            'Down' : ["LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0", "LHEScaleWeight_muR1p0_muF0p5"],
            'Up' : ["LHEScaleWeight_muR2p0_muF2p0", "LHEScaleWeight_muR2p0_muF1p0","LHEScaleWeight_muR1p0_muF2p0"]   
        }
        
        
    def binNumb_calculator(self,histo, axis, lowEdge) : #axis can be X,Y,Z
        binout = 0
        if axis =='X' :
            for i in range(1,histo.GetNbinsX()+1) :
                if abs(histo.GetXaxis().GetBinLowEdge(i)-lowEdge)<0.00001 :
                    binout = i
                    return binout
        if axis =='Y' :
            for i in range(1,histo.GetNbinsY()+1) :
                if abs(histo.GetYaxis().GetBinLowEdge(i)-lowEdge)<0.00001 :
                    binout = i
                    return binout
        if axis =='Z' :
            for i in range(1,histo.GetNbinsZ()+1) :
                if abs(histo.GetYaxis().GetBinLowEdge(i)-lowEdge)<0.00001 :
                    binout = i
                    return binout
        return binout 
            

    def MinuitLinearFit(self, yy, xx, invCov, p0, q0,p0Err,q0Err,s,e) :

        def linearChi2(npar, gin, f, par, istatus ):
            vec = yy-(par[0]+par[1]*(xx-26))#SLOPEOFFSETUNCORR
            chi2 = np.linalg.multi_dot( [vec.T, invCov, vec] )
            f.value = chi2
            return

        minuit = ROOT.TMinuit(2) #ROOT.TVirtualFitter.Fitter(0, 2)
        minuit.SetFCN(linearChi2)
        arglist = array('d',2*[0])
        arglist[0] =1.0
        ierflg = ctypes.c_int(0)
        minuit.SetPrintLevel(-1)
        minuit.mnexcm("SET ERR",arglist,1,ierflg)
        minuit.mnparm(0, 'offset',q0,p0Err,-1,1,ierflg)
        minuit.mnparm(1, 'slope',p0,q0Err,-0.3,0.3,ierflg)          
        arglist[0] = 500000 #call limit
        # arglist[1] = 0.01 #tolerance
        minuit.mnexcm("MIGRAD", arglist, 2,ierflg)
        # minuit.mnexcm("SEEk", arglist, 2,ierflg)


        out_cov = ROOT.TMatrixDSym(2)
        minuit.mnemat(out_cov.GetMatrixArray(),2)
        val0 = ctypes.c_double(0.)
        err0 = ctypes.c_double(0.)
        val1 = ctypes.c_double(0.)
        err1 = ctypes.c_double(0.)
        minuit.GetParameter(0,val0,err0)
        minuit.GetParameter(1,val1,err1)

        outdict = {}
        outdict[s+e+'offset'+'Minuit'] = val0.value
        outdict[s+e+'slope'+'Minuit'] = val1.value
        outdict[s+e+'offset'+'Minuit'+'Err'] = err0.value
        outdict[s+e+'slope'+'Minuit'+'Err'] = err1.value

        outdict[s+e+'offset*slope'+'Minuit'] = ROOT.TMatrixDRow(out_cov,0)(1)
        return outdict


    def correlatedFitter(self, fakedict) :
       
        DONT_REBIN = False

        #load the histos
        file_dict = {}
        histo_fake_dict = {}
        systList = []
        systDict = bkg_utils.bkg_systematics
        bin4corFit =  [26,28,30,32,34,36,38,40,42,44,46,48,50,52,55] #LoreHistos
        #bin4corFit =  [25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55] #redesign
        
        if DONT_REBIN : 
            bin4corFit = self.ptBinning

        binChange = 13 #LoreHistos
        #binChange = 17 #redesign

        bin4corFitS = ['{:.2g}'.format(x) for x in bin4corFit[:-1]]

        for sKind, sList in systDict.iteritems():
            for sName in sList :
                systList.append(sName)
                systdir = self.outdir.replace("bkg_"+self.systName,"bkg_"+sName+'/')
                file_dict[sName]=ROOT.TFile.Open(systdir+'bkg_plots'+self.nameSuff+".root")

                for s in self.signList :
                    for e in self.etaBinningS :
                        histo_fake_dict[sName+s+e] = file_dict[sName].Get('c_comparison_'+s+'_'+e).GetPrimitive('hFakes_pt_fake_'+s+'_'+e)
                        temp_histREB = ROOT.TH1F(histo_fake_dict[sName+s+e].GetName()+'_reb',histo_fake_dict[sName+s+e].GetName()+'_reb',len(bin4corFit)-1,array('d',bin4corFit))                            
                        for r in range(1,len(bin4corFit)) :
                            if DONT_REBIN :
                                valueBinRebinned = histo_fake_dict[sName+s+e].GetBinContent(r)
                                valueBinRebinnedErr = histo_fake_dict[sName+s+e].GetBinError(r) 
                            elif r<binChange :
                                valueBinRebinned = (histo_fake_dict[sName+s+e].GetBinContent(2*r)+histo_fake_dict[sName+s+e].GetBinContent(2*r-1))/2
                                valueBinRebinnedErr = (histo_fake_dict[sName+s+e].GetBinError(2*r)+histo_fake_dict[sName+s+e].GetBinError(2*r-1))/2
                            else :
                                valueBinRebinned = (histo_fake_dict[sName+s+e].GetBinContent(3*r-binChange)+histo_fake_dict[sName+s+e].GetBinContent(3*r-binChange-1)+histo_fake_dict[sName+s+e].GetBinContent(3*r-binChange-2))/3 #the not-simplified value of the bin number is: 2*N+3*(r-N), -1, -2.
                                valueBinRebinnedErr = (histo_fake_dict[sName+s+e].GetBinError(3*r-binChange)+histo_fake_dict[sName+s+e].GetBinError(3*r-binChange-1)+histo_fake_dict[sName+s+e].GetBinError(3*r-binChange-2))/3 #the not-simplified value of the bin number is: 2*N+3*(r-N), -1, -2.
                            temp_histREB.SetBinContent(r,valueBinRebinned)
                            temp_histREB.SetBinError(r,valueBinRebinnedErr)
                        temp_histREB.AddDirectory(False)
                        histo_fake_dict[sName+s+e+'reb'] = temp_histREB.Clone()
                file_dict[sName].Close()
        
        #nominal one:
        sNameDir = ''
        sName = 'nom'
        systdir = self.outdir.replace("bkg_"+self.systName,"bkg_"+sNameDir+'/')
        file_dict[sName]=ROOT.TFile.Open(systdir+'bkg_plots'+self.nameSuff+".root")
        for s in self.signList :
            for e in self.etaBinningS :
                histo_fake_dict[sName+s+e] = file_dict[sName].Get('c_comparison_'+s+'_'+e).GetPrimitive('hFakes_pt_fake_'+s+'_'+e)
                temp_histREB = ROOT.TH1F(histo_fake_dict[sName+s+e].GetName()+'_reb',histo_fake_dict[sName+s+e].GetName()+'_reb',len(bin4corFit)-1,array('d',bin4corFit))
                for r in range(1,len(bin4corFit)) :
                    if DONT_REBIN :
                        valueBinRebinned = histo_fake_dict[sName+s+e].GetBinContent(r)
                        valueBinRebinnedErr = histo_fake_dict[sName+s+e].GetBinError(r) 
                    elif r<binChange :
                        valueBinRebinned = (histo_fake_dict[sName+s+e].GetBinContent(2*r)+histo_fake_dict[sName+s+e].GetBinContent(2*r-1))/2
                        valueBinRebinnedErr = (histo_fake_dict[sName+s+e].GetBinError(2*r)+histo_fake_dict[sName+s+e].GetBinError(2*r-1))/2
                    else :
                        valueBinRebinned = (histo_fake_dict[sName+s+e].GetBinContent(3*r-binChange)+histo_fake_dict[sName+s+e].GetBinContent(3*r-binChange-1)+histo_fake_dict[sName+s+e].GetBinContent(3*r-binChange-2))/3 #the not-simplified value of the bin number is: 2*N+3*(r-N), -1, -2.
                        valueBinRebinnedErr = (histo_fake_dict[sName+s+e].GetBinError(3*r-binChange)+histo_fake_dict[sName+s+e].GetBinError(3*r-binChange-1)+histo_fake_dict[sName+s+e].GetBinError(3*r-binChange-2))/3 #the not-simplified value of the bin number is: 2*N+3*(r-N), -1, -2.

                    temp_histREB.SetBinContent(r,valueBinRebinned)
                    temp_histREB.SetBinError(r,valueBinRebinnedErr)
                temp_histREB.AddDirectory(False)
                histo_fake_dict[sName+s+e+'reb'] = temp_histREB.Clone()
        file_dict[sName].Close()

        #prepare the fit        
        correlatedFitterDict = {}
        for s in self.signList :
            for e in self.etaBinningS :
                np.set_printoptions(threshold=np.inf)

                dimFit = len(bin4corFitS)
                xx_ = np.zeros(dimFit)
                yy_ = np.zeros(dimFit)
                cov_ = np.zeros(( len(systList)+1,dimFit,dimFit))

                #debug lines -----------------------
                # covdict = {}
                # for syst in  systList :
                #     covdict[syst] = np.zeros((dimFit,dimFit))
                # covdict['nom'] = np.zeros((dimFit,dimFit))
                # #end of debug lines -------------------

                histo_fake_dict['current'+s+e+'reb'] = ROOT.TH1F(fakedict[s+e].GetName()+'_reb',fakedict[s+e].GetName()+'_reb',len(bin4corFit)-1,array('d',bin4corFit))
                for r in range(1,len(bin4corFit)) :
                    if DONT_REBIN :
                        valueBinRebinned = fakedict[s+e].GetBinContent(r)
                        valueBinRebinnedErr = fakedict[s+e].GetBinError(r) 
                    elif r<binChange :
                        valueBinRebinned = (fakedict[s+e].GetBinContent(2*r)+fakedict[s+e].GetBinContent(2*r-1))/2
                        valueBinRebinnedErr = (fakedict[s+e].GetBinError(2*r)+fakedict[s+e].GetBinError(2*r-1))/2
                    else :
                        valueBinRebinned = (fakedict[s+e].GetBinContent(3*r-binChange)+fakedict[s+e].GetBinContent(3*r-binChange-1)+fakedict[s+e].GetBinContent(3*r-binChange-2))/3 #the not-simplified value of the bin number is: 2*N+3*(r-N), -1, -2.
                        valueBinRebinnedErr = (fakedict[s+e].GetBinError(3*r-binChange)+fakedict[s+e].GetBinError(3*r-binChange-1)+fakedict[s+e].GetBinError(3*r-binChange-2))/3 #the not-simplified value of the bin number is: 2*N+3*(r-N), -1, -2.
                    histo_fake_dict['current'+s+e+'reb'].SetBinContent(r,valueBinRebinned)
                    histo_fake_dict['current'+s+e+'reb'].SetBinError(r,valueBinRebinnedErr)


                for pp in range(xx_.size) :

                    xx_[pp] = histo_fake_dict['current'+s+e+'reb'].GetBinCenter(pp+1)
                    yy_[pp] = histo_fake_dict['current'+s+e+'reb'].GetBinContent(pp+1)

                for pp in range(xx_.size) : #separate loop because needed xx,yy fully filled
                    for p2 in range(xx_.size) :
                        for syst in range(len(systList)+1) :
                            if pp==p2 and syst==len(systList):

                                erv = histo_fake_dict['nom'+s+e+'reb'].GetBinError(pp+1)**2
                                cov_[syst][pp][p2] = erv
                                # covdict['nom'][pp][p2] = erv
                            elif syst<len(systList):
                                if 'Up' in systList[syst] or  systList[syst] in self.LHEdict['Up']:#do not use down syst, will be symmetrized with up later
                                        systUp =systList[syst]
                                        if 'Up' in systUp :
                                            systDown =  systUp.replace("Up","Down")
                                        else :
                                            for ilhe in range(len(self.LHEdict['Up'])) :
                                                if systUp == self.LHEdict['Up'][ilhe] :
                                                        systDown = self.LHEdict['Down'][ilhe]   

                                        deltaPP = (abs(histo_fake_dict['nom'+s+e+'reb'].GetBinContent(pp+1)-histo_fake_dict[systUp+s+e+'reb'].GetBinContent(pp+1))+ abs(histo_fake_dict['nom'+s+e+'reb'].GetBinContent(pp+1)-histo_fake_dict[systDown+s+e+'reb'].GetBinContent(pp+1)))/2
                                        deltaP2 = (abs(histo_fake_dict['nom'+s+e+'reb'].GetBinContent(p2+1)-histo_fake_dict[systUp+s+e+'reb'].GetBinContent(p2+1))+ abs(histo_fake_dict['nom'+s+e+'reb'].GetBinContent(p2+1)-histo_fake_dict[systDown+s+e+'reb'].GetBinContent(p2+1)))/2
                                        # erv = (histo_fake_dict['nom'+s+e+'reb'].GetBinContent(pp+1)-histo_fake_dict[systList[syst]+s+e].GetBinContent(pp+1))*(histo_fake_dict['nom'+s+e+'reb'].GetBinContent(p2+1)-histo_fake_dict[systList[syst]+s+e].GetBinContent(p2+1))
                                        if deltaPP !=0 and deltaP2!=0 :
                                            signPP = (histo_fake_dict['nom'+s+e+'reb'].GetBinContent(pp+1)-histo_fake_dict[systUp+s+e+'reb'].GetBinContent(pp+1))/abs(histo_fake_dict['nom'+s+e+'reb'].GetBinContent(pp+1)-histo_fake_dict[systUp+s+e+'reb'].GetBinContent(pp+1))#Chosen the UP sign!
                                            signP2 = (histo_fake_dict['nom'+s+e+'reb'].GetBinContent(p2+1)-histo_fake_dict[systUp+s+e+'reb'].GetBinContent(p2+1))/abs(histo_fake_dict['nom'+s+e+'reb'].GetBinContent(p2+1)-histo_fake_dict[systUp+s+e+'reb'].GetBinContent(p2+1))#Chosen the UP sign!
                                        else :
                                            signPP =1
                                            signP2 =1
                                        erv = deltaPP*deltaP2*signPP*signP2
                                        
                                        cov_[syst][pp][p2] = erv
                                        # covdict[systList[syst]][pp][p2] = erv

                q0_ = fakedict[s+e+'offset']
                p0_ = fakedict[s+e+'slope']
                q0Err_ = fakedict[s+e+'offsetErr']
                p0Err_ = fakedict[s+e+'slopeErr']
                cov_proj = cov_.sum(axis=0) #square sum of syst

                #uncomment here -----------------------------------
                # matt = np.zeros((dimFit,dimFit)) 
                # for syst in  systList :
                #     print "rank of", syst,  np.linalg.matrix_rank(covdict[syst]), "det=", np.linalg.det(covdict[syst])
                #     matt = matt+covdict[syst]
                # matt = matt + covdict['nom']
                # print "rank of all",  np.linalg.matrix_rank(cov_proj), np.linalg.slogdet(cov_proj)
                # print "rank of matt", np.linalg.matrix_rank(matt), np.linalg.slogdet(matt)
                # # print "eigval and egvec of matt=", np.linalg.eig(matt)
                # matt = matt- cov_proj
                # print "debugg", matt
                #end of good debug ----------------------------

                invCov_ = np.linalg.inv(cov_proj)

                #do the fit
                minuitDict = self.MinuitLinearFit( yy=yy_, xx=xx_, invCov=invCov_, p0=p0_, q0=q0_,p0Err=p0Err_, q0Err=q0Err_,s=s,e=e)
                
                if self.statAna :
                    histoToy = ROOT.TH2F("histoToy"+s+e,"histoToy"+s+e,100,minuitDict[s+e+'offset'+'Minuit']-5*minuitDict[s+e+'offset'+'Minuit'+'Err'],minuitDict[s+e+'offset'+'Minuit']+5*minuitDict[s+e+'offset'+'Minuit'+'Err'],100,minuitDict[s+e+'slope'+'Minuit']-5*minuitDict[s+e+'slope'+'Minuit'+'Err'],minuitDict[s+e+'slope'+'Minuit']+5*minuitDict[s+e+'slope'+'Minuit'+'Err'])
                    Ntoys = 100
                    for toy in range(Ntoys) :
                        rrand = ROOT.TRandom3()  
                        for pp in range(xx_.size) :
                             xx_[pp] = rrand.Gaus(xx_[pp],histo_fake_dict['current'+s+e+'reb'].GetBinError(pp+1))
                        toyDict = self.MinuitLinearFit( yy=yy_, xx=xx_, invCov=invCov_, p0=p0_, q0=q0_,p0Err=p0Err_, q0Err=q0Err_,s=s,e=e)
                        histoToy.Fill(toyDict[s+e+'offset'+'Minuit'],toyDict[s+e+'slope'+'Minuit'])
                    minuitDict[s+e+'offset'+'Minuit'+'Err'] = histoToy.GetStdDev(1)
                    minuitDict[s+e+'slope'+'Minuit'+'Err'] = histoToy.GetStdDev(2)
                    minuitDict[s+e+'offset*slope'+'Minuit'] = histoToy.GetCovariance()    
                correlatedFitterDict.update(minuitDict)            
        return correlatedFitterDict
    
    
    def evaluate_erf_error(self,ERFfitResult, ERFp0,ERFp1, ERFp2) :
        ntoys = 1000
        covtoy ={}
        xvec = np.zeros(len(self.ptBinning)-1)
        for xx in range(xvec.size) :
            xvec[xx] = float(self.ptBinning[xx])+float((self.ptBinning[xx+1]-self.ptBinning[xx]))/2
        yvec = np.zeros((xvec.size,ntoys))
        parvec=np.zeros(3)
        parvec[0] = ERFp0
        parvec[1] = ERFp1
        parvec[2] = ERFp2
        covvec =np.zeros((3,3))

        for xx in range (3) :
            for yy in range (3) :
                covvec[xx][yy] = ROOT.TMatrixDRow(ERFfitResult,xx)(yy)

        def my_erf(x,par) :
            val = np.zeros(x.size)
            val = par[0]*erf(par[1]*(x)+par[2])
            return val

        for itoy in range(ntoys) :
            covtoy[itoy] = np.random.multivariate_normal(parvec, covvec)
            yvec[:,itoy] = my_erf(xvec,covtoy[itoy])
        sigmavec = np.std(yvec,axis=1)

        return sigmavec


    def dict2histConverter(self,fakedict,promptdict,hdict ,correlatedFit=True) :

        outdict = {}
        
        nameDict = {
            'prompt' : ['offset','slope','2deg', 'offset*slope', 'offset*2deg','slope*2deg', 'chi2red'],
            'fake'   : ['offset','slope','offset*slope']
            }

        dictDict = {
            'prompt' :  promptdict,
            'fake' : fakedict
            }

        for kind, par in nameDict.iteritems() :
            for pp in par :
                if correlatedFit and kind=='fake' :
                    pint = pp+'Minuit'
                else :
                    pint = pp
                outdict[kind+pp]= ROOT.TH2D(kind+'_'+pp,kind+'_'+pp,len(self.signList),array('f',[0,1,2]),len(self.etaBinning)-1,array('f',self.etaBinning))
                for s in self.signList :
                     for e in self.etaBinningS :
                         outdict[kind+pp].SetBinContent(self.signList.index(s)+1,self.etaBinningS.index(e)+1,dictDict[kind][s+e+pint])
                         if pp == 'offset' or pp=='slope' or pp=='2deg' :
                            outdict[kind+pp].SetBinError(self.signList.index(s)+1,self.etaBinningS.index(e)+1,dictDict[kind][s+e+pint+'Err'])
                            
        #erf error histogram (sign,eta,pt) to avoid to repeat the toys
        outdict['prompt'+'sigmaERFvec'] = ROOT.TH3D('prompt_sigmaERFvec','prompt_sigmaERFvec',len(self.signList),array('f',[0,1,2]),len(self.etaBinning)-1,array('f',self.etaBinning),len(self.ptBinning)-1,array('f',self.ptBinning))
        for s in self.signList :
                for e in self.etaBinningS :
                    for p in self.ptBinningS :
                        outdict['prompt'+'sigmaERFvec'].SetBinContent(self.signList.index(s)+1,self.etaBinningS.index(e)+1,self.ptBinningS.index(p)+1, promptdict[s+e+'sigmaERFvec'][self.ptBinningS.index(p)])

        return outdict


    def hist2dictConverter(self,kind, outDir, syst='nom') :
        # self.kind = kind #available kind = fake, prompt, 
        outdict = {}
        if syst == 'nom' :
            suffPar = self.corrFitSuff
        else :
            suffPar = ''
        inputfile=ROOT.TFile.Open(outDir+"/bkg_"+syst+"/bkg_parameters_file"+suffPar+self.nameSuff+".root")

        nameDict = {
            'prompt' : ['offset','slope','2deg', 'offset*slope', 'offset*2deg','slope*2deg', 'chi2red'],
            'fake'   : ['offset','slope','offset*slope']
            }
        
        for par in nameDict[kind] :
            histotemp = inputfile.Get(kind+'_'+par)
            for s in self.signList :
                    for e in self.etaBinningS :
                        outdict[s+e+par] = histotemp.GetBinContent(self.signList.index(s)+1,self.etaBinningS.index(e)+1)
                        if par == 'offset' or par=='slope' or par=='2deg' :
                            outdict[s+e+par+'Err'] = histotemp.GetBinError(self.signList.index(s)+1,self.etaBinningS.index(e)+1)

            if kind=='prompt':
                histotemp = inputfile.Get(kind+'_'+'sigmaERFvec')
                for s in self.signList :
                        for e in self.etaBinningS :
                            sigmaERFvec = []
                            for p in self.ptBinningS :
                                sigmaERFvec.append(histotemp.GetBinContent(self.signList.index(s)+1,self.etaBinningS.index(e)+1,self.ptBinningS.index(p)+1))
                            outdict[s+e+'sigmaERFvec'] = sigmaERFvec

        return outdict


    def bkg_plots(self) :
        print "> plotting..."

        inputFile = ROOT.TFile.Open(self.outdir+"/bkg_differential_fakerate"+self.corrFitSuff+self.nameSuff+".root")

        canvasList = []
        legDict = {}
        stackDict ={}

        for s in self.signList :
            for e in self.etaBinningS :

                    #---------------------------------------------COMPARISON PLOTS ---------------------------------------------#
                    c_comparison = ROOT.TCanvas("c_comparison_{sign}_{eta}".format(sign=s,eta=e),"c_comparison_{sign}_{eta}".format(sign=s,eta=e),800,600)
                    c_comparison.cd()
                    c_comparison.SetGridx()
                    c_comparison.SetGridy()
                    h_fakeMC = inputFile.Get("Fakerate/hFakes_pt_"+'fake'+"_"+s+"_"+e)
                    h_prompt = inputFile.Get("Fakerate/hFakes_pt_prompt_"+s+"_"+e)
                    h_fakeFit = inputFile.Get("Template/htempl_fake_pt_"+'fake'+"_"+s+"_"+e)
                    h_promptFit = inputFile.Get("Template/htempl_prompt_pt_"+'fake'+"_"+s+"_"+e)
                    h_fakeFit.SetName("hFakes_pt_fakeFit_"+s+"_"+e)
                    h_promptFit.SetName("hFakes_pt_promptFit_"+s+"_"+e)
                    h_fakeMC.SetLineWidth(3)
                    h_prompt.SetLineWidth(3)
                    h_fakeFit.SetLineWidth(3)
                    h_promptFit.SetLineWidth(3)
                    h_fakeMC.SetLineColor(632+2) #red
                    h_prompt.SetLineColor(600-4) #blue
                    h_fakeFit.SetLineColor(632-7)
                    h_promptFit.SetLineColor(600-7)
                    h_fakeFit.SetLineStyle(2)
                    h_promptFit.SetLineStyle(2)
                    h_fakeMC.Draw()
                    h_prompt.Draw("SAME")
                    h_fakeFit.Draw("SAME")
                    h_promptFit.Draw("SAME")
                    h_fakeMC.GetYaxis().SetRangeUser(0,1.1)
                    h_fakeMC.GetYaxis().SetTitle("Isolation Cut Efficiency")
                    h_fakeMC.GetYaxis().SetTitleOffset(1)
                    h_fakeMC.GetXaxis().SetTitle("p_{T}^{#mu} [GeV]")
                    h_fakeMC.SetTitle("Fake Rates, {min}<#eta<{max}, W {sign}".format(min=e, max=self.etaBinning[self.etaBinning.index(float(e))+1], sign=s))

                    legDict[e+s] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                    legDict[e+s].AddEntry(h_fakeMC,"Fake Rate")
                    legDict[e+s].AddEntry(h_prompt, "Prompt Rate")
                    legDict[e+s].Draw("SAME")

                    canvasList.append(c_comparison)                    

                    #---------------------------------------------TEMPLATE PLOTS ---------------------------------------------#
                    c_template = ROOT.TCanvas("c_template_{sign}_{eta}".format(sign=s,eta=e),"c_template_{sign}_{eta}".format(sign=s,eta=e),800,600)
                    c_template.cd()
                    c_template.SetGridx()
                    c_template.SetGridy()
                    h_template = inputFile.Get("Template/htempl_pt_"+'fake'+"_"+s+"_"+e)
                    h_W = inputFile.Get("Template/htempl_pt_prompt_"+s+"_"+e)
                    h_template.SetLineWidth(3)
                    h_W.SetLineWidth(3)
                    h_template.SetLineColor(632+2) #red
                    h_W.SetLineColor(600-4) #blue
                    h_W.Draw()
                    h_template.Draw("SAME")
                    h_W.GetYaxis().SetRangeUser(0,1500000)
                    h_W.GetYaxis().SetTitle("Events")
                    h_W.GetYaxis().SetTitleOffset(1)
                    h_W.GetXaxis().SetTitle("p_{T}^{#mu} [GeV]")
                    h_W.SetTitle("Compared yields, {min}<#eta<{max}, W {sign}".format(min=e, max=self.etaBinning[self.etaBinning.index(float(e))+1], sign=s))

                    legDict[e+s+"templ"] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                    legDict[e+s+"templ"].AddEntry(h_template,"Template bkg")
                    legDict[e+s+"templ"].AddEntry(h_W, "EWK+Template MC")
                    legDict[e+s].Draw("SAME")

                    canvasList.append(c_template)
                    

        outputFake = ROOT.TFile(self.outdir+"/bkg_plots"+self.corrFitSuff+self.nameSuff+".root","recreate")

        for h in range(len(canvasList)) :
                canvasList[h].Write()

        for h in range(len(canvasList)) :
            canvasList[h].IsA().Destructor(canvasList[h])


    def syst_comparison(self, outDir,systDict =bkg_utils.bkg_systematics, SymBands=True, noratio=False,statAna=False) : #see the description below
        
        # noratio=True --> In the ratio plots are plotted fake and prompt rate of systematics and nominals
        # symBans=True --> symmetric bands around nominal: 1/2*sqrt[sum_syst (up-down)^2]
        
        if statAna :
            statAnaSuff = 'statAna'
        else :
            statAnaSuff = ''

        histoNameDict = {
        'comparison' : {
            'Fakes' : ['fake', 'prompt', 'fakeFit', 'promptFit']
            },
        'template' : {
            'templ' : ['fake', 'prompt']
            }
        }
            
        #getting canvas and histoss
        finalCanvasDict = {}
        finalPlotFileDict = {}
        finalHistoDict = {}
        finalLegDict = {}

        for sKind, sList in systDict.iteritems():
            for sName in sList :
                suff_4corrFitNomOnly = self.corrFitSuff
                finalPlotFileDict[sName]=ROOT.TFile.Open(outDir+"/bkg_"+sName+"/bkg_plots"+suff_4corrFitNomOnly+self.nameSuff+".root")
                for s in self.signList :
                    for e in self.etaBinningS :
                        for canvas in histoNameDict :
                            for name in histoNameDict[canvas] :
                                for histo in histoNameDict[canvas][name] :
                                    finalHistoDict[sName+canvas+histo+s+e] =   finalPlotFileDict[sName].Get('c_'+canvas+'_'+s+'_'+e).GetPrimitive('h'+name+'_pt_'+histo+'_'+s+'_'+e)                       
        #nominal
        finalPlotFileDict['nom']=ROOT.TFile.Open(outDir+"/bkg_"+self.systName+"/bkg_plots"+self.corrFitSuff+statAnaSuff+self.nameSuff+".root")
        for s in self.signList :
            for e in self.etaBinningS :
                for canvas in histoNameDict : #comparison, ABCD...,
                    for name in histoNameDict[canvas] : #Fakes,Fakes, templ....
                        for histo in histoNameDict[canvas][name] : #fake, prompt,....
                            finalHistoDict['nom'+canvas+histo+s+e] = finalPlotFileDict['nom'].Get('c_'+canvas+'_'+s+'_'+e).GetPrimitive('h'+name+'_pt_'+histo+'_'+s+'_'+e)
        
        #building of error bands on the "nom" hisotgrams
        for s in self.signList :
            for e in self.etaBinningS :
                for canvas in histoNameDict : #comparison, ABCD...,                  
                    finalCanvasDict['nom'+canvas+s+e] = finalPlotFileDict['nom'].Get('c_'+canvas+'_'+s+'_'+e).Clone()
                    finalCanvasDict['nom'+canvas+s+e].cd()
                    fillSyle =0
                    for name in histoNameDict[canvas] : #Fakes,Fakes, templ....
                        for histo in histoNameDict[canvas][name] : #fake, prompt,....
                            finalHistoDict['nom'+canvas+histo+s+e+'error'] = ROOT.TGraphAsymmErrors()#finalCanvasDict['nom'+canvas+s+e]
                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetName(finalHistoDict['nom'+canvas+histo+s+e].GetName()+'_error')

                            for p in self.ptBinning :
                                if(p==self.ptBinning[-1]) : continue
                                varErrSum2Up = 0
                                varErrSum2Down = 0.
                                
                                #--------------------------OLD METHOD STARTS HERE---------------------------#
                                # if SymBands : #symmetric bands sum2 of the largest scarto between nom e syst(up,down)
                                #     for sKind, sList in systDict.iteritems():
                                #         for sName in sList :
                                #             if 'Up' in sName or sName in self.LHEdict['Up']:
                                #                 deltaSystNomUp = (finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))**2
                                #                 if 'Up' in sName :
                                #                      sNameDown =  sName.replace("Up","Down")
                                #                 else :
                                #                     # for i in range(len(LHEUpList)) :
                                #                     #     if sName == LHEUpList[i] :         
                                #                     #         sNameDown = LHEDownList[i]
                                #                     for i in range(len(self.LHEdict['Up'])) :
                                #                         if sName == self.LHEdict['Up'][i] :
                                #                             sNameDown = self.LHEdict['Down'][i] 
                                #                 deltaSystNomDown = (finalHistoDict[sNameDown+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))**2
                                #                 if deltaSystNomUp>deltaSystNomDown :
                                #                     deltaSystNom = deltaSystNomUp
                                #                 else :
                                #                     deltaSystNom = deltaSystNomDown
                                #                 varErrSum2Down = varErrSum2Down+deltaSystNom
                                #                 varErrSum2Up = varErrSum2Up+deltaSystNom
                                #       errHigh = math.sqrt(varErrSum2Up)
                                #       errLow = math.sqrt(varErrSum2Down)
                                #--------------------------OLD METHOD FINISHES HERE---------------------------#
                                
                                if SymBands :
                                    deltaSyst = 0 
                                    for sKind, sList in systDict.iteritems():
                                        for sName in sList :
                                            if 'Down' in sName : continue
                                            if sName in self.LHEdict['Down']: continue
                                            if 'Up' in sName :
                                                sNameDown =  sName.replace("Up","Down")
                                            else :
                                                for i in range(len(self.LHEdict['Up'])) :
                                                    if sName == self.LHEdict['Up'][i] :
                                                        sNameDown = self.LHEdict['Down'][i]
                                            deltaSyst += (finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-finalHistoDict[sNameDown+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))**2
                                            nomVal = finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)
                                            if (finalHistoDict[sNameDown+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)<nomVal and finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)<nomVal) or (finalHistoDict[sNameDown+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)>nomVal and finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)>nomVal) : #nominal not in between systs
                                                # if 'Fit' in histo :
                                                    # diffUp = abs(finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))
                                                    # diffDown = abs(finalHistoDict[sNameDown+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))
                                                    # if finalHistoDict['nom'+canvas+histo+s+e].GetBinError(self.ptBinning.index(float(p))+1) < diffUp or finalHistoDict['nom'+canvas+histo+s+e].GetBinError(self.ptBinning.index(float(p))+1) < diffDown :
                                                        print "WARNING: systematic", sName, canvas,histo," up/down not around nominal in bin", p,e,s, ">>> down,nom,up=",finalHistoDict[sNameDown+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1),nomVal,finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1), 'statErr=', finalHistoDict['nom'+canvas+histo+s+e].GetBinError(self.ptBinning.index(float(p))+1)
                                    deltaSyst = 0.5*math.sqrt(deltaSyst) 
                                    errHigh = deltaSyst
                                    errLow = deltaSyst                                          
                                            
                                
                                else : #asymmetric bands                                
                                    for sKind, sList in systDict.iteritems():
                                        for sName in sList :
                                            # if "Up" in sName :
                                                # varErrSum2Up = varErrSum2Up+ (finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))**2
                                            # if "Down" in sName :
                                                # varErrSum2Down = varErrSum2Down+ (finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))**2                                       
                                            deltaSystNom = finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)
                                            if deltaSystNom < 0 : #syst<Nom
                                                varErrSum2Down = varErrSum2Down+deltaSystNom**2
                                            else : #syst>nom
                                                varErrSum2Up = varErrSum2Up+deltaSystNom**2
                                                
                                    errHigh = math.sqrt(varErrSum2Up)
                                    errLow = math.sqrt(varErrSum2Down)
                            
                                finalHistoDict['nom'+canvas+histo+s+e+'error'].SetPoint(self.ptBinning.index(float(p)),finalHistoDict['nom'+canvas+histo+s+e].GetBinCenter(self.ptBinning.index(float(p))+1),finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))
                                finalHistoDict['nom'+canvas+histo+s+e+'error'].SetPointEYhigh(self.ptBinning.index(float(p)),errHigh)
                                finalHistoDict['nom'+canvas+histo+s+e+'error'].SetPointEYlow(self.ptBinning.index(float(p)),errLow)
                                finalHistoDict['nom'+canvas+histo+s+e+'error'].SetPointEXhigh(self.ptBinning.index(float(p)),finalHistoDict['nom'+canvas+histo+s+e].GetBinWidth(self.ptBinning.index(float(p))+1)/2)
                                finalHistoDict['nom'+canvas+histo+s+e+'error'].SetPointEXlow(self.ptBinning.index(float(p)),finalHistoDict['nom'+canvas+histo+s+e].GetBinWidth(self.ptBinning.index(float(p))+1)/2)

                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetLineColor(finalHistoDict['nom'+canvas+histo+s+e].GetLineColor()-3)
                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetMarkerColor(finalHistoDict['nom'+canvas+histo+s+e].GetLineColor()-3)
                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetFillColor(finalHistoDict['nom'+canvas+histo+s+e].GetLineColor()-3)
                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetFillStyle(3003+fillSyle)
                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetLineWidth(1)

                            finalHistoDict['nom'+canvas+histo+s+e+'error'].Draw("SAME 0P5")
                            # if(correlatedFit and canvas=='comparison' and histo=='fake' and 0) : #DISABLED (see if with the same boolean )
                            if(0) : #DISABLED (see if with the same boolean )
                                keyList = ['def','offsePlus','offsetMinus','slopePlus','slopeMinus']
                                for key in keyList :
                                    finalHistoDict['nom'+canvas+s+e+'functionG'+key].Draw("SAME")
                                    finalHistoDict['nom'+canvas+s+e+'functionG'+key].SetLineWidth(3)
                                    finalHistoDict['nom'+canvas+s+e+'functionG'+key].SetLineColor(2+keyList.index(key))


                            finalCanvasDict['nom'+canvas+s+e].Update()
                            finalCanvasDict['nom'+canvas+s+e].Modified()
                            fillSyle =fillSyle+1

        #ratio plot creation
        for s in self.signList :
            for e in self.etaBinningS :
                for canvas in histoNameDict :
                    for name in histoNameDict[canvas] :
                        for histo in histoNameDict[canvas][name] :
                            c_ratioSyst = ROOT.TCanvas("c_ratioSyst_{sign}_{eta}_{canvas}_{histo}".format(sign=s,eta=e,canvas=canvas,histo=histo),"c_ratioSyst_{sign}_{eta}_{canvas}_{histo}".format(sign=s,eta=e,canvas=canvas,histo=histo),800,600)
                            c_ratioSyst.cd()
                            c_ratioSyst.SetGridx()
                            c_ratioSyst.SetGridy()
                            finalLegDict[e+s+canvas+histo+"ratioSyst"] = ROOT.TLegend(0.1,0.7,0.48,0.9)

                            sameFlag = True
                            colorNumber = 1
                            colorList = [600,616,416,632,432,800,900]
                            colorCounter = 0

                            for sKind, sList in systDict.iteritems():
                                colorNumber = colorList[colorCounter]
                                colorCounter = colorCounter+1
                                for sName in sList :
                                    colorNumber = colorNumber-2
                                    if colorNumber < colorList[colorCounter-1]-10 :
                                        colorNumber = colorList[colorCounter]+2
                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'] = ROOT.TH1F(canvas+'_'+finalHistoDict[sName+canvas+histo+s+e].GetName()+'_'+sName+'_ratio',canvas+'_'+finalHistoDict[sName+canvas+histo+s+e].GetName()+'_'+sName+'_ratio',len(self.ptBinning)-1, array('f',self.ptBinning))
                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'].Divide(finalHistoDict[sName+canvas+histo+s+e],finalHistoDict['nom'+canvas+histo+s+e],1,1)
                                    if noratio :
                                        finalHistoDict[sName+canvas+histo+s+e+'ratio'] = finalHistoDict[sName+canvas+histo+s+e].Clone()

                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'].SetLineColor(colorNumber)
                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'].SetTitle("Ratio Syst/Nom: "+canvas+', '+histo+', '+s+', #eta='+e)
                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetXaxis().SetTitle("p_{T} [GeV]")
                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetYaxis().SetTitle("Syst/Nom")
                                    for b in range(1,finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetNbinsX()+1) :
                                        finalHistoDict[sName+canvas+histo+s+e+'ratio'].SetBinError(b,0)
                                    c_ratioSyst.cd()
                                    if sameFlag :
                                        finalHistoDict[sName+canvas+histo+s+e+'ratio'].Draw()
                                        finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetYaxis().SetRangeUser(0.6,1.7)
                                    else :
                                        finalHistoDict[sName+canvas+histo+s+e+'ratio'].Draw("SAME")
                                        finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetYaxis().SetRangeUser(0,3)
                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'].SetLineWidth(3)
                                    sameFlag=False
                                    finalLegDict[e+s+canvas+histo+"ratioSyst"].AddEntry(finalHistoDict[sName+canvas+histo+s+e+'ratio'], sName)

                            #nominal
                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'] = ROOT.TH1F(canvas+'_'+finalHistoDict['nom'+canvas+histo+s+e].GetName()+'_'+'nom'+'_ratio',canvas+'_'+finalHistoDict['nom'+canvas+histo+s+e].GetName()+'_'+'nom'+'_ratio',len(self.ptBinning)-1, array('f',self.ptBinning))
                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'].Divide(finalHistoDict['nom'+canvas+histo+s+e],finalHistoDict['nom'+canvas+histo+s+e],1,1)
                            if noratio :
                                finalHistoDict['nom'+canvas+histo+s+e+'ratio'] =finalHistoDict['nom'+canvas+histo+s+e].Clone()

                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'].SetLineColor(1)
                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'].SetTitle('nom'+', '+canvas+', '+histo+', '+s+', #eta='+e)
                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'].GetXaxis().SetTitle("p_{T} [GeV]")
                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'].GetYaxis().SetTitle("Syst/Nom")
                            for b in range(1,finalHistoDict['nom'+canvas+histo+s+e+'ratio'].GetNbinsX()+1) :
                                # finalHistoDict['nom'+canvas+histo+s+e+'ratio'].SetBinError(b,math.sqrt(2)*finalHistoDict['nom'+canvas+histo+s+e].GetBinError(b)/finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(b))
                                if finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(b)!=0 :
                                    finalHistoDict['nom'+canvas+histo+s+e+'ratio'].SetBinError(b,finalHistoDict['nom'+canvas+histo+s+e].GetBinError(b)/finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(b))
                                else :
                                     print "WARNING: bin content of ", canvas, histo, "is 0, bin (s,e,p)=", s,e,b, ". error of ratio plot not correctly normalized."
                                     finalHistoDict['nom'+canvas+histo+s+e+'ratio'].SetBinError(b,finalHistoDict['nom'+canvas+histo+s+e].GetBinError(b))
                            c_ratioSyst.cd()
                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'].Draw("SAME E2")
                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'].GetYaxis().SetRangeUser(0,3)
                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'].SetLineWidth(0)
                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'].SetFillColor(1)
                            finalHistoDict['nom'+canvas+histo+s+e+'ratio'].SetFillStyle(3002)
                            finalLegDict[e+s+canvas+histo+"ratioSyst"].AddEntry(finalHistoDict['nom'+canvas+histo+s+e+'ratio'], 'Nominal')


                            #sum2 errorband in the ratios (DISABLED, set=1 se vuoi vedere la banda dei sum2 sui ratio)
                            sum2Flag = False
                            if sum2Flag :
                                finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'] = ROOT.TGraphAsymmErrors()
                                finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].SetName(canvas+'_'+finalHistoDict['nom'+canvas+histo+s+e].GetName()+'_'+'nom'+'_ratio_sum2')
                                finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].SetLineColor(2)
                                for p in self.ptBinning :
                                    if(p==self.ptBinning[-1]) : continue
                                    bin_p = self.ptBinning.index(float(p))
                                    finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].SetPoint(bin_p,finalHistoDict['nom'+canvas+histo+s+e].GetBinCenter(bin_p+1),1)
                                    finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].SetPointEYhigh(bin_p,finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorYhigh(bin_p)/finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(bin_p+1))
                                    finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].SetPointEYlow(bin_p,finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorYlow(bin_p)/finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(bin_p+1))
                                    finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].SetPointEXhigh(bin_p,finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorXhigh(bin_p))
                                    finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].SetPointEXlow(bin_p,finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorXlow(bin_p))
                                c_ratioSyst.cd()
                                finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].Draw("SAME OP5")
                                finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].SetLineWidth(0)
                                finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].SetFillColor(2)
                                finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].SetFillStyle(3144)
                                finalLegDict[e+s+canvas+histo+"ratioSyst"].AddEntry(finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'], 'Square Sum')

                            finalLegDict[e+s+canvas+histo+"ratioSyst"].Draw("SAME")
                            finalCanvasDict['ratio'+canvas+histo+s+e] = c_ratioSyst

        #total error comparison
        for s in self.signList :
            for e in self.etaBinningS :
                c_errCompare = ROOT.TCanvas("c_errCompare_{sign}_{eta}".format(sign=s,eta=e),"c_errCompare_{sign}_{eta}".format(sign=s,eta=e),800,600)
                c_errCompare.cd()
                c_errCompare.SetGridx()
                c_errCompare.SetGridy()
                finalLegDict[e+s+'errCompare'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                canvas = 'template'
                for histo in ['fake','prompt'] :
                            finalHistoDict['nom'+canvas+histo+s+e+'sum2'] = ROOT.TGraphAsymmErrors()
                            finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetName(canvas+'_'+finalHistoDict['nom'+canvas+histo+s+e].GetName()+'_'+'nom'+'_ratio_sum2')

                            for p in self.ptBinning :
                                if(p==self.ptBinning[-1]) : continue
                                bin_p = self.ptBinning.index(float(p))
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetPoint(bin_p,finalHistoDict['nom'+canvas+histo+s+e].GetBinCenter(bin_p+1),1)
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetPointEYhigh(bin_p,math.sqrt(finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorYhigh(bin_p)**2+finalHistoDict['nom'+canvas+histo+s+e].GetBinError(bin_p+1)**2))
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetPointEYlow(bin_p,math.sqrt(finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorYlow(bin_p)**2+finalHistoDict['nom'+canvas+histo+s+e].GetBinError(bin_p+1)**2))
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetPointEXhigh(bin_p,finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorXhigh(bin_p))
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetPointEXlow(bin_p,finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorXlow(bin_p))

                            c_errCompare.cd()
                            if histo == 'fake' :
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].Draw()
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].Draw("SAME 5")
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetLineColor(632+2) #red
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetFillColor(632+2) #red
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetFillStyle(3002)
                                finalLegDict[e+s+'errCompare'].AddEntry(finalHistoDict['nom'+canvas+histo+s+e+'sum2'], 'QCD bkg')
                            else:
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].Draw("SAME 5")
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetLineColor(600-4) #blue
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetFillColor(600-4) #blue
                                finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetFillStyle(3005)
                                finalLegDict[e+s+'errCompare'].AddEntry(finalHistoDict['nom'+canvas+histo+s+e+'sum2'], 'W MC')
                            finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetLineWidth(1)
                            finalHistoDict['nom'+canvas+histo+s+e+'sum2'].SetTitle("Square Sum of Errors: "+canvas+', '+s+', #eta='+e)
                            finalHistoDict['nom'+canvas+histo+s+e+'sum2'].GetXaxis().SetTitle("p_{T} [GeV]")
                            finalHistoDict['nom'+canvas+histo+s+e+'sum2'].GetYaxis().SetTitle("Events/1 GeV^{-1}")

                            # finalLegDict[e+s+'errCompare'].AddEntry(finalHistoDict['nom'+canvas+histo+s+e+'sum2'], histo)
                finalLegDict[e+s+'errCompare'].Draw("SAME")
                finalCanvasDict['sum2Error'+s+e] = c_errCompare


        #unrolled plot creation
        #binning
        unrolledPtEta= list(self.ptBinning)
        intervalPtBin = []
        for p in self.ptBinning[:-1] :
            intervalPtBin.append(self.ptBinning[self.ptBinning.index(p)+1]-self.ptBinning[self.ptBinning.index(p)])

        for e in range(len(self.etaBinning)-2) :
            for p in intervalPtBin :
                unrolledPtEta.append(unrolledPtEta[-1]+p)
        # print "unrolledPtEta", unrolledPtEta

        #final plot and syst
        for s in self.signList :
            for canvas in histoNameDict :
                c_unrolled = ROOT.TCanvas("c_unrolled_{canvas}_{sign}".format(canvas=canvas,sign=s),"c_unrolled_{canvas}_{sign}".format(canvas=canvas,sign=s),800,600)
                c_unrolled.cd()
                c_unrolled.SetGridx()
                c_unrolled.SetGridy()
                finalLegDict[s+canvas+'unrolled'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                sameFlagUNR=True
                for name in histoNameDict[canvas] :
                    for histo in histoNameDict[canvas][name] :
                        finalHistoDict[s+canvas+histo+'unrolled'] = ROOT.TH1F(canvas+'_'+finalHistoDict['nom'+canvas+histo+s+'0'].GetName()+'_unrolled',canvas+'_'+finalHistoDict['nom'+canvas+histo+s+'0'].GetName()+'_unrolled',len(unrolledPtEta)-1, array('f',unrolledPtEta))
                        finalHistoDict[s+canvas+histo+'unrolled'].SetLineColor(finalHistoDict['nom'+canvas+histo+s+'0'].GetLineColor())
                        finalHistoDict[s+canvas+histo+'unrolled'].SetLineWidth(finalHistoDict['nom'+canvas+histo+s+'0'].GetLineWidth())
                        finalHistoDict[s+canvas+histo+'unrolled'].GetXaxis().SetTitle('Unrolled #eta, p_{T}')
                        finalHistoDict[s+canvas+histo+'unrolled'].GetYaxis().SetTitle(finalHistoDict['nom'+canvas+histo+s+'0'].GetYaxis().GetTitle())
                        finalHistoDict[s+canvas+histo+'unrolled'].SetTitle('Unorlled '+canvas+', W '+s)

                        finalHistoDict[s+canvas+histo+'unrolled'+'error'] = ROOT.TGraphAsymmErrors()
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetName(finalHistoDict[s+canvas+histo+'unrolled'].GetName()+'_error')
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetLineColor(finalHistoDict['nom'+canvas+histo+s+'0'+'error'].GetLineColor())
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetMarkerColor(finalHistoDict['nom'+canvas+histo+s+'0'+'error'].GetLineColor())
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetFillColor(finalHistoDict['nom'+canvas+histo+s+'0'+'error'].GetLineColor())
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetFillStyle(finalHistoDict['nom'+canvas+histo+s+'0'+'error'].GetFillStyle())
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetLineWidth(finalHistoDict['nom'+canvas+histo+s+'0'+'error'].GetLineWidth())



                        for e in self.etaBinningS :
                            for p in self.ptBinningS :
                                indexUNR = self.etaBinning.index(float(e))*len(self.ptBinningS)+self.ptBinning.index(float(p))
                                finalHistoDict[s+canvas+histo+'unrolled'].SetBinContent(indexUNR+1,finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))
                                finalHistoDict[s+canvas+histo+'unrolled'].SetBinError(indexUNR+1,finalHistoDict['nom'+canvas+histo+s+e].GetBinError(self.ptBinning.index(float(p))+1))

                                finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetPoint(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinCenter(indexUNR+1),finalHistoDict[s+canvas+histo+'unrolled'].GetBinContent(indexUNR+1))
                                finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetPointEYhigh(indexUNR,finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorYhigh(self.ptBinning.index(float(p))))
                                finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetPointEYlow(indexUNR,finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorYlow(self.ptBinning.index(float(p))))
                                finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetPointEXhigh(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinWidth(indexUNR+1)/2)
                                finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetPointEXlow(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinWidth(indexUNR+1)/2)

                        if sameFlagUNR :
                            finalHistoDict[s+canvas+histo+'unrolled'].Draw()
                            if(canvas!="template") :
                                finalHistoDict[s+canvas+histo+'unrolled'].GetYaxis().SetRangeUser(0,1.1)
                        else :
                            finalHistoDict[s+canvas+histo+'unrolled'].Draw("SAME")
                        sameFlagUNR=False
                        finalLegDict[s+canvas+'unrolled'].AddEntry(finalHistoDict[s+canvas+histo+'unrolled'], histo)
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].Draw("SAME 0P5")
                finalLegDict[s+canvas+'unrolled'].Draw("SAME")
                finalCanvasDict['unrolled'+canvas+histo+s] = c_unrolled


        #ratios unrolled
        for s in self.signList :
            for canvas in histoNameDict :
                for name in histoNameDict[canvas] :
                    for histo in histoNameDict[canvas][name] :
                        c_ratioSyst_unrolled = ROOT.TCanvas("c_ratioSyst_unrolled_{sign}_{canvas}_{histo}".format(sign=s,canvas=canvas,histo=histo),"c_ratioSyst_unrolled_{sign}_{canvas}_{histo}".format(sign=s,canvas=canvas,histo=histo),800,600)
                        c_ratioSyst_unrolled.cd()
                        c_ratioSyst_unrolled.SetGridx()
                        c_ratioSyst_unrolled.SetGridy()
                        finalLegDict[s+canvas+histo+"ratioSyst_unrolled"] = ROOT.TLegend(0.1,0.7,0.48,0.9)

                        sameFlagUNR = True

                        for sKind, sList in systDict.iteritems():
                            for sName in sList :
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'] = ROOT.TH1F(finalHistoDict[sName+canvas+histo+s+'0ratio'].GetName()+'_ratio_unrolled',finalHistoDict[sName+canvas+histo+s+'0ratio'].GetName()+'_ratio_unrolled',len(unrolledPtEta)-1, array('f',unrolledPtEta))
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetLineColor(finalHistoDict[sName+canvas+histo+s+'0ratio'].GetLineColor())
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetLineWidth(finalHistoDict[sName+canvas+histo+s+'0ratio'].GetLineWidth())
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].GetXaxis().SetTitle('Unrolled #eta, p_{T}')
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].GetYaxis().SetTitle(finalHistoDict[sName+canvas+histo+s+'0ratio'].GetYaxis().GetTitle())
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetTitle('Ratio Syst/Nom: Unorlled '+canvas+', W'+s)

                                for e in self.etaBinningS :
                                    for p in self.ptBinningS :
                                        indexUNR = self.etaBinning.index(float(e))*len(self.ptBinningS)+self.ptBinning.index(float(p))
                                        finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetBinContent(indexUNR+1,finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetBinContent(self.ptBinning.index(float(p))+1))
                                        finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetBinError(indexUNR+1,finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetBinError(self.ptBinning.index(float(p))+1))
                                c_ratioSyst_unrolled.cd()
                                if sameFlagUNR :
                                    finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].Draw()
                                    finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].GetYaxis().SetRangeUser(0.6,1.7)
                                else :
                                    finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].Draw("SAME")
                                    finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].GetYaxis().SetRangeUser(0,3)
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetLineWidth(3)
                                sameFlagUNR=False
                                finalLegDict[s+canvas+histo+"ratioSyst_unrolled"].AddEntry(finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'], sName)

                        #nominal
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'] = ROOT.TH1F(finalHistoDict['nom'+canvas+histo+s+'0ratio'].GetName()+'_ratio_unrolled',finalHistoDict['nom'+canvas+histo+s+'0ratio'].GetName()+'_ratio_unrolled',len(unrolledPtEta)-1, array('f',unrolledPtEta))
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].SetLineColor(finalHistoDict['nom'+canvas+histo+s+'0ratio'].GetLineColor())
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].SetLineWidth(finalHistoDict['nom'+canvas+histo+s+'0ratio'].GetLineWidth())
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].GetXaxis().SetTitle('Unrolled #eta, p_{T}')
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].GetYaxis().SetTitle(finalHistoDict['nom'+canvas+histo+s+'0ratio'].GetYaxis().GetTitle())
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].SetTitle('Ratio Syst/Nom: Unorlled '+canvas+', W '+s)

                        for e in self.etaBinningS :
                            for p in self.ptBinningS :
                                indexUNR = self.etaBinning.index(float(e))*len(self.ptBinningS)+self.ptBinning.index(float(p))
                                # print "index=", indexUNR, ", set=",indexUNR+1, "get=", self.ptBinning.index(float(p))+1
                                finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].SetBinContent(indexUNR+1,finalHistoDict['nom'+canvas+histo+s+e+'ratio'].GetBinContent(self.ptBinning.index(float(p))+1))
                                finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].SetBinError(indexUNR+1,finalHistoDict['nom'+canvas+histo+s+e+'ratio'].GetBinError(self.ptBinning.index(float(p))+1))
                        c_ratioSyst_unrolled.cd()
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].Draw("SAME E2")
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].GetYaxis().SetRangeUser(0,3)
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].SetLineWidth(0)
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].SetFillColor(1)
                        finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'].SetFillStyle(3002)
                        finalLegDict[s+canvas+histo+"ratioSyst_unrolled"].AddEntry(finalHistoDict['nom'+s+canvas+histo+'ratio_unrolled'], 'Nominal')

                        #sum2 errorband in the ratios (DISABLED!)
                        sum2Flag_unr = False
                        if sum2Flag_unr :
                            finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'] = ROOT.TGraphAsymmErrors()
                            finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetName(finalHistoDict[s+canvas+histo+'unrolled'].GetName()+'_ratio_sum2')
                            finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetLineColor(finalHistoDict['nom'+canvas+histo+s+'0'+'ratio_sum2'].GetLineColor())
                            finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetMarkerColor(finalHistoDict['nom'+canvas+histo+s+'0'+'ratio_sum2'].GetLineColor())
                            finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetFillColor(finalHistoDict['nom'+canvas+histo+s+'0'+'ratio_sum2'].GetLineColor())
                            finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetFillStyle(finalHistoDict['nom'+canvas+histo+s+'0'+'ratio_sum2'].GetFillStyle())
                            finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetLineWidth(finalHistoDict['nom'+canvas+histo+s+'0'+'ratio_sum2'].GetLineWidth())
                            for e in self.etaBinningS :
                                for p in self.ptBinningS :
                                    indexUNR = self.etaBinning.index(float(e))*len(self.ptBinningS)+self.ptBinning.index(float(p))
                                    finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetPoint(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinCenter(indexUNR+1),1)
                                    finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetPointEYhigh(indexUNR,finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].GetErrorYhigh(self.ptBinning.index(float(p))))
                                    finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetPointEYlow(indexUNR,finalHistoDict['nom'+canvas+histo+s+e+'ratio_sum2'].GetErrorYlow(self.ptBinning.index(float(p))))
                                    finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetPointEXhigh(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinWidth(indexUNR+1)/2)
                                    finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].SetPointEXlow(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinWidth(indexUNR+1)/2)
                            finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'].Draw("SAME 0P5")
                            finalLegDict[s+canvas+histo+"ratioSyst_unrolled"].AddEntry(finalHistoDict[s+canvas+histo+'unrolled'+'ratio_sum2'], 'Square Sum')

                        finalLegDict[s+canvas+histo+"ratioSyst_unrolled"].Draw("SAME")
                        finalCanvasDict['ratio_unrolled'+canvas+histo+s] = c_ratioSyst_unrolled


        #unrolled sum2 errorbands canvas
        for s in self.signList :
            c_errCompare_unrolled = ROOT.TCanvas("c_errCompare_unrolled_{sign}".format(sign=s),"c_errCompare_unrolled_{sign}".format(sign=s),800,600)
            c_errCompare_unrolled.cd()
            c_errCompare_unrolled.SetGridx()
            c_errCompare_unrolled.SetGridy()
            finalLegDict[s+'errCompare_unrolled'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
            canvas = 'template'
            for histo in ['fake','prompt'] :
                            finalHistoDict[s+canvas+histo+'unrolled'+'sum2'] = ROOT.TGraphAsymmErrors()
                            finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetName(finalHistoDict[s+canvas+histo+'unrolled'].GetName()+'_sum2')
                            finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetLineColor(finalHistoDict['nom'+canvas+histo+s+'0'+'sum2'].GetLineColor())
                            finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetMarkerColor(finalHistoDict['nom'+canvas+histo+s+'0'+'sum2'].GetLineColor())
                            finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetFillColor(finalHistoDict['nom'+canvas+histo+s+'0'+'sum2'].GetLineColor())
                            finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetFillStyle(finalHistoDict['nom'+canvas+histo+s+'0'+'sum2'].GetFillStyle())
                            finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetLineWidth(finalHistoDict['nom'+canvas+histo+s+'0'+'sum2'].GetLineWidth())
                            finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].GetXaxis().SetTitle('Unrolled #eta, p_{T}')
                            finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].GetYaxis().SetTitle(finalHistoDict['nom'+canvas+histo+s+'0'+'sum2'].GetYaxis().GetTitle())
                            finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetTitle('Square Sum of Errors: Unrolled, '+ canvas+', W '+s)
                            
                            for e in self.etaBinningS :
                                for p in self.ptBinningS :
                                    indexUNR = self.etaBinning.index(float(e))*len(self.ptBinningS)+self.ptBinning.index(float(p))
                                    finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetPoint(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinCenter(indexUNR+1),1)
                                    finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetPointEYhigh(indexUNR,finalHistoDict['nom'+canvas+histo+s+e+'sum2'].GetErrorYhigh(self.ptBinning.index(float(p))))
                                    finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetPointEYlow(indexUNR,finalHistoDict['nom'+canvas+histo+s+e+'sum2'].GetErrorYlow(self.ptBinning.index(float(p))))
                                    finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetPointEXhigh(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinWidth(indexUNR+1)/2)
                                    finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].SetPointEXlow(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinWidth(indexUNR+1)/2)
                            if histo=='fake' :
                                finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].Draw("")
                                finalLegDict[s+'errCompare_unrolled'].AddEntry(finalHistoDict[s+canvas+histo+'unrolled'+'sum2'], 'QCD bkg')
                            else :
                                finalHistoDict[s+canvas+histo+'unrolled'+'sum2'].Draw("SAME 5")
                                finalLegDict[s+'errCompare_unrolled'].AddEntry(finalHistoDict[s+canvas+histo+'unrolled'+'sum2'], 'W MC')
                            # finalLegDict[s+"errCompare_unrolled"].AddEntry(finalHistoDict[s+canvas+histo+'unrolled'+'sum2'], histo)
            # finalHistoDict[s+canvas+'fake'+'unrolled'+'sum2'].Draw("5")
            finalLegDict[s+'errCompare_unrolled'].Draw("SAME")
            finalCanvasDict['sum2Error_unrolled'+s] = c_errCompare_unrolled




        #fit parameters i.e. the real templates -------
        
        ParnameDict = {
            'prompt' : ['offset','slope','2deg', 'offset*slope', 'offset*2deg','slope*2deg', 'chi2red'],
            'fake'   : ['offset','slope','offset*slope']
            }
        for sKind, sList in systDict.iteritems():
            for sName in sList :
                suffPar = ''
                suff_4corrFitNomOnly = self.corrFitSuff
                finalPlotFileDict[sName+'parameters'] = ROOT.TFile.Open(outDir+"/bkg_"+sName+"/bkg_parameters_file"+suff_4corrFitNomOnly+self.nameSuff+".root")
                for s in self.signList :
                    for kind in ParnameDict :
                        for par in ParnameDict[kind] :
                                temp2DHist =   finalPlotFileDict[sName+'parameters'].Get(kind+'_'+par)
                                if s== 'Plus' : sCount=1
                                else : sCount =2
                                finalHistoDict[sName+kind+par+s] = temp2DHist.ProjectionY('parameters_'+kind+'_'+par+'_'+s,sCount,sCount,"e")
        #nominal
        sName=self.systName
        finalPlotFileDict['nom'+'parameters'] = ROOT.TFile.Open(outDir+"/bkg_"+sName+"/bkg_parameters_file"+self.corrFitSuff+statAnaSuff+self.nameSuff+".root")
        for s in self.signList :
            for kind in ParnameDict :
                for par in ParnameDict[kind] :
                    temp2DHist =   finalPlotFileDict['nom'+'parameters'].Get(kind+'_'+par)
                    if s== 'Plus' : sCount=1
                    else : sCount =2
                    #NB: if i run with correlatedFit=1 in bkg_parameters_file.root there are the "Minuit" (aka correlatedFit) parameters. 
                    finalHistoDict['nom'+kind+par+s] = temp2DHist.ProjectionY('parameters_'+kind+'_'+par+'_'+s,sCount,sCount,"e")     
       
        #building bands
        for s in self.signList :
            for kind in ParnameDict : 
                for par in ParnameDict[kind] :              
                    finalCanvasDict['parameters'+'nom'+kind+par+s] = ROOT.TCanvas("c_parameters_{sign}_{kind}_{par}".format(sign=s,kind=kind,par=par),"c_parameters_{sign}_{kind}_{par}".format(sign=s,kind=kind,par=par),800,600)#finalPlotFileDict['nom'].Get('c_'+canvas+'_'+s+'_'+e).Clone()
                    finalCanvasDict['parameters'+'nom'+kind+par+s].cd()
                    fillSyle =0
                    finalHistoDict['nom'+kind+par+s+'error'] = ROOT.TGraphAsymmErrors()
                    finalHistoDict['nom'+kind+par+s+'error'].SetName(finalHistoDict['nom'+kind+par+s].GetName()+'_error')

                    for e in self.etaBinning :
                        if(e==self.etaBinning[-1]) : continue
                        varErrSum2Up = 0
                        varErrSum2Down = 0.
                        
                        #--------------------------OLD METHOD STARTS HERE---------------------------#
                        # if SymBands : #symmetric bands sum2 of the largest scarto between nom e syst(up,down)
                        #     for sKind, sList in systDict.iteritems():
                        #         for sName in sList :
                        #             if 'Up' in sName :
                        #                 deltaSystNomUp = (finalHistoDict[sName+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)-finalHistoDict['nom'+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1))**2
                        #                 sNameDown = sName.replace('Up','Down')
                        #                 deltaSystNomDown = (finalHistoDict[sNameDown+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)-finalHistoDict['nom'+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1))**2
                        #                 if deltaSystNomUp>deltaSystNomDown :
                        #                     deltaSystNom = deltaSystNomUp
                        #                 else :
                        #                     deltaSystNom = deltaSystNomDown
                        #                 varErrSum2Down = varErrSum2Down+deltaSystNom
                        #                 varErrSum2Up = varErrSum2Up+deltaSystNom
                        #     errHigh = math.sqrt(varErrSum2Up)
                        #     errLow = math.sqrt(varErrSum2Down)
                        #--------------------------OLD METHOD FINISHES HERE---------------------------#

                        
                        if SymBands :
                            deltaSyst = 0 
                            for sKind, sList in systDict.iteritems():
                                for sName in sList :
                                    if 'Down' in sName : continue
                                    if sName in self.LHEdict['Down']: continue
                                    if 'Up' in sName :
                                        sNameDown =  sName.replace("Up","Down")
                                    else :
                                        for i in range(len(self.LHEdict['Up'])) :
                                            if sName == self.LHEdict['Up'][i] :
                                                sNameDown = self.LHEdict['Down'][i]
                                    deltaSyst += (finalHistoDict[sName+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)-finalHistoDict[sNameDown+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1))**2
                                    nomVal = finalHistoDict['nom'+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)
                                    if (finalHistoDict[sNameDown+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)<nomVal and finalHistoDict[sName+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)<nomVal) or (finalHistoDict[sNameDown+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)>nomVal and finalHistoDict[sName+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)>nomVal) : #nominal not in between systs
                                        # if par in ['offset','slope','2deg'] :
                                            # diffUp = abs(finalHistoDict[sName+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)-finalHistoDict['nom'+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1))
                                            # diffDown = abs(finalHistoDict[sNameDown+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)-finalHistoDict['nom'+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1))
                                            # if finalHistoDict['nom'+kind+par+s].GetBinError(self.etaBinning.index(float(e))+1) < diffUp or finalHistoDict['nom'+kind+par+s].GetBinError(self.etaBinning.index(float(e))+1) < diffDown :
                                                print "WARNING: systematic (parameters)", sName," up/down not around nominal in bin", p,e,s, ">>> down,nom,up=", finalHistoDict[sNameDown+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1), nomVal, finalHistoDict[sName+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)
                            deltaSyst = 0.5*math.sqrt(deltaSyst) 
                            errHigh = deltaSyst
                            errLow = deltaSyst
                                                
                        else : #asymmetric bands 
                            for sKind, sList in systDict.iteritems():
                                for sName in sList :
                                    # if "Up" in sName :
                                        # varErrSum2Up = varErrSum2Up+ (finalHistoDict[sName+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)-finalHistoDict['nom'+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1))**2
                                    # if "Down" in sName :
                                        # varErrSum2Down = varErrSum2Down+ (finalHistoDict[sName+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)-finalHistoDict['nom'+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1))**2
                                    deltaSystNom = finalHistoDict[sName+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)-finalHistoDict['nom'+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1)
                                    if deltaSystNom < 0 : #syst<Nom
                                        varErrSum2Down = varErrSum2Down+deltaSystNom**2
                                    else : #syst>nom
                                        varErrSum2Up = varErrSum2Up+deltaSystNom**2
                            
                            errHigh = math.sqrt(varErrSum2Up)
                            errLow = math.sqrt(varErrSum2Down)
                        finalHistoDict['nom'+kind+par+s+'error'].SetPoint(self.etaBinning.index(float(e)),finalHistoDict['nom'+kind+par+s].GetBinCenter(self.etaBinning.index(float(e))+1),finalHistoDict['nom'+kind+par+s].GetBinContent(self.etaBinning.index(float(e))+1))
                        finalHistoDict['nom'+kind+par+s+'error'].SetPointEYhigh(self.etaBinning.index(float(e)),errHigh)
                        finalHistoDict['nom'+kind+par+s+'error'].SetPointEYlow(self.etaBinning.index(float(e)),errLow)
                        finalHistoDict['nom'+kind+par+s+'error'].SetPointEXhigh(self.etaBinning.index(float(e)),finalHistoDict['nom'+kind+par+s].GetBinWidth(self.etaBinning.index(float(e))+1)/2)
                        finalHistoDict['nom'+kind+par+s+'error'].SetPointEXlow(self.etaBinning.index(float(e)),finalHistoDict['nom'+kind+par+s].GetBinWidth(self.etaBinning.index(float(e))+1)/2)

                    finalHistoDict['nom'+kind+par+s+'error'].SetLineColor(8)
                    finalHistoDict['nom'+kind+par+s+'error'].SetMarkerColor(8)
                    finalHistoDict['nom'+kind+par+s+'error'].SetFillColor(8)
                    finalHistoDict['nom'+kind+par+s+'error'].SetFillStyle(3003+fillSyle)
                    finalHistoDict['nom'+kind+par+s+'error'].SetLineWidth(1)
                    finalHistoDict['nom'+kind+par+s].Draw()
                    finalHistoDict['nom'+kind+par+s].SetLineWidth(3)
                    finalHistoDict['nom'+kind+par+s].GetXaxis().SetTitle('#eta')
                    finalHistoDict['nom'+kind+par+s].GetYaxis().SetTitle('Value of parameter')
                    finalHistoDict['nom'+kind+par+s+'error'].Draw("SAME 0P5")

                    finalCanvasDict['parameters'+'nom'+kind+par+s].Update()
                    finalCanvasDict['parameters'+'nom'+kind+par+s].Modified()
                    fillSyle =fillSyle+1

        #ratio plot creation
        for s in self.signList :
            for kind in ParnameDict : 
                for par in ParnameDict[kind] : 
                            c_ratioSyst = ROOT.TCanvas("c_parameters_ratioSyst_{sign}_{kind}_{par}".format(sign=s,kind=kind,par=par),"c_parameters_ratioSyst_{sign}_{kind}_{par}".format(sign=s,kind=kind,par=par),800,600)
                            c_ratioSyst.cd()
                            c_ratioSyst.SetGridx()
                            c_ratioSyst.SetGridy()
                            finalLegDict[s+kind+par+"ratioSyst"] = ROOT.TLegend(0.1,0.7,0.48,0.9)

                            sameFlag = True
                            colorNumber = 1
                            colorList = [600,616,416,632,432,800,900]
                            colorCounter = 0

                            for sKind, sList in systDict.iteritems():
                                colorNumber = colorList[colorCounter]
                                colorCounter = colorCounter+1
                                for sName in sList :
                                    colorNumber = colorNumber-2
                                    if colorNumber < colorList[colorCounter-1]-10 :
                                        colorNumber = colorList[colorCounter]+2
                                    finalHistoDict[sName+kind+par+s+'ratio'] = ROOT.TH1F(canvas+'_'+finalHistoDict[sName+kind+par+s].GetName()+'_'+sName+'_ratio',canvas+'_'+finalHistoDict[sName+kind+par+s].GetName()+'_'+sName+'_ratio',len(self.etaBinning)-1, array('f',self.etaBinning))
                                    finalHistoDict[sName+kind+par+s+'ratio'].Divide(finalHistoDict[sName+kind+par+s],finalHistoDict['nom'+kind+par+s],1,1)
                                    if noratio :
                                        finalHistoDict[sName+kind+par+s+'ratio'] = finalHistoDict[sName+kind+par+s].Clone()

                                    finalHistoDict[sName+kind+par+s+'ratio'].SetLineColor(colorNumber)
                                    for b in range(1,finalHistoDict[sName+kind+par+s+'ratio'].GetNbinsX()+1) :
                                        finalHistoDict[sName+kind+par+s+'ratio'].SetBinError(b,0)
                                    c_ratioSyst.cd()
                                    if sameFlag :
                                        finalHistoDict[sName+kind+par+s+'ratio'].Draw()
                                        finalHistoDict[sName+kind+par+s+'ratio'].GetYaxis().SetRangeUser(0.6,1.7)
                                    else :
                                        finalHistoDict[sName+kind+par+s+'ratio'].Draw("SAME")
                                        finalHistoDict[sName+kind+par+s+'ratio'].GetYaxis().SetRangeUser(0,3)
                                    finalHistoDict[sName+kind+par+s+'ratio'].GetXaxis().SetTitle('#eta')
                                    finalHistoDict[sName+kind+par+s+'ratio'].GetYaxis().SetTitle('Syst/Nom')
                                    finalHistoDict[sName+kind+par+s+'ratio'].SetLineWidth(3)
                                    sameFlag=False
                                    finalLegDict[s+kind+par+"ratioSyst"].AddEntry(finalHistoDict[sName+kind+par+s+'ratio'], sName)

                            #nominal
                            finalHistoDict['nom'+kind+par+s+'ratio'] = ROOT.TH1F(canvas+'_'+finalHistoDict['nom'+kind+par+s].GetName()+'_'+'nom'+'_ratio',canvas+'_'+finalHistoDict['nom'+kind+par+s].GetName()+'_'+'nom'+'_ratio',len(self.etaBinning)-1, array('f',self.etaBinning))
                            finalHistoDict['nom'+kind+par+s+'ratio'].Divide(finalHistoDict['nom'+kind+par+s],finalHistoDict['nom'+kind+par+s],1,1)
                            if noratio :
                                finalHistoDict['nom'+kind+par+s+'ratio'] =finalHistoDict['nom'+kind+par+s].Clone()

                            finalHistoDict['nom'+kind+par+s+'ratio'].SetLineColor(1)
                            for b in range(1,finalHistoDict['nom'+kind+par+s+'ratio'].GetNbinsX()+1) :
                                # finalHistoDict['nom'+kind+par+s+'ratio'].SetBinError(b,math.sqrt(2)*finalHistoDict['nom'+kind+par+s].GetBinError(b)/finalHistoDict['nom'+kind+par+s].GetBinContent(b))
                                if finalHistoDict['nom'+kind+par+s].GetBinContent(b)!=0 :
                                    finalHistoDict['nom'+kind+par+s+'ratio'].SetBinError(b,finalHistoDict['nom'+kind+par+s].GetBinError(b)/finalHistoDict['nom'+kind+par+s].GetBinContent(b))
                                else :
                                     print "WARNING: bin content of ", canvas, histo, "is 0, bin (s,e,p)=", s,e,b, ". error of ratio plot not correctly normalized."
                                     finalHistoDict['nom'+kind+par+s+'ratio'].SetBinError(b,finalHistoDict['nom'+kind+par+s].GetBinError(b))
                            c_ratioSyst.cd()
                            finalHistoDict['nom'+kind+par+s+'ratio'].Draw("SAME E2")
                            finalHistoDict['nom'+kind+par+s+'ratio'].GetYaxis().SetRangeUser(0,3)
                            finalHistoDict['nom'+kind+par+s+'ratio'].SetLineWidth(0)
                            finalHistoDict['nom'+kind+par+s+'ratio'].SetFillColor(1)
                            finalHistoDict['nom'+kind+par+s+'ratio'].SetFillStyle(3002)
                            finalHistoDict['nom'+kind+par+s+'ratio'].GetXaxis().SetTitle('#eta')
                            finalHistoDict['nom'+kind+par+s+'ratio'].GetXaxis().SetTitle('Syst/Nom')
                            finalLegDict[s+kind+par+"ratioSyst"].AddEntry(finalHistoDict['nom'+kind+par+s+'ratio'], 'Nominal')
                            
                            finalLegDict[s+kind+par+"ratioSyst"].Draw("SAME")
                            finalCanvasDict['parameters'+'ratio'+kind+par+s] = c_ratioSyst

        # outputFinal.Close()
        print "writing comparison syst plots..."

            
        outputFinal = ROOT.TFile(outDir+"/final_plots"+self.corrFitSuff+statAnaSuff+self.nameSuff+".root","recreate")

        dirFinalDict = {}
        for s in self.signList :
            for e in self.etaBinningS :
                dirFinalDict[s+e] =    outputFinal.mkdir(s+"_eta"+e)
                # dirFinalDict[s+e].cd()   #DECOMMENTA QUESTO SE NON FUNZIONANO LE DIRECTORY
            dirFinalDict[s+'unrolled'] =    outputFinal.mkdir(s+'_unrolled')
            dirFinalDict[s+'parameters'] = outputFinal.mkdir(s+'_Fit_parameters')
            
        for ind, obj in finalCanvasDict.iteritems():
                for s in self.signList :
                    
                    if ind.startswith('parameters') :
                        dirFinalDict[s+'parameters'].cd()
                        obj.Write()                        
                    elif 'unrolled' in ind and s in ind:
                        dirFinalDict[s+'unrolled'].cd()
                        obj.Write()
                    else :
                        for e in self.etaBinningS :
                            if ind.endswith(s+e) :
                                dirFinalDict[s+e].cd()
                                obj.Write()


    def main_analysis(self, correlatedFit=False,template=False,output4Plots=False, produce_ext_output=True,extrapSuff='') :
        
        ''' 
        workflow:
        - get the histogram 3D in a dictionary  [sample+region]  sample=(WToMuNu,Data)  Region=(A,B,C,D)
        - COMMENTED: get the Iso SF
        - differential_fakerate doing C/(C+A) or D/(B+D) for fake and prompt. EWKSF=1. 
        - fakerateFitter :
            - slicing to have TH1 to fit
            - prompt: fit erf and error prop. with toys
            - fake: linear fit 
            - if correlated: additional correlatedFitter
        - template produce the template (required only for validation and intermediate plots)
        - output using dict2histConverter 
        '''
        
        print "> getting histograms..."
        regionList = ['A','B','C','D']
        regionList_name = ['Sideband_aiso','Signal_aiso','Sideband','Signal']
        var_inside_histo = 'template'
        histo3DDict = {}
        sName = self.systName
        if self.systName !='' : 
            sName = '_'+sName
        for f in self.sampleList :
            for r, rl in map(None,regionList,regionList_name) :
                if rl=='Sideband' or rl=='Sideband_aiso' : #add extrapolation suffix (mapped in looseCutDict)
                    rl = rl+extrapSuff
                if f=='Data' :
                    histo3DDict[f+r] = self.rootFiles[self.sampleList.index(f)].Get('templates_'+rl+'/Nominal/'+var_inside_histo)
                else :
                    histo3DDict[f+r] = self.rootFiles[self.sampleList.index(f)].Get('templates_'+rl+'/'+self.systKind+'/'+var_inside_histo+sName)
                
        print "> evaluating fakerate..."
        hfakes = self.differential_fakerate(kind='fake',histo3D=histo3DDict)
        hprompt = self.differential_fakerate(kind='prompt',histo3D=histo3DDict)
        
        print "> fitting fakerate..."
        fakeDict = self.fakerateFitter(kind='fake',fake3D=hfakes,correlatedFit=correlatedFit)
        promptDict = self.fakerateFitter(kind='prompt',fake3D=hprompt)
                               
        if template:
            print "> evaluating template QCD..."
            templQCD = self.template(kind='fake',fakeDict=fakeDict,promptDict=promptDict, histo3D=histo3DDict,correlatedFit=correlatedFit)
            templW = self.template(kind='prompt',fakeDict=fakeDict,promptDict=promptDict, histo3D=histo3DDict)
        
        if output4Plots :
            print "> saving fakerate..."
            output = ROOT.TFile(self.outdir+"/bkg_differential_fakerate"+self.corrFitSuff+self.nameSuff+".root","recreate")
            fakerate_dir = output.mkdir("Fakerate")
            fakerate_dir.cd()        
            for s in self.signList :
                for e in self.etaBinningS :
                    fakeDict[s+e].Write()
                    promptDict[s+e].Write() 
                    
            if template:
                print "> saving template..."
                template_dir = output.mkdir("Template")
                template_dir.cd()
                for s in self.signList :
                    for e in self.etaBinningS :
                        templQCD[s+e].Write()
                        templW[s+e].Write()
                        templQCD[s+e+'fake'].Write()
                        templQCD[s+e+'prompt'].Write()   
                                    
            output.Save()
            
        if produce_ext_output :
            external_histoDict = self.dict2histConverter(fakedict=fakeDict,promptdict=promptDict,hdict=histo3DDict, correlatedFit=correlatedFit)
            external_output = ROOT.TFile(self.outdir+"/bkg_parameters_file"+self.corrFitSuff+self.nameSuff+".root","recreate")
            external_output.cd()
            for i in external_histoDict :
                external_histoDict[i].Write()
            external_output.Close()   
    
      
    def differential_fakerate(self, kind, histo3D) :    
       
        kindDict = {
            'fake' : 'Data',
            'prompt' : 'WToMuNu'
        }
        
        signBinning = [-2,0,2]
        if kind =='prompt' :
            h3Num =  histo3D[kindDict[kind]+'D'].Clone(kind+'Num') #D
            h3Den = histo3D[kindDict[kind]+'D'].Clone(kind+'Den') #D
            h3Den.Add(histo3D[kindDict[kind]+'B']) #D+B
            
        if kind == 'fake' :
            h3Num =  histo3D[kindDict[kind]+'C'].Clone(kind+'Num') #C
            h3NumW = histo3D['WToMuNu'+'C'].Clone('WToMuNu_diff_Num') #C_w
            h3Num.Add(h3NumW,-1) #C-C_w
            h3Den = histo3D[kindDict[kind]+'C'].Clone(kind+'Den') #C
            h3Den.Add(histo3D[kindDict[kind]+'A']) #C+A
            h3DenW = histo3D['WToMuNu'+'C'].Clone('WToMuNu_diff_Den') #C_w
            h3DenW.Add(histo3D['WToMuNu'+'A']) #(C+A)_w
            
            h3Den.Add(h3DenW,-1) #C+A-(C+A)_w
        h3Fakes = h3Num.Clone("h3Fakes_{kind}") 
        h3Fakes.Divide(h3Num,h3Den,1,1,'B') #num/den
        return h3Fakes  
    
    
    def fakerateFitter(self, kind,fake3D,correlatedFit=False) :
        h1Dict = {}
        parDict = {}
        for s in self.signList :
            #slicing
            if s=='Minus' :
                binS = 1
            else :
                binS = 2
            for e in self.etaBinningS :
                etaBinN= self.binNumb_calculator(fake3D,'X',self.etaBinning[self.etaBinningS.index(e)])
                h1Dict[s+e] = fake3D.ProjectionY('hFakes_pt_'+kind+'_'+s+'_'+e,etaBinN,etaBinN,binS,binS,"e")
                    
                #do the fit
                if kind == 'prompt' :
                    fitFake = ROOT.TF1("fitFake", "[0]*erf([1]*(x)+[2])",0,100,3)
                    fitFake.SetParameters(1,0.1,-3)
                    fitFake.SetParNames("offset","slope",'2deg')              
                if kind == 'fake' :
                    fitFake = ROOT.TF1("fitFake", '[0]+[1]*(x-26)',0,100,2) 
                    fitFake.SetParameters(0.5,0.1)
                    fitFake.SetParNames("offset","slope")
                
                fit_result = h1Dict[s+e].Fit(fitFake,"QS","",26,self.maxPt_linearFit)
                
                #assign the fit results
                cov = fit_result.GetCovarianceMatrix()
                parDict[s+e+'offset']=fitFake.GetParameter(0)
                parDict[s+e+'slope']=fitFake.GetParameter(1)
                parDict[s+e+'offsetErr']=fitFake.GetParError(0)
                parDict[s+e+'slopeErr']=fitFake.GetParError(1)
                parDict[s+e+'offset*slope'] = ROOT.TMatrixDRow(cov,0)(1) #covarinace
                parDict[s+e+'chi2red'] =fitFake.GetChisquare()/fitFake.GetNDF()

                if kind =='prompt' :
                    parDict[s+e+'2deg']=fitFake.GetParameter(2)
                    parDict[s+e+'2degErr']=fitFake.GetParError(2)
                    parDict[s+e+'offset*2deg'] = ROOT.TMatrixDRow(cov,0)(2)
                    parDict[s+e+'slope*2deg'] = ROOT.TMatrixDRow(cov,1)(2)
                    parDict[s+e+'fitRes'] = cov
                    parDict[s+e+'sigmaERFvec'] = self.evaluate_erf_error(ERFfitResult=cov,ERFp0=parDict[s+e+'offset'],ERFp1=parDict[s+e+'slope'],ERFp2=parDict[s+e+'2deg'])
                else:  #only to have a coherent filling
                    parDict[s+e+'2deg']=0
                    parDict[s+e+'2degErr']=0
                    parDict[s+e+'offset*2deg'] = 0
                    parDict[s+e+'slope*2deg'] = 0
                    parDict[s+e+'fitRes'] = 0
                    parDict[s+e+'sigmaERFvec'] = 0

        #for compatibility add the 1Dhistos to the parDict
        parDict.update(h1Dict)        
        
        if correlatedFit :
            parCorFit = self.correlatedFitter(parDict)
            parDict.update(parCorFit)    
        
        
        return parDict
    
    
    def template(self,kind,fakeDict,promptDict, histo3D,correlatedFit=False) :
        
        kindDict = {
            'fake' : 'Data',
            'prompt' : 'WToMuNu',
        }
        shiftX=26
        
        htempl = {}
        h1Dict = {}
        for s in self.signList :
            if s=='Minus' :
                binS = 1
            else :
                binS = 2
            for e in self.etaBinningS :
                if kind =='prompt' :
                    etaBinN= self.binNumb_calculator( histo3D[kindDict[kind]+'D'],'X',self.etaBinning[self.etaBinningS.index(e)])    
                    h1Dict[s+e] = histo3D[kindDict[kind]+'D'].ProjectionY('htempl_pt_'+kind+'_'+s+'_'+e,etaBinN,etaBinN,binS,binS,"e")

                if kind=='fake' :
                    etaBinN= self.binNumb_calculator(histo3D[kindDict[kind]+'D'],'X',self.etaBinning[self.etaBinningS.index(e)]) 
                    hDregion = histo3D[kindDict[kind]+'D'].ProjectionY('D_pt_'+kind+'_'+s+'_'+e,etaBinN,etaBinN,binS,binS,"e")
                    hBregion = histo3D[kindDict[kind]+'B'].ProjectionY('B_pt_'+kind+'_'+s+'_'+e,etaBinN,etaBinN,binS,binS,"e")

                    htempl_pt = ROOT.TH1F("htempl_pt_{kind}_{sign}_{eta}".format(kind=kind,sign=s,eta=e),"htempl_pt_{kind}_{sign}_{eta}".format(kind=kind,sign=s,eta=e), len(self.ptBinning)-1, array('f',self.ptBinning) )
                    htempl_fake_pt = ROOT.TH1F("htempl_fake_pt_{kind}_{sign}_{eta}".format(kind=kind,sign=s,eta=e),"htempl_fake_pt_{kind}_{sign}_{eta}".format(kind=kind,sign=s,eta=e), len(self.ptBinning)-1, array('f',self.ptBinning) )
                    htempl_prompt_pt = ROOT.TH1F("htempl_prompt_pt_{kind}_{sign}_{eta}".format(kind=kind,sign=s,eta=e),"htempl_prompt_pt_{kind}_{sign}_{eta}".format(kind=kind,sign=s,eta=e), len(self.ptBinning)-1, array('f',self.ptBinning) )

                    for p in self.ptBinningS :
        
                        #fake rate eval:
                        if correlatedFit :
                            CF_string = 'Minuit' 
                        else :
                            CF_string = ''
                        xf = fakeDict[s+e].GetBinCenter(self.ptBinningS.index(p)+1)
                        fr = fakeDict[s+e+'offset'+CF_string]+(xf-shiftX)*fakeDict[s+e+'slope'+CF_string] #SLOPEOFFSETUNCORR
                        dx_f =0
                        dfr = math.sqrt(fakeDict[s+e+'offset'+CF_string+'Err']**2+(dx_f**2)*(fakeDict[s+e+'slope'+CF_string]**2)+((xf-shiftX)**2)*(fakeDict[s+e+'slope'+CF_string+'Err']**2)+2*(xf-shiftX)*fakeDict[s+e+'offset*slope'+CF_string])

                        #prompt rate eval:
                        xp = promptDict[s+e].GetBinCenter(self.ptBinningS.index(p)+1)
                        pr = promptDict[s+e+'offset']*erf(promptDict[s+e+'slope']*xp+promptDict[s+e+'2deg'])
                        dpr = promptDict[s+e+'sigmaERFvec'][self.ptBinningS.index(p)]#*math.sqrt(promptDict[s+e+'chi2red'])
                        # print "WARNING: not applied x1.5 to the ERF fit error."
                        
                        #warning messages
                        if pr>1:
                            print "WARNING!!!!!!, pr>1", pr, fr, "datakind=",kind, ", with eta, pt, sign =", e, p, s,
                            print ", not fitted value=", promptDict[s+e].GetBinContent(self.ptBinningS.index(p)+1)
                            pr = 1
                        if fr>1 :
                            print "WARNING!!!!!!, fr>1,", pr, fr, "datakind=",kind, ", with eta, pt, sign =", e, p, s,
                            print ", not fitted value=", fakeDict[s+e].GetBinContent(self.ptBinningS.index(p)+1)
                            fr = 0
                        
                        #template evaluation
                        NTight= hDregion.GetBinContent(self.ptBinningS.index(p)+1)
                        NNotTight= hBregion.GetBinContent(self.ptBinningS.index(p)+1)
                        NTightErr= hDregion.GetBinError(self.ptBinningS.index(p)+1)
                        NNotTightErr= hBregion.GetBinError(self.ptBinningS.index(p)+1)
                        
                        scaleTight = -fr*(1-pr)/(pr-fr)
                        scaleNotTight = fr*pr/(pr-fr)
                        NQCD=NTight*scaleTight+NNotTight*scaleNotTight
                        NQCDErr = math.sqrt((fr**4*(dpr**2*NTight**2+2*NNotTight*dpr**2*NTight+NNotTight**2*dpr**2+NNotTightErr**2*pr**2+NTightErr**2*(pr-1)**2)-2*fr**(3)*(dpr**2*NTight**2+NNotTight*dpr**2*NTight+(NNotTightErr**2*pr**2+NTightErr**2*(pr-1)**2)*pr)+fr**2*(dpr**2*NTight**2+(NNotTightErr**2*pr**2+NTightErr**2*(pr-1)**2)*pr**2)+dfr**2*pr**2*((pr-1)*NTight+NNotTight*pr)**2)/((fr-pr)**4))
                        htempl_pt.SetBinContent(self.ptBinningS.index(p)+1,NQCD)
                        htempl_pt.SetBinError(self.ptBinningS.index(p)+1,NQCDErr)
                        
                        htempl_fake_pt.SetBinContent(self.ptBinningS.index(p)+1,fr)
                        htempl_fake_pt.SetBinError(self.ptBinningS.index(p)+1,dfr)
                        htempl_prompt_pt.SetBinContent(self.ptBinningS.index(p)+1,pr)
                        htempl_prompt_pt.SetBinError(self.ptBinningS.index(p)+1,dpr)

                        htempl[s+e] = htempl_pt
                        htempl[s+e+'fake'] = htempl_fake_pt
                        htempl[s+e+'prompt'] = htempl_prompt_pt       
            
        if kind =='prompt' :
            htempl.update( h1Dict)
            
        return htempl
        
            
    def extrapolationSyst(self,extrapDict, linearFit=True) :
        #linearFit = trend in Mt, if False-->constant trend in Mt assumed
        
        nomMeanMt=67
        kindDict = {
            'paramExtrap' : True #extrapolation of two separated parameteters
        }
        looseCutBinning = [] 
        fileDict = {}
        for lcut, lbin in looseCutDict.iteritems() :
            fileDict[lcut] = ROOT.TFile.Open(self.outdir+'/bkg_'+lcut+'/bkg_differential_fakerate.root')
            fileDict[lcut+'par'] = ROOT.TFile.Open(self.outdir+'/bkg_'+lcut+'/bkg_parameters_file.root')
            looseCutBinning.append(lbin[0])
        fileDict[self.systName] = ROOT.TFile.Open(self.outdir+'/bkg_'+self.systName+'/bkg_differential_fakerate.root')
        fileDict[self.systName+'par'] = ROOT.TFile.Open(self.outdir+'/bkg_'+self.systName+'/bkg_parameters_file.root')
        looseCutBinning.append(40)
        looseCutBinning.append(2*nomMeanMt-40)
        looseCutBinning.sort() #because from dict 
        
        looseCutDict.update({self.systName:[40,2*nomMeanMt-40]})
        
        histoDict = {}
        discrepancyMapDict = {}
        paramMapDict = {}
        for s in self.signList :        
            discrepancyMapDict[s] = ROOT.TH2F("discrepancyMap_{sign}".format(sign=s),"discrepancyMap_{sign}".format(sign=s),len(self.etaBinning)-1, array('f',self.etaBinning), len(self.ptBinning)-1, array('f',self.ptBinning) )
            discrepancyMapDict[s].GetXaxis().SetTitle("#eta")
            discrepancyMapDict[s].GetYaxis().SetTitle("p_{T} [GeV]")
            discrepancyMapDict[s].GetZaxis().SetTitle("(extrap-nom)/nom")
            
            discrepancyMapDict[s+'Err'] = ROOT.TH2F("discrepancyMapErr_{sign}".format(sign=s),"discrepancyMapErr_{sign}".format(sign=s),len(self.etaBinning)-1, array('f',self.etaBinning), len(self.ptBinning)-1, array('f',self.ptBinning) )
            discrepancyMapDict[s+'Err'].GetXaxis().SetTitle("#eta")
            discrepancyMapDict[s+'Err'].GetYaxis().SetTitle("p_{T} [GeV]")
            discrepancyMapDict[s+'Err'].GetZaxis().SetTitle("#sigma (extrap-nom)/nom")
           
            if s=='Plus' :
                chi2_extrap = ROOT.TH2F("chi2_extrap","chi2_extrap",len(self.etaBinning)-1, array('f',self.etaBinning), len(self.signList), array('f',[0,1,2]))
                chi2_extrap.GetXaxis().SetTitle("#eta")
                chi2_extrap.GetYaxis().SetTitle("sign (1st=plus,2nd=minus)")
                chi2_extrap.GetZaxis().SetTitle("reduced #chi^{2}") 
           
            for par in ['offset','slope','offset*slope'] :
                paramMapDict[s+par] = ROOT.TH1F("paramMap_{sign}_{par}".format(sign=s,par=par),"paramMap_{sign}_{par}".format(sign=s,par=par),len(self.etaBinning)-1, array('f',self.etaBinning))
                paramMapDict[s+par].GetXaxis().SetTitle("#eta")
                paramMapDict[s+par].GetYaxis().SetTitle("Parameter value")
                        
            for e in self.etaBinningS :
                dirString = 'Fakerate/hFakes_pt_fake_'
                histoDict[s+e] = ROOT.TH2F("paramExtrap_{sign}_{eta}".format(sign=s,eta=e),"paramExtrap_{sign}_{eta}".format(sign=s,eta=e),len(looseCutBinning)-1, array('f',looseCutBinning), len(self.ptBinning)-1, array('f',self.ptBinning) )
                histoDict[s+e].GetXaxis().SetTitle("M_{T} [GeV]")
                histoDict[s+e].GetYaxis().SetTitle("p_{T} [GeV]")
                histoDict[s+e].GetZaxis().SetTitle("f (iso. cut eff.)")
                for p in self.ptBinningS :
                    for lcut, lbin in looseCutDict.iteritems() :
                        if lcut==self.systName :
                            dirString = 'Template/htempl_fake_pt_fake_'
                        val = fileDict[lcut].Get(dirString+s+'_'+e).GetBinContent(self.ptBinningS.index(p)+1) 
                        err = fileDict[lcut].Get(dirString+s+'_'+e).GetBinError(self.ptBinningS.index(p)+1)
                        histoDict[s+e].SetBinContent( histoDict[s+e].FindBin(lbin[0]),self.ptBinningS.index(p)+1,val)
                        histoDict[s+e].SetBinError( histoDict[s+e].FindBin(lbin[0]),self.ptBinningS.index(p)+1,err)
                if linearFit :
                    fitFunc = ROOT.TF2("fitFunc", "[0]+[1]*x+[2]*y+[3]*x*y",0.,40.,26,55) 
                else :
                    fitFunc = fitFunc = ROOT.TF2("fitFunc", "[0]+[1]*y",0.,40.,26,55) 
                    
                fitRes = histoDict[s+e].Fit(fitFunc,"QSR","")#,0,40,26,55)
                
                if linearFit :
                    cov = fitRes.GetCovarianceMatrix()
                    a = fitFunc.GetParameter(0)
                    b = fitFunc.GetParameter(1)
                    c = fitFunc.GetParameter(2)
                    d = fitFunc.GetParameter(3)
                    wa = fitFunc.GetParError(0)
                    wb = fitFunc.GetParError(1)
                    wc = fitFunc.GetParError(2)
                    wd = fitFunc.GetParError(3)
                    wab = ROOT.TMatrixDRow(cov,0)(1)
                    wac = ROOT.TMatrixDRow(cov,0)(2)
                    wad = ROOT.TMatrixDRow(cov,0)(3)
                    wbc = ROOT.TMatrixDRow(cov,1)(2)
                    wbd = ROOT.TMatrixDRow(cov,1)(3)
                    wcd = ROOT.TMatrixDRow(cov,2)(3)
                    xf=nomMeanMt
                    for p in self.ptBinningS :
                        yf = histoDict[s+e].GetYaxis().GetBinCenter(self.ptBinningS.index(p)+1)
                        fr = fitFunc.Eval(xf)
                        dfr = math.sqrt(wa**2+(wb*xf)**2+(wc*yf)**2+(wd*xf*yf)**2+2*(wab*xf+wac*yf+wad*xf*yf+wbc*xf*yf+wbd*xf**2*yf+wcd*xf*yf**2))
                        fr_nom =histoDict[s+e].GetBinContent(histoDict[s+e].GetNbinsX(),self.ptBinningS.index(p)+1)  
                        err_nom = histoDict[s+e].GetBinError(histoDict[s+e].GetNbinsX(),self.ptBinningS.index(p)+1)
                        discr = (fr_nom - fr)/fr_nom
                        discrepancyMapDict[s].SetBinContent(self.etaBinningS.index(e)+1,self.ptBinningS.index(p)+1,discr) 
                        err = 1/(fr_nom**2)*math.sqrt((dfr*fr_nom)**2+(fr*err_nom)**2)
                        discrepancyMapDict[s].SetBinError(self.etaBinningS.index(e)+1,self.ptBinningS.index(p)+1,err)
                        discrepancyMapDict[s+'Err'].SetBinContent(self.etaBinningS.index(e)+1,self.ptBinningS.index(p)+1,err)
                    
                    paramMapDict[s+'slope'].SetBinContent(self.etaBinningS.index(e)+1,b)    
                    paramMapDict[s+'slope'].SetBinError(self.etaBinningS.index(e)+1,wb)
                    paramMapDict[s+'offset'].SetBinContent(self.etaBinningS.index(e)+1,a)    
                    paramMapDict[s+'offset'].SetBinError(self.etaBinningS.index(e)+1,wa)
                    paramMapDict[s+'offset*slope'].SetBinContent(self.etaBinningS.index(e)+1,wab)    
                    paramMapDict[s+'offset*slope'].SetBinError(self.etaBinningS.index(e)+1,0)        
                    
                else : #contant in mt fit
                    cov = fitRes.GetCovarianceMatrix()
                    dmq2 = ROOT.TMatrixDRow(cov,0)(1)
                    q = fitFunc.GetParameter(0)
                    m = fitFunc.GetParameter(1)
                    dq = fitFunc.GetParError(0)
                    dm = fitFunc.GetParError(1)             
                    xf=nomMeanMt
                    for p in self.ptBinningS :
                        yf = histoDict[s+e].GetYaxis().GetBinCenter(self.ptBinningS.index(p)+1)
                        fr = fitFunc.Eval(xf)
                        dfr = math.sqrt(dq**2+(yf*dm)**2+2*yf*dmq2)
                        fr_nom =histoDict[s+e].GetBinContent(histoDict[s+e].GetNbinsX(),self.ptBinningS.index(p)+1)  
                        err_nom = histoDict[s+e].GetBinError(histoDict[s+e].GetNbinsX(),self.ptBinningS.index(p)+1)
                        discr = (fr_nom - fr)/fr_nom
                        discrepancyMapDict[s].SetBinContent(self.etaBinningS.index(e)+1,self.ptBinningS.index(p)+1,discr) 
                        err = 1/(fr_nom**2)*math.sqrt((dfr*fr_nom)**2+(fr*err_nom)**2)
                        discrepancyMapDict[s].SetBinError(self.etaBinningS.index(e)+1,self.ptBinningS.index(p)+1,err)
                        discrepancyMapDict[s+'Err'].SetBinContent(self.etaBinningS.index(e)+1,self.ptBinningS.index(p)+1,err)
                        
                chi2_extrap.SetBinContent(self.etaBinningS.index(e)+1,self.signList.index(s)+1,fitFunc.GetChisquare()/fitFunc.GetNDF())
                chi2_extrap.SetBinError(self.etaBinningS.index(e)+1,self.signList.index(s)+1,math.sqrt(2*fitFunc.GetNDF())/fitFunc.GetNDF())
        
        if linearFit :
            linString = "_linearFit" 
        else :
            linString = "_constFit"                
        output = ROOT.TFile(self.outdir+"/extrapolation_plots_"+str(nomMeanMt)+'GeV_'+linString+".root","recreate")
        for s in self.signList :
            for e in self.etaBinningS : 
                histoDict[s+e].Write()  
            discrepancyMapDict[s].Write()
            discrepancyMapDict[s+'Err'].Write()
            if s=='Plus' : 
                chi2_extrap.Write()
            for par in ['offset','slope','offset*slope'] :
                paramMapDict[s+par].Write()


    
    def buildOutput(self,outputDir,statAna=True) :
        
        CFsuff=''
        statAnaSuff=''
        if self.correlatedFit :
            CFsuff = '_CF'
        if statAna :
            statAnaSuff='statAna'
        
        output = ROOT.TFile(outputDir+"/bkg_parameters"+CFsuff+statAnaSuff+self.nameSuff+".root","recreate")
        
        fileDict = {}
        dirDict = {}
        
        systDict = copy.deepcopy(bkg_utils.bkg_systematics)
        systDict['Nominal'] = ['']
        
        for sKind, sList in systDict.iteritems():
            if sKind =='Nominal' : 
                SAsuff = statAnaSuff
            else :
                SAsuff = ''
            for sName in sList :
                dirDict[sKind+'_'+sName] = output.mkdir(sKind+'_'+sName)
                fileDict[sKind+'_'+sName] = ROOT.TFile.Open(outputDir+"/bkg_"+sName+"/bkg_parameters_file"+CFsuff+SAsuff+self.nameSuff+".root")
                for key in fileDict[sKind+'_'+sName].GetListOfKeys():
                    histo = fileDict[sKind+'_'+sName].Get(key.GetName())
                    dirDict[sKind+'_'+sName].cd()
                    histo.Write()
