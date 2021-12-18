from __future__ import print_function 
import ROOT
from array import array
import math
import copy
import argparse
import os
import numpy as np 
from obsminimization import pmin
import root_numpy as rootnp
import ctypes
import scipy
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)
np.set_printoptions(linewidth=np.inf,precision=4, threshold=np.inf)

import jax
import jax.numpy as jnp


def getFitRes(inFile, coeffList) :
    outDict = {}
    openFile = ROOT.TFile.Open(inFile)
    for li in coeffList :
        if 'unpol' in li[0] :
            suff = '_4apo'
        else :
            suff = ''
        outDict[li[0]] = openFile.Get('coeff2D_reco/recoFitAC'+li[0]+suff)
        outDict[li[0]+'qt'] = openFile.Get('coeff2D_reco/recoFitACqt'+li[0]+suff)
        outDict[li[0]+'y'] = openFile.Get('coeff2D_reco/recoFitACy'+li[0]+suff)
    
    ext = 'unpolarizedxsec'
    outDict[ext] = openFile.Get('coeff2D_reco/recoFitAC'+ext+'_4apo')
    outDict[ext+'qt'] = openFile.Get('coeff2D_reco/recoFitACqt'+ext+'_4apo')
    outDict[ext+'y'] = openFile.Get('coeff2D_reco/recoFitACy'+ext+'_4apo')
    outDict['cov']=openFile.Get('matrices_reco/covariance_matrix_channelhelpois') #unrolled order: y0qt0a0.....y0qtNa0, y1qt0a0...y1qt0a0, ...
    outDict['cov1D']=openFile.Get('matrices_reco/covariance_matrix_channelhelmetapois') #unrolled order: y0qt0a0.....y0qtNa0, y1qt0a0...y1qt0a0, ...
    outDict['mass'] = openFile.Get('nuisance_reco/massreco')
    # outDict['mass'] = openFile.Get('matrices_reco/massreco')
    print("input W mass=", outDict['mass'].GetBinContent(1), "+/-", outDict['mass'].GetBinError(1))
    outDict['nui'+'LHEPdf'] = openFile.Get('nuisance_reco/reco_NuiConstrLHEPdfWeight')
    outDict['nui'+'WHSFSyst'] = openFile.Get('nuisance_reco/reco_NuiConstrWHSFStat')
    outDict['nui'+'LHEScale'] = openFile.Get('nuisance_reco/reco_NuiConstrLHEScaleWeight')
    outDict['nui'+'others'] = openFile.Get('nuisance_reco/reco_NuiConstrothers')
    return outDict 
    
def getRegFuncPar(inFile, coeffList,sign) :
    outDict = {}
    openFile = ROOT.TFile.Open(inFile)
    if sign=='plus' :
        sString = 'WtoMuP'
    if sign=='minus' :
        sString = 'WtoMuN'
    for li in coeffList :
        if li[0]!='unpolarizedxsec' : 
            outDict[li[0]] = openFile.Get('Y'+str(li[1])+'_Qt'+str(li[2])+'/hfit_'+sString+'_'+li[0]+'_'+str(li[1])+str(li[2]))
            funcReg = outDict[li[0]].GetFunction("fit_"+sString+'_'+li[0])
            outDict['par'+li[0]] = np.zeros((li[1],li[2]))
            for yy in range(0,li[1]) :
                for qq in range(0,li[2]) :
                    itot = li[1]*qq+yy
                    outDict['par'+li[0]][yy,qq] = funcReg.GetParameter(itot)
                    # print(li[0],itot, "y=", yy,"qt=",qq, outDict['par'+li[0]][yy,qq], "name=",funcReg.GetParName(itot) )
        else :
            temph = openFile.Get('unpol_Y3/hfit_'+sString+'_unpol_3')
            for qt in range(1, temph.GetNbinsY()+1) :
                outDict[li[0]+str(qt)] = openFile.Get('unpol_Y'+str(li[1])+'/hfit_'+sString+'_unpol_3_'+str(qt))
                funcReg = outDict[li[0]+str(qt)].GetFunction("fit_"+sString+'_unpol')
                outDict['par'+li[0]+str(qt)] = np.zeros((li[1]))
                for yy in range(0,li[1]) :
                    outDict['par'+li[0]+str(qt)][yy] = funcReg.GetParameter(yy)
            print("WARNING: unpolarized qt binning must be aligned to the fit qt binning, if fitted in Y!")
    return outDict

def chi2MassOnlyFit(modelPars,npFitRes,npCovMatInv) : #not used
    print("shapes",np.shape(modelPars), np.shape(npFitRes), np.shape(npCovMatInv))
    print("vals",modelPars, npFitRes, npCovMatInv, 50*math.sqrt(1/npCovMatInv)) #1/npnpCovMatInv??
    diff = npFitRes-modelPars
    # chi2 = jnp.matmul(diff.T, jnp.matmul(npCovMatInv, diff) )
    chi2 = diff*diff*npCovMatInv
    return chi2
    
def par2polyModel(modelPars, npFitRes, coeffList, npBinCenters,parNum,dimYQt,numNui) :#input = the parameters, output the polyinomial function
    fitModelList = []
    valY = npBinCenters[:,0]
    valQt = npBinCenters[:,1]
    npBinCenters_res = npBinCenters.reshape(dimYQt[0],dimYQt[1],2)
    valYunpol=npBinCenters_res[:,0,0]

    icoeff = 0
    ipar=0
    unpolMult=[1e5,1e5]    
    for li in coeffList :
        if not li[8] :
            if li[0]==5 :
                fitModelList_unpol = []
                for iqt in range(dimYQt[1]) :    
                    fitModelList_unpol.append(unpolMult[0]*modelPars[0+ipar]+unpolMult[1]*modelPars[1+ipar]*valYunpol+unpolMult[1]*modelPars[2+ipar]*np.power(valYunpol,2))#+unpolMult[1]*modelPars[2+ipar]*np.power(valYunpol,4))
                    ipar = ipar+parNum[icoeff]
                    icoeff+=1
                fitModelList_unpol_stack = jnp.stack(fitModelList_unpol,0)
                fitModelList.append(fitModelList_unpol_stack.T)
                        
            else :
                fitModelList_temp = []
                fitModelList_temp.append(np.zeros(jnp.shape(npBinCenters)[0])) # list of  (y,qt), with lenght (Ny*Nqt) 
                effY = li[1]
                effQt = li[2]
                if li[5] : effY-= 1
                if li[6] : effQt-= 1
                if li[7] : effQt-= 1
                for iqt in range(0, effQt) :
                    powQt = iqt
                    if li[6] : powQt+= 1 
                    if li[7] and li[0]!=4 : powQt+= 1 
                    if li[7] and li[0]==4 and iqt>0 : powQt+=1
                    for iy in range(0, effY) :
                        powY = iy
                        if li[5] : powY+= 1
                        # fitModelList_temp.append(fitModelList_temp[-1]+modelPars[iqt*(effY)+iy+ipar]*np.power(valY,powY)*np.power(valQt,powQt) )
                        fitModelList_temp.append(fitModelList_temp[-1]+modelPars[iy*(effQt)+iqt+ipar]*np.power(valY,powY)*np.power(valQt,powQt) )
                        # fitModelList_temp.append(fitModelList_temp[-1]+modelPars[iy*(effQt)+iqt+ipar]*np.power(valY,powY)*np.power(valQt,powQt)/(np.power(2.2,powY)*np.power(46,powQt)) )
                ipar = ipar+parNum[icoeff]
                icoeff+=1       
                fitModelList.append(fitModelList_temp[-1].reshape(dimYQt[0],dimYQt[1])) # reshape to (Ny x Nqt), like the fitRes histogram
        
        else : #li[8]=not regularize this
            fitModelList_temp = modelPars[ipar+0:ipar+dimYQt[0]*dimYQt[1]]*npFitRes[ipar+0:ipar+dimYQt[0]*dimYQt[1]]
            # fitModelList_temp = modelPars[ipar+0:ipar+dimYQt[0]*dimYQt[1]]
            ipar = ipar+parNum[icoeff]
            icoeff+=1
            fitModelList.append(fitModelList_temp.reshape(dimYQt[0],dimYQt[1]))
            
    fitModel = jnp.stack(fitModelList,0)
    fitModel = fitModel.flatten()
    # fitModel = jnp.append(fitModel,modelPars[-1])#mass
    # print("internal model",len(fitModel))
    #add nuisance
    if numNui> 0 :# = not nUL #DEBSIMP 
        fitModel = jnp.append(fitModel,modelPars[-numNui:])
    
    
    # fitModel = jnp.append(fitModel,npFitRes[-numNui:])
    # print("internal model after nuisance add",len(fitModel))
    
    return fitModel
     
def chi2PolyFit(modelPars,npFitRes, npCovMatInv, coeffList, npBinCenters,parNum,dimYQt) :
    numNui = len(npCovMatInv[len(coeffList)*dimYQt[0]*dimYQt[1]:])
    fitModel = par2polyModel(modelPars= modelPars, coeffList= coeffList, npBinCenters= npBinCenters,parNum= parNum,dimYQt= dimYQt,npFitRes=npFitRes,numNui=numNui)
    diff = npFitRes-fitModel
    chi2 = jnp.matmul(diff.T, jnp.matmul(npCovMatInv, diff) )
    # print("chi2=",chi2)

    return chi2
    

def polyFit(fitRes, coeffList,nUL,regPar,inPar='') : #inPar wins on regPar
    
    dimY = fitRes[coeffList[0][0]].GetNbinsX()
    dimQt = fitRes[coeffList[0][0]].GetNbinsY()
    dimCoeff = len(coeffList)
    dimYQt = [dimY,dimQt]

    npCovMatUncut = rootnp.hist2array(fitRes['cov'])
    # npCovMat = npCovMatUncut[0:dimCoeff*dimQt*dimY,0:dimCoeff*dimQt*dimY]
    if nUL :
        npCovMat = npCovMatUncut[dimQt*dimY:6*dimQt*dimY,dimQt*dimY:6*dimQt*dimY] #removed nuisance and UL #DEBSIMP
    else :
        npCovMat = npCovMatUncut 
        # npCovMat = npCovMatUncut
    # npCovMat = npCovMatUncut[1*dimQt*dimY:6*dimQt*dimY,1*dimQt*dimY:6*dimQt*dimY]
    
    #add mass to the convariance matrix
    # massBin = fitRes['cov'].GetXaxis().FindBin('mass')
    # massColumn = npCovMatUncut[massBin-1,0:dimCoeff*dimQt*dimY]
    # massRow = npCovMatUncut[0:dimCoeff*dimQt*dimY,massBin-1]
    # massEl = npCovMatUncut[massBin-1,massBin-1]
    # npCovMat = np.c_[npCovMat,massColumn]
    # massRow = np.append(massRow,massEl)
    # npCovMat = np.r_[npCovMat,[massRow]]
         
    # reweight the covariance matrix to have it in "GeV"
    # for i in range(0,np.shape(npCovMat)[0]) :
    #     for j in range(0,np.shape(npCovMat)[1]) :
    #         if i==np.shape(npCovMat)[0]-1 or j==np.shape(npCovMat)[0]-1 :
    #             if i==j :
    #                 npCovMat[i,j] = npCovMat[i,j]/100. #/100->GeV, *10000->MeV
    #             else :
    #                 npCovMat[i,j] = npCovMat[i,j]/10. #/10->GeV, *100->MeV
        

    npCovMatInv = np.linalg.inv(npCovMat)
    
    npFitResList = []
    for li in coeffList :     
        npFitResList.append(rootnp.hist2array(fitRes[li[0]])) 
    npFitRes = np.stack(npFitResList,0)  #dimensions: coeff, y, qt = (6,6,11). if flatten--> q0y0a0...qNy0a0, q0y1a0....qNy1a0, 
    # print(np.shape(npFitRes))
    # print(rootnp.hist2array(fitRes['A0']).flatten())
    # print(npFitRes)
    # npFitMass = np.zeros(1)
    # npFitMass[0] = fitRes['mass'].GetBinContent(1)    
    # qtBinCenters = np.array([1., 3., 5., 7., 9., 11., 13., 15., 18., 22., 28., 38., 60.])#Extended
    qtBinCenters = np.array([1., 3., 5., 7., 9., 11., 14., 18., 23., 30., 46.])#Extended-reduced

    print("WARNING: hardcoded qt bin centers")
    print("WARNING: hardcoded bin limits (in the integration,not used)")
    
    npBinCenters = np.zeros((dimQt*dimY,2))
    for iy in range(0, dimY): #y
        for iqt in range(0, dimQt):#qt
            itot = iy*(dimQt)+iqt
            # npBinCenters[itot] = fitRes[li[0]].GetXaxis().GetBinCenter(iy+1), fitRes[li[0]].GetYaxis().GetBinCenter(iqt+1) 
            npBinCenters[itot] = fitRes[li[0]].GetXaxis().GetBinCenter(iy+1), qtBinCenters[iqt] 
            
    #model init
    modelParsList = []
    parNum = []
    parCounter=0
    for li in coeffList :
        if not li[8] : 
            if li[0]=='unpolarizedxsec' : 
                print("Warning: hardcoded deg3 for Y unpolarized fit (in the chi2)")
                for qt in range(1, dimQt+1) :
                    modelParsC=np.zeros(li[1])
                    for yy in range(len(modelParsC)) :
                        modelParsC[yy] = 1.
                    modelParsList.append(modelParsC.copy())
                    parNum.append(len(modelParsC))
                    # print(li[0], "qt=",qt, "pars=", parNum[-1])
            else :
                effY = li[1]
                effQt = li[2]
                if li[5] : effY-= 1
                if li[6] : effQt-= 1
                if li[7] : effQt-= 1
                modelParsC=np.zeros((effY,effQt))
                for y in range(0, effY) :
                    for qt in range(0, effQt) :
                        if inPar!='' :
                            modelParsC[y,qt] = inPar[0][parCounter]
                            parCounter+=1
                        elif regPar!='':
                            pary = y 
                            parqt = qt
                            if li[5] : pary+=1
                            if li[6] : parqt+= 1 
                            if li[7] and li[0]!=4 : parqt+= 1 
                            if li[7] and li[0]==4 and qt>0 : parqt+=1
                            modelParsC[y,qt] = regPar['par'+li[0]][pary,parqt]
                            # print(li[0], "yeff=", y,"qteff=",qt, "y=",pary,"qt=",parqt,modelParsC[y,qt] )

                        else :
                            modelParsC[y,qt] = 1.
                modelParsList.append(modelParsC.copy())
                print(li[0], "pars (y,qt)=", effY, effQt)
                parNum.append(effY*effQt)
                modelParsList[-1]=modelParsList[-1].flatten()
        else : #li[8]=not regularize this
            # modelParsC = rootnp.hist2array(fitRes[li[0]])
            # modelParsList.append(modelParsC.copy())
            if inPar!='' :
                modelParsC=np.zeros((dimY,dimQt))
                for y in range(0, dimY) :
                    for qt in range(0, dimQt) : 
                        modelParsC[y,qt] = inPar[0][parCounter]
                        parCounter+=1
                modelParsList.append(modelParsC.copy())        
            else :
                modelParsList.append(np.full((dimY,dimQt),1.))
                print(li[0], "pars (y,qt)=", np.shape(modelParsList[-1]))
            parNum.append( dimY*dimQt)
            modelParsList[-1]=modelParsList[-1].flatten()
            
            
    modelPars = np.concatenate(modelParsList) #already flattened, otherwise YOU HAVE TO FLATTEN before pass it to the fit
    
    # modelMass= np.zeros(1)
    # modelMass[0] = 1.
    # modelPars = np.append(modelPars,modelMass) #added the mass to the fit
    
    coeffListJAX = copy.deepcopy(coeffList)
    for li in coeffListJAX :
        if li[0] != 'unpolarizedxsec' : 
            li[0] = coeffListJAX.index(li)-1
        else :
            li[0] = 5
    
    npFitRes = npFitRes.flatten() #to add the nuisances
    # npFitRes = np.append(npFitRes,npFitMass) #mass
    
    if not nUL  : #DEBSIMP 
        #add the nuisance parameters
        # print("modelPars=", len(modelPars), "FitRes=", len(npFitRes))
        massPosition = 0
        for nbin in range(dimCoeff*dimQt*dimY+1,fitRes['cov'].GetNbinsX()+1) :
            nName = fitRes['cov'].GetXaxis().GetBinLabel(nbin) 
            if 'LHEPdf' in nName :
                nui = fitRes['nui'+'LHEPdf'].GetBinContent(fitRes['nui'+'LHEPdf'].GetXaxis().FindBin(nName))
            elif 'WHSFSyst' in nName :
                nui = fitRes['nui'+'WHSFSyst'].GetBinContent(fitRes['nui'+'WHSFSyst'].GetXaxis().FindBin(nName))
            elif 'LHEScale' in nName :
                nui = fitRes['nui'+'LHEScale'].GetBinContent(fitRes['nui'+'LHEScale'].GetXaxis().FindBin(nName)) 
            else :
                nui = fitRes['nui'+'others'].GetBinContent(fitRes['nui'+'LHEScale'].GetXaxis().FindBin(nName)) 
                if 'mass' in nName :
                    massPosition = fitRes['cov'].GetNbinsX()+1-nbin
                    # print("mass:", massPosition,nbin,dimCoeff*dimQt*dimY+1,fitRes['cov'].GetNbinsX()+1)
            # print(nName,nui)
            npFitRes = np.append(npFitRes,nui)
        
        numNui = len(npCovMat[dimCoeff*dimQt*dimY:])
        modelNuis = np.zeros(numNui)
        modelPars = np.append(modelPars,modelNuis)
        
        
        print("modelPars=", len(modelPars), "FitRes=", len(npFitRes), "Nnuis=", numNui)
        print("modelPars-nuisance", len(modelPars[:-numNui]))
    else : #placeholders
        numNui = 0 
        massPosition=0

    print("par num", parNum)
    print("everything initialized, minimizer call...")
    print("prefit pars:")
    print(modelPars[:-numNui])  
    modelPars = pmin(chi2PolyFit, modelPars, args=(npFitRes, npCovMatInv, coeffListJAX, npBinCenters,parNum,dimYQt), doParallel=False)
    # modelPars = pmin(chi2MassOnlyFit, modelPars[-1], args=(npFitRes[-1], npCovMatInv[-1][-1]), doParallel=False)
    print ("fit ended, post fit operation...")
        
    # after fit results
    static_argnumsPF = (1,2,3,4,5,6)
    modelChi2Grad_eval = jax.jit(jax.value_and_grad(chi2PolyFit), static_argnums=static_argnumsPF)
    modelHess_eval = jax.jit(jax.hessian(chi2PolyFit), static_argnums=static_argnumsPF)
    modelJac_eval = jax.jit(jax.jacfwd(par2polyModel), static_argnums=(1,2,3,4,5,6))
    
    modelChi2,modelGrad = modelChi2Grad_eval(modelPars, npFitRes, npCovMatInv, coeffListJAX, npBinCenters,parNum,dimYQt)
    modelHess = modelHess_eval(modelPars, npFitRes, npCovMatInv, coeffListJAX, npBinCenters,parNum,dimYQt)
    modelJac = modelJac_eval(modelPars, npFitRes, coeffListJAX, npBinCenters,parNum,dimYQt,numNui)
    
    modelCov = np.linalg.inv(0.5*modelHess)
    modelErr = np.sqrt(np.diag(modelCov))

    modelEDM = 0.5*np.matmul(np.matmul(modelGrad.T,modelCov),modelGrad)
    modelNDOF = dimCoeff*dimQt*dimY+1-np.size(modelPars) #+1= mass 
    
    
    
    
    print("*-------------------------------------------*")
    print("FIT RESULTS")
    print("edm=",modelEDM) 
    print("chi2/dof=", modelChi2,"/",modelNDOF, "=", modelChi2/float(modelNDOF), ", N pars=",np.size(modelPars)) 
    # print("parameters=", modelPars)
    # print("errors=", modelErr)
    # print("relative uncertaintiy" , np.true_divide(modelErr,modelPars))
    
    print("check covariance matrix:")
    print("is covariance semi positive definite?", np.all(np.linalg.eigvals(modelCov) >= 0))
    print("is covariance symmetric?", np.allclose(np.transpose(modelCov), modelCov))
    print("is covariance approx symmetric?", np.allclose(np.transpose(modelCov), modelCov, rtol=1e-05, atol=1e-06))
    # print("mass uncertainity if everything fixed=", math.sqrt(npCovMatInv[-1][-1]))
    print("*-------------------------------------------*")
    
    outFit = (modelPars,modelCov,dimYQt, npBinCenters,parNum,modelJac,npFitRes,numNui,massPosition)
    
    print("postfit pars:")
    print(modelPars[:-numNui])        
    return outFit
    
    
    
    
    
    








def par2polyModel1D(modelPars, npFitRes, coeffList, npBinCenters,parNum,dimYQt,numNui) :#input = the parameters, output the polyinomial function #not used
    fitModelList = []
    valY = npBinCenters[:dimYQt[0]]
    valQt = npBinCenters[dimYQt[0]:]
    # npBinCenters_res = npBinCenters.reshape(dimYQt[0],dimYQt[1],2)
    # valYunpol=npBinCenters_res[:,0,0]

    icoeff = 0
    ipar=0
    unpolMult=[1e5,1e5]    
    for li in coeffList :
        if not li[8] :
            # if li[0]==5 : #unpolarized
            #     fitModelList_unpol = []
            #     for iqt in range(dimYQt[1]) :    
            #         fitModelList_unpol.append(unpolMult[0]*modelPars[0+ipar]+unpolMult[1]*modelPars[1+ipar]*valYunpol+unpolMult[1]*modelPars[2+ipar]*np.power(valYunpol,2))#+unpolMult[1]*modelPars[2+ipar]*np.power(valYunpol,4))
            #         ipar = ipar+parNum[icoeff]
            #         icoeff+=1
            #     fitModelList_unpol_stack = jnp.stack(fitModelList_unpol,0)
            #     fitModelList.append(fitModelList_unpol_stack.T)
                        
            # else :
               #Y-BLOCK  
                fitModelList_y = []
                fitModelList_y.append(np.zeros(jnp.shape(valY))) #list of y
                effY = li[1]
                if li[5] : effY-= 1
                for iy in range(0, effY) :
                    powY = iy
                    if li[5] : powY+= 1
                    fitModelList_y.append(fitModelList_y[-1]+modelPars[iy+ipar]*np.power(valY,powY) )
                # fitModelList.append(fitModelList_y[-1])
                ipar = ipar+parNum[icoeff][0]
                
                #QT-BLOCK  
                if li[0]!=5 : 
                    fitModelList_qt = []
                    fitModelList_qt.append(np.zeros(jnp.shape(valQt))) #list of qt
                    effQt = li[2]
                    if li[6] : effQt-= 1
                    if li[7] : effQt-= 1
                    for iqt in range(0, effQt) :
                        powQt = iqt
                        if li[6] : powQt+= 1 
                        if li[7] and li[0]!=4 : powQt+= 1 
                        if li[7] and li[0]==4 and iqt>0 : powQt+=1
                        fitModelList_qt.append(fitModelList_qt[-1]+modelPars[iqt+ipar]*np.power(valQt,powQt) )
                    # fitModelList.append(fitModelList_qt[-1])
                    ipar = ipar+parNum[icoeff][1]
                else : #unpolarized (never regularized in qt)
                    fitModelList_qt = []
                    fitModelList_qt.append(np.zeros(jnp.shape(valQt))) #list of qt
                    fitModelList_qt.append(modelPars[ipar+0:ipar+dimYQt[1]]*npFitRes[ipar+0:ipar+dimYQt[1]])
                    ipar = ipar+parNum[icoeff][1]
                    # fitModelList.append(fitModelList_qt[-1])
                fitModelList.append( jnp.append(fitModelList_y[-1], fitModelList_qt[-1]))
                    
                icoeff+=1       
        
        else : #li[8]=not regularize this
            fitModelList_temp = modelPars[ipar+0:ipar+dimYQt[0]+dimYQt[1]]*npFitRes[ipar+0:ipar+dimYQt[0]+dimYQt[1]]
            # fitModelList_temp = modelPars[ipar+0:ipar+dimYQt[0]*dimYQt[1]]
            ipar = ipar+parNum[icoeff][0]+parNum[icoeff][1]
            icoeff+=1
            fitModelList.append(fitModelList_temp)
            # print("internal model", modelPars[ipar+0:ipar+dimYQt[0]+dimYQt[1]])
            # print("internal fitres", npFitRes[ipar+0:ipar+dimYQt[0]+dimYQt[1]])
    
    fitModel = jnp.stack(fitModelList,0)
    fitModel = fitModel.flatten()
    # fitModel = jnp.append(fitModel,modelPars[-1])#mass
    
    #add nuisance
    fitModel = jnp.append(fitModel,modelPars[-numNui:])
    
    return fitModel
     
def chi2PolyFit1D(modelPars,npFitRes, npCovMatInv, coeffList, npBinCenters,parNum,dimYQt) : #not used
    numNui = len(npCovMatInv[len(coeffList)*(dimYQt[0]+dimYQt[1]):])
    # numNui=0
    fitModel = par2polyModel1D(modelPars= modelPars, coeffList= coeffList, npBinCenters= npBinCenters,parNum= parNum,dimYQt= dimYQt,npFitRes=npFitRes,numNui=numNui)
    diff = npFitRes-fitModel
    chi2 = jnp.matmul(diff.T, jnp.matmul(npCovMatInv, diff, precision=jax.lax.Precision.HIGHEST) )
    
    return chi2    
    
def polyFit1D(fitRes, coeffList) : #not used
    
    nuisanceIntegrated=True
    
    dimY = fitRes[coeffList[0][0]].GetNbinsX()
    dimQt = fitRes[coeffList[0][0]].GetNbinsY()
    dimCoeff = len(coeffList)
    dimYQt = [dimY,dimQt]
    
    if nuisanceIntegrated  :
        npCovMatUncut = rootnp.hist2array(fitRes['cov1D'])
        npCovMat = npCovMatUncut[(dimQt+dimY):,(dimQt+dimY):]

    else :
        npCovMatUncut = rootnp.hist2array(fitRes['cov1D'])
        npCovMat = npCovMatUncut[:dimCoeff*(dimQt+dimY),:dimCoeff*(dimQt+dimY)]
        # npCovMat = npCovMatUncut[(dimQt+dimY):6*(dimQt+dimY),(dimQt+dimY):6*(dimQt+dimY)]

    npCovMatInv = np.linalg.inv(npCovMat)

    
    # npFitResList_y = []
    # npFitResList_qt = []
    npFitResList = []
    for li in coeffList :
        # npFitResList_qt.append(rootnp.hist2array(fitRes[li[0]+'y'])) 
        # npFitResList_y.append(rootnp.hist2array(fitRes[li[0]+'qt'])) 
        npFitResList.append(  np.append(rootnp.hist2array(fitRes[li[0]+'y']),rootnp.hist2array(fitRes[li[0]+'qt']) ))
    # npFitRes = np.stack(npFitResList_y+npFitResList_qt,0)  #dimensions: coeff, y+qt---> if unrolled a0y0....a0yN,a0qt1....a0qtN,a1...
    npFitRes = np.stack(npFitResList,0)  #dimensions: coeff, y+qt---> if unrolled a0y0....a0yN,a0qt1....a0qtN,a1...
    
    qtBinCenters = np.array([1., 3., 5., 7., 9., 11., 13., 15., 18., 22., 28., 38., 60.])#Extended
    # qtBinCenters = np.array([1., 3., 5., 7., 9., 11., 14., 18., 23., 30., 46.])#Extended-reduced

    print("WARNING: hardcoded qt bin centers")
    
    npBinCenters = np.zeros(dimQt+dimY)
    for iy in range(0, dimY): #y
        npBinCenters[iy] = fitRes[li[0]].GetXaxis().GetBinCenter(iy+1)
    for iqt in range(0, dimQt):#qt    
        npBinCenters[iqt] = qtBinCenters[iqt] 
            
    #model init
    modelParsList = []
    parNum = []
    for li in coeffList :
        if not li[8] : 
            effY = li[1]
            if li[5] : effY-= 1
            modelParsCy=np.zeros(effY)
            for yy in range(len(modelParsCy)) :
                modelParsCy[yy] = 1.
            modelParsList.append(modelParsCy.copy())
            # parNum.append(len(modelParsCy))
            # print(li[0], "y=",yy, "pars=", parNum[-1])
            
            if li[0]!='unpolarizedxsec' : 
                effQt = li[2]
                if li[6] : effQt-= 1
                if li[7] : effQt-= 1
                modelParsCqt=np.zeros(effQt)
                for qt in range(0, effQt) :
                    modelParsCqt[qt] = 1.
                modelParsList.append(modelParsCqt.copy())
                # parNum.append(len(modelParsCqt))
                # print(li[0], "qt=",qt, "pars=", parNum[-1])
                parNum.append((len(modelParsCy),len(modelParsCqt)))
            
            else : #unpolarized never fitted in qt
                modelParsList.append(np.full((dimQt),1.))
                # parNum.append(dimQt)
                # print(li[0], "pars (qt)=", np.shape(modelParsList[-1]))
                parNum.append((len(modelParsCy),len(dimQt)))
            
        else :
            modelParsList.append(np.full((dimY+dimQt),1.))
            # parNum.append( dimY+dimQt)
            parNum.append((dimY,dimQt))
            print(li[0], "pars (y+qt)=", np.shape(modelParsList[-1]))
                
    modelPars = np.concatenate(modelParsList) #already flattened, otherwise YOU HAVE TO FLATTEN before pass it to the fit

    
    coeffListJAX = copy.deepcopy(coeffList)
    for li in coeffListJAX :
        if li[0] != 'unpolarizedxsec' : 
            li[0] = coeffListJAX.index(li)-1
        else :
            li[0] = 5
    
    npFitRes = npFitRes.flatten() #to add the nuisances
    print("preniisance",len(npFitRes))
    
    #add the nuisance parameters
    # print("modelPars=", len(modelPars), "FitRes=", len(npFitRes))
    if nuisanceIntegrated :
        massPosition = 0
        for nbin in range((dimCoeff+1)*(dimQt+dimY)+1,fitRes['cov1D'].GetNbinsX()+1) :
            nName = fitRes['cov1D'].GetXaxis().GetBinLabel(nbin) 
            if 'LHEPdf' in nName :
                nui = fitRes['nui'+'LHEPdf'].GetBinContent(fitRes['nui'+'LHEPdf'].GetXaxis().FindBin(nName))
            elif 'WHSFSyst' in nName :
                nui = fitRes['nui'+'WHSFSyst'].GetBinContent(fitRes['nui'+'WHSFSyst'].GetXaxis().FindBin(nName))
            elif 'LHEScale' in nName :
                nui = fitRes['nui'+'LHEScale'].GetBinContent(fitRes['nui'+'LHEScale'].GetXaxis().FindBin(nName)) 
            else :
                nui = fitRes['nui'+'others'].GetBinContent(fitRes['nui'+'LHEScale'].GetXaxis().FindBin(nName)) 
                if 'mass' in nName :
                    massPosition = fitRes['cov1D'].GetNbinsX()+1-nbin
                    print("mass:", massPosition,nbin,dimCoeff*dimQt*dimY+1,fitRes['cov1D'].GetNbinsX()+1)
            # print(nName,nui)
            npFitRes = np.append(npFitRes,nui)
        
        numNui = len(npCovMat[dimCoeff*(dimQt+dimY):])
        # modelNuis = np.zeros(numNui)
        modelNuis = np.full((numNui),0.)
        modelPars = np.append(modelPars,modelNuis)
    else :
        numNui = 0
    
    
    print("modelPars=", len(modelPars), "FitRes=", len(npFitRes), "Nnuis=", numNui)
    print("modelPars-nuisance", len(modelPars[:-numNui]))
    print("MODEL PARS=", modelPars)

    print("par num", parNum)
    
    # diff = np.full((114),0.)
    # chi2 = np.matmul(diff.T, np.matmul(npCovMatInv, diff) )
    # print("chi2out=", chi2)
    
    print("everything initialized, minimizer call...")
    modelPars = pmin(chi2PolyFit1D, modelPars, args=(npFitRes, npCovMatInv, coeffListJAX, npBinCenters,parNum,dimYQt), doParallel=False)
    # modelPars = pmin(chi2MassOnlyFit, modelPars[-1], args=(npFitRes[-1], npCovMatInv[-1][-1]), doParallel=False)
    print ("fit ended, post fit operation...")
        
    # after fit results
    static_argnumsPF = (1,2,3,4,5,6)
    modelChi2Grad_eval1D = jax.jit(jax.value_and_grad(chi2PolyFit1D), static_argnums=static_argnumsPF)
    modelHess_eval1D = jax.jit(jax.hessian(chi2PolyFit1D), static_argnums=static_argnumsPF)
    modelJac_eval1D = jax.jit(jax.jacfwd(par2polyModel1D), static_argnums=(1,2,3,4,5,6))
    
    modelChi2,modelGrad = modelChi2Grad_eval1D(modelPars, npFitRes, npCovMatInv, coeffListJAX, npBinCenters,parNum,dimYQt)
    modelHess = modelHess_eval1D(modelPars, npFitRes, npCovMatInv, coeffListJAX, npBinCenters,parNum,dimYQt)
    modelJac = modelJac_eval1D(modelPars, npFitRes, coeffListJAX, npBinCenters,parNum,dimYQt,numNui)
    
    modelCov = np.linalg.inv(0.5*modelHess)
    modelErr = np.sqrt(np.diag(modelCov))

    modelEDM = 0.5*np.matmul(np.matmul(modelGrad.T,modelCov),modelGrad)
    modelNDOF = dimCoeff*dimQt*dimY+1-np.size(modelPars) #+1= mass 
    
    
    
    
    print("*-------------------------------------------*")
    print("FIT RESULTS")
    print("edm=",modelEDM) 
    print("chi2/dof=", modelChi2,"/",modelNDOF, "=", modelChi2/float(modelNDOF), ", N pars=",np.size(modelPars)) 
    # print("parameters=", modelPars)
    # print("errors=", modelErr)
    # print("relative uncertaintiy" , np.true_divide(modelErr,modelPars))
    
    print("check covariance matrix:")
    print("is covariance semi positive definite?", np.all(np.linalg.eigvals(modelCov) >= 0))
    print("is covariance symmetric?", np.allclose(np.transpose(modelCov), modelCov))
    print("is covariance approx symmetric?", np.allclose(np.transpose(modelCov), modelCov, rtol=1e-05, atol=1e-06))
    # print("mass uncertainity if everything fixed=", math.sqrt(npCovMatInv[-1][-1]))
    print("*-------------------------------------------*")
    # print("covariance=", modelCov.diagonal())
    
    if nuisanceIntegrated :
        outFit = (modelPars,modelCov,dimYQt, npBinCenters,parNum,modelJac,npFitRes,numNui,massPosition)   
    else : 
        outFit = (modelPars,modelCov,dimYQt, npBinCenters,parNum,modelJac,npFitRes,numNui,1)    
    return outFit



# def rebuildPoly(fitPars,dimYQt,npBinCenters,coeffList,parNum) : #---------------------> same thing of par2polyModel but without jax
        
#     # output = np.zeros((len(coeffList)+1,dimYQt[0],dimYQt[1])) #+1= mass block 
#     output = np.zeros((len(coeffList),dimYQt[0],dimYQt[1]))
#     valY = npBinCenters[:,0]
#     valQt = npBinCenters[:,1]
#     npBinCenters_res = npBinCenters.reshape(dimYQt[0],dimYQt[1],2)
#     valYunpol=npBinCenters_res[:,0,0]
#     # print valYunpol

#     icoeff = 0
#     ipar=0
#     unpolMult=[1e5,1e5]
#     for li in coeffList :
#         if not li[8] :
#             if li[0]=='unpolarizedxsec' :
#                 for iqt in range(dimYQt[1]) :
#                     output[0,:,iqt] = unpolMult[0]*fitPars[0+ipar]+unpolMult[1]*fitPars[1+ipar]*valYunpol+unpolMult[1]*fitPars[2+ipar]*np.power(valYunpol,2)#+unpolMult*fitPars[2+ipar]*np.power(valYunpol,4))
#                     ipar = ipar+parNum[icoeff]
#                     icoeff+=1
#             else :
#                 tempCoeff =np.zeros(np.shape(npBinCenters)[0])
#                 effY = li[1]
#                 effQt = li[2]
#                 if li[5] : effY-= 1
#                 if li[6] : effQt-= 1
#                 if li[7] : effQt-= 1
#                 for iqt in range(0, effQt) :
#                     powQt = iqt
#                     if li[6] : powQt+= 1 
#                     if li[7] and li[0]!=4 : powQt+= 1 
#                     if li[7] and li[0]==4 and iqt>0 : powQt+=1
#                     for iy in range(0, effY) :
#                         powY = iy
#                         if li[5] : powY+= 1
#                         tempCoeff = tempCoeff+fitPars[iqt*(effY)+iy+ipar]*np.power(valY,powY)*np.power(valQt,powQt) 
#                 ipar = ipar+parNum[icoeff]
#                 icoeff+=1       
#                 output[coeffList.index(li)] = tempCoeff.reshape((dimYQt[0],dimYQt[1])) # reshape to (Ny x Nqt), like the fitRes histogram
#         else : #li[8]=not regularize this
#             tempCoeff = fitPars[ipar+0:ipar+dimYQt[0]*dimYQt[1]]
#             ipar = ipar+parNum[icoeff]
#             icoeff+=1
#             output[coeffList.index(li)] = tempCoeff.reshape((dimYQt[0],dimYQt[1])) # reshape to (Ny x Nqt), like the fitRes histogram
            
#     # output[-1] = np.full((dimYQt[0],dimYQt[1]),fitPars[-1]) #mass, entire block dimY,dimQt of same value to be able to align with the rest
#     output.flatten()
#     output = np.append(output,fitPars[-1])
      
#     return output
         
        
# def errorPoly(fitPars,cov,dimYQt,npBinCenters,coeffList,parNum) : ----> with yoys
#     ntoys = 10000
#     covtoy ={}
    
#     # toyVec = np.zeros((len(coeffList)+1,dimYQt[0],dimYQt[1],ntoys))  #+1= mass block
#     toyEl = np.zeros((len(coeffList),dimYQt[0],dimYQt[1]))  
#     toyEl.flatten()
#     toyEl = np.append(toyEl, 0.)#mass
#     toyVec = np.zeros((np.shape(toyEl)[0],ntoys))

#     for itoy in range(ntoys) :
#         covtoy[itoy] = np.random.multivariate_normal(fitPars, cov,tol=1e-6)
#         # toyVec[...,itoy] = rebuildPoly(fitPars=covtoy[itoy],dimYQt=dimYQt,npBinCenters=npBinCenters, coeffList=coeffList,parNum=parNum)
#         toyVec[...,itoy] = par2polyModel(modelPars=covtoy[itoy],coeffList=coeffList,npBinCenters=npBinCenters,parNum=parNum, dimYQt=dimYQt)
        
#     sigmavec = np.std(toyVec,axis=-1)

#     return sigmavec

def poly2QTint(modelPars, npFitRes, coeffList, npBinCenters,parNum,dimYQt,numNui,binWidth) :#integral of poly2D in Y, analythic formula #not used
    # for each coefficient f_int = sum_ij * c_ij * q^i (maxY)^{j+1}/(j+1)
    fitModelList = []
    # valY = npBinCenters[:,0]
    # valQt = npBinCenters[:,1]
    npBinCenters_res = npBinCenters.reshape(dimYQt[0],dimYQt[1],2)
    # valYunpol=npBinCenters_res[:,0,0]
    Ycent = npBinCenters_res[:,0,0]
    QTcent = npBinCenters_res[0,:,1]
    
    # print(valQt[:dimYQt[1]])
    # print(npBinCenters_res[:,0,0])
    # print(npBinCenters_res[0,:,1])

    icoeff = 0
    ipar=0
    unpolMult=[1e5,1e5]    
    for li in coeffList :
        if not li[8] :
            if li[0]==5 :
                print("integral on regularized unpol not implemented")
                FIXME
                # fitModelList_unpol = []
                # for iqt in range(dimYQt[1]) :    
                #     fitModelList_unpol.append(unpolMult[0]*modelPars[0+ipar]+unpolMult[1]*modelPars[1+ipar]*valYunpol+unpolMult[1]*modelPars[2+ipar]*np.power(valYunpol,2))#+unpolMult[1]*modelPars[2+ipar]*np.power(valYunpol,4))
                #     ipar = ipar+parNum[icoeff]
                #     icoeff+=1
                # fitModelList_unpol_stack = jnp.stack(fitModelList_unpol,0)
                # fitModelList.append(fitModelList_unpol_stack.T)
                        
            else :
                fitModelList_temp = []
                fitModelList_temp.append(np.zeros(jnp.shape(dimYQt[1]))) # list of Qt
                effY = li[1]
                effQt = li[2]
                if li[5] : effY-= 1
                if li[6] : effQt-= 1
                if li[7] : effQt-= 1
                for iqt in range(0, effQt) :
                    powQt = iqt
                    if li[6] : powQt+= 1 
                    if li[7] and li[0]!=4 : powQt+= 1 
                    if li[7] and li[0]==4 and iqt>0 : powQt+=1
                    for iy in range(0, effY) :
                        powY = iy
                        if li[5] : powY+= 1
                        # fitModelList_temp.append(fitModelList_temp[-1]+modelPars[iqt*(effY)+iy+ipar]*np.power(2.4,powY+1)*np.power(QTcent,powQt)/(powY+1) )
                        fitModelList_temp.append(fitModelList_temp[-1]+modelPars[iy*(effQt)+iqt+ipar]*np.power(2.4,powY+1)*np.power(QTcent,powQt)/(powY+1) )
                        if(li[0]==2) : print("iy=",iy, " iqt=",iqt, "powQT=", powQt, " powY=",powY, " ipar=",ipar, " par=", iy*(effQt)+iqt+ipar, "effY=", effY, " effQt=",effQt, modelPars[iy*(effQt)+iqt+ipar])#  iy*(effQt)+iqt+ipar,modelPars[iy*(effQt)+iqt+ipar])
                ipar = ipar+parNum[icoeff]
                icoeff+=1       
                fitModelList.append(fitModelList_temp[-1])#.reshape(dimYQt[0],dimYQt[1])) # reshape to (Ny x Nqt), like the fitRes histogram
        
        else : #li[8]=not regularize this
            # fitModelList_temp = []
            # fitModelList_temp.append(np.zeros(jnp.shape(valQt))) # list of Qt
            listOfQt = []
            for iqt in range(0, dimYQt[1]) :
                signleQt = []
                signleQt.append(np.zeros(1)) # list of Qt
                for iy in range(0, dimYQt[0]) : 
                    # signleQt.append(signleQt[-1] + modelPars[iqt*(dimYQt[1])+iy+ipar] * npFitRes[iqt*(dimYQt[1])+iy+ipar]*binWidth[iy])
                    signleQt.append(signleQt[-1] + modelPars[iy*(dimYQt[1])+iqt+ipar] * npFitRes[iy*(dimYQt[1])+iqt+ipar]*binWidth[iy])
                    # print("QTint", iy,iqt,modelPars[iy*(dimYQt[0])+iqt+ipar],npFitRes[iy*(dimYQt[0])+iqt+ipar],binWidth[iy], ipar)

                listOfQt.append(signleQt[-1])
                # print("final qt list,",listOfQt )
            fitModelList.append(jnp.array(listOfQt).reshape(dimYQt[1]))
        
            # fitModelList_temp = modelPars[ipar+0:ipar+dimYQt[0]*dimYQt[1]] * npFitRes[ipar+0:ipar+dimYQt[0]*dimYQt[1]]
            # fitModelList_temp = modelPars[ipar+0:ipar+dimYQt[0]*dimYQt[1]]
            ipar = ipar+parNum[icoeff]
            icoeff+=1
            # fitModelList.append(fitModelList_temp.reshape(dimYQt[0],dimYQt[1]))
    
    fitModel = jnp.stack(fitModelList,0)
    fitModel = fitModel.flatten()
    # fitModel = jnp.append(fitModel,modelPars[-1])#mass
    # print("internal model",len(fitModel))
    #add nuisance
    fitModel = jnp.append(fitModel,modelPars[-numNui:])
    # print("internal model after nuisance add",len(fitModel))
    print("OUT")
    return fitModel

def poly2Yint(modelPars, npFitRes, coeffList, npBinCenters,parNum,dimYQt,numNui,binWidth) :#integral of poly2D in qt, analythic formula #not used
    # for each coefficient f_int = sum_ij * c_ij * y^i (maxQt)^{j+1}/(j+1)
    fitModelList = []
    # valY = npBinCenters[:,0]
    # valQt = npBinCenters[:,1]
    npBinCenters_res = npBinCenters.reshape(dimYQt[0],dimYQt[1],2)
    # valYunpol=npBinCenters_res[:,0,0]
    Ycent = npBinCenters_res[:,0,0]
    QTcent = npBinCenters_res[0,:,1]

    icoeff = 0
    ipar=0
    unpolMult=[1e5,1e5]    
    for li in coeffList :
        if not li[8] :
            if li[0]==5 :
                print("integral on regularized unpol not implemented")
                FIXME
                # fitModelList_unpol = []
                # for iqt in range(dimYQt[1]) :    
                #     fitModelList_unpol.append(unpolMult[0]*modelPars[0+ipar]+unpolMult[1]*modelPars[1+ipar]*valYunpol+unpolMult[1]*modelPars[2+ipar]*np.power(valYunpol,2))#+unpolMult[1]*modelPars[2+ipar]*np.power(valYunpol,4))
                #     ipar = ipar+parNum[icoeff]
                #     icoeff+=1
                # fitModelList_unpol_stack = jnp.stack(fitModelList_unpol,0)
                # fitModelList.append(fitModelList_unpol_stack.T)
                        
            else :
                fitModelList_temp = []
                fitModelList_temp.append(np.zeros(jnp.shape(dimYQt[0]))) # list of Qt
                effY = li[1]
                effQt = li[2]
                if li[5] : effY-= 1
                if li[6] : effQt-= 1
                if li[7] : effQt-= 1
                for iqt in range(0, effQt) :
                    powQt = iqt
                    if li[6] : powQt+= 1 
                    if li[7] and li[0]!=4 : powQt+= 1 
                    if li[7] and li[0]==4 and iqt>0 : powQt+=1
                    for iy in range(0, effY) :
                        powY = iy
                        if li[5] : powY+= 1
                        # fitModelList_temp.append(fitModelList_temp[-1]+modelPars[iqt*(effY)+iy+ipar]*np.power(Ycent,powY+1)*np.power(80.,powQt+1)/(powQt+1) )
                        fitModelList_temp.append(fitModelList_temp[-1]+modelPars[iy*(effQt)+iqt+ipar]*np.power(Ycent,powY)*np.power(80.,powQt+1)/(powQt+1) )
                ipar = ipar+parNum[icoeff]
                icoeff+=1       
                fitModelList.append(fitModelList_temp[-1])#.reshape(dimYQt[0],dimYQt[1])) # reshape to (Ny x Nqt), like the fitRes histogram
        
        else : #li[8]=not regularize this
            # fitModelList_temp = []
            # fitModelList_temp.append(np.zeros(jnp.shape(valQt))) # list of Qt
            listOfY = []
            for iy in range(0, dimYQt[0]) : 
                singleY = []
                singleY.append(np.zeros(1)) # list of Qt
                for iqt in range(0, dimYQt[1]) :
                    singleY.append(singleY[-1] + modelPars[iy*(dimYQt[1])+iqt+ipar] * npFitRes[iy*(dimYQt[1])+iqt+ipar]*binWidth[iqt])
                    # print("YINT", iy,iqt,modelPars[iy*(dimYQt[1])+iqt+ipar],npFitRes[iy*(dimYQt[1])+iqt+ipar],binWidth[iqt])
                    # singleY.append(singleY[-1] + modelPars[iqt*(dimYQt[1])+iy+ipar] * npFitRes[iqt*(dimYQt[1])+iy+ipar]*binWidth[iqt])
                    # print(iy,iqt,modelPars[iqt*(dimYQt[1])+iy+ipar],npFitRes[iqt*(dimYQt[1])+iy+ipar],binWidth[iqt])
                listOfY.append(singleY[-1])
                # print("FINALIST Y",iy, listOfY)
            fitModelList.append(jnp.array(listOfY).reshape(dimYQt[0]))
        
            # fitModelList_temp = modelPars[ipar+0:ipar+dimYQt[0]*dimYQt[1]] * npFitRes[ipar+0:ipar+dimYQt[0]*dimYQt[1]]
            # fitModelList_temp = modelPars[ipar+0:ipar+dimYQt[0]*dimYQt[1]]
            ipar = ipar+parNum[icoeff]
            icoeff+=1
            # fitModelList.append(fitModelList_temp.reshape(dimYQt[0],dimYQt[1]))
            
    fitModel = jnp.stack(fitModelList,0)
    fitModel = fitModel.flatten()
    # fitModel = jnp.append(fitModel,modelPars[-1])#mass
    # print("internal model",len(fitModel))
    #add nuisance
    fitModel = jnp.append(fitModel,modelPars[-numNui:])
    # print("internal model after nuisance add",len(fitModel))
    
    return fitModel   
    
    
    
    
    
def errorPoly(cov, jacobian) :
    covariancePoly = np.matmul(jacobian,np.matmul(cov,jacobian.T))
    variances = covariancePoly.diagonal()
    stddev = np.sqrt(variances)
    return stddev, covariancePoly
    
def integrateModel(modelPars, coeffList,histo,kind) : #integrated (without jax) #not used
        
    if kind=='y' : 
        out = np.zeros(( len(coeffList), histo['A4'].GetNbinsX()  ))        
        UL = np.zeros(( len(coeffList), histo['A4'].GetNbinsX()  ))        
        helxsec = np.zeros(( len(coeffList), histo['A4'].GetNbinsX()  ))        
        for li in coeffList :
            #qt integration
            # yDist = np.zeros(histo[li[0]].GetNbinsX())
            for yy in range(0,histo[li[0]].GetNbinsX()) :
                yWid =  histo[li[0]].GetXaxis().GetBinWidth(yy+1)
                for qq in range(0, histo[li[0]].GetNbinsY()) :
                    qtVal = modelPars[coeffList.index(li)][yy][qq] 
                    # if 'unpol' in li[0]:
                    qtWid =  histo[li[0]].GetYaxis().GetBinWidth(qq+1)
                    
                    #     out[coeffList.index(li)][yy] += qtVal*qtWid
                    # else :
                    ULval = modelPars[coeffList.index(coeffList[0])][yy][qq] #unpolarized xsec 
                    helxsec[coeffList.index(li)][yy] += ULval*qtVal*(yWid*qtWid)
                    UL[coeffList.index(li)][yy] += ULval*(yWid*qtWid)
                if 'unpol' in li[0]:
                     out[coeffList.index(li)][yy] = UL[coeffList.index(li)][yy]/yWid
                else : 
                    out[coeffList.index(li)][yy] = helxsec[coeffList.index(li)][yy]/UL[coeffList.index(li)][yy]
    
    if kind =='qt' :
        out = np.zeros(( len(coeffList), histo['A4'].GetNbinsY()  ))
        UL = np.zeros(( len(coeffList), histo['A4'].GetNbinsY()  ))        
        helxsec = np.zeros(( len(coeffList), histo['A4'].GetNbinsY()  ))       
        for li in coeffList : 
            #y integration
            for qq in range(0,histo[li[0]].GetNbinsY()) :
                qtWid =  histo[li[0]].GetYaxis().GetBinWidth(qq+1)
                for yy in range(0, histo[li[0]].GetNbinsX()) :
                    yVal = modelPars[coeffList.index(li)][yy][qq] 
                    yWid =  histo[li[0]].GetXaxis().GetBinWidth(yy+1)
                    ULval = modelPars[coeffList.index(coeffList[0])][yy][qq] #unpolarized xsec 
                    helxsec[coeffList.index(li)][qq] += ULval*yVal*(yWid*qtWid)
                    UL[coeffList.index(li)][qq] += ULval*(yWid*qtWid)
                if 'unpol' in li[0]:
                    out[coeffList.index(li)][qq] = UL[coeffList.index(li)][qq]/qtWid
                else :
                    out[coeffList.index(li)][qq] = helxsec[coeffList.index(li)][qq]/UL[coeffList.index(li)][qq]
            
    # if kind =='qt' :
    #     out = np.zeros(( len(coeffList), histo['A4'].GetNbinsY()  ))
    #     for li in coeffList : 
    #         #y integration
    #         # qtDist = np.zeros(histo[li[0]].GetNbinsY())
    #         for qq in range(0,histo[li[0]].GetNbinsY()) :
    #             for yy in range(0, histo[li[0]].GetNbinsX()) :
    #                 yVal = modelPars[coeffList.index(li)][yy][qq] 
    #                 if 'unpol' in li[0]:
    #                     yWid =  histo[li[0]].GetXaxis().GetBinWidth(yy+1)
    #                 # qtDist[qq]+= yVal*yWid
    #                 # print(coeffList.index(li), qq, yy, np.shape(out))
    #                     out[coeffList.index(li)][qq] += yVal*yWid
    #                 else :#only because same bin width 
    #                     out[coeffList.index(li)][qq] += yVal
    #         # outDict[li[0]+'qt'] = qtDist.copy()
    
    return out
    
def errorIntToy(fitPars,cov,dimYQt,npBinCenters,coeffList,parNum,histo,npFitRes, numNui) : #not used
    ntoys = 10000
    covtoy ={}
    
    # npCovMat = rootnp.hist2array(histo['cov'])
    # numNui = len(npCovMat[dimCoeff*dimQt*dimY:])
    
    # toyVec = np.zeros((len(coeffList)+1,dimYQt[0],dimYQt[1],ntoys))  #+1= mass block
    toyEl = np.zeros((len(coeffList),dimYQt[0],dimYQt[1]))  
    toyEl.flatten()
    modelNuis = np.zeros(numNui)
    toyEl = np.append(toyEl,modelNuis)
    # toyEl = np.append(toyEl, 0.)#mass
    toyVec = np.zeros((np.shape(toyEl)[0],ntoys))
    
    # rng = np.random.default_rng()
    for itoy in range(ntoys) :
        print("toy=", itoy)
        covtoy[itoy] = np.random.multivariate_normal(fitPars, cov,tol=1e-8)
        # covtoy[itoy] = rng.random.multivariate_normal(fitPars, cov)
        # toyVec[...,itoy] = rebuildPoly(fitPars=covtoy[itoy],dimYQt=dimYQt,npBinCenters=npBinCenters, coeffList=coeffList,parNum=parNum)
        toyVec[...,itoy] = par2polyModel(modelPars=covtoy[itoy],coeffList=coeffList,npBinCenters=npBinCenters,parNum=parNum, dimYQt=dimYQt, npFitRes=npFitRes, numNui=numNui)
    
    toyVec = toyVec[:-numNui,:] #exclude nuisance
    toyVec= np.reshape(toyVec,(len(coeffList),dimYQt[0],dimYQt[1],ntoys))
    
    toyIntQt = out = np.zeros(( len(coeffList), histo['A4'].GetNbinsY(), ntoys))   
    toyIntY = out = np.zeros(( len(coeffList), histo['A4'].GetNbinsX(), ntoys))   
    
    for itoy in range(ntoys) :
        toyIntQt[...,itoy] = integrateModel(modelPars = toyVec[...,itoy], coeffList=coeffList,histo=histo,kind='qt')
        toyIntY[...,itoy] = integrateModel(modelPars = toyVec[...,itoy], coeffList=coeffList,histo=histo,kind='y')

    # sigmavec = np.std(toyVec,axis=-1)
    sigmaQt = np.std(toyIntQt,axis=-1)
    sigmaY = np.std(toyIntY,axis=-1)
    
    return sigmaY,sigmaQt
    
    
def poly2YintH(modelPars, coeffList,dimYQt,numNui,binWidthY,binWidthQT,nUL, extUL ) : #input (modelPars) = the double diff function, output the polyinomial function integrated
    #qt integration
    out_list = [] 
 
    ipar=0
    for li in coeffList :
      
        for yy in range(0,dimYQt[0]) :
            helxsec_list_y = []
            UL_list_y = []
            helxsec_list_y.append(np.zeros(1))
            UL_list_y.append(np.zeros(1))
            
            for qq in range(0, dimYQt[1]) :
                            
                qtVal = modelPars[yy*(dimYQt[1])+qq+ipar]
                if not nUL :
                    ULval = modelPars[yy*(dimYQt[1])+qq] #unpolarized xsec, as ipar=0 (first coeff)
                else :
                    ULval = extUL[yy][qq]
                
                # helxsec_list_y.append(helxsec_list_y[-1]+ULval*qtVal*(binWidthY[yy]*binWidthQT[qq]))
                # UL_list_y.append(UL_list_y[-1]+ULval*(binWidthY[yy]*binWidthQT[qq]))
                helxsec_list_y.append(helxsec_list_y[-1]+ULval*qtVal)
                UL_list_y.append(UL_list_y[-1]+ULval)
                
                # print(li[0], "qt,y=",qq,yy, "valqt=",qtVal, "binw qt,y=",binWidthQT[qq], binWidthY[yy], UL_list_y[-1] )
            
            if li[0]==5:
                # out_list.append(UL_list_y[-1]/binWidthY[yy])
                out_list.append(UL_list_y[-1])
            else : 
                out_list.append(helxsec_list_y[-1]/UL_list_y[-1])
        
        ipar = ipar+dimYQt[1]*dimYQt[0]
    
    model = jnp.array(out_list)
    if not nUL : #DEBSIMP
        model = jnp.append(model,modelPars[-numNui:])
    else :
        model = model.flatten()
    return model 

    
def poly2QTintH(modelPars, coeffList,dimYQt,numNui,binWidthY,binWidthQT,nUL,extUL ) :#input (modelPars) = the double diff function, output the polyinomial function integrated
    #Y integration
    out_list = [] 
    ipar=0
    for li in coeffList :
        
        for qq in range(0,dimYQt[1]) :
            helxsec_list_qt = []
            UL_list_qt = []
            helxsec_list_qt.append(np.zeros(1))
            UL_list_qt.append(np.zeros(1))
            
            for yy in range(0, dimYQt[0]) :
                yVal = modelPars[yy*(dimYQt[1])+qq+ipar]
                if not nUL : 
                    ULval =modelPars[yy*(dimYQt[1])+qq] #unpolarized xsec, as ipar=0 (first coeff)
                else :
                    ULval = extUL[yy][qq]
                
                # helxsec_list_qt.append(helxsec_list_qt[-1]+ULval*yVal*(binWidthY[yy]*binWidthQT[qq]))
                # UL_list_qt.append(UL_list_qt[-1]+ULval*(binWidthY[yy]*binWidthQT[qq]))
                helxsec_list_qt.append(helxsec_list_qt[-1]+ULval*yVal)
                UL_list_qt.append(UL_list_qt[-1]+ULval)
            
            if li[0]==5:
                # out_list.append(UL_list_qt[-1]/binWidthQT[qq])
                out_list.append(UL_list_qt[-1])
            else : 
                out_list.append(helxsec_list_qt[-1]/UL_list_qt[-1])
        
        ipar = ipar+dimYQt[1]*dimYQt[0]
    
    model = jnp.array(out_list)
    if not nUL : #DEBSIMP
        model = jnp.append(model,modelPars[-numNui:])
    else :
        model = model.flatten()
    return model 

    
    
      
# def plotterPostReg(fitResult, fitResult1D, output, save, coeffList,coeffList1D, histo) :
def plotterPostReg(fitResult, output, save, coeffList, histo,nUL) :
    print("plotting...")
    
    coeffListJAX = copy.deepcopy(coeffList)
    for li in coeffListJAX :
        if li[0] != 'unpolarizedxsec' : 
            li[0] = coeffListJAX.index(li)-1
        else :
            li[0] = 5
    
    ########################################## double differential #############################################
    #build the error vector
    # funcVec = rebuildPoly(fitPars=fitResult[0],dimYQt=fitResult[2],npBinCenters=fitResult[3], coeffList=coeffList,parNum=fitResult[4])
    funcVec = par2polyModel(modelPars=fitResult[0],coeffList=coeffListJAX,npBinCenters=fitResult[3],parNum=fitResult[4], dimYQt=fitResult[2],npFitRes=fitResult[6],numNui=fitResult[7]) 
    funcVecFull = funcVec.copy() #to integrated
    #build the function vector
    errorVec, covVec = errorPoly(cov=fitResult[1],jacobian =fitResult[5])
    
    #exclude nuisances
    if not nUL : #or 0
        nuiVec = funcVec[-fitResult[7]:].copy()
        nuiErr = errorVec[-fitResult[7]:].copy()
    
        massVal = nuiVec[-fitResult[8]]
        massErr = nuiErr[-fitResult[8]]
        
        funcVec = funcVec[:-fitResult[7]] #DEBLORE
        errorVec = errorVec[:-fitResult[7]] #DEBLORE
    
    # massVal = funcVec[-1].copy()
    # massErr = errorVec[-1].copy()
    # funcVec = np.reshape(funcVec[:-1],(len(coeffList),fitResult[2][0],fitResult[2][1]))
    # errorVec = np.reshape(errorVec[:-1],(len(coeffList),fitResult[2][0],fitResult[2][1]))
    funcVec = np.reshape(funcVec,(len(coeffList),fitResult[2][0],fitResult[2][1]))
    errorVec = np.reshape(errorVec,(len(coeffList),fitResult[2][0],fitResult[2][1]))
    
    
    outFile = ROOT.TFile(output+".root", "recreate")
    for li in coeffList :
        h = histo[li[0]].Clone("post-fit-regularization_"+li[0])
        for y in range(1, h.GetNbinsX()+1) :
            for qt in range(1, h.GetNbinsY()+1) : 
                # print("difference", li[0], y,qt, ", original=", h.GetBinContent(y,qt), ", apofitted=", funcVec[coeffList.index(li)][y-1][qt-1], ", ratio=", h.GetBinContent(y,qt)/funcVec[coeffList.index(li)][y-1][qt-1])
                h.SetBinContent(y,qt,funcVec[coeffList.index(li)][y-1][qt-1])  
                h.SetBinError(y,qt,errorVec[coeffList.index(li)][y-1][qt-1])  
        
        if 'unpol' in li[0] :
            normm = (3./16./math.pi)*35.9
            h.Scale(1/normm)
            for i in range(1, h.GetNbinsX()+1):
                for j in range(1,  h.GetNbinsY()+1):
                    h.SetBinContent(i,j, h.GetBinContent(i,j)/h.GetXaxis().GetBinWidth(i)/h.GetYaxis().GetBinWidth(j))
                    h.SetBinError(i,j, h.GetBinError(i,j)/h.GetXaxis().GetBinWidth(i)/h.GetYaxis().GetBinWidth(j)) 

        
        outFile.cd()
        h.Write()
    if nUL : #placeholder
        h = histo['unpolarizedxsec'].Clone("post-fit-regularization_unpolarizedxsec")
        h.Write()
        
    if not nUL : #DEBSIMP 
        print("WARNING: values muliplied x 50 to obtain MeV! good only if upd/down 50 MeV") # 1_nuisanceLimit : 50MeV_Given = massVal : x 
        print("predicted W mass = ", 50*massVal, "+/-", 50*massErr) #, ", only diag=",math.sqrt(fitResult[1][-1][-1])) 
    
    #save the correlations and covariance matrices  
    dim2DFit = np.shape(fitResult[1])[0] #Npar_poly2D+nuisance = 108+242
    h2DFitCorrMatrix = ROOT.TH2F("h2DFitCorrMatrix","h2DFitCorrMatrix", dim2DFit,0,dim2DFit+1,dim2DFit,0,dim2DFit+1)
    h2DFitCovMatrix = ROOT.TH2F("h2DFitCovMatrix","h2DFitCovMatrix", dim2DFit,0,dim2DFit+1,dim2DFit,0,dim2DFit+1)
    for xx in range(1, dim2DFit+1) :
        for yy in range(1, dim2DFit+1) :
            val = fitResult[1][xx-1][yy-1]/math.sqrt(fitResult[1][xx-1][xx-1]*fitResult[1][yy-1][yy-1])
            h2DFitCorrMatrix.SetBinContent(xx,yy,val)
            h2DFitCovMatrix.SetBinContent(xx,yy,fitResult[1][xx-1][yy-1])
    h2DFitCorrMatrix.Write()    
    h2DFitCovMatrix.Write()    
    
    dim2DPoly = np.shape(covVec)[0] #Nai*Nqt*Ny+nuisance =6*11*6+nuisance = 396+242 
    h2DPolyCorrMatrix = ROOT.TH2F("h2DPolyCorrMatrix","h2DPolyCorrMatrix", dim2DPoly,0,dim2DPoly+1,dim2DPoly,0,dim2DPoly+1)
    h2DPolyCovMatrix = ROOT.TH2F("h2DPolyCovMatrix","h2DPolyCovMatrix", dim2DPoly,0,dim2DPoly+1,dim2DPoly,0,dim2DPoly+1)
    for xx in range(1, dim2DPoly+1) :
        for yy in range(1, dim2DPoly+1) :
            val = covVec[xx-1][yy-1]/math.sqrt(covVec[xx-1][xx-1]*covVec[yy-1][yy-1])
            h2DPolyCorrMatrix.SetBinContent(xx,yy,val)
            h2DPolyCovMatrix.SetBinContent(xx,yy,covVec[xx-1][yy-1])
    h2DPolyCorrMatrix.Write()  
    h2DPolyCovMatrix.Write()  
    
    #save model parameters (fitted)
    dimPars = np.shape(fitResult[0])[0]
    hPars = ROOT.TH1F("hPars","hPars",dimPars,0,dimPars+1)
    for xx in range(1,dimPars+1) :
        hPars.SetBinContent(xx,fitResult[0][xx-1])
        hPars.SetBinError(xx, math.sqrt(fitResult[1][xx-1][xx-1]))
    hPars.Write()
    

    
    ############################################ Integrated (propagation with histogram with jax functionality) #############################################
    
    binWidthY = np.zeros(fitResult[2][0])
    for yy in range(1,histo[coeffList[0][0]+'y'].GetNbinsX()+1) :
        binWidthY[yy-1] = histo[coeffList[0][0]+'y'].GetXaxis().GetBinWidth(yy)
    
    binWidthQT = np.zeros(fitResult[2][1])
    for qq in range(1,histo[coeffList[0][0]+'qt'].GetNbinsX()+1) :
        binWidthQT[qq-1] = histo[coeffList[0][0]+'qt'].GetXaxis().GetBinWidth(qq)
    
    extUL = rootnp.hist2array(histo['unpolarizedxsec'])
        
    intFuncVecQt = poly2QTintH(modelPars=funcVecFull,coeffList=coeffListJAX, dimYQt=fitResult[2],numNui=fitResult[7], binWidthY = binWidthY, binWidthQT = binWidthQT, nUL = nUL, extUL=extUL)
    intFuncVecY = poly2YintH(modelPars=funcVecFull,coeffList=coeffListJAX, dimYQt=fitResult[2],numNui=fitResult[7], binWidthY = binWidthY, binWidthQT = binWidthQT, nUL = nUL, extUL=extUL)
        
    intQtJac_eval = jax.jit(jax.jacfwd(poly2QTintH), static_argnums=(1,2,3,4,5,6,7))
    intYJac_eval = jax.jit(jax.jacfwd(poly2YintH), static_argnums=(1,2,3,4,5,6,7))
    
    if nUL :
        funcVecFull = np.reshape(funcVecFull, np.shape(funcVecFull)[0])
    intQtJac = intQtJac_eval(funcVecFull,  coeffListJAX, fitResult[2],  fitResult[7], binWidthY, binWidthQT, nUL,extUL)
    intYJac = intYJac_eval(funcVecFull,  coeffListJAX,  fitResult[2], fitResult[7], binWidthY, binWidthQT, nUL,extUL)
    
    
    intFuncErrQt, covQt = errorPoly(cov=covVec,jacobian =intQtJac)
    intFuncErrY, covY = errorPoly(cov=covVec,jacobian =intYJac)

    if not nUL : #DEBSIMP
    #esclude nuisance
        intFuncVecQt = intFuncVecQt[:-fitResult[7]]
        intFuncVecY = intFuncVecY[:-fitResult[7]]
        
        intFuncErrQt = intFuncErrQt[:-fitResult[7]]
        intFuncErrY = intFuncErrY[:-fitResult[7]]
        
    intFuncVecQt = np.reshape(intFuncVecQt,(len(coeffList),fitResult[2][1]))
    intFuncVecY = np.reshape(intFuncVecY,(len(coeffList),fitResult[2][0]))
        
    intFuncErrQt = np.reshape(intFuncErrQt,(len(coeffList),fitResult[2][1]))
    intFuncErrY = np.reshape(intFuncErrY,(len(coeffList),fitResult[2][0]))
    
    for li in coeffList :
        hy = histo[li[0]+'y'].Clone("post-fit-regularization_"+li[0]+'_y')
        for y in range(1, hy.GetNbinsX()+1) :
            print("difference y", li[0], y, ", original=", hy.GetBinContent(y), ", apofitted=", intFuncVecY[coeffList.index(li)][y-1], ", ratio=", hy.GetBinContent(y)/intFuncVecY[coeffList.index(li)][y-1])
            hy.SetBinContent(y,intFuncVecY[coeffList.index(li)][y-1])  
            hy.SetBinError(y,intFuncErrY[coeffList.index(li)][y-1]) 
        
        if 'unpol' in li[0] :
            normm = (3./16./math.pi)*35.9
            hy.Scale(1/normm)
            for i in range(1, hy.GetNbinsX()+1):
                    hy.SetBinContent(i, hy.GetBinContent(i)/hy.GetXaxis().GetBinWidth(i))
                    hy.SetBinError(i, hy.GetBinError(i)/hy.GetXaxis().GetBinWidth(i)) 
                    
        hy.Write()
        
        
        
        hqt = histo[li[0]+'qt'].Clone("post-fit-regularization_"+li[0]+'_qt')     
        for qt in range(1, hqt.GetNbinsX()+1) : 
            print("difference qt", li[0], qt, ", original=", hqt.GetBinContent(qt), ", apofitted=", intFuncVecQt[coeffList.index(li)][qt-1], ", ratio=", hqt.GetBinContent(qt)/intFuncVecQt[coeffList.index(li)][qt-1])
            hqt.SetBinContent(qt,intFuncVecQt[coeffList.index(li)][qt-1])  
            hqt.SetBinError(qt,intFuncErrQt[coeffList.index(li)][qt-1])  
        
        if 'unpol' in li[0] :
            normm = (3./16./math.pi)*35.9
            hqt.Scale(1/normm)
            for i in range(1, hqt.GetNbinsX()+1):
                hqt.SetBinContent(i, hqt.GetBinContent(i)/hqt.GetXaxis().GetBinWidth(i))
                hqt.SetBinError(i, hqt.GetBinError(i)/hqt.GetXaxis().GetBinWidth(i)) 
                print("unpolarized relative error",i, hqt.GetBinError(i)/ hqt.GetBinContent(i))
            
        hqt.Write()
        
    if nUL : #placeholder
        h = histo['unpolarizedxsecqt'].Clone("post-fit-regularization_unpolarizedxsec_qt")
        h.Write()
        h = histo['unpolarizedxsecy'].Clone("post-fit-regularization_unpolarizedxsec_y")
        h.Write()
    
    
    


    
    
    ############################################ Integrated (propagation with explicit function) #############################################
    
    # binWidthY = np.zeros(fitResult[2][1])
    # for yy in range(1,histo['A4'+'y'].GetNbinsX()+1) :
    #     binWidthY[yy-1] = histo['A4'+'y'].GetXaxis().GetBinWidth(yy)
    
    # binWidthQT = np.zeros(fitResult[2][1])
    # for qq in range(1,histo['A4'+'qt'].GetNbinsX()+1) :
    #     binWidthQT[qq-1] = histo['A4'+'qt'].GetXaxis().GetBinWidth(qq)
        
    # intFuncVecQt = poly2QTint(modelPars=fitResult[0],coeffList=coeffListJAX,npBinCenters=fitResult[3],parNum=fitResult[4], dimYQt=fitResult[2],npFitRes=fitResult[6],numNui=fitResult[7], binWidth = binWidthY)
    # intFuncVecY = poly2Yint(modelPars=fitResult[0],coeffList=coeffListJAX,npBinCenters=fitResult[3],parNum=fitResult[4], dimYQt=fitResult[2],npFitRes=fitResult[6],numNui=fitResult[7], binWidth = binWidthQT)
        
    # intQtJac_eval = jax.jit(jax.jacfwd(poly2QTint), static_argnums=(1,2,3,4,5,6,7))
    # intYJac_eval = jax.jit(jax.jacfwd(poly2Yint), static_argnums=(1,2,3,4,5,6,7))
    
    # intQtJac = intQtJac_eval(fitResult[0],  fitResult[6], coeffListJAX, fitResult[3], fitResult[4], fitResult[2], fitResult[7], binWidthY)
    # intYJac = intYJac_eval(fitResult[0],  fitResult[6], coeffListJAX, fitResult[3], fitResult[4], fitResult[2], fitResult[7], binWidthQT)
    
    # intFuncErrQt = errorPoly(cov=fitResult[1],jacobian =intQtJac)
    # intFuncErrY = errorPoly(cov=fitResult[1],jacobian =intYJac)


    # #esclude nuisance
    # intFuncVecQt = intFuncVecQt[:-fitResult[7]]
    # intFuncVecY = intFuncVecY[:-fitResult[7]]
    
    # intFuncErrQt = intFuncErrQt[:-fitResult[7]]
    # intFuncErrY = intFuncErrY[:-fitResult[7]]
    
    # intFuncVecQt = np.reshape(intFuncVecQt,(len(coeffList),fitResult[2][1]))
    # intFuncVecY = np.reshape(intFuncVecY,(len(coeffList),fitResult[2][0]))
    
    # intFuncErrQt = np.reshape(intFuncErrQt,(len(coeffList),fitResult[2][1]))
    # intFuncErrY = np.reshape(intFuncErrY,(len(coeffList),fitResult[2][0]))
    
    # for li in coeffList :
    #     hy = histo[li[0]+'y'].Clone("post-fit-regularization_"+li[0]+'_y')
    #     for y in range(1, hy.GetNbinsX()+1) :
    #         # print("difference y", li[0], y,qt, ", original=", hy.GetBinContent(y,qt), ", apofitted=", funcVec[coeffList.index(li)][y-1][qt-1], ", ratio=", hy.GetBinContent(y,qt)/funcVec[coeffList.index(li)][y-1])
    #         hy.SetBinContent(y,intFuncVecY[coeffList.index(li)][y-1])  
    #         hy.SetBinError(y,intFuncErrY[coeffList.index(li)][y-1])  
    #     hy.Write()
        
    #     hqt = histo[li[0]+'qt'].Clone("post-fit-regularization_"+li[0]+'_qt')     
    #     for qt in range(1, hqt.GetNbinsX()+1) : 
    #         # print("difference qt", li[0], y,qt, ", original=", hqt.GetBinContent(y,qt), ", apofitted=", funcVec[coeffList.index(li)][y-1][qt-1], ", ratio=", hqt.GetBinContent(y,qt)/funcVec[coeffList.index(li)][y-1][hy.GetNbinsX()+qt-1])
    #         hqt.SetBinContent(qt,intFuncVecQt[coeffList.index(li)][qt-1])  
    #         hqt.SetBinError(qt,intFuncErrQt[coeffList.index(li)][qt-1])  
    #     hqt.Write()
        
        
        
    
    
    # ############################################ Integrated (propagation with the histogram, toy needed) #############################################
    # intFuncVecY = integrateModel(modelPars = funcVec, coeffList=coeffList,histo=histo,kind='y')
    # intFuncVecQt = integrateModel(modelPars = funcVec, coeffList=coeffList,histo=histo, kind='qt')

    # intFuncErrY, intFuncErrQt = errorIntToy(fitPars=fitResult[0],cov=fitResult[1],dimYQt=fitResult[2],npBinCenters=fitResult[3],coeffList=coeffList,parNum=fitResult[4],histo=histo, npFitRes=fitResult[6], numNui = fitResult[7])
    
    
    
    # intFuncVecY = np.reshape(intFuncVecY,(len(coeffList),fitResult[2][0]))
    # intFuncErrY = np.reshape(intFuncErrY,(len(coeffList),fitResult[2][0]))
    
    # intFuncVecQt = np.reshape(intFuncVecQt,(len(coeffList),fitResult[2][1]))
    # intFuncErrQt = np.reshape(intFuncErrQt,(len(coeffList),fitResult[2][1]))
    
    # for li in coeffList :
    #     hy = histo[li[0]+'y'].Clone("post-fit-regularization_"+li[0]+'_y')
    #     for y in range(1, hy.GetNbinsX()+1) :
    #         # print("difference y", li[0], y,qt, ", original=", hy.GetBinContent(y,qt), ", apofitted=", funcVec[coeffList.index(li)][y-1][qt-1], ", ratio=", hy.GetBinContent(y,qt)/funcVec[coeffList.index(li)][y-1])
    #         hy.SetBinContent(y,intFuncVecY[coeffList.index(li)][y-1])  
    #         hy.SetBinError(y,intFuncErrY[coeffList.index(li)][y-1])  
    #     hy.Write()
        
    #     hqt = histo[li[0]+'qt'].Clone("post-fit-regularization_"+li[0]+'_qt')     
    #     for qt in range(1, hqt.GetNbinsX()+1) : 
    #         # print("difference qt", li[0], y,qt, ", original=", hqt.GetBinContent(y,qt), ", apofitted=", funcVec[coeffList.index(li)][y-1][qt-1], ", ratio=", hqt.GetBinContent(y,qt)/funcVec[coeffList.index(li)][y-1][hy.GetNbinsX()+qt-1])
    #         hqt.SetBinContent(qt,intFuncVecQt[coeffList.index(li)][qt-1])  
    #         hqt.SetBinError(qt,intFuncErrQt[coeffList.index(li)][qt-1])  
    #         if li[0]=='A0' : print(qt, intFuncVecQt[coeffList.index(li)][qt-1])
    #     hqt.Write()
    
    
    
    
    ########################################## integrated (with a 1D fit, disabled) ##########################################################################################
    # #build the error vector
    # funcVec = par2polyModel1D(modelPars=fitResult1D[0],coeffList=coeffList1D,npBinCenters=fitResult1D[3],parNum=fitResult1D[4], dimYQt=fitResult1D[2],npFitRes=fitResult1D[6],numNui=fitResult1D[7]) 
    
    # #build the function vector
    # errorVec = errorPoly(cov=fitResult1D[1],jacobian =fitResult1D[5])
    
    # #exclude nuisances
    # nuiVec = funcVec[-fitResult1D[7]:].copy()
    # nuiErr = errorVec[-fitResult1D[7]:].copy()
   
    # massVal = nuiVec[-fitResult1D[8]]
    # massErr = nuiErr[-fitResult1D[8]]
    
    # funcVec = funcVec[:-fitResult1D[7]]
    # errorVec = errorVec[:-fitResult1D[7]]
    
    # funcVec = np.reshape(funcVec,(len(coeffList),fitResult1D[2][0]+fitResult1D[2][1]))
    # errorVec = np.reshape(errorVec,(len(coeffList),fitResult1D[2][0]+fitResult1D[2][1]))
    
    # outFile.cd()
    # for li in coeffList :
    #     hy = histo[li[0]+'y'].Clone("post-fit-regularization_"+li[0]+'y')
    #     for y in range(1, hy.GetNbinsX()+1) :
    #         print("difference y", li[0], y,qt, ", original=", hy.GetBinContent(y,qt), ", apofitted=", funcVec[coeffList.index(li)][y-1][qt-1], ", ratio=", hy.GetBinContent(y,qt)/funcVec[coeffList.index(li)][y-1])
    #         hy.SetBinContent(y,funcVec[coeffList.index(li)][y-1])  
    #         hy.SetBinError(y,errorVec[coeffList.index(li)][y-1])  
    #     hy.Write()
        
    #     hy = histo[li[0]+'qt'].Clone("post-fit-regularization_"+li[0]+'qt')     
    #     for qt in range(1, hqt.GetNbinsY()+1) : 
    #         print("difference qt", li[0], y,qt, ", original=", hqt.GetBinContent(y,qt), ", apofitted=", funcVec[coeffList.index(li)][y-1][qt-1], ", ratio=", hqt.GetBinContent(y,qt)/funcVec[coeffList.index(li)][y-1][hy.GetNbinsX()+qt-1])
    #         hqt.SetBinContent(qt,funcVec[coeffList.index(li)][y-1][hy.GetNbinsX()+qt-1])  
    #         hqt.SetBinError(qt,errorVec[coeffList.index(li)][y-1][hy.GetNbinsX()+qt-1])  
    #     hqt.Write()
    # # print("predicted W mass = ", funcVec[-1][0][0], "+/-", errorVec[-1][0][0], ", only diag=",math.sqrt(fitResult1D[1][-1][-1])) 
    # print("WARNING: values muliplied x 50 to obtain MeV! good only if upd/down 50 MeV") # 1_nuisanceLimit : 50MeV_Given = massVal : x 
    # print("predicted W mass = ", 50*massVal, "+/-", 50*massErr)
    ####################################################################################################################################################################################    
    


parser = argparse.ArgumentParser("")

parser.add_argument('-o','--output', type=str, default='aposterioriFit_Wplus',help="name of the output file")
parser.add_argument('-f','--fitInput', type=str, default='fitPlots_Wplus.root',help="name of the fit result root file, after plotter_fitResult")
parser.add_argument('-r','--regInput', type=str, default='../../regularization/OUTPUT_1Apr/regularizationFit_range11____nom_nom_qt7_y5_limqt60_noC1.root',help="name of the regularization study result root file")#input file which must be used here
# parser.add_argument('-r','--regInput', type=str, default='',help="name of the regularization study result root file") #to have empty default
parser.add_argument('-s','--save', type=int, default=False,help="save .png and .pdf canvas")
parser.add_argument('-nUL','--notUnpol', type=int, default=False,help="exclude unpol and nuisances from the fit")
parser.add_argument('-af','--asimovFile', type=str, default='',help="name of the asimov file, if not empty run first the fit on the asimov, use the pars. as initial par and then on the real file.if empty initialized the pars=1 ")

args = parser.parse_args()
OUTPUT = args.output
FITINPUT = args.fitInput
REGINPUT = args.regInput
SAVE= args.save
NUL= args.notUnpol
ASIMOV= args.asimovFile

print("have you done cmsenv in CMSSW_11_2_0_pre8 instead of nightlies? (jax incomatibility)")
print("WARNING: the degree written in the dictionary is the first excluded--> N-1 used")
coeffList = []#y plus, qt plus, y minus, qt minus, constraint y, constraint qt,constrain C1, noRegul
# coeffList.append(['unpolarizedxsec', 3,3,-999,-999,  0,0,0, 0]) #ORIGINAL SET
# coeffList.append(['A0', 3,4,3,4, 0,1,1, 0])
# coeffList.append(['A1', 4,4,3,4, 1,1,0, 0])
# coeffList.append(['A2', 2,5,2,4, 0,1,1, 0])
# coeffList.append(['A3', 3,4,3,3, 0,1,0, 0])
# coeffList.append(['A4', 4,5,4,6, 1,0,0, 0]) #5-->7 but no convergence! and last number=1

coeffList.append(['unpolarizedxsec', 3,3,-999,-999,  0,0,0, 1]) #------------------------------------------------------------> this
# coeffList.append(['A0', 3,4,3,4, 0,1,1, 0])
# coeffList.append(['A1', 4,4,3,4, 1,1,0, 0])
# coeffList.append(['A2', 2,5,2,4, 0,1,1, 0])
# coeffList.append(['A3', 3,4,3,3, 0,1,0, 0])
# coeffList.append(['A4', 4,5,4,6, 1,0,0, 0]) 

#optimal 80 gev
# coeffList.append(['A0', 3,5,3,5, 0,1,0, 0])
# coeffList.append(['A1', 3,6,3,6, 1,1,0, 0])
# coeffList.append(['A2', 2,6,2,5, 0,1,0, 0])
# coeffList.append(['A3', 3,6,3,6, 0,1,0, 0])
# coeffList.append(['A4', 4,6,4,6, 1,0,0, 0]) 

#optimal 60
# coeffList.append(['A0', 3,4,3,5, 0,1,0, 0])
# coeffList.append(['A1', 3,6,3,6, 1,1,0, 0])
# coeffList.append(['A2', 2,5,2,4, 0,1,0, 0])
# coeffList.append(['A3', 3,5,3,5, 0,1,0, 0])
# coeffList.append(['A4', 4,6,4,5, 1,0,0, 0]) 

#working - no unpol
# coeffList.append(['A0', 3,4,3,5, 0,1,0, 0])
# coeffList.append(['A1', 3,5,3,6, 1,1,0, 0])
# coeffList.append(['A2', 2,5,2,4, 0,1,0, 0])
# coeffList.append(['A3', 3,5,3,5, 0,1,0, 0])
# coeffList.append(['A4', 4,5,4,5, 1,0,0, 0]) 

#working
coeffList.append(['A0', 3,4,3,4, 0,1,0, 0])
coeffList.append(['A1', 3,4,3,4, 1,1,0, 0])
coeffList.append(['A2', 3,4,3,4, 0,1,0, 0])
coeffList.append(['A3', 3,4,3,4, 0,1,0, 0])
coeffList.append(['A4', 4,4,4,4, 1,0,0, 0]) 

# coeffList1D = copy.deepcopy(coeffList)
if NUL : 
    del coeffList[0]

plusOnly=False
if plusOnly :
    signList = ['plus']
else :
    signList = ['plus', 'minus']
# print("WARNING: used wplus parameterization also for wminus, and not the dictionary above (ok if in the fit has been done the same)")

for s in signList : 
    FITINPUT_s = FITINPUT.replace('plus',s)
    OUTPUT_s = OUTPUT.replace('plus',s)
    ASIMOV_s = ASIMOV.replace('plus',s)
    
    coeffList_s = copy.deepcopy(coeffList)
    if s=='minus' :
        for li in coeffList_s : #use the element 3,4 of coeffList as 1,2 for minus, to access 1,2 everywhere also for minus.
            li[1] = li[3]
            li[2] = li[4]
        
    fitPostRegResult_asimov=''
    
    if ASIMOV!='' :
        fitResDict_asimov = getFitRes(inFile=ASIMOV_s, coeffList=coeffList_s)
        fitPostRegResult_asimov = polyFit(fitRes=fitResDict_asimov, coeffList=coeffList_s, nUL =NUL)
    
    fitResDict = getFitRes(inFile=FITINPUT_s, coeffList=coeffList_s)
    if REGINPUT!='' :
        print("Used regul. study fit results as input from file:", REGINPUT )
        regFuncParDict = getRegFuncPar(inFile=REGINPUT, coeffList=coeffList_s, sign=s)
    else :
        print("Used as starting parameter 1 everywhere")
        regFuncParDict=''
    # fitPostRegResult = polyFit(fitRes=fitResDict, regFunc=regFuncDict, coeffList=coeffList_s)
    fitPostRegResult = polyFit(fitRes=fitResDict, coeffList=coeffList_s, nUL =NUL, inPar=fitPostRegResult_asimov,regPar=regFuncParDict)
    # fitPostRegResult1D = polyFit1D(fitRes=fitResDict, coeffList=coeffList1D) #works but is without unpolarized and with lower degrees of polynomals
    
    # plotterPostReg(fitResult = fitPostRegResult, fitReult1D= fitPostRegResult1D, output=OUTPUT_s, save=SAVE, coeffList=coeffList_s,coeffList1D=coeffList1D, histo=fitResDict)
    plotterPostReg(fitResult = fitPostRegResult, output=OUTPUT_s, save=SAVE, coeffList=coeffList_s, histo=fitResDict, nUL = NUL)
    
