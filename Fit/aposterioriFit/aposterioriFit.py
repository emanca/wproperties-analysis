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
        outDict[li[0]] = openFile.Get('coeff2D_reco/recoFitAC'+li[0])
    outDict['cov']=openFile.Get('matrices_reco/covariance_matrix_channelhelpois') #unrolled order: y0qt0a0.....y0qtNa0, y1qt0a0...y1qt0a0, ...
    outDict['mass'] = openFile.Get('matrices_reco/massreco')
    print("input W mass=", outDict['mass'].GetBinContent(1), "+/-", outDict['mass'].GetBinError(1))
    return outDict 
    
# def getRegFunc(inFile, coeffList) :
#     outDict = {}
#     openFile = ROOT.TFile.Open(inFile)
#     for li in coeffList :
#         if li[0]!='unpolarizedxsec' : 
#             outDict[li[0]] = openFile.Get('Y'+str(li[1])+'_Qt'+str(li[2])+'/hfit_WtoMuP_'+li[0]+'_'+str(li[1])+str(li[2]))
#         else :
#             temph = openFile.Get('unpol_Y3/hfit_WtoMuP_unpol_3')
#             for qt in range(1, temph.GetNbinsY()+1) :
#                 outDict[li[0]+str(qt)] = openFile.Get('unpol_Y'+str(li[1])+'/hfit_WtoMuP_unpol_3_'+str(qt))
#             print("WARNING: unpolarized qt binning must be aligned to the fit qt binning!")
#     return outDict

def chi2MassOnlyFit(modelPars,npFitRes,npCovMatInv) :
    print("shapes",np.shape(modelPars), np.shape(npFitRes), np.shape(npCovMatInv))
    print("vals",modelPars, npFitRes, npCovMatInv, 100*math.sqrt(1/npCovMatInv)) #1/npnpCovMatInv??
    diff = npFitRes-modelPars
    # chi2 = jnp.matmul(diff.T, jnp.matmul(npCovMatInv, diff) )
    chi2 = diff*diff*npCovMatInv
    return chi2
    
def par2polyModel(modelPars, npFitRes, coeffList, npBinCenters,parNum,dimYQt) :#input = the parameters, output the polyinomial function
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
                        fitModelList_temp.append(fitModelList_temp[-1]+modelPars[iqt*(effY)+iy+ipar]*np.power(valY,powY)*np.power(valQt,powQt) )
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
    fitModel = jnp.append(fitModel,modelPars[-1])
    
    return fitModel
     
def chi2PolyFit(modelPars,npFitRes, npCovMatInv, coeffList, npBinCenters,parNum,dimYQt) :
    
    fitModel = par2polyModel(modelPars= modelPars, coeffList= coeffList, npBinCenters= npBinCenters,parNum= parNum,dimYQt= dimYQt,npFitRes=npFitRes)
    diff = npFitRes-fitModel
    chi2 = jnp.matmul(diff.T, jnp.matmul(npCovMatInv, diff) )

    return chi2
    

# def polyFit(fitRes,regFunc, coeffList) :
def polyFit(fitRes, coeffList) :
    
    dimY = fitRes[coeffList[0][0]].GetNbinsX()
    dimQt = fitRes[coeffList[0][0]].GetNbinsY()
    dimCoeff = len(coeffList)
    dimYQt = [dimY,dimQt]

    npCovMatUncut = rootnp.hist2array(fitRes['cov'])
    npCovMat = npCovMatUncut[0:dimCoeff*dimQt*dimY,0:dimCoeff*dimQt*dimY]
    # npCovMat = npCovMatUncut[1*dimQt*dimY:6*dimQt*dimY,1*dimQt*dimY:6*dimQt*dimY]
    
    #add mass to the convariance matrix
    massBin = fitRes['cov'].GetXaxis().FindBin('mass')
    massColumn = npCovMatUncut[massBin-1,0:dimCoeff*dimQt*dimY]
    massRow = npCovMatUncut[0:dimCoeff*dimQt*dimY,massBin-1]
    massEl = npCovMatUncut[massBin-1,massBin-1]
    npCovMat = np.c_[npCovMat,massColumn]
    massRow = np.append(massRow,massEl)
    npCovMat = np.r_[npCovMat,[massRow]]
     
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
    npFitRes = np.stack(npFitResList,0)  #dimensions: coeff, y, qt = (6,6,8). if flatten--> q0y0a0...qNy0a0, q0y1a0....qNy1a0, ....
    npFitMass = np.zeros(1)
    npFitMass[0] = fitRes['mass'].GetBinContent(1)
    
    npBinCenters = np.zeros((dimQt*dimY,2))
    for iy in range(0, dimY): #y
        for iqt in range(0, dimQt):#qt
            itot = iy*(dimQt)+iqt
            npBinCenters[itot] = fitRes[li[0]].GetXaxis().GetBinCenter(iy+1), fitRes[li[0]].GetYaxis().GetBinCenter(iqt+1) 
            
    #model init
    modelParsList = []
    parNum = []
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
                    print(li[0], "qt=",qt, "pars=", parNum[-1])
            else :
                effY = li[1]
                effQt = li[2]
                if li[5] : effY-= 1
                if li[6] : effQt-= 1
                if li[7] : effQt-= 1
                modelParsC=np.zeros((effY,effQt))
                for qt in range(0, effQt) :
                    for y in range(0, effY) :
                        modelParsC[y,qt] = 1. 
                modelParsList.append(modelParsC.copy())
                print(li[0], "pars (y,qt)=", effY, effQt)
                parNum.append(effY*effQt)
                modelParsList[-1]=modelParsList[-1].flatten()
        else : #li[8]=not regularize this
            # modelParsC = rootnp.hist2array(fitRes[li[0]])
            # modelParsList.append(modelParsC.copy())
            modelParsList.append(np.full((dimY,dimQt),1.))
            print(li[0], "pars (y,qt)=", np.shape(modelParsList[-1]))
            parNum.append( dimY*dimQt)
            modelParsList[-1]=modelParsList[-1].flatten()
            
            
    modelPars = np.concatenate(modelParsList) #already flattened, otherwise YOU HAVE TO FLATTEN before pass it to the fit
    
    modelMass= np.zeros(1)
    modelMass[0] = 1.
    modelPars = np.append(modelPars,modelMass) #added the mass to the fit
    
    coeffListJAX = copy.deepcopy(coeffList)
    for li in coeffListJAX :
        if li[0] != 'unpolarizedxsec' : 
            li[0] = coeffListJAX.index(li)-1
        else :
            li[0] = 5
    
    npFitRes = npFitRes.flatten() #to add the mass
    npFitRes = np.append(npFitRes,npFitMass) #mass
    print("par num", parNum)

    print("everything initialized, minimizer call...")
    modelPars = pmin(chi2PolyFit, modelPars, args=(npFitRes, npCovMatInv, coeffListJAX, npBinCenters,parNum,dimYQt), doParallel=False)
    # modelPars = pmin(chi2MassOnlyFit, modelPars[-1], args=(npFitRes[-1], npCovMatInv[-1][-1]), doParallel=False)
    print ("fit ended, post fit operation...")
        
    # after fit results
    static_argnumsPF = (1,2,3,4,5,6)
    modelChi2Grad_eval = jax.jit(jax.value_and_grad(chi2PolyFit), static_argnums=static_argnumsPF)
    modelHess_eval = jax.jit(jax.hessian(chi2PolyFit), static_argnums=static_argnumsPF)
    modelJac_eval = jax.jit(jax.jacfwd(par2polyModel), static_argnums=(1,2,3,4,5))
    
    modelChi2,modelGrad = modelChi2Grad_eval(modelPars, npFitRes, npCovMatInv, coeffListJAX, npBinCenters,parNum,dimYQt)
    modelHess = modelHess_eval(modelPars, npFitRes, npCovMatInv, coeffListJAX, npBinCenters,parNum,dimYQt)
    modelJac = modelJac_eval(modelPars, npFitRes, coeffListJAX, npBinCenters,parNum,dimYQt)
    
    modelCov = np.linalg.inv(0.5*modelHess)
    modelErr = np.sqrt(np.diag(modelCov))

    modelEDM = 0.5*np.matmul(np.matmul(modelGrad.T,modelCov),modelGrad)
    modelNDOF = dimCoeff*dimQt*dimY+1-np.size(modelPars) #+1= mass 
    
    
    
    
    print("*-------------------------------------------*")
    print("FIT RESULTS")
    print("edm=",modelEDM) 
    print("chi2/dof=", modelChi2,"/",modelNDOF, "=", modelChi2/float(modelNDOF), "pars=",np.size(modelPars)) 
    print("parameters=", modelPars)
    print("errors=", modelErr)
    print("relative uncertaintiy" , np.true_divide(modelErr,modelPars))
    
    print("check covariance matrix:")
    print("is covariance semi positive definite?", np.all(np.linalg.eigvals(modelCov) >= 0))
    print("is covariance symmetric?", np.allclose(np.transpose(modelCov), modelCov))
    print("is covariance approx symmetric?", np.allclose(np.transpose(modelCov), modelCov, rtol=1e-05, atol=1e-06))
    # print("mass uncertainity if everything fixed=", math.sqrt(npCovMatInv[-1][-1]))
    print("*-------------------------------------------*")
    # print("covariance=", modelCov.diagonal())
    
    outFit = (modelPars,modelCov,dimYQt, npBinCenters,parNum,modelJac,npFitRes)    
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
    # ntoys = 10000
    # covtoy ={}
    
    # # toyVec = np.zeros((len(coeffList)+1,dimYQt[0],dimYQt[1],ntoys))  #+1= mass block
    # toyEl = np.zeros((len(coeffList),dimYQt[0],dimYQt[1]))  
    # toyEl.flatten()
    # toyEl = np.append(toyEl, 0.)#mass
    # toyVec = np.zeros((np.shape(toyEl)[0],ntoys))

    # for itoy in range(ntoys) :
    #     covtoy[itoy] = np.random.multivariate_normal(fitPars, cov,tol=1e-6)
    #     # toyVec[...,itoy] = rebuildPoly(fitPars=covtoy[itoy],dimYQt=dimYQt,npBinCenters=npBinCenters, coeffList=coeffList,parNum=parNum)
    #     toyVec[...,itoy] = par2polyModel(modelPars=covtoy[itoy],coeffList=coeffList,npBinCenters=npBinCenters,parNum=parNum, dimYQt=dimYQt)
        
    # sigmavec = np.std(toyVec,axis=-1)

    # return sigmavec
    
def errorPoly(cov, jacobian) :
    covariancePoly = np.matmul(jacobian,np.matmul(cov,jacobian.T))
    variances = covariancePoly.diagonal()
    stddev = np.sqrt(variances)
    return stddev
            
def plotterPostReg(fitResult, output, save, coeffList,histo) :
    print("plotting...")
    
    #build the error vector
    # funcVec = rebuildPoly(fitPars=fitResult[0],dimYQt=fitResult[2],npBinCenters=fitResult[3], coeffList=coeffList,parNum=fitResult[4])
    funcVec = par2polyModel(modelPars=fitResult[0],coeffList=coeffList,npBinCenters=fitResult[3],parNum=fitResult[4], dimYQt=fitResult[2],npFitRes=fitResult[6]) 
    
    #build the function vector
    errorVec = errorPoly(cov=fitResult[1],jacobian =fitResult[5])
    
    massVal = funcVec[-1].copy()
    massErr = errorVec[-1].copy()
    funcVec = np.reshape(funcVec[:-1],(len(coeffList),fitResult[2][0],fitResult[2][1]))
    errorVec = np.reshape(errorVec[:-1],(len(coeffList),fitResult[2][0],fitResult[2][1]))
    
    outFile = ROOT.TFile(output+".root", "recreate")
    for li in coeffList :
        h = histo[li[0]].Clone("post-fit-regularization_"+li[0])
        for y in range(1, h.GetNbinsX()+1) :
            for qt in range(1, h.GetNbinsY()+1) : 
                h.SetBinContent(y,qt,funcVec[coeffList.index(li)][y-1][qt-1])  
                h.SetBinError(y,qt,errorVec[coeffList.index(li)][y-1][qt-1])  
        outFile.cd()
        h.Write()
    # print("predicted W mass = ", funcVec[-1][0][0], "+/-", errorVec[-1][0][0], ", only diag=",math.sqrt(fitResult[1][-1][-1])) 
    print("WARNING: values muliplied x 100 to obtain MeV! good only if upd/down 100 MeV")
    print("predicted W mass = ", 100*massVal, "+/-", 100*massErr, ", only diag=",math.sqrt(fitResult[1][-1][-1])) 
        
        
    


parser = argparse.ArgumentParser("")

parser.add_argument('-o','--output', type=str, default='aposterioriFit',help="name of the output file")
parser.add_argument('-f','--fitInput', type=str, default='fitResult.root',help="name of the fit result root file, after plotter_fitResult")
parser.add_argument('-r','--regInput', type=str, default='../../regularization/OUTPUT_poly/regularizationFit_range11_rebuild____nom_nom.root',help="name of the regularization study result root file")
parser.add_argument('-s','--save', type=int, default=False,help="save .png and .pdf canvas")

args = parser.parse_args()
OUTPUT = args.output
FITINPUT = args.fitInput
REGINPUT = args.regInput
SAVE= args.save

print("have you done cmsenv in CMSSW_11_2_0_pre8 instead of nightlies? (jax incomatibility)")
coeffList = []#y plus, qt plus, y minus, qt minus, constraint y, constraint qt,constrain C1, noRegul
# coeffList.append(['unpolarizedxsec', 3,3,-999,-999,  0,0,0, 0]) #ORIGINAL SET
# coeffList.append(['A0', 3,4,3,4, 0,1,1, 0])
# coeffList.append(['A1', 4,4,3,4, 1,1,0, 0])
# coeffList.append(['A2', 2,5,2,4, 0,1,1, 0])
# coeffList.append(['A3', 3,4,3,3, 0,1,0, 0])
# coeffList.append(['A4', 4,5,4,6, 1,0,0, 0]) #5-->7 but no convergence! and last number=1

coeffList.append(['unpolarizedxsec', 3,3,-999,-999,  0,0,0, 1])
coeffList.append(['A0', 3,4,3,4, 0,1,1, 0])
coeffList.append(['A1', 4,4,3,4, 1,1,0, 0])
coeffList.append(['A2', 2,5,2,4, 0,1,1, 0])
coeffList.append(['A3', 3,4,3,3, 0,1,0, 0])
coeffList.append(['A4', 4,5,4,6, 1,0,0, 0]) #5-->7 but no convergence! and last number=1

fitResDict = getFitRes(inFile=FITINPUT, coeffList=coeffList)
# regFuncDict = getRegFunc(inFile=REGINPUT, coeffList=coeffList)
# fitPostRegResult = polyFit(fitRes=fitResDict, regFunc=regFuncDict, coeffList=coeffList)
fitPostRegResult = polyFit(fitRes=fitResDict, coeffList=coeffList)
plotterPostReg(fitResult = fitPostRegResult, output=OUTPUT,save=SAVE, coeffList=coeffList,histo=fitResDict)