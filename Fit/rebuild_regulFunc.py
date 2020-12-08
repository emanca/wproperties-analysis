import ROOT
from array import array
import math
import copy
import argparse
import os
import numpy as np 
import root_numpy as rootnp
import ctypes
import scipy
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)


def getFitRes(inFile,coeffList, bins) :
    infile = ROOT.TFile.Open(inFile)
    inTree = infile.Get('fitresults')
    
    branchList = []
    for li in coeffList :
        if 'unpol' in li[0] : continue 
        for iy in range(0,bins[0]) :
            for iqt in range(0,bins[1]) :
                branchList.append('poly2d'+li[0]+'_polycoeff2d_'+str(iy)+'_'+str(iqt))
             
    npTree = rootnp.tree2array(inTree,branches=branchList, selection='1', start=0,stop=2)
    npTreeNP = np.zeros((len(coeffList)-1)*bins[0]*bins[1])
    for obj in npTree :
        counter=0
        for boh in obj :
            npTreeNP[counter] = boh
            counter+=1
    npTree = np.reshape(npTreeNP, ((len(coeffList)-1),bins[0],bins[1]))
    return npTree
    

def rebuildFunc(pars, coeffList, href) :
    
    dimY = href.GetNbinsX()
    dimQt = href.GetNbinsY()
    
    npBinCenters = np.zeros((dimQt*dimY,2))
    for iy in range(0, dimY): 
        for iqt in range(0, dimQt):
            itot = iy*(dimQt)+iqt
            npBinCenters[itot] = href.GetXaxis().GetBinCenter(iy+1), href.GetYaxis().GetBinCenter(iqt+1) 
    
    output = np.zeros((len(coeffList),dimY,dimQt))  
    valY = npBinCenters[:,0]
    valQt = npBinCenters[:,1]
    npBinCenters_res = npBinCenters.reshape(dimY,dimQt,2)
    valYunpol=npBinCenters_res[:,0,0]

    for li in coeffList :
        if li[0]=='unpolarizedxsec' :
            for iqt in range(dimQt) :
                output[0,:,iqt] = 0
        else :
            tempCoeff =np.zeros(np.shape(npBinCenters)[0])
            for iqt in range(0, dimQt) :
                for iy in range(0, dimY) :
                   tempCoeff = tempCoeff+pars[coeffList.index(li)-1][iy][iqt]*np.power(valY,iy)*np.power(valQt,iqt) 
            output[coeffList.index(li)] = tempCoeff.reshape((dimY,dimQt))
            
    return output    
    
    

# def rebuildError(pars, cov) :
#     ...


parser = argparse.ArgumentParser("")

parser.add_argument('-o','--output', type=str, default='fit_regularizationFunc',help="name of the output file")
parser.add_argument('-f','--fitInput', type=str,default='fit_Wplus_reco.root' ,help="name of the fit output root file")
parser.add_argument('-i','--extraInput', type=str,default='fitResult_asimovPlus_2nov_legacy.root' ,help="after-plotting input file")

args = parser.parse_args()
OUTPUT = args.output
INPUT = args.fitInput
EXTRAIN = args.extraInput


coeffList = []#y plus, qt plus, y minus, qt minus, constraint y, constraint qt
coeffList.append(['unpolarizedxsec', 3,3,-999,-999,0,0,0])
coeffList.append(['A0', 3,4,3,4, 0,1,1])
coeffList.append(['A1', 4,4,3,4, 1,1,0])
coeffList.append(['A2', 2,5,2,4, 0,1,1])
coeffList.append(['A3', 3,4,3,3, 0,1,0])
coeffList.append(['A4', 4,5,4,6, 1,0,0]) #5-->7 but no convergence!

extraInFile = ROOT.TFile.Open(EXTRAIN)
href = extraInFile.Get('coeff2D_reco/recoFitACA0')

parVec = getFitRes(inFile=INPUT,coeffList=coeffList,bins = [href.GetNbinsX(),href.GetNbinsY()])
funcVec = rebuildFunc(pars=parVec,coeffList=coeffList,href=href)

outFile = ROOT.TFile(OUTPUT+".root", "recreate")
for li in coeffList :
    h = href.Clone("post-fit-regularization_"+li[0])
    for y in range(1, h.GetNbinsX()+1) :
        for qt in range(1, h.GetNbinsY()+1) : 
            h.SetBinContent(y,qt,funcVec[coeffList.index(li)][y-1][qt-1])  
            h.SetBinError(y,qt,0)  
    outFile.cd()
    h.Write()
