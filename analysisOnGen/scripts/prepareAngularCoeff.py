import os
import ROOT
import copy
import sys
import argparse
import math


ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)


parser = argparse.ArgumentParser("")

parser.add_argument('-o','--output', type=str, default='genInput',help="name of the output file")
parser.add_argument('-i','--input', type=str, default='./GenInfo/genInfo.root',help="name of the input root file")
parser.add_argument('-u','--uncorrelate', type=int, default=True,help="if true uncorrelate num and den of Angular Coeff in MC scale variation")


args = parser.parse_args()
OUTPUT = args.output
INPUT = args.input
UNCORR = args.uncorrelate

coeffDict = {
    'A0' : 1.,
    'A1' : 5.,
    'A2' : 20.,
    'A3' : 4.,
    'A4' : 4.,
    'A5' : 5.,
    'A6' : 5.,
    'A7' : 4.,
    'AUL' : 1.,
}

systDict = {
    "_LHEScaleWeight" : ["_LHEScaleWeight_muR0p5_muF0p5", "_LHEScaleWeight_muR0p5_muF1p0","_LHEScaleWeight_muR1p0_muF0p5","_LHEScaleWeight_muR1p0_muF2p0","_LHEScaleWeight_muR2p0_muF1p0", "_LHEScaleWeight_muR2p0_muF2p0"],
    "_LHEPdfWeight" : ["_LHEPdfWeightHess" + str(i)  for i in range(1, 61)],
    "" : [""]
}

Wcharge = ["_Wplus","_Wminus"]
for charge in Wcharge:
    #GET HISTOS
    hDict = {}
    inFile = ROOT.TFile.Open(INPUT)
    for sKind, sList in systDict.items():
        for sName in sList :
            hDict[sName+'mapTot'] =  inFile.Get('angularCoefficients'+charge+sKind+'/mapTot'+sName)
            hDict[sName+'Y'] =  inFile.Get('angularCoefficients'+charge+sKind+'/Y'+sName)
            hDict[sName+'Pt'] =  inFile.Get('angularCoefficients'+charge+sKind+'/Pt'+sName)
            for coeff,div in coeffDict.items() :
                hDict[sName+coeff] =  inFile.Get('angularCoefficients'+charge+sKind+'/harmonics'+coeff+sName)
                hDict[sName+coeff+'Err'] =  inFile.Get('angularCoefficients'+charge+sKind+'/harmonicsSq'+coeff+sName)
                hDict[sName+coeff+'Y'] =  inFile.Get('angularCoefficients'+charge+sKind+'/harmonicsY'+coeff+sName)
                hDict[sName+coeff+'Pt'] =  inFile.Get('angularCoefficients'+charge+sKind+'/harmonicsPt'+coeff+sName)
    # get maps
    mapTot = inFile.Get('angularCoefficients'+charge+'/mapTot')
    mapAccEta = inFile.Get('angularCoefficients'+charge+'/mapAccEta')
    mapAcc = inFile.Get('angularCoefficients'+charge+'/mapAcc')
    sumw = inFile.Get('angularCoefficients'+charge+'/sumw')

    #BUILD OUTPUT
    outFile = ROOT.TFile(OUTPUT+charge+'.root', "RECREATE")
    outFile.cd()
    for sKind, sList in systDict.items():    
        outFile.mkdir('angularCoefficients'+sKind)
        outFile.cd('angularCoefficients'+sKind)

        sListMod = copy.deepcopy(sList)
        if sKind=='_LHEScaleWeight' and UNCORR :
            sListMod.append("") #add nominal variation

        for sName in sListMod :
                hDict[sName+'mapTot'].Write()
                hDict[sName+'Y'].Write()
                hDict[sName+'Pt'].Write()
                for sNameDen in sListMod :
                    if sNameDen!=sName and not (sKind=='_LHEScaleWeight' and UNCORR) : #PDF or correlated Scale
                        continue 
                    for coeff,div in coeffDict.items() :
                        hist = hDict[sName+coeff].Clone()
                        histY = hDict[sName+coeff+'Y'].Clone()
                        histPt = hDict[sName+coeff+'Pt'].Clone()
                        if coeff!='AUL' : 
                            hist.Divide(hDict[sNameDen+'mapTot'])
                            histY.Divide(hDict[sNameDen+'Y'])
                            histPt.Divide(hDict[sNameDen+'Pt'])
                            hist.Scale(div)
                            histY.Scale(div)
                            histPt.Scale(div)
                        if coeff=='A0' :
                            for xx in range(1,hist.GetNbinsX()+1) :
                                for yy in range(1,hist.GetNbinsY()+1) :
                                    content = hist.GetBinContent(xx,yy)
                                    hist.SetBinContent(xx,yy,20./3*(content+1./10.))
                            for xx in range(1,histY.GetNbinsX()+1):
                                    content = histY.GetBinContent(xx)
                                    histY.SetBinContent(xx,20./3*(content+1./10.))
                            for xx in range(1,histPt.GetNbinsX()+1):
                                    content = histPt.GetBinContent(xx)
                                    histPt.SetBinContent(xx,20./3*(content+1./10.))


                        #error assignment                    
                        if coeff!='AUL' : 
                            for xx in range(1,hist.GetNbinsX()+1) :
                                for yy in range(1,hist.GetNbinsY()+1) :
                                    N_eff = hDict[sNameDen+'mapTot'].GetBinContent(xx,yy)**2/(hDict[sNameDen+'mapTot'].GetBinError(xx,yy)**2)
                                    f2w = hDict[sName+coeff+'Err'].GetBinContent(xx,yy)/hDict[sNameDen+'mapTot'].GetBinContent(xx,yy)
                                    fw2 = (hDict[sName+coeff].GetBinContent(xx,yy)/hDict[sName+'mapTot'].GetBinContent(xx,yy))**2
                                    A_err = math.sqrt(f2w - fw2)/math.sqrt(N_eff)
                                    A_err = A_err*div
                                    if coeff=='A0' :
                                        A_err = 20./3*(A_err)
                                    hist.SetBinError(xx,yy, A_err )

                        # if sKind=="" :
                        #     for xx in range(1,histPt.GetNbinsX()+1) :
                        #         N_eff = hDict[sNameDen+'mapTot'].GetBinContent(xx)**2/(hDict[sNameDen+'mapTot'].GetBinError(xx)**2)
                        #         f2w = hDict[sName+coeff+'Err'].GetBinContent(xx)/hDict[sNameDen+'mapTot'].GetBinContent(xx)
                        #         fw2 = (histPt.GetBinContent(xx))**2
                        #         A_err = math.sqrt(f2w - fw2)/math.sqrt(N_eff)
                        #         A_err = A_err/div
                        #         if coeff=='A0' :
                        #             A_err = 20./3*(A_err+1./10.)
                        #         histPt.SetBinError(i,j, A_err )
                        #     for xx in range(1,histY.GetNbinsX()+1):
                        #         N_eff = hDict[sNameDen+'mapTot'].GetBinContent(xx)**2/(hDict[sNameDen+'mapTot'].GetBinError(xx)**2)
                        #         f2w = hDict[sName+coeff+'Err'].GetBinContent(xx)/hDict[sNameDen+'mapTot'].GetBinContent(xx)
                        #         fw2 = (histY.GetBinContent(xx))**2
                        #         A_err = math.sqrt(f2w - fw2)/math.sqrt(N_eff)
                        #         A_err = A_err/div
                        #         if coeff=='A0' :
                        #             A_err = 20./3*(A_err+1./10.)
                        #         histY.SetBinError(i,j, A_err )

                        suff = ''
                        suffDen = ''
                        if sName == "" : suff = '_nom'
                        if sNameDen == '' : suffDen = '_nom'                    
                        hist.SetName(hist.GetName()+suff+sNameDen+suffDen)
                        hist.SetTitle(hist.GetTitle()+suff+sNameDen+suffDen)
                        hist.Write()
                        histY.SetName(histY.GetName()+suff+sNameDen+suffDen)
                        histY.SetTitle(histY.GetTitle()+suff+sNameDen+suffDen)
                        histY.Write()
                        histPt.SetName(histPt.GetName()+suff+sNameDen+suffDen)
                        histPt.SetTitle(histPt.GetTitle()+suff+sNameDen+suffDen)
                        histPt.Write()

    outFile.mkdir('accMaps')
    outFile.cd('accMaps')
    mapTot.Write()
    mapAccEta.Write()
    mapAcc.Write()
    sumw.Write()


    outFile.Close()

