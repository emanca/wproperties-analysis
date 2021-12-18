import itertools
from termcolor import colored
from time import time
import ROOT
import copy
import sys
import argparse
from collections import OrderedDict
sys.path.append('../data')
from systematics import systematics
import math
import copy
import multiprocessing

ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)


def uncorrelate1PtVar(t):
    h = t[0]
    hvar = t[1]
    i = t[2]
    updown = t[3]
    uncorrh = []
    for j in range(1, 4+1):  # eta macro bins
        #create one histogram per macro bin
        haux = h.Clone()
        haux.SetName(h.GetName() +'_stateig{}Eta{}{}'.format(i, j, updown))
        #print haux.GetName()
        for t in range(1, 12+1):  # 48 eta bins /4 macro bins
            for k in range(1, hvar.GetNbinsY()+1):  # loop over pt bins
                bin1D = hvar.GetBin((j-1)*12+t, k)
                varcont = hvar.GetBinContent(bin1D)
                haux.SetBinContent(bin1D, varcont)
        uncorrh.append(haux)
    return [uncorrh[0], uncorrh[1], uncorrh[2], uncorrh[3]]

class plotter:
    
    def __init__(self, outDir, inDir = '',simult=False):
        self.indir = inDir # indir containig the various outputs
        self.outdir = outDir
        self.simult = simult
        self.sampleFile = 'WToMu_AC_plots.root'

        self.extSyst = copy.deepcopy(systematics)
        # self.extSyst['Nominal'] =  (['', '_massUp', '_massDown'],"Nom")
        self.extSyst['Nominal'] =  ([''],"Nom")
        self.extSyst['jme'] = (['_jesTotalUp', '_jesTotalDown', '_unclustEnUp', '_unclustEnDown'],"jme")
        self.extSyst["ptScale"] =  (["_correctedUp", "_correctedDown"],"ptScale") 

        
        del self.extSyst['LHEScaleWeight_WQTlow'] #remove wqt-shape variation (because we are fitting the wqt)
        del self.extSyst['LHEScaleWeight_WQTmid'] 
        del self.extSyst['LHEScaleWeight_WQThigh'] 
        del self.extSyst['Nom_WQTlow'] 
        del self.extSyst['Nom_WQTmid'] 
        del self.extSyst['Nom_WQThigh'] 

        self.charges = ["Wplus", "Wminus"]
        self.yields = {}
        self.yields["Wplus"] = {}
        self.yields["Wminus"] = {}
        self.bins = []
        self.nBinsQt = 11#13#9#29#8
        self.nBinsY = 6#7
        
        for iY in range(1, self.nBinsY+1):
            for iQt in range(1, self.nBinsQt+1):
                self.yields["Wplus"][(iY, iQt)] = 0.
                self.yields["Wminus"][(iY, iQt)] = 0.
                self.bins.append("y_{}_qt_{}".format(iY,iQt))
        self.clist = ["L", "I", "T", "A", "P", "7", "8", "9", "UL"]
        
        self.histoDict ={}
        for charge in self.charges:
            self.histoDict[charge] = {}
            for c in self.clist:
                self.histoDict[charge][c] = {}
                for bin in self.bins:
                    self.histoDict[charge][c][bin] = {}
                    for sKind in self.extSyst:
                        self.histoDict[charge][c][bin][sKind] = []

        self.varName = 'helXsecs'
        self.inFile = ROOT.TFile.Open(self.indir+'/'+self.sampleFile)
        if not self.inFile :
            print(self.inFile, ' does not exist')
            sys.exit(1)
        self.fileACplus = ROOT.TFile.Open("../../analysisOnGen/genInput_Wplus.root")
        self.fileACminus = ROOT.TFile.Open("../../analysisOnGen/genInput_Wminus.root")
        # self.fileACplus = ROOT.TFile.Open("../../analysisOnGen/genInput_4GeVBinning_udShift_Wplus.root")
        # self.fileACminus = ROOT.TFile.Open("../../analysisOnGen/genInput_4GeVBinning_udShift_Wminus.root")
        self.imapPlus = self.fileACplus.Get("angularCoefficients/mapTot")
        self.imapMinus = self.fileACminus.Get("angularCoefficients/mapTot")
        self.closPlus = copy.deepcopy(ROOT.TH2D("closPlus", "closPlus", self.imapPlus.GetXaxis().GetNbins(),self.imapPlus.GetXaxis().GetXbins().GetArray(),self.imapPlus.GetYaxis().GetNbins(),self.imapPlus.GetYaxis().GetXbins().GetArray()))
        self.closMinus = copy.deepcopy(ROOT.TH2D("closMinus", "closMinus", self.imapPlus.GetXaxis().GetNbins(), self.imapPlus.GetXaxis().GetXbins().GetArray(), self.imapPlus.GetYaxis().GetNbins(), self.imapPlus.GetYaxis().GetXbins().GetArray()))
        
        
        self.helXsecs = OrderedDict()
        self.helXsecs["L"] = "A0"
        self.helXsecs["I"] = "A1" 
        self.helXsecs["T"] = "A2" 
        self.helXsecs["A"] = "A3" 
        self.helXsecs["P"] = "A4" 
        self.helXsecs["7"] = "A5" 
        self.helXsecs["8"] = "A6" 
        self.helXsecs["9"] = "A7" 
        self.helXsecs["UL"] = "AUL"
        
        self.factors = {}
        self.factors["A0"]= 2.
        self.factors["A1"]=2.*math.sqrt(2)
        self.factors["A2"]=4.
        self.factors["A3"]=4.*math.sqrt(2)
        self.factors["A4"]=2.
        self.factors["A5"]=2.
        self.factors["A6"]=2.*math.sqrt(2)
        self.factors["A7"]=4.*math.sqrt(2)
        self.factors["AUL"]=1.
    def unroll2D(self, th2):

        nbins = th2.GetNbinsX()*th2.GetNbinsY()
        #th2.Sumw2() #don't think it's necessary
        new = th2.GetName()
        old = new + '_roll'
        th2.SetName(old)
        unrolledth2 = ROOT.TH1F(new, '', nbins, 1., nbins+1)
        for ibin in range(1, th2.GetNbinsX()+1):
            for jbin in range(1, th2.GetNbinsY()+1):
                bin1D = th2.GetBin(ibin, jbin)
                unrolledth2.SetBinContent(bin1D, th2.GetBinContent(ibin, jbin))
                unrolledth2.SetBinError(bin1D, th2.GetBinError(ibin, jbin))
        return unrolledth2
    def normaliseYields(self, th2, coeff, iY, iQt, imap, fAC, sKind,sName):
        #normalise templates to its helicity xsec
        if sKind in ['LHEScaleWeight','LHEPdfWeight','alphaS']:
            imap = fAC.Get("angularCoefficients_"+sKind+"/mapTot"+sName)
        
        nsum = (3./16./math.pi)*imap.GetBinContent(iY, iQt)
        if not 'UL' in coeff:
            if sKind in ['LHEScaleWeight','LHEPdfWeight','alphaS']:
                hAC = fAC.Get("angularCoefficients_"+sKind+"/harmonics"+str(self.helXsecs[coeff])+sName+""+sName)
            else :
                hAC = fAC.Get("angularCoefficients/harmonics{}_nom_nom".format(self.helXsecs[coeff]))
            nsum = nsum * hAC.GetBinContent(iY, iQt)/self.factors[self.helXsecs[coeff]]
        th2.Scale(nsum)
        return
    def makeTH3slices(self, th3, systname):
        
        hname=th3.GetName()
        charge = hname.split('_')[0]
        iQt = int(hname.split('_')[2])
        coeff = hname.split('_')[4]
        
        try:
            # syst = "_" + hname.split('_')[5]
            syst = "_" + '_'.join(hname.split('_')[5:])
        except IndexError:
            syst = ""
        if syst == '_' : syst=''        
        for iY in range(1, th3.GetNbinsZ()+1):
            
            if not self.simult :
                slicename = 'helXsecs'+ coeff + '_y_{}'.format(iY)+'_qt_{}'.format(iQt) + syst
            else :
                if 'Prefire' in syst or 'WHSF' in syst or 'jesTotal' in syst or 'unclustEn' in syst or 'corrected' in syst : 
                    chargeSting = '_'+charge
                elif syst!='' :
                    chargeSting = '_'
                else :
                    chargeSting = '' #nominal
                slicename = charge+'helXsecs'+ coeff + '_y_{}'.format(iY)+'_qt_{}'.format(iQt) +chargeSting + syst[1:]
                # print("syst=",syst,slicename )

            th3.GetZaxis().SetRange(iY, iY)
            
            th2slice=th3.Project3D("y_{iY}_yxe".format(iY=iY))
            th2slice.SetName(slicename)
            if charge == "Wplus":
                self.normaliseYields(th2slice, coeff, iY, iQt, self.imapPlus, self.fileACplus,systname,syst)
            else:
                self.normaliseYields(th2slice, coeff, iY, iQt, self.imapMinus, self.fileACminus,systname,syst)
            th2slice.SetDirectory(0)
            self.histoDict[charge][coeff]["y_{}_qt_{}".format(iY, iQt)][systname].append(th2slice)

            if syst=="":
                self.yields[charge][(iY,iQt)]+=th2slice.Integral(0,th2slice.GetNbinsX()+2,0,th2slice.GetNbinsY()+2)
    def closureMap(self):
        for iY in range(1, self.nBinsY+1):
            for iQt in range(1, self.nBinsQt+1):
                self.closPlus.SetBinContent(iY,iQt, self.yields["Wplus"][(iY,iQt)])
                self.closMinus.SetBinContent(iY,iQt, self.yields["Wminus"][(iY, iQt)])
        self.sumOfTemplatePlus = self.closPlus.Clone("sumOfTemplatePlus")
        self.sumOfTemplateMinus = self.closPlus.Clone("sumOfTemplateMinus")
        self.closPlus.Divide(self.imapPlus)
        self.closMinus.Divide(self.imapMinus)
        for iY in range(1, self.nBinsY+1):
            for iQt in range(1, self.nBinsQt+1):
                print(colored(self.closPlus.GetBinContent(iY,iQt),'magenta'))
        fout = ROOT.TFile("accMap.root","recreate")
        fout.cd()
        self.closPlus.Write()
        self.closMinus.Write()
        self.imapPlus.Write()
        self.imapMinus.Write()
        self.sumOfTemplatePlus.Write()
        self.sumOfTemplateMinus.Write()
        fout.Close()
    def getHistos(self) :
        basepath="templatesAC_Signal/"
        for charge in self.charges:
            for c in self.clist:
                for iQt in range(1, self.nBinsQt+1):
                    for sKind, sList in self.extSyst.items():
                        for sname in sList[0]:  # variations of each sKind
                            fpath = basepath + sKind + '/'+charge+ '_qt_{}_helXsecs_'.format(iQt) + c + sname
                            print("Histo read:", fpath)
                            th3 = self.inFile.Get(fpath)
                            if not th3: print(colored('fpath not found', 'red'))
                            self.makeTH3slices(th3, sKind)
                            #assert(0)
        self.uncorrelateEff()
        self.symmetrisePDF()
        self.symmetriseScale()
        # self.uncorrelateJME_eta()
        # self.symmetrisePDF_shift()
        # self.alphaVariations()
        # self.uncorrelatePtVars() ALWAYS COMMENTED, ROCHESTER CORR. NOT RUN
        self.closureMap()
        # self.rescaleTempl() #1/6 UL, +1/6UL for the others
        self.writeHistos()
    def uncorrelatePtVars(self):
        for charge in self.charges:
            for c in self.clist:
                for bin in self.bins:
                    hlist = []
                    hnom = [h for h in self.histoDict[charge][c][bin]['Nominal'] if not 'mass' in h.GetName()][0]
                    for hvar in self.histoDict[charge][c][bin]['ptScale']:
                        for j in range(1, 4+1):  # eta macro bins
                            #create one histogram per macro bin
                            haux = hnom.Clone()
                            if 'stateig' in hvar.GetName():
                                haux.SetName(hvar.GetName().replace('stateig', 'Eta{}stateig'.format(j)))
                            elif 'zptsyst' in hvar.GetName():
                                haux.SetName(hvar.GetName().replace('zptsyst', 'Eta{}zptsyst'.format(j)))
                            elif 'Ewksyst' in hvar.GetName():
                                haux.SetName(hvar.GetName().replace('Ewksyst', 'Eta{}Ewksyst'.format(j)))
                            elif 'deltaM' in hvar.GetName():
                                haux.SetName(hvar.GetName().replace('deltaM', 'Eta{}deltaM'.format(j)))
                            #print haux.GetName()
                            for t in range(1, 12+1):  # 48 eta bins /4 macro bins
                                for k in range(1, hvar.GetNbinsY()+1):  # loop over pt bins
                                    bin1D = hvar.GetBin((j-1)*12+t, k)
                                    varcont = hvar.GetBinContent(bin1D)
                                    haux.SetBinContent(bin1D, varcont)
                            hlist.append(haux)
                    self.histoDict[charge][c][bin]['ptScale'] = hlist
    def uncorrelateEff(self):
        for charge in self.charges:
            for c in self.clist:
                for bin in self.bins:
                    hlist = []
                    hnom = [h for h in self.histoDict[charge][c][bin]['Nominal'] if not 'mass' in h.GetName()][0]
                    for hvar in self.histoDict[charge][c][bin]['WHSF']:
                        if not 'Flat' in hvar.GetName():
                            for j in range(1, hvar.GetNbinsX()+1):  # loop over eta bins
                                #create one histogram per eta bin
                                haux = hnom.Clone()
                                if 'Up' in hvar.GetName():
                                    haux.SetName(hvar.GetName().replace('Up', 'Eta{}Up'.format(j)))
                                else:
                                    haux.SetName(hvar.GetName().replace('Down', 'Eta{}Down'.format(j)))
                                #print haux.GetName()
                                for k in range(1, hvar.GetNbinsY()+1):  # loop over pt bins
                                    bin1D = hvar.GetBin(j, k)
                                    varcont = hvar.GetBinContent(bin1D)
                                    haux.SetBinContent(bin1D, varcont)
                                hlist.append(haux)
                        else:
                            hlist.append(hvar)
                    self.histoDict[charge][c][bin]['WHSF'] = hlist
    def uncorrelateJME(self):
        for charge in self.charges:
            for c in self.clist:
                for bin in self.bins:
                    hlist = []
                    hnom = [h for h in self.histoDict[charge][c][bin]['Nominal'] if not 'mass' in h.GetName()][0]
                    for hvar in self.histoDict[charge][c][bin]['jme']:
                        # hlist.append(hvar) #fully correlated variation
                        for j in range(1, hvar.GetNbinsY()+1):  # loop over pt bins
                            #create one histogram per eta bin
                            haux = hnom.Clone()
                            if 'Up' in hvar.GetName():
                                haux.SetName(hvar.GetName().replace('Up', 'pt{}Up'.format(j)))
                            else:
                                haux.SetName(hvar.GetName().replace('Down', 'pt{}Down'.format(j)))
                            #print haux.GetName()
                            for k in range(1, hvar.GetNbinsX()+1):  # loop over eta bins
                                bin1D = hvar.GetBin(k, j)
                                varcont = hvar.GetBinContent(bin1D)
                                haux.SetBinContent(bin1D, varcont)
                            hlist.append(haux)
                    self.histoDict[charge][c][bin]['jme'] = hlist
    def uncorrelateJME_eta(self):
        for charge in self.charges:
            for c in self.clist:
                for bin in self.bins:
                    hlist = []
                    hnom = [h for h in self.histoDict[charge][c][bin]['Nominal'] if not 'mass' in h.GetName()][0]
                    for hvar in self.histoDict[charge][c][bin]['jme']:
                        # hlist.append(hvar) #fully correlated variation
                        for j in range(1, hvar.GetNbinsX()+1):  # loop over eta bins
                            #create one histogram per eta bin
                            haux = hnom.Clone()
                            if 'Up' in hvar.GetName():
                                haux.SetName(hvar.GetName().replace('Up', 'eta{}Up'.format(j)))
                            else:
                                haux.SetName(hvar.GetName().replace('Down', 'eta{}Down'.format(j)))
                            #print haux.GetName()
                            for k in range(1, hvar.GetNbinsY()+1):  # loop over pt bins
                                bin1D = hvar.GetBin(j, k)
                                varcont = hvar.GetBinContent(bin1D)
                                haux.SetBinContent(bin1D, varcont)
                            hlist.append(haux)
                    self.histoDict[charge][c][bin]['jme'] = hlist
    def symmetrisePDF(self):
        for charge in self.charges:
            for c in self.clist:
                for bin in self.bins:
                    hlist = []
                    hnom = [h for h in self.histoDict[charge][c][bin]['Nominal'] if not 'mass' in h.GetName()][0]
                    for hvar in self.histoDict[charge][c][bin]['LHEPdfWeight']:
                        print(hvar.GetName())
                        th2c = hnom.Clone()
                        th2varD = hvar.Clone()
                        hvar.Divide(th2c)
                        th2c.Divide(th2varD)
                        nbinsX = hnom.GetXaxis().GetNbins()
                        nbinsY = hnom.GetYaxis().GetNbins()
                        th2Up = ROOT.TH2D("up", "up", nbinsX, hnom.GetXaxis().GetBinLowEdge(1), hnom.GetXaxis().GetBinUpEdge(nbinsX), nbinsY, hnom.GetYaxis().GetBinLowEdge(1), hnom.GetYaxis().GetBinUpEdge(nbinsY))
                        th2Up.Sumw2()
                        th2Down = ROOT.TH2D("down", "down", nbinsX, hnom.GetXaxis().GetBinLowEdge(1), hnom.GetXaxis().GetBinUpEdge(nbinsX), nbinsY, hnom.GetYaxis().GetBinLowEdge(1), hnom.GetYaxis().GetBinUpEdge(nbinsY))
                        th2Down.Sumw2()
                        for j in range(1, hnom.GetNbinsX()+1):
                            for k in range(1, hnom.GetNbinsY()+1):
                                th2Up.SetBinContent(j, k, hnom.GetBinContent(j, k)*hvar.GetBinContent(j, k))
                                th2Down.SetBinContent(j, k, hnom.GetBinContent(j, k)*th2c.GetBinContent(j, k))
                        th2Up.SetName(hvar.GetName() + 'Up')
                        th2Down.SetName(hvar.GetName() + 'Down')
                        hlist.append(th2Up)
                        hlist.append(th2Down)
                    self.histoDict[charge][c][bin]['LHEPdfWeight'] = hlist
    def symmetriseScale(self):
        for charge in self.charges:
            for c in self.clist:
                for bin in self.bins:
                    hlist = []
                    hnom = [h for h in self.histoDict[charge][c][bin]['Nominal'] if not 'mass' in h.GetName()][0]
                    for hvar in self.histoDict[charge][c][bin]['LHEScaleWeight']:
                        print(hvar.GetName())
                        th2c = hnom.Clone()
                        th2varD = hvar.Clone()
                        hvar.Divide(th2c)
                        th2c.Divide(th2varD)
                        nbinsX = hnom.GetXaxis().GetNbins()
                        nbinsY = hnom.GetYaxis().GetNbins()
                        th2Up = ROOT.TH2D("up", "up", nbinsX, hnom.GetXaxis().GetBinLowEdge(1), hnom.GetXaxis().GetBinUpEdge(nbinsX), nbinsY, hnom.GetYaxis().GetBinLowEdge(1), hnom.GetYaxis().GetBinUpEdge(nbinsY))
                        th2Up.Sumw2()
                        th2Down = ROOT.TH2D("down", "down", nbinsX, hnom.GetXaxis().GetBinLowEdge(1), hnom.GetXaxis().GetBinUpEdge(nbinsX), nbinsY, hnom.GetYaxis().GetBinLowEdge(1), hnom.GetYaxis().GetBinUpEdge(nbinsY))
                        th2Down.Sumw2()
                        for j in range(1, hnom.GetNbinsX()+1):
                            for k in range(1, hnom.GetNbinsY()+1):
                                th2Up.SetBinContent(j, k, hnom.GetBinContent(j, k)*hvar.GetBinContent(j, k))
                                th2Down.SetBinContent(j, k, hnom.GetBinContent(j, k)*th2c.GetBinContent(j, k))
                        th2Up.SetName(hvar.GetName() + 'Up')
                        th2Down.SetName(hvar.GetName() + 'Down')
                        hlist.append(th2Up)
                        hlist.append(th2Down)
                    self.histoDict[charge][c][bin]['LHEScaleWeight'] = hlist
    def symmetrisePDF_shift(self): #nom+-var (insted of nom*(nom/var), nom*(var/nom))
        for charge in self.charges:
            for c in self.clist:
                for bin in self.bins:
                    hlist = []
                    hnom = [h for h in self.histoDict[charge][c][bin]['Nominal'] if not 'mass' in h.GetName()][0]
                    for hvar in self.histoDict[charge][c][bin]['LHEPdfWeight']:
                        print(hvar.GetName())
                        hdiff = hvar.Clone()
                        hdiff.Add(hnom,-1)
                        nbinsX = hnom.GetXaxis().GetNbins()
                        nbinsY = hnom.GetYaxis().GetNbins()
                        th2Up = ROOT.TH2D("up", "up", nbinsX, hnom.GetXaxis().GetBinLowEdge(1), hnom.GetXaxis().GetBinUpEdge(nbinsX), nbinsY, hnom.GetYaxis().GetBinLowEdge(1), hnom.GetYaxis().GetBinUpEdge(nbinsY))
                        th2Up.Sumw2()
                        th2Down = ROOT.TH2D("down", "down", nbinsX, hnom.GetXaxis().GetBinLowEdge(1), hnom.GetXaxis().GetBinUpEdge(nbinsX), nbinsY, hnom.GetYaxis().GetBinLowEdge(1), hnom.GetYaxis().GetBinUpEdge(nbinsY))
                        th2Down.Sumw2()
                        for j in range(1, hnom.GetNbinsX()+1):
                            for k in range(1, hnom.GetNbinsY()+1):
                                th2Up.SetBinContent(j, k, hnom.GetBinContent(j, k)+hdiff.GetBinContent(j, k))
                                th2Down.SetBinContent(j, k, hnom.GetBinContent(j, k)-hdiff.GetBinContent(j, k))
                        th2Up.SetName(hvar.GetName() + 'Up')
                        th2Down.SetName(hvar.GetName() + 'Down')
                        hlist.append(th2Up)
                        hlist.append(th2Down)
                        self.histoDict[charge][c][bin]['LHEPdfWeight'] = hlist
    def alphaVariations(self):
        for charge in self.charges:
            for c in self.clist:
                for bin in self.bins:
                    hlist = []
                    hnom = [h for h in self.histoDict[charge][c][bin]['Nominal'] if not 'mass' in h.GetName()][0]
                    for hvar in self.histoDict[charge][c][bin]['alphaS']:
                        hvar.Divide(hnom)
                        hlist.append(hvar)
                    self.histoDict[charge][c][bin]['alphaS'] = hlist
    def writeHistos(self):
        print('writing histograms')
        for charge in self.charges:
            foutName = charge+'_reco'
            foutName += '.root'
            fout = ROOT.TFile.Open(self.outdir + '/' + foutName, "recreate")
            fout.cd()
            for c in self.clist:
                for bin in self.bins:
                    for sKind in self.extSyst:
                        for h in self.histoDict[charge][c][bin][sKind]:
                            fout.cd()
                            h.Write()
            fout.cd()
            fout.Save()
            fout.Close()
    # def rescaleTempl(self) :
    #     for charge in self.charges:
    #         for bin in self.bins: 
    #             for sKind, sList in self.extSyst.items() :
    #                 for hind in range(0, len(self.histoDict[charge]['UL'][bin][sKind])) :
    #                     # print("UL main=", self.histoDict[charge]['UL'][bin][sKind][hind].GetName())
    #                     self.histoDict[charge]['UL'][bin][sKind][hind].Scale(1/6.)
    #                     for c in self.clist:
    #                         if 'UL' not in c :
    #                             # print(self.histoDict[charge][c][bin][sKind][hind].GetName(),self.histoDict[charge]['UL'][bin][sKind][hind].GetName())
    #                             self.histoDict[charge][c][bin][sKind][hind].Add(self.histoDict[charge]['UL'][bin][sKind][hind])


parser = argparse.ArgumentParser("")
parser.add_argument('-o','--output', type=str, default='./',help="name of the output directory")
parser.add_argument('-i','--input', type=str, default='./',help="name of the input directory root file")
parser.add_argument('-s', '--simult',      type=int, default=False, help="simultaneous fit of two W charges")

args = parser.parse_args()
OUTPUT = args.output
INPUT = args.input
SIMULT = args.simult
p=plotter(outDir=OUTPUT, inDir = INPUT, simult=SIMULT)
p.getHistos()
