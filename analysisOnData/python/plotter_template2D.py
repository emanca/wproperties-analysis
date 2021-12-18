import os
import ROOT
import copy
import sys
sys.path.append('../data')
from systematics import systematics
import argparse
import math

ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)

class plotter:
    
    def __init__(self, outDir, inDir = '',simult=False):

        self.indir = inDir # indir containig the various outputs
        self.outdir = outDir
        self.simult= simult
        
        self.sampleDict = { "data_obs"      :  ('WToMu_plots.root', 2 ),
                            "DYJets"      : ('DYJets_plots.root', 2 ),
                            "WtoTau"      : ('WToTau_plots.root', 2 ),
                            "Top"         : ('Top_plots.root', 1),
                            "DiBoson"     : ('Diboson_plots.root', 1),
                            "Fake" : ('FakeFromData_plots.root', 0), 
                            "Data"        : ('Data_plots.root', 0),
                            "LowAcc": ('WToMu_AC_plots.root', 2)
                          }

        self.selections =["Signal"] 
        self.selectionFake = ['fakes']       
        self.extSyst = copy.deepcopy(systematics)
        self.extSyst['Nominal'] = ['']
        self.extSyst['jme'] = (['_jesTotalUp', '_jesTotalDown', '_unclustEnUp', '_unclustEnDown'],"jme")
        self.extSyst["ptScale"] =  (["_correctedUp", "_correctedDown"],"ptScale") 
        self.extSyst["lumi"] =  (['_lumiUp', '_lumiDown'], 'lumi'),
        del self.extSyst['Nom_WQTlow'] #remove wqt-shape "nom" variation (already included in fake vars)
        del self.extSyst['Nom_WQTmid'] 
        del self.extSyst['Nom_WQThigh'] 
        self.histoDict ={} 

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
    def getHistoforSample(self, sample, infile, chargeBin) :
        # if not 'Fake' in sample:
            for sKind in self.extSyst:
                if 'lumi' in sKind and 'Fake' not in sample : continue 
                if 'WQT' in sKind and sample not in ['Fake','LowAcc','WtoTau']: continue
                self.histoDict[sample][sKind] = []
                basepath = 'templates_Signal/' + sKind
                if 'LowAcc' in sample: basepath = 'templatesLowAcc_Signal/' + sKind
                if 'Fake' in sample: basepath = 'templates_fakes/' + sKind
                # if 'Fake' in sample and sKind!='Nominal' and sKind!='LHEScaleWeight' and sKind!='jme' : continue  
                if infile.GetDirectory(basepath):
                    for key in infile.Get(basepath).GetListOfKeys():
                        th3=infile.Get(basepath+'/'+key.GetName())
                        if sample!='Fake' :
                            if sKind=='LHEScaleWeight_WQTlow' :
                                print(basepath.replace(sKind,'Nom_WQTmid')+'/'+key.GetName().split('_')[0], sample)
                                th3.Add(infile.Get(basepath.replace(sKind,'Nom_WQTmid')+'/'+key.GetName().split('_')[0]))
                                th3.Add(infile.Get(basepath.replace(sKind,'Nom_WQThigh')+'/'+key.GetName().split('_')[0]))
                            if sKind=='LHEScaleWeight_WQTmid' :
                                th3.Add(infile.Get(basepath.replace(sKind,'Nom_WQTlow')+'/'+key.GetName().split('_')[0]))
                                th3.Add(infile.Get(basepath.replace(sKind,'Nom_WQThigh')+'/'+key.GetName().split('_')[0]))
                            if sKind=='LHEScaleWeight_WQThigh' :
                                th3.Add(infile.Get(basepath.replace(sKind,'Nom_WQTlow')+'/'+key.GetName().split('_')[0]))
                                th3.Add(infile.Get(basepath.replace(sKind,'Nom_WQTmid')+'/'+key.GetName().split('_')[0]))
                            if 'WQT' in sKind : 
                                th3.SetName(th3.GetName()+sKind.replace('LHEScaleWeight','') )# add WQTlow/mid/high at the end of the name
                        
                        #plus charge bin 2, minus bin 1
                        th3.GetZaxis().SetRange(chargeBin,chargeBin)
                        th2=th3.Project3D("yx")
                        th2.SetDirectory(0)
                        th2.SetName(th3.GetName().replace('templates',sample))
                        self.histoDict[sample][sKind].append(th2)
                        
        # else:
        #     for sKind in self.extSyst:
        #         self.histoDict[sample][sKind] = []
        #         basepath = 'templates_fakes/' + sKind
        #         if infile.GetDirectory(basepath):
        #             print(basepath)
        #             for key in infile.Get(basepath).GetListOfKeys():
        #                 print(key.GetName())
        #                 th3=infile.Get(basepath+'/'+key.GetName())
        #                 #plus charge bin 2, minus bin 1
        #                 th3.GetZaxis().SetRange(chargeBin,chargeBin)
        #                 th2=th3.Project3D("yx")
        #                 print(th2.GetName())
        #                 th2.SetDirectory(0)
        #                 th2.SetName(th3.GetName())
        #                 self.histoDict[sample][sKind].append(th2)
    def uncorrelateEff(self, sample):
        aux = {}
        aux['WHSF'] = []
        for h in self.histoDict[sample]['Nominal']:
            if 'mass' in h.GetName():
                continue
            for i in range(3):
                for hvar in self.histoDict[sample]['WHSF']:
                    for updown in ['Up', 'Down']:
                        #print h.GetName() + 'WHSFSyst{}{}'.format(i, updown), "match"
                        if hvar.GetName() == h.GetName() + '_WHSFSyst{}{}'.format(i, updown):
                            for j in range(1, hvar.GetNbinsX()+1):  # loop over eta bins
                                #create one histogram per eta bin
                                haux = h.Clone()
                                haux.SetName(
                                    h.GetName() + '_WHSFSyst{}Eta{}{}'.format(i, j, updown))
                                #print haux.GetName()
                                for k in range(1, hvar.GetNbinsY()+1):  # loop over pt bins
                                    bin1D = hvar.GetBin(j, k)
                                    varcont = hvar.GetBinContent(bin1D)
                                    haux.SetBinContent(bin1D, varcont)
                                aux['WHSF'].append(haux)
        for hvar in self.histoDict[sample]['WHSF']:
            if 'Flat' in hvar.GetName():  # leave syst uncertainty as it is
                aux['WHSF'].append(hvar)
        self.histoDict[sample].update(aux)
    def uncorrelateJME(self, sample):
        aux = {}
        aux['jme'] = []
        for h in self.histoDict[sample]['Nominal']:
            if 'mass' in h.GetName():
                continue
            for i in ['jesTotal','unclustEn']:
                for hvar in self.histoDict[sample]['jme']:
                    for updown in ['Up', 'Down']:
                        #print h.GetName() + 'WHSFSyst{}{}'.format(i, updown), "match"
                        if hvar.GetName() == h.GetName() + '_{}{}'.format(i, updown):
                            for j in range(1, hvar.GetNbinsY()+1):  # loop over pt bins
                                #create one histogram per eta bin
                                haux = h.Clone()
                                haux.SetName(
                                    h.GetName() + '_{}pt{}{}'.format(i, j, updown))
                                #print haux.GetName()
                                for k in range(1, hvar.GetNbinsX()+1):  # loop over eta bins
                                    bin1D = hvar.GetBin(k, j)
                                    varcont = hvar.GetBinContent(bin1D)
                                    haux.SetBinContent(bin1D, varcont)
                                aux['jme'].append(haux)
        self.histoDict[sample].update(aux)
    
    def uncorrelateJME_eta(self, sample):
        aux = {}
        aux['jme'] = []
        for h in self.histoDict[sample]['Nominal']:
            if 'mass' in h.GetName():
                continue
            for i in ['jesTotal','unclustEn']:
                for hvar in self.histoDict[sample]['jme']:
                    for updown in ['Up', 'Down']:
                        #print h.GetName() + 'WHSFSyst{}{}'.format(i, updown), "match"
                        if hvar.GetName() == h.GetName() + '_{}{}'.format(i, updown):
                            for j in range(1, hvar.GetNbinsX()+1):  # loop over eta bins
                                #create one histogram per eta bin
                                haux = h.Clone()
                                haux.SetName(
                                    h.GetName() + '_{}eta{}{}'.format(i, j, updown))
                                #print haux.GetName()
                                for k in range(1, hvar.GetNbinsY()+1):  # loop over pt bins
                                    bin1D = hvar.GetBin(j, k)
                                    varcont = hvar.GetBinContent(bin1D)
                                    haux.SetBinContent(bin1D, varcont)
                                aux['jme'].append(haux)
        self.histoDict[sample].update(aux)

    def symmetrisePDF(self,sample):
        if not self.histoDict[sample]['LHEPdfWeight']==[]:
            aux = {}
            aux['LHEPdfWeight']=[]
            for h in self.histoDict[sample]['Nominal']:
                if 'mass' in h.GetName(): continue
                for i in range(60):
                    for hvar in self.histoDict[sample]['LHEPdfWeight']:
                        if hvar.GetName() == h.GetName()+ '_LHEPdfWeightHess{}'.format(i+1):
                            th2var = hvar
                            break
                    th2c = h.Clone()
                    th2varD = th2var.Clone()
                    th2var.Divide(th2c)
                    th2c.Divide(th2varD)
                    nbinsX = h.GetXaxis().GetNbins()
                    nbinsY = h.GetYaxis().GetNbins()
                    th2Up = ROOT.TH2D("up","up",nbinsX,h.GetXaxis().GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(nbinsX),nbinsY,h.GetYaxis().GetBinLowEdge(1),h.GetYaxis().GetBinUpEdge(nbinsY))
                    th2Up.Sumw2()
                    th2Down =ROOT.TH2D("down","down",nbinsX,h.GetXaxis().GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(nbinsX),nbinsY,h.GetYaxis().GetBinLowEdge(1),h.GetYaxis().GetBinUpEdge(nbinsY))
                    th2Down.Sumw2()
                    for j in range(1,h.GetNbinsX()+1):
                        for k in range(1,h.GetNbinsY()+1):
                        
                            th2Up.SetBinContent(j,k,h.GetBinContent(j,k)*th2var.GetBinContent(j,k))
                            th2Down.SetBinContent(j,k,h.GetBinContent(j,k)*th2c.GetBinContent(j,k))
                    th2Up.SetName(h.GetName()+ '_LHEPdfWeightHess{}Up'.format(i+1))
                    th2Down.SetName(h.GetName()+ '_LHEPdfWeightHess{}Down'.format(i+1))
                    aux['LHEPdfWeight'].append(th2Up)
                    aux['LHEPdfWeight'].append(th2Down)
            self.histoDict[sample].update(aux)
    def symmetriseSyst(self,sample):
        for sKind,sList in self.extSyst.items() :
            if not "LHEScale" in sKind : continue #oly pDF, MCscale need symmetrization
            # if 'WQT' in sKind and sample!='Fake' : continue
            if 'WQT' in sKind and sample not in ['Fake','LowAcc','WtoTau']: continue
            if sKind == 'LHEScaleWeight' and sample!='DYJets' : continue # only DY has flat scale
            if not self.histoDict[sample][sKind]==[]:
                aux = {}
                aux[sKind]=[]
                for h in self.histoDict[sample]['Nominal']:
                    if 'mass' in h.GetName(): continue
                    for sName in sList[0]:
                        for hvar in self.histoDict[sample][sKind]:
                            if hvar.GetName() == h.GetName()+ sName:
                                th2var = hvar
                                break
                        th2c = h.Clone()
                        th2varD = th2var.Clone()
                        th2var.Divide(th2c)
                        th2c.Divide(th2varD)
                        nbinsX = h.GetXaxis().GetNbins()
                        nbinsY = h.GetYaxis().GetNbins()
                        th2Up = ROOT.TH2D("up","up",nbinsX,h.GetXaxis().GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(nbinsX),nbinsY,h.GetYaxis().GetBinLowEdge(1),h.GetYaxis().GetBinUpEdge(nbinsY))
                        th2Up.Sumw2()
                        th2Down =ROOT.TH2D("down","down",nbinsX,h.GetXaxis().GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(nbinsX),nbinsY,h.GetYaxis().GetBinLowEdge(1),h.GetYaxis().GetBinUpEdge(nbinsY))
                        th2Down.Sumw2()
                        for j in range(1,h.GetNbinsX()+1):
                            for k in range(1,h.GetNbinsY()+1):
                                th2Up.SetBinContent(j,k,h.GetBinContent(j,k)*th2var.GetBinContent(j,k))
                                th2Down.SetBinContent(j,k,h.GetBinContent(j,k)*th2c.GetBinContent(j,k))
                        th2Up.SetName(h.GetName()+ sName+'Up')
                        th2Down.SetName(h.GetName()+ sName+'Down')
                        aux[sKind].append(th2Up)
                        aux[sKind].append(th2Down)
                self.histoDict[sample].update(aux)
    def symmetriseSyst_shift(self,sample): #nom+-var (insted of nom*(nom/var), nom*(var/nom))
        for sKind,sList in self.extSyst.items() :
            if not "LHE" in sKind : continue #oly pDF, MCscale need symmetrization
            # if 'WQT' in sKind and sample!='Fake' : continue
            if 'WQT' in sKind and sample not in ['Fake','LowAcc','WtoTau']: continue
            if sKind == 'LHEScaleWeight' and sample!='DYJets' : continue # only DY has flat scale
            if not self.histoDict[sample][sKind]==[]:
                aux = {}
                aux[sKind]=[]
                for h in self.histoDict[sample]['Nominal']:
                    if 'mass' in h.GetName(): continue
                    for sName in sList[0]:
                        for hvar in self.histoDict[sample][sKind]:
                            if hvar.GetName() == h.GetName()+ sName:
                                th2var = hvar
                                break
                        hdiff = hvar.Clone()
                        hdiff.Add(h,-1)
                        nbinsX = h.GetXaxis().GetNbins()
                        nbinsY = h.GetYaxis().GetNbins()
                        th2Up = ROOT.TH2D("up","up",nbinsX,h.GetXaxis().GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(nbinsX),nbinsY,h.GetYaxis().GetBinLowEdge(1),h.GetYaxis().GetBinUpEdge(nbinsY))
                        th2Up.Sumw2()
                        th2Down =ROOT.TH2D("down","down",nbinsX,h.GetXaxis().GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(nbinsX),nbinsY,h.GetYaxis().GetBinLowEdge(1),h.GetYaxis().GetBinUpEdge(nbinsY))
                        th2Down.Sumw2()
                        for j in range(1,h.GetNbinsX()+1):
                            for k in range(1,h.GetNbinsY()+1):
                                th2Up.SetBinContent(j, k, h.GetBinContent(j, k)+hdiff.GetBinContent(j, k))
                                th2Down.SetBinContent(j, k, h.GetBinContent(j, k)-hdiff.GetBinContent(j, k))
                        th2Up.SetName(h.GetName()+ sName+'Up')
                        th2Down.SetName(h.GetName()+ sName+'Down')
                        aux[sKind].append(th2Up)
                        aux[sKind].append(th2Down)
                self.histoDict[sample].update(aux)
    def alphaVariations(self,sample):
        aux = {}
        aux['alphaS'] = []
        for h in self.histoDict[sample]['Nominal']:
            if 'mass' in h.GetName():
                continue
            for hvar in self.histoDict[sample]['alphaS']:
                if hvar.GetName() == h.GetName() + '_alphaSUp' or hvar.GetName() == h.GetName() + '_alphaSDown':
                    hvar.Divide(h)
                    aux['alphaS'].append(hvar)
        self.histoDict[sample].update(aux)
    def getHistos(self, chargeBin) :
        print('writing histograms')
        foutName = 'Wplus_reco' if chargeBin == 2 else 'Wminus_reco'
        foutName += '.root'
        fout = ROOT.TFile.Open(self.outdir + '/' + foutName, "UPDATE")
        for sample,fname in self.sampleDict.items():
            # if not 'LowAcc' in sample and not 'data_obs' in sample:
            #     continue
            print("Processing sample:", sample)
            infile = ROOT.TFile(self.indir + '/' + fname[0])
            systType = fname[1]
            if not infile:
                print(infile,' does not exist')
                continue
            self.histoDict[sample] = {}
            self.getHistoforSample(sample,infile, chargeBin)
            if sample not in  self.histoDict : 
                print("No histo dict for sample:", sample, " What have you done??!!!!")
                continue
            self.symmetrisePDF(sample)
            self.symmetriseSyst(sample)
            # self.uncorrelateJME_eta(sample)
            # self.symmetriseSyst_shift(sample)
            self.uncorrelateEff(sample)
            # self.alphaVariations(sample)
            for syst, hlist in self.histoDict[sample].items():
                #fout.mkdir(syst)
                #fout.cd(syst)
                fout.cd()
                for h in hlist:
                    #th1=self.unroll2D(h)
                    if self.simult :
                        if 'data_obs' not in sample :
                            chargeS = foutName.replace('_reco.root','')
                            h.SetName(chargeS+h.GetName())
                            if 'Prefire' in h.GetName() or 'WHSF' in h.GetName() or 'jesTotal' in h.GetName() or 'unclustEn' in h.GetName() or 'corrected' in h.GetName() : 
                                # print(h.GetName(), chargeS+sample+'_'+chargeS+h.GetName().replace(chargeS+sample+'_',''))
                                h.SetName(chargeS+sample+'_'+chargeS+h.GetName().replace(chargeS+sample+'_',''))
                    h.Write()
                #fout.cd()
        fout.Save()
        fout.Write()

parser = argparse.ArgumentParser("")
parser.add_argument('-o','--output', type=str, default='./',help="name of the output directory")
parser.add_argument('-i','--input', type=str, default='./',help="name of the input directory root file")
parser.add_argument('-s', '--simult',      type=int, default=False, help="simultaneous fit of two W charges")

args = parser.parse_args()
OUTPUT = args.output
INPUT = args.input
SIMULT = args.simult

p=plotter(outDir=OUTPUT, inDir = INPUT,simult = SIMULT)
p.getHistos(1)
p.getHistos(2)
