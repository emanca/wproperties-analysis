import os
import sys
import json
import argparse
import ROOT
import pprint 
import math
import numpy as np
from array import array
import scipy
import copy
from scipy import special
from scipy import stats
# from ROOT import *
from ROOT import gStyle
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)
import multiprocessing


pp = pprint.PrettyPrinter(indent=2)

parser = argparse.ArgumentParser("")
parser.add_argument('-o', '--output_dir',type=str, default='TEST', help="")
parser.add_argument('-i','--input_gen', type=str, default='../analysisOnGen/genInput_fineBinned_4regularization_fix_Wplus.root ../analysisOnGen/genInput_fineBinned_4regularization_fix_Wminus.root',nargs='*', help="list of two gen input files (plus, minus), separated by space")
parser.add_argument('-validation_only', '--validation_only',type=int, default=False, help="False: skip fit and do validation only")
parser.add_argument('-input_reg', '--input_reg',type=str, default='regularizationFit', help="Name of the input file for validation_only (without .root)")
parser.add_argument('-map01', '--map01',type=int, default=False, help="if true map the polynomial between 0 and 1, otherwise between -1,1")
parser.add_argument('-suffix', '--suffix',type=str, default='', help="suffix for output file")
parser.add_argument('-syst_kind', '--syst_kind',type=str, default='', help="name of the kind of syst empty=nominal")
parser.add_argument('-syst_name', '--syst_name',type=str, default='_nom_nom', help="name of the particular syst empty=nominal")
parser.add_argument('-syst_val', '--syst_val',type=int, default=False, help="activate the systs validation, using all the syst of the dictionary")
parser.add_argument('-ncores', '--ncores',type=int, default=1, help="number of cores used")

args = parser.parse_args()
output_dir = args.output_dir
VALONLY = args.validation_only
IN_REG = args.input_reg
IN_GEN = args.input_gen
SUFFIX = args.suffix
MAP01 = args.map01
SYST_kind= args.syst_kind
SYST_name= args.syst_name
SYST_VALIDATION = args.syst_val
NCORES=args.ncores
if NCORES>1 :
    MULTICORE=True
else :
    MULTICORE=False


if VALONLY :
    SUFFIX = 'rebuild_'+SUFFIX
    
if MAP01 :
    SUFFIX = 'range01_'+SUFFIX
else :
    SUFFIX = 'range11_'+SUFFIX

# if not SYST_VALIDATION :
    SUFFIX = SUFFIX+'_'+SYST_kind+'_'+SYST_name


def buildSystDict() :
    outDict = {}
        
    outDict[''] = ['_nom_nom'] 
    outDict["_LHEPdfWeight"] = ["_LHEPdfWeightHess" + str(i) +'_LHEPdfWeightHess'+str(i) for i in range(1, 61)]
    scaleList = ["_nom", "muR0p5_muF0p5", "muR0p5_muF1p0", "muR1p0_muF0p5","muR1p0_muF2p0","muR2p0_muF1p0","muR2p0_muF2p0"]
    uncorrList = []
    for sName in scaleList :
        if sName!='_nom' : sName = '_LHEScaleWeight_'+sName
        for sNameDen in scaleList :    
            if sName==sNameDen and sName=='_nom' : continue
            if sNameDen!='_nom' : sNameDen = '_LHEScaleWeight_'+sNameDen
            uncorrList.append(sName+sNameDen)
    outDict["_LHEScaleWeight"] = uncorrList 
    return outDict
    
    
systDict  = buildSystDict()

signDict = {
        'Plus' : 'WtoMuP',
        'Minus' : 'WtoMuN'
        }
        
# coeffDict = {
#         'A0' : 'qt->0',
#         'A1' : 'qt->0, y->0',
#         'A2' : 'qt->0',
#         'A3' : 'qt->0, y->0',
#         'A4' : 'y->0'
#     }
# coeffDict = {
#         'A0' : 'qt->0, C1',
#         'A1' : 'qt->0, y->0, ',
#         'A2' : 'qt->0, C1',
#         'A3' : 'qt->0',
#         'A4' : 'y->0, C1'
#     }
coeffDict = {
        'A0' : 'qt->0',
        'A1' : 'qt->0, y->0, ',
        'A2' : 'qt->0',
        'A3' : 'qt->0',
        'A4' : 'y->0'
    }
# coeffDict = {
#         'A0' : 'qt->0',
#         'A1' : 'qt->0',
#         'A2' : 'qt->0',
#         'A3' : 'qt->0',
#         'A4' : ''
#     }
        
#oder of Chebyshev polinomials:
# degreeYList = [2,3,4,5,6,7]
degreeQtList = [2,3,4,5,6,7]
# degreeYList = [2,3]
# degreeQtList = [2,3]
degreeYList = [2,3,4,5]
# degreeQtList = [2,3,4,5]




# ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")

groupedSystColors = {
        # "WHSFVars"  : [ROOT.kGreen+1, 'Scale Factors'],
        "_LHEScaleWeight" : [ROOT.kViolet-2, 'MC Scale'],
        # "ptScaleVars" : [ROOT.kBlue-4, 'pT Scale'],
        # "jmeVars" : [ROOT.kAzure+10, 'MET'],
        "_LHEPdfWeight" : [ROOT.kRed+1, 'PDF'],
        "" : [1, 'Stat Unc.']
        }









class fit_func01:

    def __init__(self, id_y,id_qt, constraints='qt->0'):
        self.id_qt = id_qt
        self.len_qt = int(len(id_qt)-1)
        self.id_y = id_y
        self.len_y = int(len(id_y)-1)
        self.nparams =  len(id_qt)*len(id_y)
        self.limqt = 'qt->0' in constraints
        self.limy = 'y->0' in constraints
        
        # even-odd dependent members
        if self.len_y%2==0 : 
            self.yeven = True
        else : 
            self.yeven = False
        
        if self.len_qt%2==0 : 
            self.qteven = True
        else : 
            self.qteven = False
            

    def __call__(self, x, parameters):
        y,qt = x[0],x[1]
        val = 0.0
        p = [0.]*self.nparams                
        
        for iqt,vqt in enumerate(self.id_qt):
            for iy,vy in enumerate(self.id_y):
                itot = iqt*(len(self.id_y))+iy
                p[itot] = parameters[itot]

        if self.limqt:
            p1 = copy.deepcopy(p)
            for iqt,vqt in enumerate(self.id_qt):
                if self.qteven :
                    if iqt!=self.len_qt: continue
                    lastId=self.len_qt
                else : #N odd
                    if iqt!=(self.len_qt-1): continue
                    lastId=self.len_qt-1
                for iy,vy in enumerate(self.id_y):
                    itot = iqt*(len(self.id_y))+iy
                    constr = 0.0
                    for k in range(0, lastId): 
                        # constr += np.power(-1., self.len_qt+k+1)*p[ k*(len(self.id_y))+iy]
                        if k%2!=0 : continue
                        constr += np.power(-1.,0.5*(k-lastId)-1) *p[ k*(len(self.id_y))+iy]
                    p1[itot] = constr
            p = p1

        if self.limy:
            p1 = copy.deepcopy(p)
            for iy,vy in enumerate(self.id_y):
                if self.yeven :
                    if iy!=self.len_y: continue
                    lastId=self.len_y
                else : #N odd
                    if iy!=(self.len_y-1): continue
                    lastId=self.len_y-1
                for iqt,vqt in enumerate(self.id_qt):
                    itot = iqt*(len(self.id_y))+iy
                    constr = 0.0
                    for k in range(0, lastId): 
                        # constr += np.power(-1., self.len_y+k+1)*p[ vqt*(len(self.id_y))+k]
                        if k%2!=0 : continue #odd terms auto-cancel at y=0.
                        constr += np.power(-1.,0.5*(k-lastId)-1) *p[ vqt*(len(self.id_y))+k]
                    p1[itot] = constr
            p = p1

        for iqt,vqt in enumerate(self.id_qt):
            for iy,vy in enumerate(self.id_y):
                itot = iqt*(len(self.id_y))+iy
                val += p[itot]*scipy.special.eval_chebyt(vqt, qt)*scipy.special.eval_chebyt(vy, y) 
        return val
    
    
    
class fit_func:

    def __init__(self, id_y,id_qt, constraints='qt->0,y-even'):
        self.id_qt = id_qt
        self.len_qt = int(len(id_qt)-1)
        self.id_y = id_y
        self.len_y = int(len(id_y)-1)
        self.nparams =  len(id_qt)*len(id_y)
        self.limqt = 'qt->0' in constraints
        self.limy = 'y->0' in constraints
        
    def __call__(self, x, parameters):
        y,qt = x[0],x[1]
        val = 0.0
        p = [0.]*self.nparams                
        
        for iqt,vqt in enumerate(self.id_qt):
            for iy,vy in enumerate(self.id_y):
                itot = iqt*(len(self.id_y))+iy
                p[itot] = parameters[itot]

        if self.limqt:
            p1 = copy.deepcopy(p)
            for iqt,vqt in enumerate(self.id_qt):
                if iqt!=self.len_qt: continue
                for iy,vy in enumerate(self.id_y):
                    itot = iqt*(len(self.id_y))+iy
                    constr = 0.0
                    for k in range(0, self.len_qt): 
                        constr += np.power(-1., self.len_qt+k+1)*p[ k*(len(self.id_y))+iy]
                    p1[itot] = constr
            p = p1

        if self.limy:
            p1 = copy.deepcopy(p)
            for iy,vy in enumerate(self.id_y):
                if iy!=self.len_y: continue
                for iqt,vqt in enumerate(self.id_qt):
                    itot = iqt*(len(self.id_y))+iy
                    constr = 0.0
                    for k in range(0, self.len_y): 
                        constr += np.power(-1., self.len_y+k+1)*p[ vqt*(len(self.id_y))+k]
                    p1[itot] = constr
            p = p1

        for iqt,vqt in enumerate(self.id_qt):
            for iy,vy in enumerate(self.id_y):
                itot = iqt*(len(self.id_y))+iy
                val += p[itot]*scipy.special.eval_chebyt(vqt, qt)*scipy.special.eval_chebyt(vy, y) 
        return val






def fit_coeff(fname=["genInput_Wplus.root","genInput_Wminus.root"] , charge="WtoMuP", coeff="A4", constraint='y->0', qt_max = 24., y_max = 2.5, degreeY=4,degreeQt=3,Map01=False): #fname="WJets_LHE.root"

    print "Coeff=",coeff, "charge=",charge," (ydeg, qtdeg)= ",degreeY,degreeQt
    # -------------- prepare the histos for the fit ----------------------- #
    if charge == 'WtoMuP' : fnameUsed = fname[0]
    if charge == 'WtoMuN' : fnameUsed = fname[1]
    inputFile   = ROOT.TFile.Open(fnameUsed)
    hA    = ROOT.gDirectory.Get('angularCoefficients'+SYST_kind+'/harmonics'+coeff+SYST_name)
        
    #settings for the output 
    hA.SetTitle('coeff_'+charge+'_'+coeff)
    hA.SetName('coeff_'+charge+'_'+coeff)
    
    hA.GetXaxis().SetTitle('Y')
    hA.GetYaxis().SetTitle('q_{T} [GeV]')
    
    y_max_bin  = hA.GetXaxis().FindBin( y_max+0.0001 )
    y_min_bin  = hA.GetXaxis().FindBin( 0.0 )
    y_max_cut  = hA.GetXaxis().GetBinLowEdge(y_max_bin+1)        
    y_min_cut  = 0.0   #hMC.GetXaxis().GetBinCenter(y_min_bin)        

    qt_max_bin = hA.GetYaxis().FindBin( qt_max )
    qt_min_bin = hA.GetYaxis().FindBin( 0.0 )
    qt_max_cut = hA.GetYaxis().GetBinLowEdge(qt_max_bin+1)        
    qt_min_cut = 0.0    #hMC.GetYaxis().GetBinCenter(qt_min_bin)        
     
    x,y = [], []
    if Map01 :
        for i in range(y_min_bin, y_max_bin+1):   x.append( (hA.GetXaxis().GetBinCenter(i)- y_min_cut)/(y_max_cut-  y_min_cut)) #range between 0,1
        for i in range(qt_min_bin, qt_max_bin+1): y.append( (hA.GetYaxis().GetBinCenter(i)-qt_min_cut)/(qt_max_cut-qt_min_cut)) #range between 0,1
    else :
        for i in range(y_min_bin, y_max_bin+1):   x.append( 2*(hA.GetXaxis().GetBinCenter(i)- y_min_cut)/(y_max_cut-  y_min_cut) - 1.) #range between -1,1
        for i in range(qt_min_bin, qt_max_bin+1): y.append( 2*(hA.GetYaxis().GetBinCenter(i)-qt_min_cut)/(qt_max_cut-qt_min_cut) - 1.) #range between -1,1 

    xx = array('f', x)
    yy = array('f', y)
    
    hA_cut = ROOT.TH2F("hfit_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt), "hfit_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt), len(xx)-1, xx, len(yy)-1, yy)
    hA_cut.GetXaxis().SetTitle('Y, 1='+str(y_max))
    hA_cut.GetYaxis().SetTitle('q_{T} [GeV], 1='+str(qt_max))
    
    for i in range(1, hA_cut.GetXaxis().GetNbins()+1):
        for j in range(1, hA_cut.GetYaxis().GetNbins()+1):
            hA_cut.SetBinContent(i,j, hA.GetBinContent(i,j) )
            hA_cut.SetBinError(i,j, hA.GetBinError(i,j) )
            

    id_y  = [] #poly in y
    id_qt = [] #poly in qt
    for yy in range(0,degreeY) :
        id_y.append(yy)
    for qq in range(0,degreeQt) :
        id_qt.append(qq)    
    

    if Map01 :
        func = fit_func01(id_y,id_qt, constraint)
    else :
        func = fit_func(id_y,id_qt, constraint)
    nparams = len(id_qt)*len(id_y)
    
    if Map01 :
        fit = ROOT.TF2("fit_"+charge+"_"+coeff, func, 0, 1.0, 0, 1.0, nparams)
        # fit = ROOT.TF2("fit_"+charge+"_"+coeff, func, 0, 1.0, 0.2, 1.0, nparams) #exclude first bins in pt to avoid negative values, doesn't work
    else :
        fit = ROOT.TF2("fit_"+charge+"_"+coeff, func, -1, 1.0, -1, 1.0, nparams)    
        # fit = ROOT.TF2("fit_"+charge+"_"+coeff, func, -1, 1.0, -0.6, 1.0, nparams) #exclude first bins in pt to avoid negative values, doesn't work
    
    
    for iqt,vqt in enumerate(id_qt):
        for iy,vy in enumerate(id_y):
            itot = iqt*(len(id_y))+iy
            fit.SetParName(itot, "T"+str(iqt)+"_"+"T"+str(iy))
            fit.SetParError(itot, 0.002)

    if 'qt->0' in constraint:
        for iqt,vqt in enumerate(id_qt):
            if iqt!=int(len(id_qt)-1): continue
            for iy,vy in enumerate(id_y):
                itot = iqt*(len(id_y))+iy
                fit.SetParName(itot, "T"+str(iqt)+"_"+"T"+str(iy))
                fit.SetParError(itot, 0.002)
                # print "qt->0: Fix parameter", itot
                fit.FixParameter(itot, 0.0)

    if 'y->0' in constraint:
        for iy,vy in enumerate(id_y):
            if iy!=int(len(id_y)-1): continue
            for iqt,vqt in enumerate(id_qt):
                itot = iqt*(len(id_y))+iy
                fit.SetParName(itot, "T"+str(iqt)+"_"+"T"+str(iy))
                fit.SetParError(itot, 0.002)
                # print "y->0: Fix parameter", itot
                fit.FixParameter(itot, 0.0)
    
    res = hA_cut.Fit(fit, "RSQ") 

    fit.SetNpx(5000)
    
    #build the errorbands for the fit
    h_fitError = hA_cut.Clone("hfitErr_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt))
    ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_fitError,0.68)
    
    

    #DEBUG block------------------------------------#
    # print "-------------------- Post fit debug:"
    # res.Print()
    print coeff, "(qt,y)=", degreeQt, degreeY, charge, ", NDF=",fit.GetNDF()
    for i in range(nparams) :
        print "par", i, fit.GetParName(i), "=",fit.GetParameter(i), "+/-",fit.GetParError(i), 
        if fit.GetParameter(i)!=0 : print ", deltaPar/par=",  fit.GetParError(i)/fit.GetParameter(i)
        else : print ""
    # print res.GetCorrelationMatrix()
    # print res.GetCovarianceMatrix()
    print "chi2 reduced:", fit.GetChisquare()/fit.GetNDF(), "+/-", math.sqrt(2*fit.GetNDF())/fit.GetNDF()
    print "-------------------- "
    #------------------------------------------------#

    #addittional debug ---------------#
    # for y in [-1.0, 0.0, +1.0]:
    #     print 'f(qt=0, y=', y, ')=', fit.Eval(y, -1.0 )
    # for qt in [-1.0, 0.0, +1.0]:
    #     print 'f(qt=', qt, ', y=0)=', fit.Eval(-1.0, qt )
    #--------------------------------#
    
    
    #-------------------------pulls building ------------------------------#
    
    hPulls2D = hA_cut.Clone("hPulls2D_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt))
    hPulls2D.SetTitle("hPulls2D_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt))
    hPulls2D.SetName("hPulls2D_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt))
    hPulls2D.Reset()
    for vy in range(1, hPulls2D.GetXaxis().GetNbins()+1):
        for iq in range(1, hPulls2D.GetYaxis().GetNbins()+1):
            hPulls2D.SetBinContent( vy, iq, (fit.Eval( hPulls2D.GetXaxis().GetBinCenter(vy), hPulls2D.GetYaxis().GetBinCenter(iq) ) - hA_cut.GetBinContent(vy,iq))/hA_cut.GetBinError(vy,iq) )
    hPulls2D.SetMinimum(-8)
    hPulls2D.SetMaximum(+8)

    hPulls1D = ROOT.TH1F("hPulls1D_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt), "hPulls1D_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt), 25,-4,4)
    for vy in range(1, hPulls2D.GetXaxis().GetNbins()+1):
        for iq in range(1, hPulls2D.GetYaxis().GetNbins()+1):
            hPulls1D.Fill( hPulls2D.GetBinContent( vy, iq) )
    # hPulls1D.Sumw2()

    c_pull = ROOT.TCanvas("c_pull_"+charge+'_'+coeff+'_'+str(degreeY)+str(degreeQt), "c_pull_"+charge+'_'+coeff+'_'+str(degreeY)+str(degreeQt), 1200, 600)
    c_pull.cd()
    c_pull.Divide(2,1)

    c_pull.cd(1)
    hPulls2D.SetStats(0)
    hPulls2D.DrawCopy("colz")
    c_pull.cd(2)
    hPulls1D.Fit("gaus", "Q", "", -4,4)
    gaus = hPulls1D.GetFunction("gaus")
    hPulls1D.DrawCopy("hist")
    # if hPulls1D.GetBinContent(0)+hPulls1D.GetBinContent(hPulls1D.GetNbinsX()+1)>=hPulls1D.GetEntries() : #there are entry instide the histogram
    try :
        gaus.Draw("same")
    except :
        print "empty pulls"
    # hPulls1D.SetStats()
    # stats = hPulls1D.FindObject("stats")
    ROOT.gStyle.SetOptFit()
    

    c_pull.Update()
    
    c_scratch = ROOT.TCanvas("c_scratch"+charge+coeff+str(degreeY)+str(degreeQt), "c_scratch"+charge+coeff+str(degreeY)+str(degreeQt), 1800, 600)
    c_scratch.cd()
    # c_pull.cd()
    # c_pull.Draw()

    
    #---------------prepare the output ----------------#
    outDict ={}
    outDict['hh'+charge+coeff] = hA
    outDict['fit'+charge+coeff] = hA_cut
    outDict['fit'+'Res'+charge+coeff] = res
    outDict['fit'+'Err'+charge+coeff] = h_fitError
    outDict['pull'+'c'+charge+coeff] = c_pull
    outDict['pull'+'1D'+charge+coeff] = hPulls1D
    outDict['pull'+'2D'+charge+coeff] = hPulls2D
    
    return outDict



class fit_func_unpol:

    def __init__(self, id_y):
        self.id_yy= id_y
        self.nparamss =  len(id_y)
       
    def __call__(self, x, parameters):
        y = x[0]
        val = 0.0
        p = [0.]*self.nparamss    
        for iy,vy in enumerate(self.id_yy):
            if vy%2!=0 : continue 
            val+= parameters[iy]*scipy.special.eval_chebyt(vy, y) 
        return val

def fit_unpol(fname=["genInput_Wplus.root","genInput_Wminus.root"], charge="WtoMuP", qt_max = 24., y_max = 2.5, degreeY=4,Map01=True,coeff='AUL') : #fname="WJets_LHE.root
    
    print "unpol cross sec, charge=", charge, " ydeg=, ",degreeY
    
    if charge == 'WtoMuP' : fnameUsed = fname[0]
    if charge == 'WtoMuN' : fnameUsed = fname[1]
    inputFile   = ROOT.TFile.Open(fnameUsed)
    hMC   = ROOT.gDirectory.Get('angularCoefficients'+SYST_kind+'/harmonics'+coeff+SYST_name)


    hMC.SetTitle('MC_'+charge+'_unpol')
    hMC.SetName('MC_'+charge+'_unpol')
    hMC.GetXaxis().SetTitle('Y')
    hMC.GetYaxis().SetTitle('q_{T} [GeV]')
    
    #------ rebin input histo ------ #
    y_max_bin  = hMC.GetXaxis().FindBin( y_max+0.0001 )
    y_min_bin  = hMC.GetXaxis().FindBin( 0.0 )
    y_max_cut  = hMC.GetXaxis().GetBinLowEdge(y_max_bin+1)        
    y_min_cut  = 0.0   #hMC.GetXaxis().GetBinCenter(y_min_bin)        

    qt_max_bin = hMC.GetYaxis().FindBin( qt_max )
    qt_min_bin = hMC.GetYaxis().FindBin( 0.0 )
    qt_max_cut = hMC.GetYaxis().GetBinLowEdge(qt_max_bin+1)        
    qt_min_cut = 0.0    #hMC.GetYaxis().GetBinCenter(qt_min_bin)
    
    x,y = [], []
    if Map01 :
        for i in range(y_min_bin, y_max_bin+1):   x.append( (hMC.GetXaxis().GetBinCenter(i)- y_min_cut)/(y_max_cut-  y_min_cut)) #range between 0,1
        for i in range(qt_min_bin, qt_max_bin+1): y.append( (hMC.GetYaxis().GetBinCenter(i)-qt_min_cut)/(qt_max_cut-qt_min_cut)) #range between 0,1
    else :
        for i in range(y_min_bin, y_max_bin+1):   x.append( 2*(hMC.GetXaxis().GetBinCenter(i)- y_min_cut)/(y_max_cut-  y_min_cut) - 1.) #range between -1,1
        for i in range(qt_min_bin, qt_max_bin+1): y.append( 2*(hMC.GetYaxis().GetBinCenter(i)-qt_min_cut)/(qt_max_cut-qt_min_cut) - 1.) #range between -1,1 

    xx = array('f', x)
    yy = array('f', y)
    
    hMC_rebin = ROOT.TH2F("hfit_"+charge+'_unpol_'+str(degreeY), "hfit_"+charge+'_unpol_'+str(degreeY), len(xx)-1, xx, len(yy)-1, yy)
    hMC_rebin.GetXaxis().SetTitle('Y, 1='+str(y_max))
    hMC_rebin.GetYaxis().SetTitle('q_{T} [GeV], 1='+str(qt_max))
    for i in range(1, hMC_rebin.GetXaxis().GetNbins()+1):
        for j in range(1, hMC_rebin.GetYaxis().GetNbins()+1):
            hMC_rebin.SetBinContent(i,j, hMC.GetBinContent(i,j))
            hMC_rebin.SetBinError(i,j, hMC.GetBinError(i,j))
    
    #----- slice the histogram and define fit functions-----#
    hList = []
    funcList = []
    resList = []
    fitErrList = []
    
    id_y  = [] #poly in y
    for yy in range(0,degreeY) :
        id_y.append(yy)
    nparams = len(id_y)
    
    for i in range(1, hMC_rebin.GetYaxis().GetNbins()+1):
        hList.append(hMC_rebin.ProjectionX(hMC_rebin.GetName()+'_'+str(i),i,i))
        func = fit_func_unpol(id_y)
        if Map01 :
            funcList.append(ROOT.TF1("fit_"+charge+"_unpol", func, 0, 1.0, nparams))    
        else :
            funcList.append(ROOT.TF1("fit_"+charge+"_unpol", func, -1, 1.0, nparams))
        for yy in range(0,degreeY) :
            funcList[-1].SetParName(yy,'T'+str(yy))
            if yy%2!=0 :     
                funcList[-1].FixParameter(yy,0.0)
    #---fit----#
    for h,f in zip(hList,funcList) :
        resList.append(h.Fit(f, "RSQ"))
        fitErrList.append(h.Clone("hfitErr_"+charge+'_unpol_'+str(degreeY)+'_'+str(hList.index(h)+1)))
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(fitErrList[-1],0.68)

    # fit.SetNpx(500)    
    
    #-------------------------pulls building ------------------------------#
    
    hPulls2D = hMC_rebin.Clone("hPulls2D_"+charge+"_"+'unpol'+'_'+str(degreeY))
    hPulls2D.SetTitle("hPulls2D_"+charge+"_"+'unpol'+'_'+str(degreeY))
    hPulls2D.SetName("hPulls2D_"+charge+"_"+'unpol'+'_'+str(degreeY))
    hPulls2D.Reset()
    for vy in range(1, hPulls2D.GetXaxis().GetNbins()+1):
        for iq in range(1, hPulls2D.GetYaxis().GetNbins()+1):
            hPulls2D.SetBinContent( vy, iq, (funcList[iq-1].Eval( hPulls2D.GetXaxis().GetBinCenter(vy)) - hMC_rebin.GetBinContent(vy,iq))/hMC_rebin.GetBinError(vy,iq) )
    hPulls2D.SetMinimum(-4)
    hPulls2D.SetMaximum(+4)

    hPulls1D = ROOT.TH1F("hPulls1D_"+charge+"_"+'unpol'+'_'+str(degreeY), "hPulls1D_"+charge+"_"+'unpol'+'_'+str(degreeY), 25,-4,4)
    for vy in range(1, hPulls2D.GetXaxis().GetNbins()+1):
        for iq in range(1, hPulls2D.GetYaxis().GetNbins()+1):
            hPulls1D.Fill( hPulls2D.GetBinContent( vy, iq) )
    # hPulls1D.Sumw2()

    c_pull = ROOT.TCanvas("c_pull_"+charge+'_'+'unpol'+'_'+str(degreeY), "c_pull_"+charge+'_'+'unpol'+'_'+str(degreeY), 1200, 600)
    c_pull.cd()
    c_pull.Divide(2,1)

    c_pull.cd(1)
    hPulls2D.SetStats(0)
    hPulls2D.DrawCopy("colz")
    c_pull.cd(2)
    hPulls1D.Fit("gaus", "Q", "", -4,4)
    gaus = hPulls1D.GetFunction("gaus")
    hPulls1D.DrawCopy("hist")
    gaus.Draw("same")
    # hPulls1D.SetStats()
    # stats = hPulls1D.FindObject("stats")
    ROOT.gStyle.SetOptFit()
    

    c_pull.Update()
    
    # c_scratch = ROOT.TCanvas("c_scratch"+charge+'unpol'+str(degreeY), "c_scratch"+charge+'unpol'+str(degreeY), 1800, 600)
    # c_scratch.cd()
    
    c_scratchBIS = ROOT.TCanvas("c_scratchBIS"+charge+'unpol'+str(degreeY), "c_scratchBIS"+charge+'unpol'+str(degreeY), 1800, 600)
    c_scratchBIS.cd()
    
    
    
    #--- prepare the output ---#
    outDict = {}
    # for h,r in zip(hList,resList) :
    for i in range(1, hMC_rebin.GetYaxis().GetNbins()+1): 
        outDict["fit"+charge+'unpol'+str(i)] = hList[i-1]
        outDict["fitRes"+charge+'unpol'+str(i)] = resList[i-1]
        outDict["fitErr"+charge+'unpol'+str(i)] = fitErrList[i-1]
    outDict['pull'+'c'+charge+'unpol'] = c_pull
    outDict['pull'+'1D'+charge+'unpol'] = hPulls1D
    outDict['pull'+'2D'+charge+'unpol'] = hPulls2D
    outDict['hh'+charge+'unpol'+'MC'] = hMC_rebin


    
    return outDict
    

















class fit_poly:

    def __init__(self, id_y,id_qt, constraints='qt->0'):
        self.id_qt = id_qt
        self.id_y = id_y
        self.nparams =  len(id_qt)*len(id_y)
        self.limqt = 'qt->0' in constraints
        self.limy = 'y->0' in constraints
        self.limC1 = 'C1' in constraints

    def __call__(self, x, parameters):
        y,qt = x[0],x[1]
        val = 0.0
        p = [0.]*self.nparams                
        
        for iqt,vqt in enumerate(self.id_qt):
            for iy,vy in enumerate(self.id_y):
                itot = iqt*(len(self.id_y))+iy
                p[itot] = parameters[itot]
                if self.limy and vy==0 : 
                    p[itot]=0
                    # print "itot=", itot
                if self.limqt and vqt==0 : 
                    p[itot]=0
                    # print "itot=",itot
                if self.limC1 and vqt==1 :
                    p[itot]=0
                # if self.limy and vy%2==0 : 
                #     p[itot]=0
                # if ((not self.limy) and vy%2!=0 ): 
                #     p[itot]=0
                val += p[itot]*pow(qt,vqt)*pow(y,vy) 
        return val



def fit_coeff_poly(fname=["genInput_Wplus.root","genInput_Wminus.root"] , charge="WtoMuP", coeff="A4", constraint='y->0', qt_max = 24., y_max = 2.5, degreeY=4,degreeQt=3,Map01=False): #fname="WJets_LHE.root"

    print "Coeff=",coeff, "charge=",charge," (ydeg, qtdeg)= ",degreeY,degreeQt
    # -------------- prepare the histos for the fit ----------------------- #
    if charge == 'WtoMuP' : fnameUsed = fname[0]
    if charge == 'WtoMuN' : fnameUsed = fname[1]
    inputFile   = ROOT.TFile.Open(fnameUsed)
    hA    = ROOT.gDirectory.Get('angularCoefficients'+SYST_kind+'/harmonics'+coeff+SYST_name)
 
    #settings for the output 
    hA.SetTitle('coeff_'+charge+'_'+coeff)
    hA.SetName('coeff_'+charge+'_'+coeff)
    hA.GetXaxis().SetTitle('Y')
    hA.GetYaxis().SetTitle('q_{T} [GeV]')
    
    y_max_bin  = hA.GetXaxis().FindBin( y_max+0.0001 )
    qt_max_bin = hA.GetYaxis().FindBin( qt_max )
        
    x,y = [], []
    for i in range(1, y_max_bin+1): x.append( hA.GetXaxis().GetBinLowEdge(i))
    for i in range(1, qt_max_bin+1): y.append( hA.GetYaxis().GetBinLowEdge(i))    

    xx = array('f', x)
    yy = array('f', y)
    print "bin y=",xx
    print "bin qt=", yy
    
    hA_cut = ROOT.TH2F("hfit_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt), "hfit_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt), len(xx)-1, xx, len(yy)-1, yy)
    hA_cut.GetXaxis().SetTitle('Y, 1='+str(y_max))
    hA_cut.GetYaxis().SetTitle('q_{T} [GeV], 1='+str(qt_max))
    
    for i in range(1, hA_cut.GetXaxis().GetNbins()+1):
        for j in range(1, hA_cut.GetYaxis().GetNbins()+1):
            hA_cut.SetBinContent(i,j, hA.GetBinContent(i,j) )
            hA_cut.SetBinError(i,j, hA.GetBinError(i,j) )

    id_y  = [] #poly in y
    id_qt = [] #poly in qt
    for yy in range(0,degreeY) :
        id_y.append(yy)
    for qq in range(0,degreeQt) :
        id_qt.append(qq)    
    
    func = fit_poly(id_y,id_qt, constraint)
    nparams = len(id_qt)*len(id_y)
    
    fit = ROOT.TF2("fit_"+charge+"_"+coeff, func, 0, y_max, 0, qt_max, nparams)
    
    for iqt,vqt in enumerate(id_qt):
        for iy,vy in enumerate(id_y):
            itot = iqt*(len(id_y))+iy
            fit.SetParName(itot, "qt"+str(iqt)+"_"+"y"+str(iy))
            fit.SetParError(itot, 0.002)
            if 'qt->0' in constraint and iqt==0 : 
                fit.FixParameter(itot, 0.0)
                print "fixed:", fit.GetParName(itot)
            if 'y->0' in constraint and iy==0 : 
                fit.FixParameter(itot, 0.0)
                print "fixed:", fit.GetParName(itot)
            if 'C1' in constraint and vqt==1 :
                fit.FixParameter(itot, 0.0)
                print "fixed:", fit.GetParName(itot)
            # if 'y->0' in constraint and vy%2==0 : #y=0-->dispari
            #     fit.FixParameter(itot, 0.0)
            #     print "fixed:", fit.GetParName(itot)
            # if ((not 'y->0' in constraint) and vy%2!=0) : #y!=0-->pari
            #     fit.FixParameter(itot, 0.0)
            #     print "fixed:", fit.GetParName(itot)    

    # res = hA_cut.Fit(fit, "RSV") #range, result, verbose
    res = hA_cut.Fit(fit, "RSQ") 
    fit.SetNpx(5000)
    
    #build the errorbands for the fit
    h_fitError = hA_cut.Clone("hfitErr_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt))
    ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_fitError,0.68)


    #DEBUG block------------------------------------#
    # print "-------------------- Post fit debug:"
    # res.Print()
    print coeff, "(qt,y)=", degreeQt, degreeY, charge, ", NDF=",fit.GetNDF()
    for i in range(nparams) :
        print "par", i, fit.GetParName(i), "=",fit.GetParameter(i), "+/-",fit.GetParError(i), 
        if fit.GetParameter(i)!=0 : print ", deltaPar/par=",  fit.GetParError(i)/fit.GetParameter(i)
        else : print ""
    print "chi2 reduced:", fit.GetChisquare()/fit.GetNDF(), "+/-", math.sqrt(2*fit.GetNDF())/fit.GetNDF()
    print "-------------------- "
    #------------------------------------------------#

    #addittional debug ---------------#
    # for y in [-1.0, 0.0, +1.0]:
    #     print 'f(qt=0, y=', y, ')=', fit.Eval(y, -1.0 )
    # for qt in [-1.0, 0.0, +1.0]:
    #     print 'f(qt=', qt, ', y=0)=', fit.Eval(-1.0, qt )
    #--------------------------------#
    
    
    #-------------------------pulls building ------------------------------#
    
    hPulls2D = hA_cut.Clone("hPulls2D_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt))
    hPulls2D.SetTitle("hPulls2D_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt))
    hPulls2D.SetName("hPulls2D_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt))
    hPulls2D.Reset()
    for vy in range(1, hPulls2D.GetXaxis().GetNbins()+1):
        for iq in range(1, hPulls2D.GetYaxis().GetNbins()+1):
            hPulls2D.SetBinContent( vy, iq, (fit.Eval( hPulls2D.GetXaxis().GetBinCenter(vy), hPulls2D.GetYaxis().GetBinCenter(iq) ) - hA_cut.GetBinContent(vy,iq))/hA_cut.GetBinError(vy,iq) )
    hPulls2D.SetMinimum(-8)
    hPulls2D.SetMaximum(+8)

    hPulls1D = ROOT.TH1F("hPulls1D_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt), "hPulls1D_"+charge+"_"+coeff+'_'+str(degreeY)+str(degreeQt), 25,-4,4)
    for vy in range(1, hPulls2D.GetXaxis().GetNbins()+1):
        for iq in range(1, hPulls2D.GetYaxis().GetNbins()+1):
            hPulls1D.Fill( hPulls2D.GetBinContent( vy, iq) )
    # hPulls1D.Sumw2()

    c_pull = ROOT.TCanvas("c_pull_"+charge+'_'+coeff+'_'+str(degreeY)+str(degreeQt), "c_pull_"+charge+'_'+coeff+'_'+str(degreeY)+str(degreeQt), 1200, 600)
    c_pull.cd()
    c_pull.Divide(2,1)

    c_pull.cd(1)
    hPulls2D.SetStats(0)
    hPulls2D.DrawCopy("colz")
    c_pull.cd(2)
    hPulls1D.Fit("gaus", "Q", "", -4,4)
    gaus = hPulls1D.GetFunction("gaus")
    hPulls1D.DrawCopy("hist")
    # if hPulls1D.GetBinContent(0)+hPulls1D.GetBinContent(hPulls1D.GetNbinsX()+1)>=hPulls1D.GetEntries() : #there are entry instide the histogram
    try :
        gaus.Draw("same")
    except :
        print "empty pulls"
    # hPulls1D.SetStats()
    # stats = hPulls1D.FindObject("stats")
    ROOT.gStyle.SetOptFit()
    
    c_pull.Update()
    
    c_scratch = ROOT.TCanvas("c_scratch"+charge+coeff+str(degreeY)+str(degreeQt), "c_scratch"+charge+coeff+str(degreeY)+str(degreeQt), 1800, 600)
    c_scratch.cd()
    # c_pull.cd()
    # c_pull.Draw()

    
    #---------------prepare the output ----------------#
    outDict ={}
    outDict['hh'+charge+coeff] = hA
    outDict['fit'+charge+coeff] = hA_cut
    outDict['fit'+'Res'+charge+coeff] = res
    outDict['fit'+'Err'+charge+coeff] = h_fitError
    outDict['pull'+'c'+charge+coeff] = c_pull
    outDict['pull'+'1D'+charge+coeff] = hPulls1D
    outDict['pull'+'2D'+charge+coeff] = hPulls2D
    return outDict




class fit_func_unpol_poly:

    def __init__(self, id_y):
        self.id_yy= id_y
        self.nparamss =  len(id_y)
       
    def __call__(self, x, parameters):
        y = x[0]
        val = 0.0
        p = [0.]*self.nparamss    
        for iy,vy in enumerate(self.id_yy):
            # if vy%2!=0 : continue 
            val += parameters[iy]*pow(y,vy)
        return val

def fit_unpol_poly(fname=["genInput_Wplus.root","genInput_Wminus.root"], charge="WtoMuP", qt_max = 24., y_max = 2.5, degreeY=4,Map01=True,coeff='AUL') : #fname="WJets_LHE.root
    
    print "unpol cross sec, charge=", charge, " ydeg=, ",degreeY
    
    if charge == 'WtoMuP' : fnameUsed = fname[0]
    if charge == 'WtoMuN' : fnameUsed = fname[1]
    inputFile   = ROOT.TFile.Open(fnameUsed)
    hMC   = ROOT.gDirectory.Get('angularCoefficients'+SYST_kind+'/harmonics'+coeff+SYST_name)

    hMC.SetTitle('MC_'+charge+'_unpol')
    hMC.SetName('MC_'+charge+'_unpol')
    hMC.GetXaxis().SetTitle('Y')
    hMC.GetYaxis().SetTitle('q_{T} [GeV]')
    
    #------ rebin input histo ------ #
    y_max_bin  = hMC.GetXaxis().FindBin( y_max+0.0001 )
    qt_max_bin = hMC.GetYaxis().FindBin( qt_max )   
    
    x,y = [], []
    for i in range(1, y_max_bin+1): x.append( hMC.GetXaxis().GetBinLowEdge(i))
    for i in range(1, qt_max_bin+1): y.append( hMC.GetYaxis().GetBinLowEdge(i))    

    xx = array('f', x)
    yy = array('f', y)
    
    hMC_rebin = ROOT.TH2F("hfit_"+charge+'_unpol_'+str(degreeY), "hfit_"+charge+'_unpol_'+str(degreeY), len(xx)-1, xx, len(yy)-1, yy)
    hMC_rebin.GetXaxis().SetTitle('Y, 1='+str(y_max))
    hMC_rebin.GetYaxis().SetTitle('q_{T} [GeV], 1='+str(qt_max))
    for i in range(1, hMC_rebin.GetXaxis().GetNbins()+1):
        for j in range(1, hMC_rebin.GetYaxis().GetNbins()+1):
            hMC_rebin.SetBinContent(i,j, hMC.GetBinContent(i,j))
            hMC_rebin.SetBinError(i,j, hMC.GetBinError(i,j))
    
    #----- slice the histogram and define fit functions-----#
    hList = []
    funcList = []
    resList = []
    fitErrList = []
    
    id_y  = [] #poly in y
    for yy in range(0,degreeY) :
        id_y.append(yy)
    nparams = len(id_y)
    
    for i in range(1, hMC_rebin.GetYaxis().GetNbins()+1):
        hList.append(hMC_rebin.ProjectionX(hMC_rebin.GetName()+'_'+str(i),i,i))
        func = fit_func_unpol_poly(id_y)
        funcList.append(ROOT.TF1("fit_"+charge+"_unpol", func, 0, y_max, nparams))    
        for yy in range(0,degreeY) :
            funcList[-1].SetParName(yy,'y'+str(yy))
            # if yy%2!=0 :     
            #     funcList[-1].FixParameter(yy,0.0)
    
    #---fit----#
    for h,f in zip(hList,funcList) :
        resList.append(h.Fit(f, "RSQ"))
        fitErrList.append(h.Clone("hfitErr_"+charge+'_unpol_'+str(degreeY)+'_'+str(hList.index(h)+1)))
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(fitErrList[-1],0.68)
    
    #-------------------------pulls building ------------------------------#
    
    hPulls2D = hMC_rebin.Clone("hPulls2D_"+charge+"_"+'unpol'+'_'+str(degreeY))
    hPulls2D.SetTitle("hPulls2D_"+charge+"_"+'unpol'+'_'+str(degreeY))
    hPulls2D.SetName("hPulls2D_"+charge+"_"+'unpol'+'_'+str(degreeY))
    hPulls2D.Reset()
    for vy in range(1, hPulls2D.GetXaxis().GetNbins()+1):
        for iq in range(1, hPulls2D.GetYaxis().GetNbins()+1):
            hPulls2D.SetBinContent( vy, iq, (funcList[iq-1].Eval( hPulls2D.GetXaxis().GetBinCenter(vy)) - hMC_rebin.GetBinContent(vy,iq))/hMC_rebin.GetBinError(vy,iq) )
    hPulls2D.SetMinimum(-4)
    hPulls2D.SetMaximum(+4)

    hPulls1D = ROOT.TH1F("hPulls1D_"+charge+"_"+'unpol'+'_'+str(degreeY), "hPulls1D_"+charge+"_"+'unpol'+'_'+str(degreeY), 25,-4,4)
    for vy in range(1, hPulls2D.GetXaxis().GetNbins()+1):
        for iq in range(1, hPulls2D.GetYaxis().GetNbins()+1):
            hPulls1D.Fill( hPulls2D.GetBinContent( vy, iq) )

    c_pull = ROOT.TCanvas("c_pull_"+charge+'_'+'unpol'+'_'+str(degreeY), "c_pull_"+charge+'_'+'unpol'+'_'+str(degreeY), 1200, 600)
    c_pull.cd()
    c_pull.Divide(2,1)

    c_pull.cd(1)
    hPulls2D.SetStats(0)
    hPulls2D.DrawCopy("colz")
    c_pull.cd(2)
    hPulls1D.Fit("gaus", "Q", "", -4,4)
    gaus = hPulls1D.GetFunction("gaus")
    hPulls1D.DrawCopy("hist")
    gaus.Draw("same")
    ROOT.gStyle.SetOptFit()
    c_pull.Update()
    
    
    c_scratchBIS = ROOT.TCanvas("c_scratchBIS"+charge+'unpol'+str(degreeY), "c_scratchBIS"+charge+'unpol'+str(degreeY), 1800, 600)
    c_scratchBIS.cd()
    
    #--- prepare the output ---#
    outDict = {}
    for i in range(1, hMC_rebin.GetYaxis().GetNbins()+1): 
        outDict["fit"+charge+'unpol'+str(i)] = hList[i-1]
        outDict["fitRes"+charge+'unpol'+str(i)] = resList[i-1]
        outDict["fitErr"+charge+'unpol'+str(i)] = fitErrList[i-1]
    outDict['pull'+'c'+charge+'unpol'] = c_pull
    outDict['pull'+'1D'+charge+'unpol'] = hPulls1D
    outDict['pull'+'2D'+charge+'unpol'] = hPulls2D
    outDict['hh'+charge+'unpol'+'MC'] = hMC_rebin

    return outDict
       
            
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
            
    
def rebuildHistoDict(charge="WtoMuP", coeff="A4", degreeY=4,degreeQt=3):
    inputFile   = ROOT.TFile.Open(output_dir+'/'+IN_REG+'.root')
    
    path = 'Y'+str(degreeY)+'_Qt'+str(degreeQt)+'/'
  
    nameDict = { #name inside : [prefix, suffix,degSuffixFlag]
        # 'MC_' : ['hh','MC',False],
        'hPulls2D_' : ['pull2D','',True], #True
        'hPulls1D_' : ['pull1D','',True], #add _
        'hfit_' : ['fit','',True],
        'hfitErr_' :['fitErr','',True],
        'coeff_' : ['hh','',False],
        # 'coeffErr_' : ['hh','Err',False],
        'c_pull_' : ['pullc','',True]
    }
    
    outDict = {}
    
    for ind,val in nameDict.iteritems() :
        val2 = ''
        if val[2] :
            val2 = '_'+str(degreeY)+str(degreeQt)
        # if val[0] == 'fitErr':
        #     val2+='_Err'
        getstring=path+ind+charge+'_'+coeff+val2
        outDict[val[0]+charge+coeff+val[1]] = inputFile.Get(getstring)
    #fit result are missing-->rebuild the chi2 form the function
    outDict['fit'+'Res'+charge+coeff] = outDict['fit'+charge+coeff].GetFunction("fit_"+charge+"_"+coeff)
    
    return outDict

def F_test(hdict,valDict,signDict, degreeYList, degreeQtList,Npt,VALONLY,   s,coeff,yMax,qtMax,yMin,qtMin) :#return the pvalue of the test or 1 if chi2max>chi2min
        
    binYmax=degreeYList.index(yMax)+1
    binYmin=degreeYList.index(yMin)+1
    
    if coeff!='unpol' :
        binQtmax=degreeQtList.index(qtMax)+1
        binQtmin=degreeQtList.index(qtMin)+1
    
        if not (VALONLY or MULTICORE) :
            NDFMax = hdict[str(yMax)+str(qtMax)+s+coeff]['fit'+'Res'+signDict[s]+coeff].Ndf()
            NDFMin = hdict[str(yMin)+str(qtMin)+s+coeff]['fit'+'Res'+signDict[s]+coeff].Ndf()
        else :
            NDFMax = hdict[str(yMax)+str(qtMax)+s+coeff]['fit'+'Res'+signDict[s]+coeff].GetNDF()
            NDFMin = hdict[str(yMin)+str(qtMin)+s+coeff]['fit'+'Res'+signDict[s]+coeff].GetNDF()
        parMax = Npt-NDFMax
        parMin = Npt-NDFMin
        
        #evaluate chi2 in the two bins
        chi2Max =  valDict['h2_chi2'+s+coeff].GetBinContent(binYmax,binQtmax)*NDFMax
        chi2Min =  valDict['h2_chi2'+s+coeff].GetBinContent(binYmin,binQtmin)*NDFMin
    
    else : #unpol analysis
        NDFMax = hdict[str(yMax)+s+coeff]['fit'+'Res'+signDict[s]+coeff+str(qtMax)].Ndf()
        NDFMin = hdict[str(yMin)+s+coeff]['fit'+'Res'+signDict[s]+coeff+str(qtMin)].Ndf()
        parMax = Npt-NDFMax
        parMin = Npt-NDFMin
        chi2Max =  valDict['h2_chi2'+s+coeff].GetBinContent(binYmax,qtMax)*NDFMax
        chi2Min =  valDict['h2_chi2'+s+coeff].GetBinContent(binYmin,qtMin)*NDFMin
    
    try : 
        Fobs = (chi2Min-chi2Max)*(Npt-parMax)/(chi2Max*(parMax-parMin))
        if Fobs>0 :
            pvalue=1-scipy.stats.f.cdf(Fobs,1,(Npt-parMax))
        else :
            pvalue=1
    except :
        print "WARNING: issue in F-test evaluation,",s,coeff, ",parMax=parMin?", parMax, parMin
        pvalue = 1
        
    return pvalue

def findBest_path(hdict,valDict,signDict, degreeYList, degreeQtList,Npt,VALONLY, s,coeff) :# look for the best combination y,qt with a search following the best path from lower combination
    pthr=0.05
    outCell = []
    
    print "new optimal path:", s, coeff,
    
    yy = degreeYList[0]
    qq = degreeQtList[0]
    contY=True
    contQt=True
    while contY or contQt :
        if yy<degreeYList[-1] :
            yMax = yy+1
            pvalY = F_test(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=degreeQtList,Npt=Npt,VALONLY=VALONLY,s=s,coeff=coeff,yMax=yMax,qtMax=qq,yMin=yy,qtMin=qq)
        else :
            contY=False
            pvalY=1
        if qq<degreeQtList[-1] :
            qtMax = qq+1
            pvalQt = F_test(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=degreeQtList,Npt=Npt,VALONLY=VALONLY,s=s,coeff=coeff,yMax=yy,qtMax=qtMax,yMin=yy,qtMin=qq)
        else :
            contQt=False
            pvalQt=1
        if pvalY<pvalQt :
            if pvalY<pthr :
                yy+=1
            else: 
                contY=False
                contQt=False
        elif pvalQt<=pvalY:
            if pvalQt<pthr :
                qq+=1
            else: 
                contY=False
                contQt=False
        # print "new step=", yy, qq 
    print "best cell=",yy,qq
    outCell.append((yy,qq))
    return outCell

def findBest_search(hdict,valDict,signDict, degreeYList, degreeQtList,Npt,VALONLY, s,coeff) :# look for the best combination y,qt with with all possible combination (optimized). Ambiguities -->chi2 comparsion
    
    print "new optimal search", s, coeff, 
    
    def FindMinDeg(goodPoints,k) :
        degMin= goodPoints[0]
        for deg in goodPoints :
            if k=='y' :
                if deg[0]<degMin[0] :
                    degMin= deg
                elif deg[0]==degMin[0] :
                    if deg[1]<degMin[1] :
                        degMin=deg
            if k=='qt' :
                if deg[1]<degMin[1] :
                    degMin= deg
                elif deg[1]==degMin[1] :
                    if deg[0]<degMin[0] :
                        degMin=deg
        return degMin
    
    def BuildGoodList(goodPoints,degMin,endPath=False) :
        goodCombDict[(degMin[0],degMin[1])] = []
        endPathList = []
        for yDeg in degreeYList :
            if yDeg<degMin[0] : continue
            for qtDeg in degreeQtList :
                if qtDeg<degMin[1] : continue
                if yDeg==degMin[0] and qtDeg==degMin[1] : continue
                if (yDeg,qtDeg) not in goodPoints : continue 
                pval = F_test(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=degreeQtList,Npt=Npt,VALONLY=VALONLY,s=s,coeff=coeff,yMax=yDeg,qtMax=qtDeg,yMin=degMin[0],qtMin=degMin[1]) 
                if pval<pthr :
                    if not endPath :
                        goodCombDict[(degMin[0],degMin[1])].append((yDeg,qtDeg)) 
                    else :
                        endPathList.append((yDeg,qtDeg)) 
        if endPath:
            return endPathList
    
    def endPathLoop(endPathPoints) :
        # print "internal",endPathPoints
        degMinY = FindMinDeg(endPathPoints,'y')
        degMinQt = FindMinDeg(endPathPoints,'qt')
        if degMinY == degMinQt :
            endPathList = BuildGoodList(goodPoints=endPathPoints,degMin=degMinY,endPath=True)
            if len(endPathList) ==0 :
                endPathList.append(degMinY)    
            return endPathList
        else :
            print "WARNING: manual comparison needed ", coeff, s, "the list of candidates is:", endPathPoints
            endPathList = []
            endPathList.append("WARNING")
            return endPathList
                                
    pthr=0.05
    endPathPoints = []
    
    goodCombDict = {}
    yy = degreeYList[0]
    qq = degreeQtList[0]
    
    #comparison with init cell
    goodCombDict[(yy,qq)] = []
    for yDeg in degreeYList :
        for qtDeg in degreeQtList :
           if yDeg==degreeYList[0] and qtDeg==degreeQtList[0] : continue
           pval = F_test(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=degreeQtList,Npt=Npt,VALONLY=VALONLY,s=s,coeff=coeff,yMax=yDeg,qtMax=qtDeg,yMin=yy,qtMin=qq) 
           if pval<pthr :
               goodCombDict[(yy,qq)].append((yDeg,qtDeg))
    # goodCombDict[(yy,qq)].append("DONE")
    
    
    # for initPoint,goodPoints in goodCombDict.iteritems():
    while len(list(goodCombDict.keys()))>0 :
        
        # print "keys:", list(goodCombDict.keys())
        for initPoint in list(goodCombDict.keys()) :
            goodPoints = goodCombDict[initPoint]
            # print "START process:", initPoint, ", N good point:", len(goodPoints)
            if len(goodPoints) == 0 :
                # print "path ended:", initPoint
                endPathPoints.append(initPoint)
                del goodCombDict[initPoint]
                continue 
            degMinY = FindMinDeg(goodPoints,'y')
            degMinQt = FindMinDeg(goodPoints,'qt')
            if degMinY == degMinQt :
                branch = False
                # print "NOT BRANCHED," "(minY,minQt)", degMinY 
                BuildGoodList(goodPoints=goodPoints,degMin=degMinY) #add the new comparison list with respect to the new "minimum deg"
                del goodCombDict[initPoint] #because it is not branched the previous goodPoints are not useful (if degMin is the only min and is better wrt init, can be the new init)
            else :
                branch = True
                # print "BRANCHING", "degMinY=", degMinY, ", degMinQt=", degMinQt
                BuildGoodList(goodPoints=goodPoints,degMin=degMinY) #add the new comparison list with respect to the new "minimum deg"
                BuildGoodList(goodPoints=goodPoints,degMin=degMinQt) #add the new comparison list with respect to the new "minimum deg"
                # goodCombDict[initPoint].append("DO NOT BRANCH ME AGAIN") 
                del goodCombDict[initPoint]
    
    while(len(endPathPoints)>1) :
        #    print "endloop!"
           endPathPoints = endPathLoop(endPathPoints)
           if "WARNING" in endPathPoints : break
           
           
    print "chosen cell=",endPathPoints
    if "WARNING" in endPathPoints :
        endPathPoints[0] = (0,0) 
    
    return endPathPoints
             
    # print "END OF DEBUG------------------------", s, coeff        
            
def findBest_unpol(hdict,valDict,signDict, degreeYList, qtBin,Npt,s) :
    pthr=0.05
    outCell = degreeYList[0]
    yy = degreeYList[0]
    finalStop = False
    contY=True 
    skip=1
    while contY :
        if (yy+skip)<=degreeYList[-1] :
            yMax = yy+skip
            # print "comparing: ", yy, yMax
            pvalY = F_test(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=qtBin,Npt=Npt,VALONLY=VALONLY,s=s,coeff='unpol',yMax=yMax,qtMax=qtBin,yMin=yy,qtMin=qtBin)
        else :
            contY=False
            pvalY=1
        if pvalY<pthr :
            yy+=skip
            skip=1
        else :
            skip+=1
    outCell = yy
    print "best unpol serarch",s, "=",yy
    return yy        
    
    
def Validation(hdict, signDict, coeffDict, degreeYList,degreeQtList) :
    
    valDict = {}
    
    # chi2 (single degree) analysis
    for yDeg in degreeYList :
        for qtDeg in degreeQtList :
            for s,sName in signDict.iteritems() :
                valDict[str(yDeg)+str(qtDeg)+'chi2'+s] = ROOT.TH1F("hChi2"+s+'_'+str(yDeg)+str(qtDeg),"hChi2"+s+'_'+str(yDeg)+str(qtDeg),len(coeffDict),-0.5,len(coeffDict)-0.5)
                valDict[str(yDeg)+str(qtDeg)+'chi2'+s].SetLineWidth(3)
                valDict[str(yDeg)+str(qtDeg)+'chi2'+s].GetXaxis().SetTitle("Angular Coefficient")
                valDict[str(yDeg)+str(qtDeg)+'chi2'+s].GetYaxis().SetTitle("#chi^{2}/NDF")
                valDict[str(yDeg)+str(qtDeg)+'chi2'+s].SetStats(0)
                valDict[str(yDeg)+str(qtDeg)+'chi2'+s].GetYaxis().SetRangeUser(0,6)
                # binN=1
                for coeff,constr in coeffDict.iteritems() :
                    binN=int(coeff.replace('A',''))+1
                    if not (VALONLY or MULTICORE) :
                        if hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Ndf()>0 :
                            valDict[str(yDeg)+str(qtDeg)+'chi2'+s].SetBinContent(binN,hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Chi2()/hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Ndf())
                        else :
                            print "WARNING: NDF=0 for", coeff, s, " y=",yDeg, " qt=",qtDeg
                            valDict[str(yDeg)+str(qtDeg)+'chi2'+s].SetBinContent(binN,0)
                        valDict[str(yDeg)+str(qtDeg)+'chi2'+s].SetBinError(binN,math.sqrt(2*hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Ndf())/hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Chi2())
                    else :
                        if hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetNDF() >0 :  
                            valDict[str(yDeg)+str(qtDeg)+'chi2'+s].SetBinContent(binN,hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetChisquare()/hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetNDF())
                            # print('DEBUG', "h1", s, yDeg, qtDeg,binN, hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetChisquare(), hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetNDF(),valDict[str(yDeg)+str(qtDeg)+'chi2'+s].GetBinContent(binN))
                        else :
                            print "WARNING: NDF=0 for", coeff, s, " y=",yDeg, " qt=",qtDeg
                            valDict[str(yDeg)+str(qtDeg)+'chi2'+s].SetBinContent(binN,0)
                        valDict[str(yDeg)+str(qtDeg)+'chi2'+s].SetBinError(binN,math.sqrt(2*hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetNDF())/hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetChisquare())

                    # binN+=1
                    
            valDict[str(yDeg)+str(qtDeg)+'chi2_c'] = ROOT.TCanvas("c_chi2"+'_'+str(yDeg)+str(qtDeg), "c_chi2"+'_'+str(yDeg)+str(qtDeg), 800, 600)
            valDict[str(yDeg)+str(qtDeg)+'chi2_c'].cd()
            valDict[str(yDeg)+str(qtDeg)+'chi2_c'].SetGridx()
            valDict[str(yDeg)+str(qtDeg)+'chi2_c'].SetGridy()
            valDict[str(yDeg)+str(qtDeg)+'chi2'+'Plus'].SetLineColor(1)
            valDict[str(yDeg)+str(qtDeg)+'chi2'+'Minus'].SetLineColor(2)
            valDict[str(yDeg)+str(qtDeg)+'chi2'+'Plus'].Draw()
            valDict[str(yDeg)+str(qtDeg)+'chi2'+'Minus'].Draw("SAME")
            valDict[str(yDeg)+str(qtDeg)+'chi2_leg'] = ROOT.TLegend(0.7,0.8,0.95,0.95)
            for s,sname in signDict.iteritems() : 
                valDict[str(yDeg)+str(qtDeg)+'chi2_leg'].AddEntry(valDict[str(yDeg)+str(qtDeg)+'chi2'+s],s)
            valDict[str(yDeg)+str(qtDeg)+'chi2_leg'].Draw("SAME") 
    
    #chi2 varying degree analysis
    degreeYBins_temp = []
    degreeQtBins_temp =  []
    for yDeg in degreeYList :
        degreeYBins_temp.append(yDeg-0.5)
    degreeYBins_temp.append(degreeYList[-1]+0.5)
    for qtDeg in degreeQtList :
        degreeQtBins_temp.append(qtDeg-0.5)
    degreeQtBins_temp.append(degreeQtList[-1]+0.5)
    

    
    degreeQtBins = array('f',degreeQtBins_temp)
    degreeYBins = array('f',degreeYBins_temp)
    for s,sName in signDict.iteritems() :
         for coeff,constr in coeffDict.iteritems() :
                valDict['h2_chi2'+s+coeff] = ROOT.TH2F("h2_chi2_"+s+"_"+coeff, "h2_chi2_"+s+"_"+coeff, len(degreeYList), degreeYBins , len(degreeQtList), degreeQtBins)
                binY = 0 
                valDict['h2_chi2'+s+coeff].GetXaxis().SetTitle("Y polinomial degree")
                valDict['h2_chi2'+s+coeff].GetYaxis().SetTitle("q_{T} polinomial degree")
                for yDeg in degreeYList :
                    binY+=1
                    binQt=0
                    for qtDeg in degreeQtList : 
                       binQt+=1
                       if not (VALONLY or MULTICORE) :
                            if hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Ndf()>0 :   
                                valDict['h2_chi2'+s+coeff].SetBinContent(binY,binQt,hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Chi2()/hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Ndf())   
                            else :
                                print "WARNING: NDF=0 for", coeff, s, " y=",yDeg, " qt=",qtDeg
                                valDict['h2_chi2'+s+coeff].SetBinContent(binY,binQt,0)
                            valDict['h2_chi2'+s+coeff].SetBinError(binY,binQt,math.sqrt(2*hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Ndf())/hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Chi2())
                       else :
                            if hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetNDF()>0 :  
                                valDict['h2_chi2'+s+coeff].SetBinContent(binY,binQt,hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetChisquare()/hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetNDF())   
                                # print('DEBUG2', "h2", s, yDeg, qtDeg,coeff, hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetChisquare(), hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetNDF(), valDict['h2_chi2'+s+coeff].GetBinContent(binY,binQt))
                            else :
                                print "WARNING: NDF=0 for", coeff, s, " y=",yDeg, " qt=",qtDeg
                                valDict['h2_chi2'+s+coeff].SetBinContent(binY,binQt,0)
                            valDict['h2_chi2'+s+coeff].SetBinError(binY,binQt,math.sqrt(2*hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetNDF())/hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetChisquare())
    
    #summary of width of the pulls
    for var in ['h2_pullsSigma', 'h2_pullsSigmaFit'] : #direct from histo o gaussian fit
        for s,sName in signDict.iteritems() :
            for coeff,constr in coeffDict.iteritems() :
                valDict[var+s+coeff] = ROOT.TH2F(var+s+"_"+coeff, var+s+"_"+coeff, len(degreeYList), degreeYBins , len(degreeQtList), degreeQtBins)
                binY = 0 
                valDict[var+s+coeff].GetXaxis().SetTitle("Y polinomial degree")
                valDict[var+s+coeff].GetYaxis().SetTitle("q_{T} polinomial degree")
                for yDeg in degreeYList :
                    binY+=1
                    binQt=0
                    for qtDeg in degreeQtList : 
                        binQt+=1
                        if 'Fit' in var :
                            try :
                                valDict[var+s+coeff].SetBinContent(binY,binQt,hdict[str(yDeg)+str(qtDeg)+s+coeff]['pull1D'+sName+coeff].GetFunction("gaus").GetParameter(2))
                                valDict[var+s+coeff].SetBinError(binY,binQt,hdict[str(yDeg)+str(qtDeg)+s+coeff]['pull1D'+sName+coeff].GetFunction("gaus").GetParError(2))  
                            except :
                                valDict[var+s+coeff].SetBinContent(binY,binQt,0)
                                valDict[var+s+coeff].SetBinError(binY,binQt,0) 
                        else :
                            valDict[var+s+coeff].SetBinContent(binY,binQt,hdict[str(yDeg)+str(qtDeg)+s+coeff]['pull1D'+sName+coeff].GetStdDev())   
                            valDict[var+s+coeff].SetBinError(binY,binQt,hdict[str(yDeg)+str(qtDeg)+s+coeff]['pull1D'+sName+coeff].GetStdDevError())
    
    
    
    #evaluating p-value for the F-test. see my logbook for details
    
    #evaluating the number of bins of 2D istograms. these loops are a single event only actually (see the breaks)
    for s,sName in signDict.iteritems() :
        for coeff,constr in coeffDict.iteritems() :
            NptX= hdict[str(degreeYList[0])+str(degreeQtList[0])+s+coeff]['fit'+sName+coeff].GetXaxis().GetNbins()
            NptY= hdict[str(degreeYList[0])+str(degreeQtList[0])+s+coeff]['fit'+sName+coeff].GetYaxis().GetNbins()
            Npt = NptX*NptY
            break
        break
    # print "WARNING: hardcoded number of point fitted: 19x19=361" 
    # Npt=361#number of point fitted
    # print "CHECK: number of point fitted=", Npt
    for kind in ['pval_alongY','pval_alongQt','pval_alongDiagonal'] :
        for s,sName in signDict.iteritems() :
            for coeff,constr in coeffDict.iteritems() :
                valDict[kind+s+coeff] = ROOT.TH2F(kind+s+"_"+coeff, kind+s+"_"+coeff, len(degreeYList), degreeYBins , len(degreeQtList), degreeQtBins)
                binY = 0 
                valDict[kind+s+coeff].GetXaxis().SetTitle("Y polinomial degree")
                valDict[kind+s+coeff].GetYaxis().SetTitle("q_{T} polinomial degree")
                for yDeg in degreeYList :
                    binY+=1
                    binQt=0
                    for qtDeg in degreeQtList : 
                        binQt+=1  
                        if binQt==1 and binY==1 : continue #no deltaChi2 for first bin

                        # #evaluate number of parameter
                        # if (kind=='pval_alongQt' and binQt!=1) or binY==1:#vertical delta
                        #     NDFstring = str(yDeg)+str(degreeQtList[degreeQtList.index(qtDeg)-1])
                        # if (kind=='pval_alongY' and binY!=1) or binQt==1: #horizontal delta
                        #     NDFstring = str(degreeYList[degreeYList.index(yDeg)-1])+str(qtDeg) 
                        # if not VALONLY :
                        #    NDFMax = hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].Ndf()
                        #    NDFMin = hdict[NDFstring+s+coeff]['fit'+'Res'+sName+coeff].Ndf()
                        # else :
                        #     NDFMax = hdict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+'Res'+sName+coeff].GetNDF()
                        #     NDFMin = hdict[NDFstring+s+coeff]['fit'+'Res'+sName+coeff].GetNDF()                  
                        
                        # parMax = Npt-NDFMax
                        # parMin = Npt-NDFMin
                        
                        # #evaluate chi2 in the bin and in the previous
                        # chi2Max =  valDict['h2_chi2'+s+coeff].GetBinContent(binY,binQt)*NDFMax
                        # # parMin=parMax-1                        
                        # if (kind=='pval_alongQt' and binQt!=1) or binY==1:#vertical delta
                        #     chi2Min =  valDict['h2_chi2'+s+coeff].GetBinContent(binY,binQt-1)*NDFMin
                        # if (kind=='pval_alongY' and binY!=1) or binQt==1: #horizontal delta
                        #     chi2Min =  valDict['h2_chi2'+s+coeff].GetBinContent(binY-1,binQt)*NDFMin
                        
                        # #evaluate F
                        # Fobs = (chi2Min-chi2Max)*(Npt-parMax)/(chi2Max*(parMax-parMin))
                        # if Fobs>0 :
                        #     pvalue=1-scipy.stats.f.cdf(Fobs,1,(Npt-parMax))
                        # else :
                        #     pvalue=1
                        if kind=='pval_alongDiagonal' and binQt!=1 and binY!=1 : # and binQt==binY:
                            valAlt = F_test(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=degreeQtList,Npt=Npt,VALONLY=VALONLY,   s=s,coeff=coeff,yMax=yDeg,qtMax=qtDeg,yMin=yDeg-1,qtMin=qtDeg-1)
                        if (kind=='pval_alongQt' and binQt!=1) or binY==1:  #or (kind=='pval_alongDiagonal' and binY<binQt and binQt!=1): #vertical delta
                            pvalAlt = F_test(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=degreeQtList,Npt=Npt,VALONLY=VALONLY,   s=s,coeff=coeff,yMax=yDeg,qtMax=qtDeg,yMin=yDeg,qtMin=qtDeg-1) 
                        if (kind=='pval_alongY' and binY!=1) or binQt==1 :# or (kind=='pval_alongDiagonal' and binY>binQt and binY!=1): #horizontal delta
                            pvalAlt = F_test(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=degreeQtList,Npt=Npt,VALONLY=VALONLY,   s=s,coeff=coeff,yMax=yDeg,qtMax=qtDeg,yMin=yDeg-1,qtMin=qtDeg) 
                        # if pvalue!=pvalAlt :
                        #     print "DEBUG", binY,binQt, NDFMin,NDFMax
                        #     print "--->pvalue, pvalueAlt=", pvalue, pvalAlt,s,coeff, kind," ydeg,qtdeg=",yDeg,qtDeg
                        valDict[kind+s+coeff].SetBinContent(binY,binQt,pvalAlt)
                        
    #best combination histogram filling
    for kind in ['chosenDeg_path','chosenDeg_search'] :
        h_chosen = ROOT.TH2F(kind, kind, len(degreeYList), degreeYBins , len(degreeQtList), degreeQtBins)
        valDict[kind] = ROOT.TCanvas(kind, kind, 800, 600)
        valDict[kind].cd()
        h_chosen.GetXaxis().SetTitle("Y polinomial degree")
        h_chosen.GetYaxis().SetTitle("q_{T} polinomial degree")
        valDict[kind].SetGridx()
        valDict[kind].SetGridy()
        h_chosen.SetStats(0)
        h_chosen.DrawCopy()
        printedList = [] #multiple list needed to avoid overlapping writing if a 2 or 3 coeff. have same combination
        printedList2 = []
        for s,sName in signDict.iteritems() :
            for coeff,constr in coeffDict.iteritems() :
                if kind == 'chosenDeg_path' : 
                    cell = findBest_path(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=degreeQtList,Npt=Npt,VALONLY=VALONLY,s=s,coeff=coeff)
                if kind == 'chosenDeg_search' :
                    cell = findBest_search(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=degreeQtList,Npt=Npt,VALONLY=VALONLY,s=s,coeff=coeff)
                cellname = ROOT.TLatex()
                cellname.SetTextSize(0.035)
                cellname.SetTextAlign(12)
                if cell not in printedList : 
                    cellname.DrawLatex(cell[0][0]-0.3,cell[0][1], coeff+'_'+s)
                elif cell not in printedList2:
                    cellname.DrawLatex(cell[0][0]-0.3,cell[0][1]-0.2, coeff+'_'+s)
                    printedList2.append(cell)
                else :
                    cellname.DrawLatex(cell[0][0]-0.3,cell[0][1]-0.4, coeff+'_'+s)
                printedList.append(cell)
                 
                # useful output lines 
                # if not VALONLY :
                #     chi2abs = hdict[str(cell[0][0])+str(cell[0][1])+s+coeff]['fit'+'Res'+sName+coeff].Chi2()
                #     ndf4chi2 = hdict[str(cell[0][0])+str(cell[0][1])+s+coeff]['fit'+'Res'+sName+coeff].Ndf()
                # else :
                #     chi2abs = hdict[str(cell[0][0])+str(cell[0][1])+s+coeff]['fit'+'Res'+sName+coeff].GetChisquare()
                #     ndf4chi2 = hdict[str(cell[0][0])+str(cell[0][1])+s+coeff]['fit'+'Res'+sName+coeff].GetNDF()
                # chi2red = chi2abs/ndf4chi2
                # chi2redErr=math.sqrt(2*ndf4chi2)/chi2abs
                # print kind, ":::", s,coeff,"(y,qt)=", cell, ", chi2_abs=", chi2abs, ", NDF=",ndf4chi2, ", chi2_reduced=", chi2red, "+/-",chi2redErr
                
        valDict[kind].Update()
        
        
    # -------- validation of unpol analysis ------------ #
    NbinsQt = hdict[str(degreeYList[0])+s+'unpol']['hh'+signDict[s]+'unpol'+'MC'].GetYaxis().GetNbins()
    NbinsY = hdict[str(degreeYList[0])+s+'unpol']['hh'+signDict[s]+'unpol'+'MC'].GetXaxis().GetNbins()
    qtBins_temp = []
    for qi in range(1,NbinsQt+2) :
        qtBins_temp.append( hdict[str(degreeYList[0])+s+'unpol']['hh'+signDict[s]+'unpol'+'MC'].GetYaxis().GetBinLowEdge(qi)) 
    qtBins = array('f', qtBins_temp)
    
    #chi2 2D plots
    for s,sName in signDict.iteritems() : 
        valDict['h2_chi2'+s+'unpol'] = ROOT.TH2F("h2_chi2_"+s+"_"+'unpol', "h2_chi2_"+s+"_"+'unpol', len(degreeYList), degreeYBins ,NbinsQt, qtBins)
        valDict['h2_chi2'+s+'unpol'].GetXaxis().SetTitle("Y polinomial degree")
        valDict['h2_chi2'+s+'unpol'].GetYaxis().SetTitle("q_{T} bins [1=24 GeV]")
        binY=0
        for yDeg in degreeYList :
            binY +=1 
            for iq in range(1, NbinsQt+1) :
                if hdict[str(yDeg)+s+'unpol']["fitRes"+sName+'unpol'+str(iq)].Ndf()>0 :
                    valDict['h2_chi2'+s+'unpol'].SetBinContent(binY,iq, hdict[str(yDeg)+s+'unpol']["fitRes"+sName+'unpol'+str(iq)].Chi2()/hdict[str(yDeg)+s+'unpol']["fitRes"+sName+'unpol'+str(iq)].Ndf())
                else :
                    print "WARNING: NDF=0 for, unpol ",s, " y=",binY, " qt=",iq
                    valDict['h2_chi2'+s+'unpol'].SetBinContent(binY,iq, 0)
                valDict['h2_chi2'+s+'unpol'].SetBinError(binY,iq,math.sqrt(2*hdict[str(yDeg)+s+'unpol']["fitRes"+sName+'unpol'+str(iq)].Ndf())/hdict[str(yDeg)+s+'unpol']["fitRes"+sName+'unpol'+str(iq)].Chi2())

    #F-Test plots
    for s,sName in signDict.iteritems() : 
        valDict['pval_alongY'+s+'unpol'] = ROOT.TH2F("pval_alongY_"+s+"_"+'unpol', "pval_alongY_"+s+"_"+'unpol', len(degreeYList), degreeYBins ,NbinsQt, qtBins)
        valDict['pval_alongY'+s+'unpol'].GetXaxis().SetTitle("Y polinomial degree")
        valDict['pval_alongY'+s+'unpol'].GetYaxis().SetTitle("q_{T} bins [1=24 GeV]")
        binY=0
        for yDeg in degreeYList :
            binY +=1 
            if yDeg==degreeYList[0] :
                continue
            for iq in range(1, NbinsQt+1) :
                pvalAlt = F_test(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, degreeQtList=degreeQtList,Npt=NbinsY,VALONLY=VALONLY,   s=s,coeff='unpol',yMax=yDeg,qtMax=iq,yMin=yDeg-1,qtMin=iq) 
                valDict['pval_alongY'+s+'unpol'].SetBinContent(binY,iq,pvalAlt)
    
    #best combination histogram
    ccol=2
    valDict['chosen_unpol'] = ROOT.TCanvas("chosen_unpol","chosen_unpol",800,600)
    valDict['chosen_unpol'].cd()
    valDict['chosen_unpol'].SetGridx()
    valDict['chosen_unpol'].SetGridy()
    sameFlag='P'
    valDict['chosen_unpol'+'leg'] = ROOT.TLegend(0.7,0.8,0.95,0.95)
    for s,sName in signDict.iteritems() :
        valDict['chosen_unpol'+s] = ROOT.TH1F("chosen_unpol_"+s,"chosen_unpol_"+s,NbinsQt, qtBins)
        valDict['chosen_unpol'+s].GetXaxis().SetTitle("q_{T} bins [1=24 GeV]")
        valDict['chosen_unpol'+s].GetYaxis().SetTitle("best Y polinomial degree")
        valDict['chosen_unpol'+s].SetLineColor(ccol)
        valDict['chosen_unpol'+s].SetMarkerColor(ccol)
        valDict['chosen_unpol'+s].SetLineWidth(3)
        valDict['chosen_unpol'+s].SetMarkerStyle(18+ccol)
        valDict['chosen_unpol'+s].GetYaxis().SetRangeUser(0,degreeYList[-1]+1)
        valDict['chosen_unpol'+s].SetStats(0)
        ccol+=2
        for iq in range(1, NbinsQt+1) :
            cell = findBest_unpol(hdict=hdict,valDict=valDict,signDict=signDict, degreeYList=degreeYList, qtBin=iq,Npt=NbinsY,s=s) 
            valDict['chosen_unpol'+s].SetBinContent(iq,cell)
        valDict['chosen_unpol'+s].DrawCopy(sameFlag)
        sameFlag = 'P SAME'
        valDict['chosen_unpol'+'leg'].AddEntry(valDict['chosen_unpol'+s],s)
    valDict['chosen_unpol'+'leg'].Draw("SAME")
            
        
        
                            
    return valDict
    
        
def figureSaver(cumulativeDict, valDict, outputname,signDict, coeffDict,degreeYList,degreeQtList, validationList,unpolValList) :
    
    if not os.path.isdir(output_dir+'/'+outputname): os.system('mkdir '+output_dir+'/'+outputname)
    
    for s,sname in signDict.iteritems() :
        for coeff,constr in coeffDict.iteritems() : 
            for var in  validationList :
                if 'chosenDeg' in var :
                    if s=='Plus' and coeff=='A0':
                        valDict[var].SaveAs(output_dir+'/'+outputname+'/'+var+'.png')  
                else :  
                    c_general = ROOT.TCanvas('c_general', 'c_general', 800, 600)
                    c_general.cd()
                    c_general.SetLogz()
                    valDict[var+s+coeff].SetStats(0)
                    valDict[var+s+coeff].SetContour(200)
                    valDict[var+s+coeff].Draw('text colz')
                    if "pval" in var :
                        valDict[var+s+coeff].GetZaxis().SetRangeUser(10**(-16),1)
                    c_general.SaveAs(output_dir+'/'+outputname+'/'+var+'_'+s+'_'+coeff+'.png')
    
    initDict = {
        'MC_' : ['hh','MC'],
        'coeff_' : ['hh',''],
        'coeffErr_' : ['hh','Err'],
    }
    
    for s,sname in signDict.iteritems() :
        for coeff,constr in coeffDict.iteritems() :
            for ind, val in initDict.iteritems() :
                    c_init = ROOT.TCanvas('c_init', 'c_init', 800, 600)
                    c_init.cd()
                    cumulativeDict[str(3)+str(3)+s+coeff][val[0]+sname+coeff+val[1]].SetStats(0)                        
                    cumulativeDict[str(3)+str(3)+s+coeff][val[0]+sname+coeff+val[1]].Draw('colz')
                    c_init.SaveAs(output_dir+'/'+outputname+'/'+s+'_'+coeff+'_'+ind+'.png')

    
    for yDeg in degreeYList :
        for qtDeg in degreeQtList :
            for s,sname in signDict.iteritems() :
                for coeff,constr in coeffDict.iteritems() :
                    
                    c_fit3D = ROOT.TCanvas('c_fit3D', 'c_fit3D', 800, 600)
                    c_fit3D.cd()
                    cumulativeDict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+sname+coeff].SetStats(0)                        
                    cumulativeDict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+sname+coeff].DrawCopy('LEGO2')
                    c_fit3D.SaveAs(output_dir+'/'+outputname+'/'+s+'_'+coeff+'_'+str(yDeg)+str(qtDeg)+'_Fit3D.png')
                    
                    c_fit = ROOT.TCanvas('c_fit', 'c_fit', 800, 600)
                    c_fit.cd()
                    cumulativeDict[str(yDeg)+str(qtDeg)+s+coeff]['fit'+sname+coeff].Draw('colz')
                    c_fit.SaveAs(output_dir+'/'+outputname+'/'+s+'_'+coeff+'_'+str(yDeg)+str(qtDeg)+'_Fit.png')
                    
                    c_scratch = ROOT.TCanvas('c_scratch', 'c_scratch', 1200, 600)
                    c_scratch.cd()
                    ROOT.gROOT.SetBatch(True)
                    cumulativeDict[str(yDeg)+str(qtDeg)+s+coeff]['pullc'+sname+coeff].Draw()
                    cumulativeDict[str(yDeg)+str(qtDeg)+s+coeff]['pullc'+sname+coeff].SetCanvasSize(1200,600)
                    cumulativeDict[str(yDeg)+str(qtDeg)+s+coeff]['pullc'+sname+coeff].SaveAs(output_dir+'/'+outputname+'/'+s+'_'+coeff+'_'+str(yDeg)+str(qtDeg)+'_Pulls.png')
            
            valDict[str(yDeg)+str(qtDeg)+'chi2_c'].SaveAs(output_dir+'/'+outputname+'/chi2_'+str(yDeg)+str(qtDeg)+'.png')
    
    #unpol figure
    for s,sname in signDict.iteritems() :
        for var in  unpolValList :
            if 'chosen' in var :
                if s=='Plus' :
                    valDict[var].SaveAs(output_dir+'/'+outputname+'/'+var+'.png')  
            else :  
                c_general = ROOT.TCanvas('c_general', 'c_general', 800, 600)
                c_general.cd()
                c_general.SetLogz()
                valDict[var+s+'unpol'].SetStats(0)
                valDict[var+s+'unpol'].SetContour(200)
                valDict[var+s+'unpol'].Draw('text colz')
                if "pval" in var :
                    valDict[var+s+'unpol'].GetZaxis().SetRangeUser(10**(-16),1)
                c_general.SaveAs(output_dir+'/'+outputname+'/'+var+'_'+s+'_'+'unpol'+'.png')
    
    for yDeg in degreeYList :
        for s,sname in signDict.iteritems() : 
                        
            c_scratch = ROOT.TCanvas('c_scratch', 'c_scratch', 1200, 600)
            c_scratch.cd()
            ROOT.gROOT.SetBatch(True)
            cumulativeDict[str(yDeg)+s+'unpol']['pullc'+sname+'unpol'].Draw()
            cumulativeDict[str(yDeg)+s+'unpol']['pullc'+sname+'unpol'].SetCanvasSize(1200,600)
            cumulativeDict[str(yDeg)+s+'unpol']['pullc'+sname+'unpol'].SaveAs(output_dir+'/'+outputname+'/'+s+'_'+'unpol'+'_'+str(yDeg)+'_Pulls.png')    
        
            NbinsQt = cumulativeDict[str(degreeYList[0])+s+'unpol']['hh'+signDict[s]+'unpol'+'MC'].GetYaxis().GetNbins()
            for iq in range(1, NbinsQt+1) :
                c_fit = ROOT.TCanvas('c_fit', 'c_fit', 800, 600)
                c_fit.cd()
                cumulativeDict[str(yDeg)+s+'unpol']['fit'+sname+'unpol'+str(iq)].Draw()
                c_fit.SaveAs(output_dir+'/'+outputname+'/'+s+'_'+'unpol'+'_'+str(yDeg)+str(iq)+'_Fit.png')
    

def systComparison(bestCombDict, saveRoot=True,saveFig=False,excludeSyst='',statMC=True) :
    #statMC=True --> use statistical error form nominalfit, False-->error= |fittedValue-unfittedMC| (for coeff_ plot only)
    
    systHistoDict = {}
    systFileDict = {}
    systCanvasDict = {}
    systLegDict = {}
    rebSuffix = ''
    #open files and get histos
    for sKind, sList in systDict.iteritems() :
        for key in sList : 
            fileSuff = SUFFIX.replace('_'+SYST_kind+'_'+SYST_name,'_'+sKind+'_'+key)
            if key == '_nom_nom' : key='nominal_'
            # if key =='nominal_' or 'Scale' in key :
            # if 'Scale' in key :
            #     rebSuffix = 'rebuild_'
            # else :
            #     rebSuffix = ''
            systFileDict[key]=ROOT.TFile.Open(output_dir+"/regularizationFit_"+fileSuff+rebSuffix+".root")
            # systFileDict[key]=ROOT.TFile.Open(output_dir+"/regularizationFit_"+SUFFIX+".root")
            
            # systFileDict[key]=ROOT.TFile.Open(output_dir+"/regularizationFit_"+SUFFIX+"rebuild__"+key+".root")
            rebSuffix = ''
            for s, sName in signDict.iteritems() :
                for coeff,constr in coeffDict.iteritems() :  
                    systHistoDict['chi2'+s+coeff+key] = systFileDict[key].Get('h2_chi2_'+s+'_'+coeff)
                    for yDeg in degreeYList :
                        for qtDeg in degreeQtList :
                            systHistoDict['coeff'+s+coeff+key+str(yDeg)+str(qtDeg)] = systFileDict[key].Get('Y'+str(yDeg)+'_Qt'+str(qtDeg)+'/hfit_'+sName+'_'+coeff+'_'+str(yDeg)+str(qtDeg))
                            systHistoDict['coeffErr'+s+coeff+key+str(yDeg)+str(qtDeg)] = systFileDict[key].Get('Y'+str(yDeg)+'_Qt'+str(qtDeg)+'/hfitErr_'+sName+'_'+coeff+'_'+str(yDeg)+str(qtDeg))
                            # systHistoDict['coeffUnfit'+s+coeff+key+str(yDeg)+str(qtDeg)] = systFileDict[key].Get('Y'+str(yDeg)+'_Qt'+str(qtDeg)+'/coeff_'+sName+'_'+coeff)
                            # systHistoDict['coeffErrUnfit'+s+coeff+key+str(yDeg)+str(qtDeg)] = systFileDict[key].Get('Y'+str(yDeg)+'_Qt'+str(qtDeg)+'/coeffErr_'+sName+'_'+coeff)
                            

                #unpol
                NbinsQt = systHistoDict['coeff'+s+coeff+key+str(yDeg)+str(qtDeg)].GetYaxis().GetNbins()
                systHistoDict['chi2'+s+'unpol'+key] = systFileDict[key].Get('h2_chi2_'+s+'_'+'unpol')
                for yDeg in degreeYList :
                    for i in range(1, NbinsQt+1) :
                        systHistoDict['coeff'+s+'unpol'+key+str(yDeg)+str(i)] = systFileDict[key].Get('unpol_Y'+str(yDeg)+'/hfit_'+sName+'_'+'unpol'+'_'+str(yDeg)+'_'+str(i))
                        systHistoDict['coeffErr'+s+'unpol'+key+str(yDeg)+str(i)] = systFileDict[key].Get('unpol_Y'+str(yDeg)+'/hfitErr_'+sName+'_'+'unpol'+'_'+str(yDeg)+'_'+str(i))

    #unrolled binning creation: ydeg-qtdeg
    degreeYBins_temp = []
    degreeQtBins_temp =  []
    for yDeg in degreeYList :
        degreeYBins_temp.append(yDeg-0.5)
    degreeYBins_temp.append(degreeYList[-1]+0.5)
    for qtDeg in degreeQtList :
        degreeQtBins_temp.append(qtDeg-0.5)
    degreeQtBins_temp.append(degreeYList[-1]+0.5)
    
    unrolledYqTdeg= list(degreeYBins_temp)
    intervalQtdegBin = []
    for y in degreeYBins_temp[:-1] :
        intervalQtdegBin.append(degreeYBins_temp[degreeYBins_temp.index(y)+1]-degreeYBins_temp[degreeYBins_temp.index(y)])
    for q in range(len(degreeQtBins_temp)-2) :
        for y in intervalQtdegBin :
            unrolledYqTdeg.append(unrolledYqTdeg[-1]+y)        
    
    #unrolled binning creation: y-qt
    histo4Binning = systFileDict['nominal_'].Get('Y'+str(degreeYList[0])+'_Qt'+str(degreeQtList[0])+'/hfit_WtoMuP_A0_'+str(degreeYList[0])+str(degreeQtList[0]))
    yBinning = []
    qtBinning = []
    NbinsQt = histo4Binning.GetYaxis().GetNbins()
    NbinsY = histo4Binning.GetXaxis().GetNbins()
    for qi in range(1,NbinsQt+2) :
        qtBinning.append( histo4Binning.GetYaxis().GetBinLowEdge(qi))     
    for yi in range(1,NbinsY+2) :
        yBinning.append(histo4Binning.GetXaxis().GetBinLowEdge(yi))  
        
    unrolledYqT= list(yBinning)
    intervalQtBin = []
    for y in yBinning[:-1] :
        intervalQtBin.append(yBinning[yBinning.index(y)+1]-yBinning[yBinning.index(y)])
    for q in range(len(qtBinning)-2) :
        for y in intervalQtBin :
            unrolledYqT.append(unrolledYqT[-1]+y) 
            
    #binning for unpol chi2: y-qt
    histo4Binning_unpol = systFileDict['nominal_'].Get('unpol_Y'+str(degreeYList[0])+'/hfit_WtoMuP_unpol_'+str(degreeYList[0])+'_1')
    yBinning_unpol = []
    qtBinning_unpol = []
    NbinsQt = histo4Binning.GetYaxis().GetNbins() #NOT A BUG!
    NbinsY_unpol = histo4Binning_unpol.GetXaxis().GetNbins()
    for qi in range(1,NbinsQt+2) :
        qtBinning_unpol.append( histo4Binning_unpol.GetYaxis().GetBinLowEdge(qi))     
    for yi in range(1,NbinsY+2) :
        yBinning_unpol.append(histo4Binning_unpol.GetXaxis().GetBinLowEdge(yi))  
        
    unrolledYqT_unpol= list(yBinning_unpol)
    intervalQtBin_unpol = []
    for y in yBinning_unpol[:-1] :
        intervalQtBin_unpol.append(yBinning_unpol[yBinning_unpol.index(y)+1]-yBinning_unpol[yBinning_unpol.index(y)])
    for q in range(len(qtBinning_unpol)-2) :
        # print unrolledYqT_unpol
        for y in intervalQtBin_unpol :
            unrolledYqT_unpol.append(unrolledYqT_unpol[-1]+y) 
    
    # print unrolledYqT_unpol
    
    #binning for unpol chi2: ydeg-qt
    unrolledYdegqT = list(degreeYBins_temp)
    for q in range(len(qtBinning_unpol)-2) :
        for y in intervalQtBin_unpol :
            unrolledYdegqT.append(unrolledYdegqT[-1]+y)
    
            
    colorList = [600,616,416,632,432,800,861,830]
    #blue,magenta,green,red,cyan,orange,azure+1,kspring+10
    
    # for i in colorList :
    #     colorList.append(i+5)
        
    # #qt deg - qt deg plot unrolled 
    # degValidationList = ['chi2_ratio_unrolled', 'chi2_unrolled', 'chi2_sumOfPulls'] #chi2_sumOfPulls=like pull_fNom_hVar, summing all the bins (in quadrature)
    # for kind in degValidationList :
    #     for s, sName in signDict.iteritems() :
    #         for coeff,constr in coeffDict.iteritems() :
    #             colorNumber = 0
    #             systCanvasDict[kind+s+coeff] = ROOT.TCanvas(kind+'_'+s+"_"+coeff, kind+'_'+s+"_"+coeff, 800,600)
    #             systCanvasDict[kind+s+coeff].cd()
    #             systCanvasDict[kind+s+coeff].SetGridx()
    #             systCanvasDict[kind+s+coeff].SetGridy()
    #             systLegDict[kind+s+coeff] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                
    #             #nominal
    #             systHistoDict[kind+s+coeff+'nominal_'] = ROOT.TH1F(kind+'_'+s+'_'+coeff+'_nominal_',kind+'_'+s+'_'+coeff+'_nominal_', len(unrolledYqTdeg)-1, array('f',unrolledYqTdeg))
    #             systHistoDict[kind+s+coeff+'nominal_'].SetLineColor(1)
    #             systHistoDict[kind+s+coeff+'nominal_'].SetTitle("chi2 "+s+" " +coeff+" nominal_")
    #             systHistoDict[kind+s+coeff+'nominal_'].GetXaxis().SetTitle("Unrolled (y=fine,q_{T}=large)")
    #             if kind=='chi2_ratio_unrolled' :
    #                 systHistoDict[kind+s+coeff+'nominal_'].GetYaxis().SetTitle("Syst/Nom")
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetLineWidth(0)
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetFillColor(1)
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetFillStyle(3002)
    #             if kind=='chi2_unrolled' :
    #                 systHistoDict[kind+s+coeff+'nominal_'].GetYaxis().SetTitle("#chi^{2}/NDF")
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetLineWidth(3)
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetFillColor(1)
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetFillStyle(3002)
    #             if kind=='chi2_sumOfPulls' :
    #                 systHistoDict[kind+s+coeff+'nominal_'].GetYaxis().SetTitle("#Sigma_{y,q_{T}}(fitted_Nom - histo_Syst)^{2} / histo_Syst_err^{2}")
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetLineWidth(3)
    #             for q in degreeQtBins_temp[:-1] :
    #                 for y in degreeYBins_temp[:-1] :
    #                     indexUNR = degreeYBins_temp.index(y)*(len(degreeQtBins_temp)-1)+degreeQtBins_temp.index(q)+1
    #                     nbinY = degreeYBins_temp.index(y)+1
    #                     nbinQt = degreeQtBins_temp.index(q)+1
    #                     if kind=='chi2_ratio_unrolled' :
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,1.)
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,systHistoDict['chi2'+s+coeff+'nominal_'].GetBinError(nbinY,nbinQt)/systHistoDict['chi2'+s+coeff+'nominal_'].GetBinContent(nbinY,nbinQt))
    #                     if kind=='chi2_unrolled' :
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,systHistoDict['chi2'+s+coeff+'nominal_'].GetBinContent(nbinY,nbinQt))
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,systHistoDict['chi2'+s+coeff+'nominal_'].GetBinError(nbinY,nbinQt))
    #                     if kind=='chi2_sumOfPulls' :
    #                         counterNDF=0
    #                         chi2_sumOfPulls=0
    #                         # if s=='Minus' and coeff=='A4':
    #                         #     print "---------------------------------------"
    #                         for qq in qtBinning[:-1] :
    #                             for yy in yBinning[:-1] :
    #                                 counterNDF+=1
    #                                 nbinYY = yBinning.index(yy)+1
    #                                 nbinQtt = qtBinning.index(qq)+1
    #                                 yVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetXaxis().GetBinCenter(nbinYY)
    #                                 qtVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetYaxis().GetBinCenter(nbinQtt)
    #                                 nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal,qtVal)
    #                                 pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetBinContent(nbinYY,nbinQtt)
    #                                 pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetBinError(nbinYY,nbinQtt)
    #                                 pull_fNom_hVar = pull_fNom_hVar**2
    #                                 chi2_sumOfPulls+=pull_fNom_hVar
    #                                 # if s=='Minus' and coeff=='A4' and qq==qtBinning[2] and yy==yBinning[2]:
    #                                 #     print qq, yy, chi2_sumOfPulls
    #                                 #     print "Y DEG=", degreeYList[degreeYBins_temp.index(y)],y, pull_fNom_hVar
    #                                 #     print "fit only value", nominalFitVal
    #                         chi2_sumOfPulls=chi2_sumOfPulls/counterNDF   
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,chi2_sumOfPulls)
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,0)
    #             systHistoDict[kind+s+coeff+'nominal_'].Draw()
    #             if kind=='chi2_ratio_unrolled' :
    #                 systHistoDict[kind+s+coeff+'nominal_'].Draw('SAME E2')
    #             if kind=='chi2_unrolled' :
    #                 systHistoDict[kind+s+coeff+'nominal_'].Draw('SAME E3 L')
    #             systLegDict[kind+s+coeff].AddEntry(systHistoDict[kind+s+coeff+'nominal_'], 'nominal_')
                
    #             #variation
    #             for key, val in systDict.iteritems() :
    #                 if key=='nominal_' : continue
    #                 systHistoDict[kind+s+coeff+key] = ROOT.TH1F(kind+'_'+s+'_'+coeff+'_'+key,kind+'_'+s+'_'+coeff+'_'+key, len(unrolledYqTdeg)-1, array('f',unrolledYqTdeg))
    #                 systHistoDict[kind+s+coeff+key].SetLineColor(colorList[colorNumber])
    #                 systHistoDict[kind+s+coeff+key].SetTitle("chi2 "+s+" " +coeff+" "+key)
    #                 systHistoDict[kind+s+coeff+key].GetXaxis().SetTitle("Unrolled (y=fine,q_{T}=large)")
    #                 if kind=='chi2_ratio_unrolled' :
    #                     systHistoDict[kind+s+coeff+key].GetYaxis().SetTitle("Syst/Nom")  
    #                 if kind=='chi2_unrolled' :
    #                     systHistoDict[kind+s+coeff+key].GetYaxis().SetTitle("#chi^{2}/NDF") 
    #                 if kind=='chi2_sumOfPulls' :
    #                     systHistoDict[kind+s+coeff+key].GetYaxis().SetTitle("#Sigma_{y,q_{T}}(fitted_Nom - histo_Syst)^{2} / histo_Syst_err^{2}")
    #                 systHistoDict[kind+s+coeff+key].SetLineWidth(3)  
    #                 for q in degreeQtBins_temp[:-1] :
    #                     for y in degreeYBins_temp[:-1] :
    #                         indexUNR = degreeYBins_temp.index(y)*(len(degreeQtBins_temp)-1)+degreeQtBins_temp.index(q)+1
    #                         nbinY = degreeYBins_temp.index(y)+1
    #                         nbinQt = degreeQtBins_temp.index(q)+1
    #                         if kind=='chi2_ratio_unrolled' :
    #                             systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,systHistoDict['chi2'+s+coeff+key].GetBinContent(nbinY,nbinQt)/systHistoDict['chi2'+s+coeff+'nominal_'].GetBinContent(nbinY,nbinQt))
    #                         if kind=='chi2_unrolled' :
    #                             systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,systHistoDict['chi2'+s+coeff+key].GetBinContent(nbinY,nbinQt))
    #                         if kind=='chi2_sumOfPulls' :
    #                             counterNDF=0
    #                             chi2_sumOfPulls=0
    #                             for qq in qtBinning[:-1] :
    #                                 for yy in yBinning[:-1] :
    #                                     counterNDF+=1
    #                                     nbinY = yBinning.index(yy)+1
    #                                     nbinQt = qtBinning.index(qq)+1
    #                                     yVal = systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetXaxis().GetBinCenter(nbinY)
    #                                     qtVal = systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetYaxis().GetBinCenter(nbinQt)
    #                                     nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal,qtVal)
    #                                     pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetBinContent(nbinY,nbinQt)
    #                                     pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetBinError(nbinY,nbinQt)
    #                                     pull_fNom_hVar = pull_fNom_hVar**2
    #                                     chi2_sumOfPulls+=pull_fNom_hVar
    #                             chi2_sumOfPulls=chi2_sumOfPulls/counterNDF
    #                             # print kind, s, coeff, key, chi2_sumOfPulls
    #                             systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,chi2_sumOfPulls)
    #                             # print "get", indexUNR, systHistoDict[kind+s+coeff+key].GetBinContent(indexUNR)
    #                 systHistoDict[kind+s+coeff+key].Draw("SAME hist")
    #                 systLegDict[kind+s+coeff].AddEntry(systHistoDict[kind+s+coeff+key], key)
    #                 colorNumber += 1
    #             systLegDict[kind+s+coeff].Draw("SAME")
    
    # # y - qt unrolled plots 
    # yqtValidationList = ['coeff_ratio_unrolled', 'coeff_unrolled', 'pull_fNom_hVar']
    # for kind in yqtValidationList :
    #     for s, sName in signDict.iteritems() :
    #         for coeff,constr in coeffDict.iteritems() :
    #             colorNumber = 0
    #             systCanvasDict[kind+s+coeff] = ROOT.TCanvas(kind+'_'+s+"_"+coeff, kind+'_'+s+"_"+coeff, 800,600)
    #             systCanvasDict[kind+s+coeff].cd()
    #             systCanvasDict[kind+s+coeff].SetGridx()
    #             systCanvasDict[kind+s+coeff].SetGridy()
    #             systLegDict[kind+s+coeff] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                
    #             #nominal
    #             systHistoDict[kind+s+coeff+'nominal_'] = ROOT.TH1F(kind+'_'+s+'_'+coeff+'_nominal_',kind+'_'+s+'_'+coeff+'_nominal_', len(unrolledYqT)-1, array('f',unrolledYqT))
    #             systHistoDict[kind+s+coeff+'nominal_'].SetLineColor(1)
    #             systHistoDict[kind+s+coeff+'nominal_'].SetTitle("coeff "+s+" " +coeff+ ", chosen deg (y, q_{T})=" +str(bestCombDict[s+coeff][0])+","+str(bestCombDict[s+coeff][1]) )
    #             systHistoDict[kind+s+coeff+'nominal_'].GetXaxis().SetTitle("Unrolled (y=fine,q_{T}=large)")
    #             if kind=='coeff_ratio_unrolled' :
    #                 systHistoDict[kind+s+coeff+'nominal_'].GetYaxis().SetTitle("Syst/Nom")
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetLineWidth(0)
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetFillColor(1)
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetFillStyle(3002)
    #             if kind=='coeff_unrolled' :
    #                 systHistoDict[kind+s+coeff+'nominal_'].GetYaxis().SetTitle("Coeff fitted value")
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetLineWidth(3)
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetFillColor(1)
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetFillStyle(3002)
    #             if kind=='pull_fNom_hVar' :
    #                 systHistoDict[kind+s+coeff+'nominal_'].GetYaxis().SetTitle("(fitted_Nom - histo_Syst) / histo_Syst_err")
    #                 systHistoDict[kind+s+coeff+'nominal_'].SetLineWidth(3)
    #             for q in qtBinning[:-1] :
    #                 for y in yBinning[:-1] :
    #                     indexUNR = yBinning.index(y)*(len(qtBinning)-1)+qtBinning.index(q)+1
    #                     nbinY = yBinning.index(y)+1
    #                     nbinQt = qtBinning.index(q)+1
    #                     if kind=='coeff_ratio_unrolled' :
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,1)
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,systHistoDict['coeffErr'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinError(nbinY,nbinQt)/systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinContent(nbinY,nbinQt))
    #                     if kind=='coeff_unrolled' :
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,systHistoDict['coeffErr'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinContent(nbinY,nbinQt))
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,systHistoDict['coeffErr'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinError(nbinY,nbinQt))
    #                     if kind=='pull_fNom_hVar' :
    #                         yVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetXaxis().GetBinCenter(nbinY)
    #                         qtVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetYaxis().GetBinCenter(nbinQt)
    #                         nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal,qtVal)
    #                         pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinContent(nbinY,nbinQt)
    #                         pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinError(nbinY,nbinQt)
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,pull_fNom_hVar)
    #                         systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,0)

    #             systHistoDict[kind+s+coeff+'nominal_'].Draw()
    #             if kind=='coeff_ratio_unrolled' :
    #                 systHistoDict[kind+s+coeff+'nominal_'].Draw('SAME E2')
    #             if kind=='coeff_unrolled' :
    #                 systHistoDict[kind+s+coeff+'nominal_'].Draw('SAME E3 L')
    #             systLegDict[kind+s+coeff].AddEntry(systHistoDict[kind+s+coeff+'nominal_'], 'nominal_')
                
    #             #variation
    #             for key, val in systDict.iteritems() :
    #                 if key=='nominal_' : continue
    #                 systHistoDict[kind+s+coeff+key] = ROOT.TH1F(kind+'_'+s+'_'+coeff+key,kind+'_'+s+'_'+coeff+key, len(unrolledYqT)-1, array('f',unrolledYqT))
    #                 systHistoDict[kind+s+coeff+key].SetLineColor(colorList[colorNumber])
    #                 systHistoDict[kind+s+coeff+key].SetTitle("coeff "+s+" " +coeff+" "+key)
    #                 systHistoDict[kind+s+coeff+key].GetXaxis().SetTitle("Unrolled (y=fine,q_{T}=large)")
    #                 if kind=='coeff_ratio_unrolled' :
    #                     systHistoDict[kind+s+coeff+key].GetYaxis().SetTitle("Syst/Nom")
    #                 if kind=='coeff_unrolled' :
    #                     systHistoDict[kind+s+coeff+'nominal_'].GetYaxis().SetTitle("Coeff fitted value")
    #                 if kind=='pull_fNom_hVar' :
    #                     systHistoDict[kind+s+coeff+'nominal_'].GetYaxis().SetTitle("(fitted_Nom - histo_Syst) / histo_Syst_err")
    #                 systHistoDict[kind+s+coeff+key].SetLineWidth(3)  
    #                 for q in qtBinning[:-1] :
    #                     for y in yBinning[:-1] :
    #                         indexUNR = yBinning.index(y)*(len(qtBinning)-1)+qtBinning.index(q)+1
    #                         nbinY = yBinning.index(y)+1
    #                         nbinQt = qtBinning.index(q)+1
    #                         yVal = systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetXaxis().GetBinCenter(nbinY)
    #                         qtVal = systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetYaxis().GetBinCenter(nbinQt)
    #                         varFitVal = systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal,qtVal)
    #                         nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal,qtVal)
    #                         if kind=='coeff_ratio_unrolled' : 
    #                             systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,varFitVal/nominalFitVal)
    #                         if kind=='coeff_unrolled' : 
    #                             systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,varFitVal)
    #                         if kind=='pull_fNom_hVar' :
    #                             pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinContent(nbinY,nbinQt)
    #                             pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinError(nbinY,nbinQt)
    #                             systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,pull_fNom_hVar)
    #                 systHistoDict[kind+s+coeff+key].Draw("SAME hist")
    #                 systLegDict[kind+s+coeff].AddEntry(systHistoDict[kind+s+coeff+key], key)
    #                 colorNumber += 1
    #             systLegDict[kind+s+coeff].Draw("SAME")
                
    validationDict = {
        'chi2_ratio_unrolled' : ['ydeg_qtdeg',"Syst/Nom",0,'SAME E2'], # kind : [loop kind, xAxis name, nominal line width, nominal draw option]
        'chi2_unrolled' : ['ydeg_qtdeg',"#chi^{2}/NDF",3,'SAME E3 L'],
        'chi2_sumOfPulls' : ['ydeg_qtdeg',"#Sigma_{y,q_{T}}(fitted_Nom - histo_Syst)^{2} / histo_Syst_err^{2}",3,''],
        'coeff_ratio_unrolled' : ['y_qt',"Syst/Nom",0,'SAME E2'],
        'coeff_unrolled' : ['y_qt',"Coeff fitted value",3,'SAME E3 L'],
        'pull_fNom_hVar' : ['y_qt',"(fitted_Nom - histo_Syst) / histo_Syst_err",3,'']
    }
    coeff_and_unpol = copy.deepcopy(coeffDict)
    coeff_and_unpol['unpol'] = ''
    
    for kind,prop in validationDict.iteritems() :
        for s, sName in signDict.iteritems() :
            for coeff,constr in coeff_and_unpol.iteritems() :
                
                #loop definition
                if prop[0] == 'ydeg_qtdeg' :
                    if coeff!='unpol' :
                        unrolledBinning =  unrolledYqTdeg
                        qtLoop = degreeQtBins_temp
                        yLoop = degreeYBins_temp
                    else : #unpol
                        unrolledBinning = unrolledYdegqT
                        qtLoop = qtBinning_unpol
                    yLoop = degreeYBins_temp
                if prop[0] == 'y_qt' :
                    if coeff!='unpol' :
                        unrolledBinning =  unrolledYqT
                        qtLoop = qtBinning
                        yLoop = yBinning
                    else : #unpol
                        unrolledBinning = unrolledYqT_unpol
                        qtLoop = qtBinning_unpol
                        yLoop = yBinning_unpol
                    
                    
                colorNumber = 0
                # replicaDrawnFlag = False #to draw only a replica
                systCanvasDict[kind+s+coeff] = ROOT.TCanvas(kind+'_'+s+"_"+coeff, kind+'_'+s+"_"+coeff, 800,600)
                systCanvasDict[kind+s+coeff].cd()
                systCanvasDict[kind+s+coeff].SetGridx()
                systCanvasDict[kind+s+coeff].SetGridy()
                systLegDict[kind+s+coeff] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                
                #nominal
                systHistoDict[kind+s+coeff+'nominal_'] = ROOT.TH1F(kind+'_'+s+'_'+coeff+'_nominal_',kind+'_'+s+'_'+coeff+'_nominal_', len(unrolledBinning)-1, array('f',unrolledBinning))
                systHistoDict[kind+s+coeff+'nominal_'].SetLineColor(1)
                systHistoDict[kind+s+coeff+'nominal_'].SetTitle(kind+" "+s+" " +coeff+" nominal_")
                systHistoDict[kind+s+coeff+'nominal_'].GetXaxis().SetTitle("Unrolled (y=fine,q_{T}=large)")
                systHistoDict[kind+s+coeff+'nominal_'].GetYaxis().SetTitle(prop[1])
                systHistoDict[kind+s+coeff+'nominal_'].SetLineWidth(prop[2])
                systHistoDict[kind+s+coeff+'nominal_'].SetFillColor(1)
                systHistoDict[kind+s+coeff+'nominal_'].SetFillStyle(3001)
                
                systHistoDict[kind+s+coeff+'nominal_'+'unfitted'] = systHistoDict[kind+s+coeff+'nominal_'].Clone(systHistoDict[kind+s+coeff+'nominal_'].GetName()+'_unfitted')
                systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetLineColor(922)#gray
                systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetFillStyle(0)
                systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetLineWidth(3)
                # if coeff=='unpol' :
                    # print qtLoop, yLoop
                for q in qtLoop[:-1] :
                    for y in yLoop[:-1] :
                        # indexUNR = yLoop.index(y)*(len(qtLoop)-1)+qtLoop.index(q)+1
                        indexUNR = qtLoop.index(q)*(len(yLoop)-1)+yLoop.index(y)+1
                        nbinY = yLoop.index(y)+1
                        nbinQt = qtLoop.index(q)+1
                        if kind=='chi2_ratio_unrolled' :
                            systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,1.)
                            systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,systHistoDict['chi2'+s+coeff+'nominal_'].GetBinError(nbinY,nbinQt)/systHistoDict['chi2'+s+coeff+'nominal_'].GetBinContent(nbinY,nbinQt))
                        if kind=='chi2_unrolled' :
                            systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,systHistoDict['chi2'+s+coeff+'nominal_'].GetBinContent(nbinY,nbinQt))
                            systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,systHistoDict['chi2'+s+coeff+'nominal_'].GetBinError(nbinY,nbinQt))
                        if kind=='chi2_sumOfPulls' :
                            if coeff != 'unpol' :
                                counterNDF=0
                                chi2_sumOfPulls=0
                                for qq in qtBinning[:-1] :
                                    for yy in yBinning[:-1] :
                                        counterNDF+=1
                                        nbinYY = yBinning.index(yy)+1
                                        nbinQtt = qtBinning.index(qq)+1
                                        yVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetXaxis().GetBinCenter(nbinYY)
                                        qtVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetYaxis().GetBinCenter(nbinQtt)
                                        nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal,qtVal)
                                        pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetBinContent(nbinYY,nbinQtt)
                                        pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetBinError(nbinYY,nbinQtt)
                                        pull_fNom_hVar = pull_fNom_hVar**2
                                        chi2_sumOfPulls+=pull_fNom_hVar
                                chi2_sumOfPulls=chi2_sumOfPulls/counterNDF   
                                systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,chi2_sumOfPulls)
                                systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,0)
                            else : #unpol
                                counterNDF=0
                                chi2_sumOfPulls=0
                                for yy in yBinning_unpol[:-1] :
                                    counterNDF+=1
                                    nbinYY = yBinning_unpol.index(yy)+1
                                    nbinQtt = qtBinning_unpol.index(q)+1
                                    yVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(nbinQtt)].GetXaxis().GetBinCenter(nbinYY)
                                    nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(nbinQtt)].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal)
                                    pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(nbinQtt)].GetBinContent(nbinYY)
                                    pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(nbinQtt)].GetBinError(nbinYY)
                                    pull_fNom_hVar = pull_fNom_hVar**2
                                    chi2_sumOfPulls+=pull_fNom_hVar
                                chi2_sumOfPulls=chi2_sumOfPulls/counterNDF   
                                systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,chi2_sumOfPulls)
                                systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,0)
                        if prop[0] == 'y_qt' :
                            if coeff != 'unpol' : 
                                yVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetXaxis().GetBinCenter(nbinY)
                                qtVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetYaxis().GetBinCenter(nbinQt)
                                nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal,qtVal) 
                            else : #unpol 
                                yVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetXaxis().GetBinCenter(nbinY)
                                nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal) 
                        if kind=='coeff_ratio_unrolled' :
                            systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,1)
                            if nominalFitVal==0 :
                                nominalFitVal = 1 
                                print "WARNING, nominalFitVal==0, set to 1", s, coeff, "qt=",q, " y=", y
                            if coeff != 'unpol' :
                                if statMC :
                                    systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,systHistoDict['coeffErr'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinError(nbinY,nbinQt)/nominalFitVal)
                                    systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetBinContent(indexUNR,systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinContent(nbinY,nbinQt)/nominalFitVal)
                                    systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetBinError(indexUNR,systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinError(nbinY,nbinQt)/nominalFitVal)
                                else :
                                    systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,abs(systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinContent(nbinY,nbinQt)-nominalFitVal)/nominalFitVal)
                            else : #unpol
                                if statMC :
                                    systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,systHistoDict['coeffErr'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinError(nbinY)/nominalFitVal)
                                    systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetBinContent(indexUNR, systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinContent(nbinY,nbinQt)/nominalFitVal)
                                    systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetBinError(indexUNR, systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinError(nbinY,nbinQt)/nominalFitVal)
                                else :
                                    systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,abs(systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinContent(nbinY)-nominalFitVal)/nominalFitVal)
                        if kind=='coeff_unrolled' :
                            systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,nominalFitVal)
                            if coeff != 'unpol' :
                                if statMC :                                
                                    systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,systHistoDict['coeffErr'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinError(nbinY,nbinQt))
                                    systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetBinContent(indexUNR,systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinContent(nbinY,nbinQt))
                                    systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetBinError(indexUNR,systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinError(nbinY,nbinQt))
                                else :
                                    systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,abs(systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinContent(nbinY,nbinQt)-nominalFitVal))
                            else : #unpol
                                if statMC :
                                    systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,systHistoDict['coeffErr'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinError(nbinY))
                                    systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetBinContent(indexUNR, systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinContent(nbinY,nbinQt))
                                    systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].SetBinError(indexUNR, systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinError(nbinY,nbinQt))
                                else :
                                    systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,abs(systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinContent(nbinY)-nominalFitVal))
                        if kind=='pull_fNom_hVar' :
                            if coeff != 'unpol' :
                                pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinContent(nbinY,nbinQt)
                                pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinError(nbinY,nbinQt)
                            else :#unpol
                                pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinContent(nbinY)
                                pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinError(nbinY)   
                            systHistoDict[kind+s+coeff+'nominal_'].SetBinContent(indexUNR,pull_fNom_hVar)
                            systHistoDict[kind+s+coeff+'nominal_'].SetBinError(indexUNR,0)        
                systHistoDict[kind+s+coeff+'nominal_'].Draw()
                if prop[3] != '' :
                    systHistoDict[kind+s+coeff+'nominal_'].Draw(prop[3])
                systLegDict[kind+s+coeff].AddEntry(systHistoDict[kind+s+coeff+'nominal_'], 'nominal_')
                if statMC :
                    if prop[0] == 'y_qt'  and "coeff" in kind:
                        systHistoDict[kind+s+coeff+'nominal_'+'unfitted'].Draw("same")
                        systLegDict[kind+s+coeff].AddEntry(systHistoDict[kind+s+coeff+'nominal_'+'unfitted'], 'unfitted')
                    
                
                #variation
                hReplicasDict = {}
                for q in qtLoop[:-1] :
                        for y in yLoop[:-1] :
                            hReplicasDict[str(q)+str(y)] = ROOT.TH1F("hReplicas"+str(q)+str(y),"hReplicas"+str(q)+str(y),10000,-10000,10000)
                for sKind, sList in systDict.iteritems() :
                    if sKind=='' : continue 
                    for key in sList : 
                        if key == '_nom_nom' : key='nominal_'
                        if key=='nominal_' : continue
                        systHistoDict[kind+s+coeff+key] = ROOT.TH1F(kind+'_'+s+'_'+coeff+'_'+key,kind+'_'+s+'_'+coeff+'_'+key, len(unrolledBinning)-1, array('f',unrolledBinning))
                        if 'Pdf' in key :
                            systHistoDict[kind+s+coeff+key].SetLineColor(30) #light green
                            # systHistoDict[kind+s+coeff+key].SetLineColorAlpha(30,0.8) #light green
                            systHistoDict[kind+s+coeff+key].SetLineWidth(1)
                        else :
                            # systHistoDict[kind+s+coeff+key].SetLineColor(colorList[colorNumber])
                            systHistoDict[kind+s+coeff+key].SetLineWidth(3)        
                        systHistoDict[kind+s+coeff+key].SetTitle(kind+" "+s+" " +coeff+" "+key)
                        systHistoDict[kind+s+coeff+key].GetXaxis().SetTitle("Unrolled (y=fine,q_{T}=large)")
                        systHistoDict[kind+s+coeff+'nominal_'].GetYaxis().SetTitle(prop[1])
                        for q in qtLoop[:-1] :
                            for y in yLoop[:-1] :
                                # indexUNR = yLoop.index(y)*(len(qtLoop)-1)+qtLoop.index(q)+1
                                indexUNR = qtLoop.index(q)*(len(yLoop)-1)+yLoop.index(y)+1
                                nbinY = yLoop.index(y)+1
                                nbinQt = qtLoop.index(q)+1
                                if kind=='chi2_ratio_unrolled' :
                                    # systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,systHistoDict['chi2'+s+coeff+key].GetBinContent(nbinY,nbinQt)/systHistoDict['chi2'+s+coeff+'nominal_'].GetBinContent(nbinY,nbinQt))
                                    BinContentValue = systHistoDict['chi2'+s+coeff+key].GetBinContent(nbinY,nbinQt)/systHistoDict['chi2'+s+coeff+'nominal_'].GetBinContent(nbinY,nbinQt)
                                if kind=='chi2_unrolled' :
                                    # systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,systHistoDict['chi2'+s+coeff+key].GetBinContent(nbinY,nbinQt))
                                    BinContentValue = systHistoDict['chi2'+s+coeff+key].GetBinContent(nbinY,nbinQt)
                                if kind=='chi2_sumOfPulls' :
                                    if coeff != 'unpol' :
                                        counterNDF=0
                                        chi2_sumOfPulls=0
                                        for qq in qtBinning[:-1] :
                                            for yy in yBinning[:-1] :
                                                counterNDF+=1
                                                nbinYY = yBinning.index(yy)+1
                                                nbinQtt = qtBinning.index(qq)+1
                                                yVal = systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetXaxis().GetBinCenter(nbinYY)
                                                qtVal = systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetYaxis().GetBinCenter(nbinQtt)
                                                nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal,qtVal)
                                                pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetBinContent(nbinYY,nbinQtt)
                                                pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(degreeQtList[degreeQtBins_temp.index(q)])].GetBinError(nbinYY,nbinQtt)
                                                pull_fNom_hVar = pull_fNom_hVar**2
                                                chi2_sumOfPulls+=pull_fNom_hVar
                                        chi2_sumOfPulls=chi2_sumOfPulls/counterNDF
                                        # systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,chi2_sumOfPulls)
                                        BinContentValue = chi2_sumOfPulls
                                    else : #unpol 
                                        counterNDF=0
                                        chi2_sumOfPulls=0
                                        for yy in yBinning_unpol[:-1] :
                                            counterNDF+=1
                                            nbinYY = yBinning_unpol.index(yy)+1
                                            nbinQtt = qtBinning_unpol.index(q)+1
                                            yVal = systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(nbinQtt)].GetXaxis().GetBinCenter(nbinYY)
                                            nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(degreeYList[degreeYBins_temp.index(y)])+str(nbinQtt)].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal)
                                            pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(nbinQtt)].GetBinContent(nbinYY)
                                            pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+key+str(degreeYList[degreeYBins_temp.index(y)])+str(nbinQtt)].GetBinError(nbinYY)
                                            pull_fNom_hVar = pull_fNom_hVar**2
                                            chi2_sumOfPulls+=pull_fNom_hVar
                                        chi2_sumOfPulls=chi2_sumOfPulls/counterNDF   
                                        # systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,chi2_sumOfPulls)
                                        BinContentValue = chi2_sumOfPulls
                                if prop[0] == 'y_qt' :
                                    if coeff != 'unpol' : 
                                        yVal = systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetXaxis().GetBinCenter(nbinY)
                                        qtVal = systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetYaxis().GetBinCenter(nbinQt)
                                        varFitVal = systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal,qtVal)
                                        nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal,qtVal) 
                                    else : #unpol 
                                        yVal = systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetXaxis().GetBinCenter(nbinY)
                                        varFitVal = systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal)
                                        nominalFitVal = systHistoDict['coeff'+s+coeff+'nominal_'+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetFunction("fit_"+sName+"_"+coeff).Eval(yVal) 
                                        # print yVal, varFitVal, nominalFitVal, nbinY, qtBinning.index(q)+1
                                if kind=='coeff_ratio_unrolled' : 
                                    if nominalFitVal==0 :
                                        nominalFitVal = 1 
                                        print "WARNING, nominalFitVal==0, set to 1", s, coeff, "qt=",q, " y=", y
                                    # print kind+s+coeff+key
                                    # systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,varFitVal/nominalFitVal)
                                    BinContentValue = varFitVal/nominalFitVal
                                if kind=='coeff_unrolled' : 
                                    # systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,varFitVal)
                                    BinContentValue = varFitVal
                                if kind=='pull_fNom_hVar' :
                                    if coeff != 'unpol' :  
                                        pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinContent(nbinY,nbinQt)
                                        pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff][0])+str(bestCombDict[s+coeff][1])].GetBinError(nbinY,nbinQt)
                                    else : #unpol 
                                        pull_fNom_hVar = nominalFitVal-systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinContent(nbinY)
                                        pull_fNom_hVar = pull_fNom_hVar/systHistoDict['coeff'+s+coeff+key+str(bestCombDict[s+coeff])+str(qtBinning_unpol.index(q)+1)].GetBinError(nbinY)
                                    # systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,pull_fNom_hVar)
                                    BinContentValue = pull_fNom_hVar
                                systHistoDict[kind+s+coeff+key].SetBinContent(indexUNR,BinContentValue)
                                if 'Pdf' in key :
                                    hReplicasDict[str(q)+str(y)].Fill(BinContentValue)
                        # if not 'Pdf' in key : 
                        #     systHistoDict[kind+s+coeff+key].Draw("SAME hist")  
                        # if 'Pdf' in key :
                        #     continue #not drawn replica lines, comment this line to add them 
                        #     systHistoDict[kind+s+coeff+key].Draw("SAME hist")
                        #     if not replicaDrawnFlag :
                        #         systLegDict[kind+s+coeff].AddEntry(systHistoDict[kind+s+coeff+key], 'PDF replicas')
                        #         replicaDrawnFlag = True
                        # else :
                        #     colorNumber += 1
                        #     systLegDict[kind+s+coeff].AddEntry(systHistoDict[kind+s+coeff+key], key)
                

                    #band hist:
                    directions = ['Up','Down']
                    if sKind=='' : continue 
                    for direc in directions :
                        systHistoDict[kind+s+coeff+'band'+sKind+direc] = ROOT.TH1F(kind+'_'+s+'_'+coeff+'_'+'band'+sKind+'_'+direc,kind+'_'+s+'_'+coeff+'_'+'band'+sKind+'_'+direc, len(unrolledBinning)-1, array('f',unrolledBinning))
                        systHistoDict[kind+s+coeff+'band'+sKind+direc].SetLineColor(groupedSystColors[sKind][0])
                        systHistoDict[kind+s+coeff+'band'+sKind+direc].SetLineWidth(3)
                        for q in qtLoop[:-1] :
                            for y in yLoop[:-1] :
                                indexUNR = qtLoop.index(q)*(len(yLoop)-1)+yLoop.index(y)+1
                                valBand = 0
                                
                                if 'Pdf' in sKind :
                                    for key in sList : 
                                        valBand+= (systHistoDict[kind+s+coeff+'nominal_'].GetBinContent(indexUNR)-systHistoDict[kind+s+coeff+key].GetBinContent(indexUNR))**2
                                if 'Scale' in sKind :
                                    for key in sList : 
                                        deltaScale = (systHistoDict[kind+s+coeff+'nominal_'].GetBinContent(indexUNR)-systHistoDict[kind+s+coeff+key].GetBinContent(indexUNR))**2
                                        if deltaScale > valBand :
                                            valBand = deltaScale
                                valBand = math.sqrt(valBand)
                        
                                if direc =='Up' :
                                    BinContentValue = systHistoDict[kind+s+coeff+'nominal_'].GetBinContent(indexUNR)+ valBand
                                if direc =='Down' :
                                    BinContentValue = systHistoDict[kind+s+coeff+'nominal_'].GetBinContent(indexUNR)- valBand
                                systHistoDict[kind+s+coeff+'band'+sKind+direc].SetBinContent(indexUNR,BinContentValue)
                        systHistoDict[kind+s+coeff+'band'+sKind+direc].Draw("SAME hist")
                    systLegDict[kind+s+coeff].AddEntry(systHistoDict[kind+s+coeff+'band'+sKind+'Up'],groupedSystColors[sKind][1]) 
                    systLegDict[kind+s+coeff].Draw("SAME")
                
                
                # #replica hist:
                # directions = ['Up','Down']
                # for direc in directions :
                #     systHistoDict[kind+s+coeff+'replicaAverage'+direc] = ROOT.TH1F(kind+'_'+s+'_'+coeff+'_'+'replicaAverage'+direc,kind+'_'+s+'_'+coeff+'_'+'replicaAverage'+direc, len(unrolledBinning)-1, array('f',unrolledBinning))
                #     systHistoDict[kind+s+coeff+'replicaAverage'+direc].SetLineColor(30)
                #     systHistoDict[kind+s+coeff+'replicaAverage'+direc].SetLineWidth(3)
                #     for q in qtLoop[:-1] :
                #         for y in yLoop[:-1] :
                #             indexUNR = qtLoop.index(q)*(len(yLoop)-1)+yLoop.index(y)+1
                #             if direc =='Up' :
                #                 BinContentValue = hReplicasDict[str(q)+str(y)].GetMean()+ hReplicasDict[str(q)+str(y)].GetStdDev()
                #             if direc =='Down' :
                #                 BinContentValue = hReplicasDict[str(q)+str(y)].GetMean()- hReplicasDict[str(q)+str(y)].GetStdDev() 
                #             BinContentErr = math.sqrt(hReplicasDict[str(q)+str(y)].GetStdDevError()**2+ hReplicasDict[str(q)+str(y)].GetMeanError()**2)
                #             systHistoDict[kind+s+coeff+'replicaAverage'+direc].SetBinContent(indexUNR,BinContentValue)
                #             systHistoDict[kind+s+coeff+'replicaAverage'+direc].SetBinError(indexUNR,BinContentErr)
                #     systHistoDict[kind+s+coeff+'replicaAverage'+direc].Draw("SAME hist")
                #     systLegDict[kind+s+coeff].AddEntry(systHistoDict[kind+s+coeff+'replicaAverage'+direc], 'PDF replicas, '+ direc) 
                
                # systLegDict[kind+s+coeff].Draw("SAME")
    
    
    if saveRoot :
        output = ROOT.TFile(output_dir+"/systComparison_"+SUFFIX+".root", "recreate")
        output.cd()
        for k,c in systCanvasDict.iteritems() :
            c.Write()
    if saveFig :
        if not os.path.isdir(output_dir+'/systComparison'): os.system('mkdir '+output_dir+'/systComparison')
        for k,c in systCanvasDict.iteritems() :
            c.SaveAs(output_dir+'/systComparison/'+c.GetName()+'.png')

  
            
            
    
def runSignleFit(par) :
    # fit_coeff_call = fit_coeff(fname=par[0], charge=par[1],coeff=par[2],constraint=par[3],degreeY=par[4],degreeQt=par[5],Map01=par[6],y_max=par[7],qt_max=par[8] )    
    fit_coeff_call = fit_coeff_poly(fname=par[0], charge=par[1],coeff=par[2],constraint=par[3],degreeY=par[4],degreeQt=par[5],Map01=par[6],y_max=par[7],qt_max=par[8] )    
                  
    ss = signDict.keys()[signDict.values().index(par[1])] #find the key relative to "charge"
    fit_coeff_call['KEYNAME'] = str(par[4])+str(par[5])+ss+par[2]
    
    del fit_coeff_call['fit'+'Res'+par[1]+par[2]]
    fit_coeff_call['fit'+'Res'+par[1]+par[2]] = fit_coeff_call['fit'+par[1]+par[2]].GetFunction("fit_"+par[1]+"_"+par[2])
    
    return fit_coeff_call


if __name__ == "__main__":
    
    print "WARNING: chosen Degree = N --> Chebyshev polynomial max = N-1 (or standard poly. max = N-1)"
    
    if not SYST_VALIDATION : #analysis part
    
        cumulativeDict = {} 
        valDict ={}   
        
        #do the analysis
        print "analysis started..."
        processesMain = []
        for yDeg in degreeYList :
            for qtDeg in degreeQtList :
                for s,sname in signDict.iteritems() :
                    for coeff,constr in coeffDict.iteritems() :
                        if not VALONLY :
                            if MULTICORE : 
                                processesMain.append((IN_GEN, sname,coeff,constr,yDeg,qtDeg,MAP01,2.4,60.))
                            else :
                                # tempDict = fit_coeff(fname=IN_GEN, charge=sname,coeff=coeff,constraint=constr,degreeY=yDeg,degreeQt=qtDeg,Map01=MAP01,y_max=2.4,qt_max=32 )
                                tempDict = fit_coeff_poly(fname=IN_GEN, charge=sname,coeff=coeff,constraint=constr,degreeY=yDeg,degreeQt=qtDeg,Map01=MAP01,y_max=2.4,qt_max=60. )
                        else : #skip analyisis, do validation only
                            tempDict = rebuildHistoDict(charge=sname,coeff=coeff,degreeY=yDeg,degreeQt=qtDeg)    
                        if not MULTICORE : cumulativeDict[str(yDeg)+str(qtDeg)+s+coeff] = tempDict
        if MULTICORE :
            poolMain = multiprocessing.Pool(NCORES)
            poolResultList = poolMain.map(runSignleFit,processesMain)
            poolMain.close()
            poolMain.join()

            for poolResult in poolResultList :
                cumulativeDict[poolResult['KEYNAME']] = copy.deepcopy(poolResult)
                del cumulativeDict[poolResult['KEYNAME']]['KEYNAME']
                                
                            
        for yDeg in degreeYList :
            for s,sname in signDict.iteritems() :
                # tempDict = fit_unpol(fname=IN_GEN,coeff='AUL',charge=sname,degreeY=yDeg,Map01=True,y_max=2.4,qt_max=32.) #Map01 should be ALWAYS true to remove odd degree easily
                tempDict = fit_unpol_poly(fname=IN_GEN,coeff='AUL',charge=sname,degreeY=yDeg,Map01=True,y_max=2.4,qt_max=60.) #Map01 should be ALWAYS true to remove odd degree easily
                cumulativeDict[str(yDeg)+s+'unpol'] = tempDict
            
        #validation of the fit
        print "validation...."
        valDict = Validation(hdict=cumulativeDict, signDict=signDict, coeffDict=coeffDict,degreeYList=degreeYList,degreeQtList=degreeQtList)

        #writing
        print "writing..."
        if not os.path.isdir(output_dir): os.system('mkdir '+output_dir)
        output = ROOT.TFile(output_dir+"/regularizationFit_"+SUFFIX+".root", "recreate")
        output.cd()
        
        #general histos
        validationList=  ['h2_chi2','h2_pullsSigma', 'h2_pullsSigmaFit', 'pval_alongY','pval_alongQt','pval_alongDiagonal','chosenDeg_path','chosenDeg_search']
        for s,sname in signDict.iteritems() :
            for coeff,constr in coeffDict.iteritems() :
                for var in validationList :
                    if 'chosenDeg' in var :
                        if s=='Plus' and coeff=='A0':
                            valDict[var].Write()   
                    else :
                        valDict[var+s+coeff].Write()
        
        unpolValList = ['h2_chi2','pval_alongY', 'chosen_unpol']
        for s,sname in signDict.iteritems() :
            for var in unpolValList :
                if 'chosen_unpol' in var :
                    if s=='Plus' :
                        valDict[var].Write()
                else :
                    valDict[var+s+'unpol'].Write()
                
        #degree dependent histos                
        for yDeg in degreeYList :
            for qtDeg in degreeQtList :
                direc = output.mkdir('Y'+str(yDeg)+'_Qt'+str(qtDeg))
                direc.cd()
                for s,sname in signDict.iteritems() :
                    for coeff,constr in coeffDict.iteritems() :
                        for ind,val in cumulativeDict[str(yDeg)+str(qtDeg)+s+coeff].iteritems() :
                            if 'Res' in ind : continue #do not write the fit result (used in validation)
                            # print "ind", ind, SYST_kind
                            if SYST_kind != '' : #nominal
                                if 'fitErr' in ind : 
                                    'fitErr skipped in', SYST_kind, SYST_name
                                    continue 
                            # print "val",  val
                            val.Write()
                for ind, val in valDict.iteritems() :
                    if 'h2_chi2' in ind : continue
                    if str(yDeg)+str(qtDeg) not in ind : continue 
                    val.Write()
                    
            #unpol histos    
            direc = output.mkdir('unpol_Y'+str(yDeg))
            direc.cd()
            for s,sname in signDict.iteritems() :    
                for ind,val in cumulativeDict[str(yDeg)+s+'unpol'].iteritems() :  
                    if 'Res' in ind : continue
                    val.Write()
        
        print "saving figures"
        # figureSaver(cumulativeDict = cumulativeDict, valDict=valDict, outputname='plots_'+SUFFIX, signDict=signDict, coeffDict=coeffDict,degreeYList=degreeYList,degreeQtList=degreeQtList, validationList = validationList, unpolValList=unpolValList)
        
    if SYST_VALIDATION : #systematic comparison
        
        # bestCombDict = { #OLD FW
        #     'PlusA0' : [3,4],
        #     'PlusA1' : [5,4],
        #     'PlusA2' : [2,4],
        #     'PlusA3' : [5,3],
        #     'PlusA4' : [7,3],
        #     'MinusA0' : [2,5],
        #     'MinusA1' : [5,4],
        #     'MinusA2' : [2,3],
        #     'MinusA3' : [7,3],
        #     'MinusA4' : [5,3],
        #     'Plusunpol' : 3,
        #     'Minusunpol' : 3
        # }
        # bestCombDict = { #[y,qt] #CHEBYC
        #     'PlusA0' : [2,4],
        #     'PlusA1' : [5,5],
        #     'PlusA2' : [2,3],
        #     'PlusA3' : [5,4],
        #     'PlusA4' : [6,4],
        #     'MinusA0' : [3,4], 
        #     'MinusA1' : [5,4], 
        #     'MinusA2' : [2,3], 
        #     'MinusA3' : [6,3], 
        #     'MinusA4' : [5,4], 
        #     'Plusunpol' : 3,
        #     'Minusunpol' : 5 
        # }
        bestCombDict = { #[y,qt] #poly
            'PlusA0' : [3,4],
            'PlusA1' : [4,4],
            'PlusA2' : [2,5],
            'PlusA3' : [3,4],
            'PlusA4' : [4,7],
            'MinusA0' : [3,4], 
            'MinusA1' : [3,4], 
            'MinusA2' : [2,4], 
            'MinusA3' : [3,3], 
            'MinusA4' : [4,6], 
            'Plusunpol' : 3,
            'Minusunpol' : 3 
        }

        
        print "WARNING: harcoded best combination for syst comparsion"
        print "WARNING: hardcoded outdated combination"
    
        systComparison(bestCombDict=bestCombDict, saveRoot=True,saveFig=False,excludeSyst='',statMC=True)