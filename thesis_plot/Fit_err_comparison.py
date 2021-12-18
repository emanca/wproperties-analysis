import ROOT
import ROOT
import math
from array import array
from ROOT import gStyle, gPad
ROOT.gROOT.SetBatch(True)
gStyle.SetOptStat(0)

#------------------------- init
# inDir = '../../../Fit/OUTPUT_1Apr_noCF_qtMax60/'
inDir = '../../../Fit/OUTPUT_25Apr_noCF_qtMax60_fixLowAcc/'

fitDict = {
    'fit' :     ['fullFit_asimov_noBBB_noFittedMass/',          'fitPlots_Wplus_angCoeff.root'        ,'FitRatioAC',     'coeff_UNRyqt_RatioAC'    , ROOT.kBlack, 'Non-regularized', 20],
    # 'gauss' :   ['full_fit_asimov_tau100_noBBB_noFittedMass/',  'fitPlots_Wplus_angCoeff.root'        ,'FitRatioAC',     'coeff_UNRyqt_RatioAC'    , ROOT.kRed-4, 'Simult. Gaussian constraint regul.',22],
    # 'gauss' :   ['full_fit_asimov_tau500_noBBB_noFittedMass/',  'fitPlots_Wplus_angCoeff.root'        ,'FitRatioAC',     'coeff_UNRyqt_RatioAC'    , ROOT.kRed-4, 'Simult. Gaussian constraint regul.',22],
    'poly' :    ['fullFit_asimov_noBBB_noFittedMass/',          'fitPlots_Wplus_APO_angCoeff.root'    , 'apoRatio' ,     'apoRatio_UNRyqt'    , ROOT.kAzure+1, 'Post-fit regularized',21],
}

coeffList = ['A0','A1','A2','A3','A4', 'unpolarizedxsec']
coeffDict = {
                'A0' :  'A_{0}',
                'A1' : 'A_{1}',
                'A2' :  'A_{2}',
                'A3' :  'A_{3}',
                'A4' : 'A_{4}',
                'unpolarizedxsec' : '#sigma^{U+L}'
            }
catList = ['qt', 'y', 'UNRyqt']

yArr = [0, 0.4, 0.8, 1.2 ,1.6, 2.0, 2.4]
qtArr =[0.,  2.,  4., 6.,  8., 10., 12., 16., 20., 26., 36., 60.,]
dimTot = (len(yArr)-1)*(len(qtArr)-1)
unrolledYQt= list(qtArr)
intervalQtBin = []
for q in qtArr[:-1] :
    intervalQtBin.append(qtArr[qtArr.index(q)+1]-qtArr[qtArr.index(q)])
for y in range(len(yArr)-2) :
    for q in intervalQtBin :
        unrolledYQt.append(unrolledYQt[-1]+q)

cmslab = "#bf{CMS} #scale[0.7]{#it{Simulation Work in progress}}"
lumilab = " #scale[0.7]{35.9 fb^{-1} (13 TeV)}"
cmsLatex = ROOT.TLatex()

#----------------

histoDict = {}
canDict = {}
legDict ={}
filedict = {}
#get histos:
# i=0
for f,fval in fitDict.items() :
    filedict[f] = ROOT.TFile.Open(inDir+fval[0]+fval[1])
    for c in coeffList :
        for cat in catList :
            if 'UNR' in cat :
                hh = filedict[f].Get('coeff2D_reco/reco_'+fval[3]+c)
                histoDict[f+c+cat] = ROOT.TH1F('UNRyqt'+c+f,'UNRyqt'+c+f,len(unrolledYQt)-1, array('f',unrolledYQt))
            else :
                hh = filedict[f].Get('coeff2D_reco/'+fval[2]+cat+c)
                histoDict[f+c+cat] = hh.Clone(f+'_'+c+'_'+cat)
            # can.GetPrimitive('Impact_unpolarizedxsec_tot_qt')
            # stack = can.GetPrimitive('reco_impact_'+cat+c)
            # hists = stack.GetHists()
            # hists.At(5)
            for x in range(1, hh.GetNbinsX()+1) : #error only, ingored bias
                histoDict[f+c+cat].SetBinContent(x,hh.GetBinError(x))
                histoDict[f+c+cat].SetBinError(x,0)
            histoDict[f+c+cat].SetLineColor(fval[4])
            histoDict[f+c+cat].SetLineWidth(3)
            histoDict[f+c+cat].SetMarkerColor(fval[4])
            histoDict[f+c+cat].SetMarkerStyle(fval[6])
            histoDict[f+c+cat].SetMarkerSize(2)
            histoDict[f+c+cat].SetFillStyle(0)
    # i=i+1

#histo of ratio   
for c in coeffList :
    for cat in catList :
        canDict[c+cat] = ROOT.TCanvas('fit_err_comp_'+c+'_'+cat,'fit_err_comp_'+c+'_'+cat,1600,1200)
        # legDict[c+cat] = ROOT.TLegend(0.25,0.7,0.75,0.9)
        legDict[c+cat] = ROOT.TLegend(0.4,0.7,0.75,0.9)
        legDict[c+cat].SetLineWidth(0)
        legDict[c+cat].SetFillStyle(0)
        canDict[c+cat].cd()
        canDict[c+cat].SetGridy()
        if 'UNR' not in cat : 
            canDict[c+cat].SetGridx()
        
        sameflag=''
        for f,fval in fitDict.items() :
            histoDict[f+c+cat].Draw("lp hist"+sameflag)
            sameflag = 'same'
            legDict[c+cat].AddEntry(histoDict[f+c+cat], fval[5])
            
            
            
            # histoDict[f+c+cat].GetXaxis().SetLabelSize(0.02)   
            # histoDict[f+c+cat].GetXaxis().SetTitleSize(0.02)   
            # histoDict[f+c+cat].GetYaxis().SetLabelSize(0.02)   
            # histoDict[f+c+cat].GetYaxis().SetTitleSize(0.02) 
            # histoDict[f+c+cat].SetMaximum(0)
            # histoDict[f+c+cat].SetMinimum(0)
            
        histoDict['fit'+c+cat].SetMinimum(0)
        histoDict['fit'+c+cat].SetMaximum(histoDict['fit'+c+cat].GetMaximum()*1.05)
        if cat=='y' and 'unpol' in c : #root is a shit
            histoDict['fit'+c+cat].SetMaximum(0.05)
        # print(c,cat, histoDict['fit'+c+cat].GetMaximum()*1.1)
        # histoDict['fit'+c+cat].SetMaximum(0.2)
        histoDict['fit'+c+cat].GetYaxis().SetRangeUser(0,histoDict['fit'+c+cat].GetMaximum()*1.05)
        histoDict['fit'+c+cat].GetYaxis().SetTitleSize(0.05)
        histoDict['fit'+c+cat].GetYaxis().SetTitleOffset(0.98)
        histoDict['fit'+c+cat].GetYaxis().SetLabelSize(0.04)
        histoDict['fit'+c+cat].GetYaxis().SetTitle('Absolute uncertainity on '+coeffDict[c])
        if 'unpol' in c :
            histoDict['fit'+c+cat].GetYaxis().SetTitle('Relative uncertainity on '+coeffDict[c])
        histoDict['fit'+c+cat].GetXaxis().SetTitleSize(0.04)
        histoDict['fit'+c+cat].GetXaxis().SetTitleOffset(1)
        # histoDict['fit'+c+cat].GetXaxis().SetLabelSize(0.04)
        # histoDict['fit'+c+cat].SetLabelSize(0.04,'x')                    
        histoDict['fit'+c+cat].GetXaxis().SetLabelSize(0.04)
        
        if 'UNR' in cat :
            canDict[c+cat].SetBottomMargin(0.15)
            histoDict['fit'+c+'UNRyqt'].SetTitle('')
            histoDict['fit'+c+'UNRyqt'].GetXaxis().SetTickSize(0)
            histoDict['fit'+c+'UNRyqt'].GetXaxis().SetTitle('fine binning: q_{T}^{W} 0#rightarrow 60 GeV')
            histoDict['fit'+c+'UNRyqt'].GetXaxis().SetTitleOffset(1.9)
            for y in yArr[:-1] :
                for q in qtArr[:-1] :
                    indexUNRyqt = yArr.index(float(y))*(len(qtArr)-1)+qtArr.index(float(q))
                    if qtArr.index(q)==0 :
                        histoDict['fit'+c+'UNRyqt'].GetXaxis().SetNdivisions(-1)
                        histoDict['fit'+c+'UNRyqt'].GetXaxis().SetBinLabel(indexUNRyqt+1,"|Y_{W}|#in[%.1f,%.1f]" % (y, yArr[yArr.index(y)+1]))
                        histoDict['fit'+c+'UNRyqt'].GetXaxis().ChangeLabel(indexUNRyqt+1,340,0.035)    
                        histoDict['fit'+c+'UNRyqt'].GetXaxis().LabelsOption("d")
                        histoDict['fit'+c+'UNRyqt'].GetXaxis().SetLabelOffset(0.008)
        
                
        # if 'UNR' in cat :
        #     # for f,fval in fitDict.items() : 
        #     for x in range(1, histoDict['fit'+c+'UNRyqt'].GetNbinsX()+1) :
        #             print(histoDict['fit'+c+'UNRyqt'].GetXaxis().GetBinLabel(x))
        #             histoDict['fit'+c+'UNRyqt'].GetXaxis().SetBinLabel(x,"")
        #             histoDict['fit'+c+'UNRyqt'].GetXaxis().SetTickSize(0)
                    # histoDict['fit'+c+'UNRyqt'].GetXaxis().ChangeLabel(x,340,0.01)
                    # histoDict['fit'+c+'UNRyqt'].GetXaxis().SetBinLabel(x,"a")
            # axx = ROOT.TGAxis(dimTot,0,dimTot)
            # axx = ROOT.TGAxis(0,0,dimTot,0)
            
                    # 
            
            # for x in range(1, histoDict['fit'+c+cat].GetNbinsX()+1) :
            #         histoDict[f+c+cat].GetXaxis().SetBinLabel(x,"1")
            # histoDict[f+c+cat].GetXaxis().SetLabelOffset(0)
            # [suff+'FitAC'+'UNRyqt'+c].GetXaxis().SetBinLabel(indexUNRyqt+1,"|Y_{W}|#in[%.1f,%.1f]" % (y, self.yArr[self.yArr.index(y)+1]))


        topY4lineErr = histoDict['fit'+c+cat].GetMaximum()*1.05
        if 'UNR' in cat :
            bottomY4lineErr = 0
            if 'A4' in c :
                topY4lineErr = 1.975
            elif 'unpol' in c :
                topY4lineErr = 0.375
            
            vErrLines = []
            for yqt in unrolledYQt :
                if yqt%(qtArr[-1])==0 :
                    if yqt==0 :continue
                    if yqt==len(unrolledYQt)-1 : continue 
                    vErrLines.append(ROOT.TLine(yqt-2,bottomY4lineErr,yqt-2,topY4lineErr))
                    vErrLines[-1].SetLineStyle(ROOT.kDotted)
                    vErrLines[-1].DrawLine(yqt-2,bottomY4lineErr,yqt-2,topY4lineErr)
        
        legDict[c+cat].Draw("same")
        
        cmsLatex.SetNDC()
        cmsLatex.SetTextFont(42)
        cmsLatex.SetTextColor(ROOT.kBlack)
        cmsLatex.SetTextAlign(31) 
        cmsLatex.DrawLatex(1-canDict[c+cat].GetRightMargin(),1-0.8*canDict[c+cat].GetTopMargin(),lumilab)
        cmsLatex.SetTextAlign(11) 
        cmsLatex.DrawLatex(canDict[c+cat].GetLeftMargin(),1-0.8*canDict[c+cat].GetTopMargin(),cmslab)

for c in coeffList :
    for cat in catList :
        if 'UNR' not in cat : continue 
        for x in range(1, histoDict['fit'+c+cat].GetNbinsX()+1) :
                    histoDict[f+c+cat].GetXaxis().SetBinLabel(x,"")
        
        
               


outfile = ROOT.TFile('Fit_err_comparison.root',"recreate")
outfile.cd()
for c in coeffList :
    for cat in catList :
       canDict[c+cat].Write()
       if 'unpol' in c or 'A4' in c : 
            canDict[c+cat].SaveAs('Fit_err_comparison_'+c+'_'+cat+'_plus_nogauss.png')
            canDict[c+cat].SaveAs('Fit_err_comparison_'+c+'_'+cat+'_plus_nogauss.pdf')
            