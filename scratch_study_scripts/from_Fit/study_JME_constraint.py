import ROOT
import math
ROOT.gStyle.SetOptStat(0)
import copy
from array import array

# path = '../OUTPUT_1Apr_noCF_qtMax60/'
path = '../OUTPUT_25Apr_noCF_qtMax60_fixLowAcc/'
# path = '../test_1Apr_noCF_qtMax80/'
inputShapeFile = ROOT.TFile.Open(path+'/Wplus_reco.root')
inputFitFile = ROOT.TFile.Open(path+'/fullFit_asimov_noBBB_noFittedMass/fit_Wplus_reco.root')
# inputFitFile = ROOT.TFile.Open(path+'/fullFit_noFittedMass_noBBB/fit_Wplus_reco.root')



coeffList = ["L","I","T","A","P","7","8","9","UL"]
# coeffList = ["L","I","T","A","P","UL"]
bkgList = ["DYJets","DiBoson","Top","Fake","WtoTau","LowAcc"] 
# bkgList = ["DYJets","DiBoson","Top","Fake","WtoTau"] 
# print("warning: no low acc")
# bkgList = ["DYJets","DiBoson","Top","WtoTau","LowAcc"] 
# bkgList = ["DYJets","DiBoson","Top","WtoTau","Fake"] 
# bkgList = ["Fake","LowAcc"]
# bkgList = ["DYJets","Fake",'DiBoson']
# bkgList = []
procList = copy.deepcopy(coeffList)
procList.extend(bkgList)
qtList = ['1','2','3','4','5','6','7','8','9','10','11']
# qtList = ['1','2','3','4','5','6','7','8','9','10','11','12','13']
yList = ['1','2','3','4','5','6']
syst = {
    'nom' : [''],
    'JME' : ['_jesTotalUp','_jesTotalDown'],
    'ptscale' : ['_correctedUp', '_correctedDown'],
    'JME2' : ['_unclustEnUp','_unclustEnDown'],
}



#fit file ...............
hdict = {}
hpostfittemp = inputFitFile.Get('expfull_postfit')

ptArr = []
etaArr = []
for ipt in range(0,61) :
    ptArr.append(25. + ipt*(55.-25.)/60.)
for ieta in range(0,49) :
    etaArr.append(-2.4 + ieta * (4.8) / 48.)

unrolledEtaPt= list(ptArr)
intervalPtBin = []
for q in ptArr[:-1] :
    intervalPtBin.append(ptArr[ptArr.index(q)+1]-ptArr[ptArr.index(q)])
for y in range(len(etaArr)-2) :
    for q in intervalPtBin :
        unrolledEtaPt.append(unrolledEtaPt[-1]+q)


# hpostfittemp = hdict['postfit'].Clone('temp')
hdict['postfit'] = ROOT.TH1F(hpostfittemp.GetName()+'_', hpostfittemp.GetTitle()+'_', len(unrolledEtaPt)-1, array('f',unrolledEtaPt))
for ietapt in range(1,hdict['postfit'].GetNbinsX()+1) :
    hdict['postfit'].SetBinContent(ietapt, hpostfittemp.GetBinContent(ietapt))
    hdict['postfit'].SetBinError(ietapt, hpostfittemp.GetBinError(ietapt))

hdict['postfit'+'noerr'] = hdict['postfit'].Clone('postfit'+'noerr')
for x in range(1, hdict['postfit'+'noerr'].GetNbinsX()+1) :
    hdict['postfit'+'noerr'].SetBinError(x,0)

hdict['postfit'+'ratio'] = hdict['postfit'].Clone(hdict['postfit'].GetName()+'ratio')
hdict['postfit'+'ratio'].Divide(hdict['postfit'+'noerr'])

print("get")
for s,sList in syst.items() :
    for sName in sList :
        for c in coeffList :
            for qt in qtList :
                for y in yList :
                    # print("helXsecs"+c+"_y_"+y+"_qt_"+qt+sName)
                    hdict[c+qt+y+sName] = inputShapeFile.Get("helXsecs"+c+"_y_"+y+"_qt_"+qt+sName)
                    # print(hdict[c+qt+y+sName].Integral())
        for b in bkgList :
            hdict[b+sName] = inputShapeFile.Get(b+sName)
            # print(hdict[b+sName].Integral())
    
print("add")
for s,sList in syst.items() :
    for sName in sList :
        # hdict['tot'+sName] = hdict['LowAcc'+sName].Clone("tot"+sName)   
        hdict['tot'+sName] = hdict['L'+'1'+'1'+sName].Clone("tot"+sName)   
        for c in coeffList :
            for qt in qtList :
                for y in yList :
                    if c=='L' and qt=='1' and y=='1' : continue
                    # print("helXsecs"+c+"_y_"+y+"_qt_"+qt+sName)
                    hdict['tot'+sName].Add(hdict[c+qt+y+sName])
                    
        for b in bkgList :
            # if b=='LowAcc' : continue 
            hdict['tot'+sName].Add(hdict[b+sName])



print("unroll")


for s,sList in syst.items() :
    for sName in sList :
        hdict['unr'+sName] = ROOT.TH1F('unr'+sName,'unr'+sName, len(unrolledEtaPt)-1, array('f',unrolledEtaPt))
        for ieta in range(1,hdict['tot'+sName].GetNbinsX()+1) :
            for ipt in range(1,hdict['tot'+sName].GetNbinsY()+1) :
                indexUNRetpt = (ieta-1)*(len(ptArr)-1)+(ipt-1)
                hdict['unr'+sName].SetBinContent(indexUNRetpt+1,hdict['tot'+sName].GetBinContent(ieta,ipt))
                hdict['unr'+sName].SetBinError(indexUNRetpt+1,hdict['tot'+sName].GetBinError(ieta,ipt))



print("ratio") 
hdict['unr'+'noerr'] = hdict['unr'].Clone('unr'+'noerr')
for x in range(1, hdict['unr'+'noerr'].GetNbinsX()+1) :
    hdict['unr'+'noerr'].SetBinError(x,0)
for s,sList in syst.items() :
    for sName in sList :
        hdict['unrratio'+sName] = hdict['unr'+sName].Clone('unr_ratio'+sName)
        hdict['unrratio'+sName].Divide(hdict['unr'+'noerr'])



#plotting

c_compSum_jes = ROOT.TCanvas("c_compSum_jes","c_compSum_jes", 800,600)
c_compSum_jes.cd()
hdict['unrratio'].SetLineColor(2)
hdict['unrratio'].Draw()
hdict['unrratio'+'_jesTotalUp'].Draw("histsame")
hdict['unrratio'+'_jesTotalDown'].SetLineColor(3)
hdict['unrratio'+'_jesTotalDown'].Draw("histsame")

c_compSum_abs_jes = ROOT.TCanvas("c_compSum_abs_jes","c_compSum_abs_jes", 800,600)
c_compSum_abs_jes.cd()
hdict['unr'].SetLineColor(2)
hdict['unr'].Draw()
hdict['unr'+'_jesTotalUp'].Draw("histsame")
hdict['unr'+'_jesTotalDown'].SetLineColor(3)
hdict['unr'+'_jesTotalDown'].Draw("histsame")

c_compPostfit_jes = ROOT.TCanvas("c_compPostfit_jes","c_compPostfit_jes", 800,600)
c_compPostfit_jes.cd()
hdict['postfitratio'].SetLineColor(2)
hdict['postfitratio'].Draw()
hdict['unrratio'+'_jesTotalUp'].Draw("histsame")
hdict['unrratio'+'_jesTotalDown'].SetLineColor(3)
hdict['unrratio'+'_jesTotalDown'].Draw("histsame")

c_compPostfit_pt = ROOT.TCanvas("c_compPostfit_pt","c_compPostfit_pt", 800,600)
c_compPostfit_pt.cd()
hdict['postfitratio'].SetLineColor(2)
hdict['postfitratio'].Draw()
hdict['unrratio'+'_correctedUp'].Draw("histsame")
hdict['unrratio'+'_correctedDown'].SetLineColor(3)
hdict['unrratio'+'_correctedDown'].Draw("histsame")


c_compPostfit_abs_jes = ROOT.TCanvas("c_compPostfit_abs_jes","c_compPostfit_abs_jes", 800,600)
c_compPostfit_abs_jes.cd()
hdict['postfit'].SetLineColor(2)
hdict['postfit'].Draw()
hdict['unr'+'_jesTotalUp'].Draw("histsame")
hdict['unr'+'_jesTotalDown'].SetLineColor(3)
hdict['unr'+'_jesTotalDown'].Draw("histsame")

c_compPostfit_abs_pt = ROOT.TCanvas("c_compPostfit_abs_pt","c_compPostfit_abs_pt", 800,600)
c_compPostfit_abs_pt.cd()
hdict['postfit'].SetLineColor(2)
hdict['postfit'].Draw()
hdict['unr'+'_correctedUp'].Draw("histsame")
hdict['unr'+'_correctedDown'].SetLineColor(3)
hdict['unr'+'_correctedDown'].Draw("histsame")



print("write")  
outfile = ROOT.TFile("study_JME_constraint.root","recreate")
outfile.cd()
for s,sList in syst.items() :
    for sName in sList :
        hdict['tot'+sName].Write()
        hdict['unr'+sName].Write()
        hdict['unrratio'+sName].Write()
hdict['postfit'].Write()
hdict['postfit'+'ratio'].Write()

c_compSum_jes.Write()
c_compSum_abs_jes.Write()
c_compPostfit_jes.Write()
c_compPostfit_pt.Write()
c_compPostfit_abs_jes.Write()
c_compPostfit_abs_pt.Write()
                



            
        

             
            

