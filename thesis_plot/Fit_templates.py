import ROOT
import math
from ROOT import gStyle, gPad
ROOT.gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptStat(0)


hDict = {}
canDict = {}

fileDict = {}
fileDict['plus'] = ROOT.TFile.Open('/scratch/bertacch/wmass/wproperties-analysis/Fit/OUTPUT_1Apr_noCF_qtMax60/Wplus_reco.root')
fileDict['minus'] = ROOT.TFile.Open('/scratch/bertacch/wmass/wproperties-analysis/Fit/OUTPUT_1Apr_noCF_qtMax60/Wminus_reco.root')

xsecList = ["UL","L","I","T","A","P"]
# coeffList = ["U+L","A_{0}","A_{1}","A_{2}","A_{3}","A_{4}"]
kindSig = ['ylow'+'qtlow', 'ylow'+'qthigh','yhigh'+'qtlow','yhigh'+'qthigh']
kindSigLabel = ['#splitline{|Y|<0.4}{q_{T}<2 GeV}', '#splitline{|Y|<0.4}{16 GeV<q_{T}<20 GeV}','#splitline{2<|Y|<2.4}{q_{T}<2 GeV}','#splitline{2<|Y|<2.4}{16 GeV<q_{T}<20 GeV}']

nameDict ={
    "signalUL" : '#sigma_{U+L}',
    "signalL" : '#sigma_{L} (A_{0})',
    "signalI" : '#sigma_{I} (A_{1})',
    "signalT" : '#sigma_{T} (A_{2})',
    "signalA" : '#sigma_{A} (A_{3})',
    "signalP" : '#sigma_{P} (A_{4})',
    "WtoTau" :  'W#rightarrow#tau#nu',
    "Top"   : 'Top',
    "DiBoson" : 'Diboson',
    "DY"        : 'Drell-Yan',
    "QCD"        : 'QCD',
    "LowAcc" : 'Low-acceptance',
    "Data" : 'Data',
    
}

coordDict = {
    "signalUL" : [0.001,0.88],
    "signalL" :  [0.001,0.72],
    "signalI" :  [0.001,0.56],
    "signalT" :  [0.001,0.40],
    "signalA" :  [0.001,0.24],
    "signalP" :  [0.001,0.08],
    'ylow'+'qtlow' : [0.18,0.97],
    'ylow'+'qthigh': [0.35,0.97],
    'yhigh'+'qtlow': [0.60,0.97],
    'yhigh'+'qthigh' : [0.78,0.97],
}


# get histograms
for s,f in fileDict.items() :
    for xs in xsecList :
        hDict[s+'signal'+xs+'ylow'+'qtlow'] =f.Get('helXsecs'+xs+'_y_1_qt_1') 
        hDict[s+'signal'+xs+'ylow'+'qthigh'] =f.Get('helXsecs'+xs+'_y_1_qt_8') 
        hDict[s+'signal'+xs+'yhigh'+'qtlow'] =f.Get('helXsecs'+xs+'_y_6_qt_1') 
        hDict[s+'signal'+xs+'yhigh'+'qthigh'] =f.Get('helXsecs'+xs+'_y_6_qt_8') 

    hDict[s+'Top'] = f.Get('Top') 
    hDict[s+'DiBoson'] = f.Get('DiBoson') 
    hDict[s+'WtoTau'] = f.Get('WtoTau') 
    hDict[s+'DY'] = f.Get('DYJets') 
    hDict[s+'QCD'] = f.Get('Fake') 
    hDict[s+'LowAcc'] = f.Get('LowAcc') 
    hDict[s+'Data'] = f.Get('Data') 
    
    

#----------------------- output -----------------------#
    
#output
outfile = ROOT.TFile('Fit_templates.root',"recreate")
# outfile.cd()
# for s,f in fileDict.items() :  
#     canDict[s+'signal'].Write()  
#     canDict[s+'signal'+'extra'].Write()
#     canDict[s+'signal'+'extra'].SaveAs("Fit_templates_signal"+s+".png")
#     canDict[s+'signal'+'extra'].SaveAs("Fit_templates_signal"+s+".pdf")
#     for n,name in nameDict.items() :
#         if 'signal' in n : continue 
#         canDict[s+n].Write()
#         canDict[s+n].SaveAs("Fit_templates_"+n+"_"+s+".png")
#         canDict[s+n].SaveAs("Fit_templates_"+n+"_"+s+".pdf")
        
    
# outfile.Close()   


#--------------------------background-----------------------------#
for hi,hval in hDict.items() :
    if 'signal' in hi : continue 
    if 'plus' in hi : ch = 'W^{+}'
    if 'minus' in hi : ch = 'W^{-}'
    for n,name in nameDict.items() :
        if n in hi : 
            hval.SetTitle(name+', '+ch)
    
    hval.SetStats(0)
    hval.GetYaxis().SetTitle('p_{T}^{#mu}')  
    hval.GetXaxis().SetTitle('#eta^{#mu}')
    hval.GetZaxis().SetTitle("Events")  
    hval.GetZaxis().SetTitleOffset(1.4)
    
    canDict[hi] = ROOT.TCanvas("c_"+hi,"c_"+hi, 1600,1200)
    canDict[hi].cd()
    canDict[hi].SetGridx()
    canDict[hi].SetGridy()
    canDict[hi].SetRightMargin(0.15)
    hval.Draw("colz")
    
outfile.cd()
for s,f in fileDict.items() :  
    for n,name in nameDict.items() :
        if 'signal' in n : continue 
        canDict[s+n].Write()
        canDict[s+n].SaveAs("Fit_templates_"+n+"_"+s+".png")
        canDict[s+n].SaveAs("Fit_templates_"+n+"_"+s+".pdf")


# --------------------- signal -------------------------------- #

#features
for hi,hval in hDict.items() :
    for n,name in nameDict.items() :
        if n in hi : 
            hval.SetTitle(name)
    
    if 'ylow' in hi : hval.SetTitle(hval.GetTitle() + ' |Y|<0.4 ')
    if 'yhigh' in hi : hval.SetTitle(hval.GetTitle() + ' 2<|Y|<2.4 ')
    if 'qtlow' in hi : hval.SetTitle(hval.GetTitle() + ' q_{T}<2 GeV ')
    if 'qthigh' in hi : 
        hval.SetTitle(hval.GetTitle() + ' 16 GeV <q_{T}<22 GeV ')
        # hval.Scale(3)
    hval.SetTitle('')
    
    # hval.SetTitleSize(1,"t")
    # gStyle.SetTitleX(0.1) 
    # gStyle.SetTitleY(0.9) 
    # gStyle.SetTitleW(1)
    # gStyle.SetTitleH(0.2) 
    hval.SetStats(0)
    # hval.SetTitleSize(0.2)
    hval.GetYaxis().SetTitle("p_{T} [GeV]")       
    hval.GetXaxis().SetTitle("#eta")       
    # hval.GetZaxis().SetTitle("Events")   
    hval.GetXaxis().SetTitleSize(0.10)      
    hval.GetYaxis().SetTitleSize(0.09) 
    hval.GetZaxis().SetTitleSize(0.07) 
    hval.GetXaxis().SetTitleOffset(0.6)      
    hval.GetYaxis().SetTitleOffset(0.8) 
    # hval.GetZaxis().SetTitleOffset(1.1) 
    hval.GetXaxis().SetLabelSize(0.06)      
    hval.GetYaxis().SetLabelSize(0.06) 
    hval.GetZaxis().SetLabelSize(0.06) 

# align all the z scales 
minSig = {}    
maxSig = {}    
minSig['plus'] = 0.
maxSig['plus'] = 0.
minSig['minus'] = 0.
maxSig['minus'] = 0.

for hi,hval in hDict.items() :
    if 'signal' not in hi : continue 
    hval.SetContour(80)    
    if 'signalUL' in hi : continue 
    maxabs = max(abs(hval.GetMaximum()), abs(hval.GetMinimum()))
    hval.GetZaxis().SetRangeUser(-1.05*maxabs,1.05*maxabs)    

#tlatex boxes
titleDict = {}
for n,name in nameDict.items() :
    if 'signal' in n :
        titleDict[n] =  ROOT.TLatex(coordDict[n][0],coordDict[n][1],name)
        titleDict[n].SetTextSize(0.035)
for k in kindSig:
        titleDict[k] =  ROOT.TLatex(coordDict[k][0],coordDict[k][1],kindSigLabel[kindSig.index(k)])
        titleDict[k].SetTextSize(0.025)
 
# canvas 
for s,f in fileDict.items() :  
    canDict[s+'signal'] = ROOT.TCanvas("c_"+s+"signal","c_"+s+"signal", 1000,1400)
    canDict[s+'signal'].cd()
    # canDict[s+'signal'].Divide(4,6,0.02,0.02)
    canDict[s+'signal'].Divide(4,6,-1,-1)
    it = 0
    for xs in xsecList :
        for kind in kindSig :
            it+=1
            canDict[s+'signal'].cd(it)
            hDict[s+'signal'+xs+kind].Draw("colz")
            
            gPad.SetBottomMargin(0.15)
            gPad.SetTopMargin(0.05)
            gPad.SetLeftMargin(0.18)
            gPad.SetRightMargin(0.18)

    canDict[s+'signal'+'extra'] = ROOT.TCanvas("c_"+s+"signal"+"extra","c_"+s+"signal"+"extra", 1000,1400)
    canDict[s+'signal'+'extra'].cd()
    canDict[s+'signal'+'extra'+'pad'] = ROOT.TPad("p_"+s+"signal"+"extra", "p_"+s+"signal"+"extra",0.1,0,1,0.96)
    canDict[s+'signal'+'extra'+'pad'].cd()
    canDict[s+'signal'].DrawClonePad()
    canDict[s+'signal'+'extra'].cd()
    canDict[s+'signal'+'extra'+'pad'].Draw()
    
    #label 
    canDict[s+'signal'+'extra'].cd()
    for t, tlat in titleDict.items() :
        tlat.Draw("same")       
    

# -------------- other background + data -------------------#

     


    
    
    
for s,f in fileDict.items() :  
    canDict[s+'signal'].Write()  
    canDict[s+'signal'+'extra'].Write()
    canDict[s+'signal'+'extra'].SaveAs("Fit_templates_signal"+s+".png")
    canDict[s+'signal'+'extra'].SaveAs("Fit_templates_signal"+s+".pdf")















#scale, absolute 
# for s,f in fileDict.items() : 
#     for hi,hval in hDict.items() :
#         if 'signal' not in hi : continue 
#         if s not in hi : continue 
#     # if 'signal' in hi :
#         mm  = hval.GetMinimum()
#         if mm < minSig[s] :
#             minSig[s] = mm
#             print(hi, minSig[s])
#         MM  = hval.GetMaximum()
#         if MM > maxSig[s] :
#             maxSig[s] = MM
#     for hi,hval in hDict.items() :
#         if 'signal' not in hi : continue     
#         if s not in hi : continue      
#         hval.GetZaxis().SetRangeUser(minSig[s]-0.1*minSig[s],maxSig[s]+0.1*maxSig[s])

# z scale, xsec based            
# for s,f in fileDict.items() : 
#     for xs in xsecList :
#         minSig[s+xs] = 0.
#         maxSig[s+xs] = 0.
#         for k in kindSig :
#             mm  = hDict[s+'signal'+xs+k].GetMinimum()
#             if mm < minSig[s+xs] :
#                 minSig[s+xs] = mm
#             MM  = hDict[s+'signal'+xs+k].GetMaximum()
#             if MM > maxSig[s+xs] :
#                 maxSig[s+xs] = MM    
#         print(s,xs,minSig[s+xs],maxSig[s+xs])  
#         for k in kindSig : 
#             hDict[s+'signal'+xs+k].GetZaxis().SetRangeUser(minSig[s+xs]-0.1*minSig[s+xs],maxSig[s+xs]+0.1*maxSig[s+xs])

#z scale, qt,y based
# for s,f in fileDict.items() : 
#     for k in kindSig :
#         minSig[s+k] = 0.
#         maxSig[s+k] = 0.
#         for xs in xsecList :
#             mm  = hDict[s+'signal'+xs+k].GetMinimum()
#             if mm < minSig[s+k] :
#                 minSig[s+k] = mm
#             MM  = hDict[s+'signal'+xs+k].GetMaximum()
#             if MM > maxSig[s+k] :
#                 maxSig[s+k] = MM    
#         print(s,k,minSig[s+k],maxSig[s+k])  
#         for xs in xsecList : 
#             hDict[s+'signal'+xs+k].GetZaxis().SetRangeUser(minSig[s+k]-0.1*minSig[s+k],maxSig[s+k]+0.1*maxSig[s+k])