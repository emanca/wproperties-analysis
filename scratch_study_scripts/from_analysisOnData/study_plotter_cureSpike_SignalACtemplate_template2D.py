import ROOT
import numpy as np
import os
import multiprocessing


# this script produce an alternative version of the shape file. 
# to avoid spikes for each shape, evaluate the alternative tempalte as nominal*(intergral_var/integral_nom)
# only for test and debug purpose
    
def runSingleCat(par) :
    singleCat_call = singleCat_cure(categ = par[0], keyList = par[1], path = par[2], sign=par[3])

def singleCat_cure(categ, keyList, path,sign) :
    tightCut = 2
    looseCut = 2#10
    print("cured category ", sign, categ)
    # outFile.cd()
    inFile = ROOT.TFile.Open(path+'W'+s+'_reco.root')
    outFile = ROOT.TFile('W'+sign+'_reco_cureSpike_'+categ+'.root',"recreate")
        
    if 'helXsecs' in categ :
        nomName = categ
        # nomName = key.GetName().split('_')
        # nomName = nomName[:-1]
        # nomName = '_'.join(nomName)
        hnom = inFile.Get(nomName)
    
    for key in keyList :
                                
        if key.GetName().endswith('Up') or key.GetName().endswith('Down') :
            
            if 'helXsecs' not in categ :
                if 'WtoTau' in key.GetName() : nomName = 'WtoTau'
                elif 'data_obs' in key.GetName() : nomName = 'data_obs'
                elif 'Fake' in key.GetName() : nomName = 'Fake'
                elif 'LowAcc' in key.GetName() : nomName = 'LowAcc'
                elif 'DYJets' in key.GetName() : nomName = 'DYJets'
                elif 'DiBoson' in key.GetName() : nomName = 'DiBoson'
                elif 'Top' in key.GetName() : nomName = 'Top'
                elif 'Data' in key.GetName() : nomName = 'Data'
                hnom = inFile.Get(nomName)
                    
            hvar = inFile.Get(key.GetName())
            hvarCured = hvar.Clone(hvar.GetName())
            for x in range(1,hvar.GetNbinsX()+1) :
                for y in range(1,hvar.GetNbinsY()+1) :
                    
                    if 'alpha' in key.GetName() or 'LHEPdf' in key.GetName():
                        if (abs(hvar.GetBinContent(x,y))>looseCut or (abs(hvar.GetBinContent(x,y))>tightCut and abs(hvar.GetBinContent(x,y)-1)<hvar.GetBinError(x,y))) : #alpha is already div/nom
                            hvarCured.SetBinContent(x,y, hvar.GetBinContent(x,y))
                        else : 
                            hvarCured.SetBinContent(x,y, hvar.GetBinContent(x,y)) 
                            
                    else : 
                        if abs(hnom.GetBinContent(x,y))==0.0 and hvar.GetBinContent(x,y)!=0.0 :
                            print("nom==0:", x,y, key.GetName(), hnom.GetName(), ", var=",hvar.GetBinContent(x,y), ", nom=", hnom.GetBinContent(x,y), ", diff=",hvar.GetBinContent(x,y)-hnom.GetBinContent(x,y), ", errVar=",hvar.GetBinError(x,y), ", errNom=",hnom.GetBinError(x,y))                            
                            hvarCured.SetBinContent(x,y, hnom.GetBinContent(x,y))
                        
                        elif hnom.GetBinContent(x,y)!=0.0 and abs(hvar.GetBinContent(x,y)/hnom.GetBinContent(x,y))>tightCut :
                            # print("fix needed", x,y, key.GetName(), hvar.GetBinContent(x,y)/hnom.GetBinContent(x,y))
                            hdiv = hvar.Clone("hdiv")
                            hdiv.Divide(hnom)
                            
                            if abs(hvar.GetBinContent(x,y)-hnom.GetBinContent(x,y))<hvar.GetBinError(x,y)+hnom.GetBinError(x,y) or abs(hdiv.GetBinContent(x,y)-1)<hdiv.GetBinError(x,y): #if difference<error or ratio compatible with 1.
                                hvarCured.SetBinContent(x,y, hnom.GetBinContent(x,y))
                            elif abs(hdiv.GetBinContent(x,y))>looseCut : #if the spike must be cured however
                                hvarCured.SetBinContent(x,y, hnom.GetBinContent(x,y))
                                print("hard strange spike:", x,y, key.GetName(), hnom.GetName(), "var/nom=",hvar.GetBinContent(x,y)/hnom.GetBinContent(x,y), ", var=",hvar.GetBinContent(x,y), ", nom=", hnom.GetBinContent(x,y), ", diff=",hvar.GetBinContent(x,y)-hnom.GetBinContent(x,y), ", errVar=",hvar.GetBinError(x,y), ", errNom=",hnom.GetBinError(x,y))                            
                            else : #the spike is not curable but is small
                                hvarCured.SetBinContent(x,y, hvar.GetBinContent(x,y))
                        else :
                            hvarCured.SetBinContent(x,y, hvar.GetBinContent(x,y))    
            hvarCured.Write()
        elif 'bkg' in categ :
            hnominal = inFile.Get(key.GetName())
            hnominal.Write()
    outFile.Close()

NCORES = 128
Nqt=9
Ny=6
print("N bins qt=",Nqt, ", N bins y=",Ny)

# path = '../test_shiftUpDown/'
# path = '../test_ACvars_fromRemovedScale_filterTemplate/'
path = '../test_4Mar_SB_reb4/'
outList = []
chargeList = ['plus','minus']
for s in chargeList :
    print(s)
    inFile = ROOT.TFile.Open(path+'W'+s+'_reco.root')
    # outFile = ROOT.TFile('W'+s+'_reco_cureSpike.root',"recreate")
    # outFile.cd()
    
    AC_list = ['L',"I","T","A","P","7","8","9","UL"]
    
    catDict = {}
    catDict["bkg"] = []
    for c in AC_list :
        for qt in range(1,Nqt+1) :
            for y in range(1,Ny+1) :
                catDict['helXsecs'+c+'_y_'+str(y)+'_qt_'+str(qt)] = []    
    # catDict = {
    #     "helXsecsL" : [], 
    #     "helXsecsI" : [],
    #     "helXsecsT" : [],
    #     "helXsecsA" : [],
    #     "helXsecsP" : [],
    #     "helXsecs7" : [],
    #     "helXsecs8" : [],
    #     "helXsecs9" : [],
    #     "helXsecsUL" : [],
    #     "bkg" : []
    # }
    
    print("nhistos=", len(inFile.GetListOfKeys()))
    print("building categories...")
    for key in inFile.GetListOfKeys() : 
        if 'helXsecs' in key.GetName() and (key.GetName().endswith('Up') or key.GetName().endswith('Down')) :
            # catName = key.GetName().split('_')[0]+'_'+key.GetName().split('_')[1]+'_'+key.GetName().split('_')[2]+'_'+key.GetName().split('_')[3]+'_'+key.GetName().split('_')[4]
            catName = '_'.join(key.GetName().split('_')[:-1])
            catDict[catName].append(key)
        else :
            catDict['bkg'].append(key)
    
    print("let's start the cure...")
    processesMain = []
    for cat,keyList in catDict.items() :
        processesMain.append((cat,keyList,path, s))
    poolMain = multiprocessing.Pool(NCORES)
    poolMain.map(runSingleCat,processesMain)
    
    print("hadd output...")
    command = 'hadd -f W'+s+'_reco_cureSpike.root '
    for cat,keyList in catDict.items() :
        command+= ' W'+s+'_reco_cureSpike_'+cat+'.root '
    os.system(command)
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        # catName = key.GetName().split('_')[0]
        # if 'hel' in catName :
        #     catDict[catName].append(key)
        # else :
        #     catDict['bkg'].append(key)
        #     print('bkg',key)
            
            
        # if 'alpha' in key.GetName() : continue
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        # if key.GetName().endswith('Up') or key.GetName().endswith('Down') :
        # # if key.GetName().endswith('Up'):
        #     hvar = inFile.Get(key.GetName())
        #     hvarCured = hvar.Clone(hvar.GetName())
        #     if 'WtoTau' in key.GetName() : nomName = 'WtoTau'
        #     elif 'data_obs' in key.GetName() : nomName = 'data_obs'
        #     elif 'Fake' in key.GetName() : nomName = 'Fake'
        #     elif 'LowAcc' in key.GetName() : nomName = 'LowAcc'
        #     elif 'DYJets' in key.GetName() : nomName = 'DYJets'
        #     else : 
        #         nomName = key.GetName().split('_')
        #         nomName = nomName[:-1]
        #         nomName = '_'.join(nomName)
        #     hnom = inFile.Get(nomName)
        #     for x in range(1,hvar.GetNbinsX()+1) :
        #         for y in range(1,hvar.GetNbinsY()+1) :
                    
        #             if 'alpha' in key.GetName() :
        #                 if (abs(hvar.GetBinContent(x,y))>10 or (abs(hvar.GetBinContent(x,y))>2 and abs(hvar.GetBinContent(x,y)-1)<hvar.GetBinError(x,y))) : #alpha is already div/nom
        #                  hvarCured.SetBinContent(x,y, hvar.GetBinContent(x,y))
        #                 elif 'alpha' in key.GetName() and abs(hvar.GetBinContent(x,y))>2 : 
        #                     print("porcodio", key.GetName())
        #                 else : 
        #                     hvarCured.SetBinContent(x,y, hvar.GetBinContent(x,y)) 
                            
        #             else : 
        #                 if abs(hnom.GetBinContent(x,y))==0.0 and hvar.GetBinContent(x,y)!=0.0 :
        #                     print("strange", x,y, key.GetName(), hvar.GetBinContent(x,y))
        #                     hvarCured.SetBinContent(x,y, hnom.GetBinContent(x,y))
                        
                        
        #                 if hnom.GetBinContent(x,y)!=0.0 and abs(hvar.GetBinContent(x,y)/hnom.GetBinContent(x,y))>2 :
        #                     # print("fix needed", x,y, key.GetName(), hvar.GetBinContent(x,y)/hnom.GetBinContent(x,y))
        #                     hdiv = hvar.Clone("hdiv")
        #                     hdiv.Divide(hnom)
                            
                            
        #                     #average between close bin method.... Problema: oscillazione selvaggia fa si che i bin accanto non è detto che siano buoni e spesso amplificano la discrepanza
        #                     # nearBin = np.array([])
        #                     # for i in range(x-1,x+2) :
        #                     #     for j in range(y-1,y+2) :
        #                     #         if not (i==x and j==y) and hvar.GetBinContent(i,j)!=0.0 :
        #                     #             if hnom.GetBinContent(i,j)!=0.0 and abs(hvar.GetBinContent(i,j)/hnom.GetBinContent(i,j))>10 and abs(hvar.GetBinContent(i,j)/hnom.GetBinContent(x,y)): continue
        #                     #             nearBin = np.append(nearBin,hvar.GetBinContent(i,j))
        #                     # if len(nearBin)>0 :
        #                     #     mean = np.mean(nearBin)
        #                     # else :
        #                     #     mean = 0.
        #                     # hvarCured.SetBinContent(x,y, mean)
        #                     # print(x,y, key.GetName(), "var/nom=",hvar.GetBinContent(x,y)/hnom.GetBinContent(x,y), ", var=",hvar.GetBinContent(x,y), ", nom=", hnom.GetBinContent(x,y), ", mu=", mean, ", N=", len(nearBin))    
        #                     # if mean/hnom.GetBinContent(x,y)>10 :
        #                         # print( x,y, key.GetName(), "mu/nom=",mean/hnom.GetBinContent(x,y), hvar.GetBinContent(x,y)/hnom.GetBinContent(x,y), ", var=",hvar.GetBinContent(x,y), ", nom=", hnom.GetBinContent(x,y), ", mu=", mean, ", N=", len(nearBin))
                            
                            
        #                     #use nominal if in error method:
        #                     if abs(hvar.GetBinContent(x,y)-hnom.GetBinContent(x,y))<hvar.GetBinError(x,y)+hnom.GetBinError(x,y) or abs(hdiv.GetBinContent(x,y)-1)<hdiv.GetBinError(x,y): #if difference<error or ratio compatible with 1.
        #                         hvarCured.SetBinContent(x,y, hnom.GetBinContent(x,y))
        #                     elif abs(hdiv.GetBinContent(x,y))>10 : #if the spike must be cured however
        #                         hvarCured.SetBinContent(x,y, hnom.GetBinContent(x,y))
        #                         print("hard strange spike:", x,y, key.GetName(), "var/nom=",hvar.GetBinContent(x,y)/hnom.GetBinContent(x,y), ", var=",hvar.GetBinContent(x,y), ", nom=", hnom.GetBinContent(x,y), ", diff=",hvar.GetBinContent(x,y)-hnom.GetBinContent(x,y), ", errVar=",hvar.GetBinError(x,y), ", errNom=",hnom.GetBinError(x,y))                            
        #                     else : #the spike is not curable but is small
        #                         hvarCured.SetBinContent(x,y, hvar.GetBinContent(x,y))
        #                 else :
        #                     hvarCured.SetBinContent(x,y, hvar.GetBinContent(x,y))    
        #     hvarCured.Write()
        # else :
        #     hnom = inFile.Get(key.GetName())
        #     hnom.Write()