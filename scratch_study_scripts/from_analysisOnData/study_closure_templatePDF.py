import ROOT
import math
ROOT.gStyle.SetOptStat(0)


# path='./test_ACvars_fromRemovedScale/plot/hadded/'
# path='./test_4GeVflatBinning/plot/hadded/'
# path = './test_ACvars_fromRemovedScale/plot/hadded/''
# inputShapeFile = ROOT.TFile.Open(path+'/WToMu_AC_plots.root')
# inputShapeFile = ROOT.TFile('../Fit/test_ACvars_fromRemovedScale/Wplus_reco.root')
# inputShapeFile = ROOT.TFile('../Fit/test_shiftUpDown/Wplus_reco.root')
# inputShapeFile = ROOT.TFile('../Fit/test_4GeVflatBinning/Wplus_reco.root')


# path = 'test_extraDEBUG_ACvars_fromRemovedScale/plot/hadded/'
path = 'test_ACvars_fromRemovedScale_filterTemplate/plot/hadded/'


inputAC = ROOT.TFile.Open(path+'WToMu_AC_plots.root')
inputSig = ROOT.TFile.Open(path+'WToMu_plots.root')
# inputShapeFile = ROOT.TFile('../Fit/test_ACvars_fromRemovedScale_filterTemplate/Wplus_reco.root')
inputShapeFile = ROOT.TFile('python/Wplus_reco.root')
inputAnaGen = ROOT.TFile.Open('../analysisOnGen/genInput_Wplus.root')
inputAccMap = ROOT.TFile.Open('python/accMap.root')




coeffList = ["L","I","T","A","P","7","8","9","UL"]
# coeffList = ["L","I","T","A","P","UL"]
# qtList = ['1','2','3','4','5','6','7','8','9']
qtList = ['1','2','3','4','5','6','7','8','9']
yList = ['1','2','3','4','5','6']

# factors = {}
# factors["L"] = 2.
# factors["I"] = 2.*math.sqrt(2)
# factors["T"] = 4.
# factors["A"] = 4.*math.sqrt(2)
# factors["P"] = 2.
# factors["7"] = 2.
# factors["8"] = 2.*math.sqrt(2)
# factors["9"] = 4.*math.sqrt(2)
# factors["UL"] = 1
# scale= (3./16./math.pi)


varDict = {
    'nom' : ['Nominal',''],
    'PDF2' : ['LHEPdfWeight','_LHEPdfWeightHess2'],
    'PDF34' : ['LHEPdfWeight','_LHEPdfWeightHess34'],
    'PDF7' : ['LHEPdfWeight','_LHEPdfWeightHess7'],
    'jesUp' : ['jme','_jesTotal'],
}

hdict = {}

outfile = ROOT.TFile("study_closure_templatePDF.root","recreate")

hAccMap = inputAccMap.Get('closPlus')
hAccMap_new = hAccMap.Clone("closPlus_new")
for x in range(1,hAccMap_new.GetNbinsX()+1) :
    for y in range(1,hAccMap_new.GetNbinsY()+1) :
        hAccMap_new.SetBinContent(x,y,0)
hmapTot = inputAnaGen.Get('angularCoefficients/mapTot')

        

for var,varList in varDict.items() :
    if var=='jesUp': #put yours hands up for Jesoo!
        htemp = inputSig.Get('templates_Signal/'+varList[0]+'/templates'+varList[1]+'Up')  
    else :
        htemp = inputSig.Get('templates_Signal/'+varList[0]+'/templates'+varList[1])
    
    # htemp = inputSig.Get('templatesAC_Signal/'+varList[0]+'/templates'+varList[1])
    htemp.GetZaxis().SetRange(2,2)#plus charge
    hdict['templ'+var] = htemp.Project3D('yx')
    hdict['templ'+var].SetName("inclusive_template_"+var)
    # if var=='nom' :
    #     hdict['mapTot'+var] = inputAnaGen.Get('angularCoefficients/mapTot'+varList[1])
    # else :
    #     hdict['mapTot'+var] = inputAnaGen.Get('angularCoefficients_'+varList[0]+'/mapTot'+varList[1])
    
for var,varList in varDict.items() :
    emptyFlag = 1
    for c in coeffList :
        # if c!='UL' :continue
        for qt in qtList :
            for y in yList :
                # if qt!='2' or y!='2' : continue
                # if qt!='2' : continue
                if var=='nom' :
                    htemp = inputShapeFile.Get("helXsecs"+c+"_y_"+y+"_qt_"+qt+varList[1])
                    htemp.Write()
                else :
                    htemp = inputShapeFile.Get("helXsecs"+c+"_y_"+y+"_qt_"+qt+varList[1]+'Up')
                
                if var=='nom' :
                    # val_hAccMap_new = hAccMap_new.GetBinContent(qtList.index(qt)+1,yList.index(y)+1)+htemp.Integral(0,htemp.GetNbinsX()+2,0,htemp.GetNbinsY()+2)
                    # val_hAccMap_new = hAccMap_new.GetBinContent(qtList.index(qt)+1,yList.index(y)+1)+htemp.Integral()
                    # hAccMap_new.SetBinContent(qtList.index(qt)+1,yList.index(y)+1,val_hAccMap_new)
                    
                    val_hAccMap_new = hAccMap_new.GetBinContent(yList.index(y)+1,qtList.index(qt)+1)+htemp.Integral()
                    hAccMap_new.SetBinContent(yList.index(y)+1,qtList.index(qt)+1,val_hAccMap_new)
                
                
                # htemp.Scale(factors[c])
                # htemp = inputShapeFile.Get('templatesAC_Signal/'+varList[0]+'/Wplus_qt_'+qt+'_helXsecs_'+c+varList[1])
                # tempProj = htemp.Project3D('xy')
                # tempProj.Scale(scale*factors[c])
                # if c!="UL" :
                #     for x in range(1, tempProj.GetNbinsX()+1) :            #    tempProj.Multiply(hdict['mapTot'+var])
                #        for y in range(1,tempProj.GetNbinsY()+1) :
                #            tempProj.SetBinContent(x,y,tempProj.GetBinContent(x,y)*hdict['mapTot'+var].GetBinContent(x,y))
                if emptyFlag :
                    hdict['sum'+var] = htemp.Clone('sum'+var)
                    emptyFlag = 0
                else :
                    # print(var,c, qt, y),
                    # print(htemp.GetBinContent(20,20))
                    hdict['sum'+var].Add(htemp)
    if var=='nom' : 
        # htempLA = inputAC.Get('templatesLowAcc_Signal/'+varList[0]+'/templates'+varList[1])
        # htempLA.GetZaxis().SetRange(2,2)#plus charge
        # hlowacc = htempLA.Project3D('yx')
        # hlowacc.SetName("lowacc_"+var)
        hdict['lowacc'+var] = inputShapeFile.Get("LowAcc"+varList[1])
    else :
        hdict['lowacc'+var] = inputShapeFile.Get("LowAcc"+varList[1]+'Up')
    hdict['lowacc'+var].Write()
    
    
    hdict['sum'+var].Add(hdict['lowacc'+var])
    # print("warning:not summed lowacc")

#change x y axis
    # sum_var =  hdict['templ'+var].Clone("sum_"+var)
    # for x in range(1, sum_var.GetNbinsX()+1) :          
    #     for y in range(1,sum_var.GetNbinsY()+1) : 
    #         sum_var.SetBinContent(x,y,hdict['sum'+var].GetBinContent(y,x))
    # hdict['sum'+var] = sum_var      
    hdict['sum'+var].SetName("sum_of_differential_templates_"+var)
    

#external, partial total template (one bin integrated in angular coeff) :
# hOneBin3D = inputAC.Get('templatesAC_Signal/Nominal/Wplus_qt_2_helXsecs_')
# hOneBin3D.GetZaxis().SetRange(2,2)
# hOneBin2D = hOneBin3D.Project3D('yx')
# hOneBin2D.Write()

#resumming histos from template builder (without harmonics weight)
# emptyFlag2 = 1
# for qt in qtList :
#     for y in yList :
#         # if qt!='2' or y!='2' : continue
#         # if qt!='2' : continue
#         hOneBin3D = inputAC.Get('templatesAC_Signal/Nominal/Wplus_qt_'+qt+'_helXsecs_')
#         hOneBin3D.GetZaxis().SetRange(yList.index(y)+1,yList.index(y)+1)
#         htemp2 = hOneBin3D.Project3D('yx')
#         htemp2.Write()
#         if emptyFlag2 :
#             hdict['RESUM'] = htemp2.Clone('RESUM_from_templ')
#             emptyFlag2 = 0
#         else :
#             hdict['RESUM'].Add(htemp2)

#global histo produced inside template builder
# hTempBuildTot = inputAC.Get('templatesAC_Signal/Nominal/templatesCUT')
# hTempBuildTot.GetZaxis().SetRange(2,2)
# hTempBuildTot2 = hTempBuildTot.Project3D('yx')
# hTempBuildTot2.SetName("templatesCUT")
# hTempBuildTot2.Write()













# print("DEBUG") 
# print(hdict['templ'+var].GetNbinsX(), hdict['templ'+var].GetNbinsY())
# print(hdict['templ'+var].GetXaxis().GetBinCenter(47), hdict['templ'+var].GetYaxis().GetBinCenter(20)) 
# print("templ", hdict['templ'+var].GetBinContent(47,20))
# print("sum", hdict['sum'+var].GetBinContent(47,20))
# print("lowacc", hlowacc.GetBinContent(47,20))
# print("ratio", hdict['templ'+var].GetBinContent(47,20)/hdict['sum'+var].GetBinContent(47,20))


for var,varList in varDict.items() :
    #  hdict['rRESUM_HARM'+var] = hdict['RESUM'].Clone('ratio_RESUM_over_sumHARM_'+var)
    #  hdict['rRESUM_HARM'+var].Divide(hdict['sum'+var])
     
    hdict['rTEMPL_HARM'+var] = hdict['templ'+var].Clone('ratio_templ_over_sumHARM_'+var)
    #  temm = hdict['sum'+var].Clone('temm')
    #  temm.Add(hlowacc)
    #  hdict['rTEMPL_HARM'+var].Divide(temm)
    hdict['rTEMPL_HARM'+var].Divide(hdict['sum'+var])
     
    #  hdict['rCUT_HARM'+var] = hTempBuildTot2.Clone("ratio_templCUT_over_sumHARM_"+var)
    #  hdict['rCUT_HARM'+var].Divide(hdict['sum'+var])
     
    #  hdict['rTEMPL_CUT'+var] = hdict['templ'+var].Clone('ratio_templ_over_templCUT_'+var)
    #  temm2 = hTempBuildTot2.Clone('temm2')
    #  temm2.Add(hlowacc)
    #  hdict['rTEMPL_CUT'+var].Divide(temm2)
     
    #  hdict['rTEMPL_RESUM'+var] = hdict['templ'+var].Clone('ratio_templ_over_RESUM_'+var)
    #  temm3 = hdict['RESUM'].Clone('temm2')
    #  temm3.Add(hlowacc)
    #  hdict['rTEMPL_RESUM'+var].Divide(temm3)

hAccMap_new.Divide(hmapTot)
hAccMap_ratio = hAccMap_new.Clone("hAccMap_ratio")
hAccMap_ratio.Divide(hAccMap)


for var,varList in varDict.items() :
    # for k in ['sum', 'templ', 'rRESUM_HARM', 'rTEMPL_HARM','rCUT_HARM', 'rTEMPL_CUT', 'rTEMPL_RESUM'] :
    for k in ['sum', 'templ',  'rTEMPL_HARM'] :
         hdict[k+var].Write()

hAccMap.Write()
hAccMap_new.Write()
hAccMap_ratio.Write()
# hdict['RESUM'].Write()
                 

#debug
# for var,varList in varDict.items() :
#      for k in ['sum', 'templ'] :
#           print(hdict[k+var])
#           for x in range(0,hdict[k+var].GetNbinsX()+1) :
#               print(hdict[k+var].GetXaxis().GetBinUpEdge(x))
#           for y in range(0,hdict[k+var].GetNbinsY()+1) :
#               print(hdict[k+var].GetYaxis().GetBinUpEdge(y)) 