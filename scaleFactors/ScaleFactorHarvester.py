import math
import ROOT
import sys
import argparse
ROOT.gROOT.SetBatch(True)

#########################  DOCUMENTATION BEFORE THE MAIN FUNCTION #######################

def SFHarvester_nom(fileDict) :
    outDict = {}
    outDict['Trigger'+'Plus'] = fileDict['Trigger'+'Plus'].Get('scaleFactor')
    outDict['Trigger'+'Minus'] = fileDict['Trigger'+'Minus'].Get('scaleFactor')
    outDict['Reco'] = fileDict['Reco'].Get('scaleFactor_etaInterpolated')
    for ind, h in outDict.iteritems() :
        if 'Reco' in ind :
            h.Rebin2D(3,5)
            rebFactor=15
        else :
            h.RebinY(5)
            rebFactor=5
        for xx in range(1,h.GetXaxis().GetNbins()+1) :
            for yy in range(1,h.GetYaxis().GetNbins()+1) :
                h.SetBinContent(xx,yy,h.GetBinContent(xx,yy)/rebFactor)
    return outDict


def SFHarvester_triggerSyst(fileDict, SFhisto, ratio=True, WHelicity=False) :
    outDict = {}
    signList = ['Plus', 'Minus']
    parList = ['0','1','2']
    systDict = {}
    
    for s in signList :
        for par in parList :
            systDict['Trigger'+s+'Syst'+par+'Ratio']= fileDict['Trigger'+s+'syst'].Get('p'+par)
            outDict['Trigger'+s+'Syst'+par+'Up'] = SFhisto['Trigger'+s].Clone('Trigger'+s+'Syst'+par+'Up')
            outDict['Trigger'+s+'Syst'+par+'Down'] = SFhisto['Trigger'+s].Clone('Trigger'+s+'Syst'+par+'Down')
            nEtaBins = outDict['Trigger'+s+'Syst'+par+'Up'].GetXaxis().GetNbins()
            nPtBins = outDict['Trigger'+s+'Syst'+par+'Up'].GetYaxis().GetNbins()
            for eta in range(1,nEtaBins+1) :
                for pt in range(1,nPtBins+1) :
                    etaBin = systDict['Trigger'+s+'Syst'+par+'Ratio'].GetXaxis().FindBin(outDict['Trigger'+s+'Syst'+par+'Up'].GetXaxis().GetBinCenter(eta))
                    ptBin = systDict['Trigger'+s+'Syst'+par+'Ratio'].GetYaxis().FindBin(outDict['Trigger'+s+'Syst'+par+'Up'].GetYaxis().GetBinCenter(pt))
                    if not WHelicity :
                        valUp = outDict['Trigger'+s+'Syst'+par+'Up'].GetBinContent(eta,pt)*(1+systDict['Trigger'+s+'Syst'+par+'Ratio'].GetBinContent(etaBin,ptBin))
                        valDown = outDict['Trigger'+s+'Syst'+par+'Down'].GetBinContent(eta,pt)*(1-systDict['Trigger'+s+'Syst'+par+'Ratio'].GetBinContent(etaBin,ptBin))
                    else :
                        valUp = outDict['Trigger'+s+'Syst'+par+'Up'].GetBinContent(eta,pt)*(1+math.sqrt(2)*systDict['Trigger'+s+'Syst'+par+'Ratio'].GetBinContent(etaBin,ptBin))
                        valDown = outDict['Trigger'+s+'Syst'+par+'Down'].GetBinContent(eta,pt)*(1-math.sqrt(2)*systDict['Trigger'+s+'Syst'+par+'Ratio'].GetBinContent(etaBin,ptBin))
                    outDict['Trigger'+s+'Syst'+par+'Up'].SetBinContent(eta,pt,valUp)
                    outDict['Trigger'+s+'Syst'+par+'Down'].SetBinContent(eta,pt,valDown)
    if ratio :
       outDict.update(systDict)
                
    return outDict

           
def SFHarvester_recoSyst(fileDict, SFhisto, era_ratios, ratio=True) :
    outDict = {}
    print "WARNING: era ratios hardcoded [BCDEF, GH]=", era_ratios
    parList = ['Syst','Stat']
    systDict = {}
    
    for par in parList :
        hBCDEF=fileDict['Iso'+'BCDEF'].Get('NUM_TightRelIso_DEN_MediumID_eta_pt_'+par.lower())
        hGH=fileDict['Iso'+'GH'].Get('NUM_TightRelIso_DEN_MediumID_eta_pt_'+par.lower())
        outDict['Reco'+par+'Up'] = SFhisto['Reco'].Clone('Reco'+par+'Up')
        outDict['Reco'+par+'Down'] = SFhisto['Reco'].Clone('Reco'+par+'Down')
        nEtaBins = outDict['Reco'+par+'Up'].GetXaxis().GetNbins()
        nPtBins = outDict['Reco'+par+'Up'].GetYaxis().GetNbins()
        
        if ratio :
           systDict['Reco'+par+'Ratio'] =  SFhisto['Reco'].Clone('Reco'+par+'Ratio') #dummy clone to build the binning
           systDict['Reco'+par+'Ratio'].GetZaxis().SetRangeUser(-0.01,0.01)
           
        for eta in range(1,nEtaBins+1) :
            for pt in range(1,nPtBins+1) :
                etaBin = hBCDEF.GetXaxis().FindBin(outDict['Reco'+par+'Up'].GetXaxis().GetBinCenter(eta))
                ptBin = hBCDEF.GetYaxis().FindBin(outDict['Reco'+par+'Up'].GetYaxis().GetBinCenter(pt))
                val = hBCDEF.GetBinError(etaBin,ptBin)*era_ratios[0]+hGH.GetBinError(etaBin,ptBin)*era_ratios[1]
                valUp = outDict['Reco'+par+'Up'].GetBinContent(eta,pt)+val
                valDown = outDict['Reco'+par+'Down'].GetBinContent(eta,pt)-val
                outDict['Reco'+par+'Up'].SetBinContent(eta,pt,valUp)
                outDict['Reco'+par+'Down'].SetBinContent(eta,pt,valDown)
                if RATIO :
                    systDict['Reco'+par+'Ratio'].SetBinContent(eta,pt,val/SFhisto['Reco'].GetBinContent(eta,pt))
                    systDict['Reco'+par+'Ratio'].SetBinError(eta,pt,0)              
    if ratio :
       outDict.update(systDict)
                
    return outDict

def SFHarvester_IDCorr(fileDict,SFhisto,era_ratios, syst=True) : 
    
    # This function produce the Iso SF alone, with the assumpion ISO*ID=(ISO*ID)_Wheilicy=(RECO)_WHelicity, ISO/ID=(ISO/ID)_POG
    # To do that: 
    # -build ISO/ID from POG
    # -rebin in WHelicity
    # -build the histo of s=sqrt((RECO)_WHelicity/(ISO/ID)_POG), solution of the assumpion system for ISO
    # - repeat all the previous steps for all the variation of ISO-POG
    
        
    outDict = {}
    print "WARNING: era ratios hardcoded [BCDEF, GH]=", era_ratios
    
    #build nom id/reco, pog binning
    h_reco_BCDEF=fileDict['Iso'+'BCDEF'].Get('NUM_TightRelIso_DEN_MediumID_eta_pt')
    h_reco_GH=fileDict['Iso'+'GH'].Get('NUM_TightRelIso_DEN_MediumID_eta_pt')
    h_id_BCDEF=fileDict['Id'+'BCDEF'].Get('NUM_MediumID_DEN_genTracks_eta_pt')
    h_id_GH=fileDict['Id'+'GH'].Get('NUM_MediumID_DEN_genTracks_eta_pt')
    h_reco_pog = h_reco_BCDEF.Clone("h_reco_pog")
    h_reco_pog.Add(h_reco_BCDEF,h_reco_GH,era_ratios[0],era_ratios[1])
    h_id_pog = h_id_BCDEF.Clone("h_id_pog")
    h_id_pog.Add(h_id_BCDEF,h_id_GH,era_ratios[0],era_ratios[1])
    IdReco_ratio = h_id_pog.Clone('IdReco_ratio')
    IdReco_ratio.Divide(h_reco_pog)
    
    #build nom iso, whelicity binning
    outDict['Corr_Iso'] = SFhisto['Reco'].Clone('Corr_Iso')
    nEtaBins = outDict['Corr_Iso'].GetXaxis().GetNbins()
    nPtBins = outDict['Corr_Iso'].GetYaxis().GetNbins()
    for eta in range(1,nEtaBins+1) :
        for pt in range(1,nPtBins+1) :
            etaBin = IdReco_ratio.GetXaxis().FindBin(outDict['Corr_Iso'].GetXaxis().GetBinCenter(eta))
            ptBin = IdReco_ratio.GetYaxis().FindBin(outDict['Corr_Iso'].GetYaxis().GetBinCenter(pt))
            IdReco_ratio_val = IdReco_ratio.GetBinContent(etaBin,ptBin)
            val = math.sqrt(SFhisto['Reco'].GetBinContent(eta,pt)/IdReco_ratio_val)
            outDict['Corr_Iso'].SetBinContent(eta,pt,val)
            # outDict['IdReco_ratioCorr'].SetBinContent(eta,pt,IdReco_ratio_val)
    
    # outDict['Corr_Iso'] = SFhisto['Reco'].Clone('IdReco_ratioCorr')
    # for eta in range(1,nEtaBins+1) :
    #     for pt in range(1,nPtBins+1) :
    #         val = math.sqrt(SFhisto['Reco'].GetBinContent(eta,pt)/outDict['IdReco_ratioCorr'].GetBinContent(eta,pt))
    #         outDict['Corr_Iso'].SetBinContent(eta,pt,val)
    
    if syst :
        parList = ['Syst','Stat']
        direcList = ['Up','Down']
        for par in parList :
            hBCDEF=fileDict['Iso'+'BCDEF'].Get('NUM_TightRelIso_DEN_MediumID_eta_pt_'+par.lower())
            hGH=fileDict['Iso'+'GH'].Get('NUM_TightRelIso_DEN_MediumID_eta_pt_'+par.lower())
            for direc in direcList :
                outDict['Corr_Iso'+par+direc] = SFhisto['Reco'].Clone('Corr_Iso'+par+direc)
                for eta in range(1,nEtaBins+1) :
                    for pt in range(1,nPtBins+1) :
                        etaBin = h_id_pog.GetXaxis().FindBin(outDict['Corr_Iso'+par+direc].GetXaxis().GetBinCenter(eta))
                        ptBin = h_id_pog.GetYaxis().FindBin(outDict['Corr_Iso'+par+direc].GetYaxis().GetBinCenter(pt))
                        valVar = hBCDEF.GetBinError(etaBin,ptBin)*era_ratios[0]+hGH.GetBinError(etaBin,ptBin)*era_ratios[1]
                        valNom = hBCDEF.GetBinContent(etaBin,ptBin)*era_ratios[0]+hGH.GetBinContent(etaBin,ptBin)*era_ratios[1]
                        if direc == 'Up' :
                            val = valNom+valVar
                        if direc == 'Down' :
                            val = valNom-valVar                        
                        IdReco_ratio_val = val/SFhisto['Reco'+par+direc].GetBinContent(eta,pt)
                        val = math.sqrt(SFhisto['Reco'].GetBinContent(eta,pt)/IdReco_ratio_val)
                        outDict['Corr_Iso'+par+direc].SetBinContent(eta,pt,val)
                        # outDict['IdReco_ratioCorr'+par+direc].SetBinContent(eta,pt,val)
        
    return outDict 
    
    
def SFHarvester_POG(fileDict,SFhisto,era_ratios) :
    
    print "WARNING: era ratios hardcoded [BCDEF, GH]=", era_ratios

    outDict = {}
    tempDict = {}
    histo4bin = fileDict['Trigger'+'Plus'].Get('scaleFactor')
    nEtaBins = histo4bin.GetXaxis().GetNbins()
    nPtBins = histo4bin.GetYaxis().GetNbins()
    
    sfDict = {
        'Trigger' : 'IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio',
        'Id' : 'NUM_MediumID_DEN_genTracks_eta_pt',
        'Iso' : 'NUM_TightRelIso_DEN_MediumID_eta_pt'
    }
    varList = ['Syst','Stat']
    direcList = ['Up','Down']

    #merge eras:
    for sf,name in sfDict.iteritems() :
        tempDict[sf+'BCDEF']=fileDict[sf+'BCDEF'].Get(name)
        tempDict[sf+'GH']=fileDict[sf+'GH'].Get(name)
        tempDict[sf+'POG'+'unreb'] =  tempDict[sf+'BCDEF'].Clone(sf+'POG'+'unreb')
        tempDict[sf+'POG'+'unreb'].Add(tempDict[sf+'BCDEF'],tempDict[sf+'GH'],era_ratios[0],era_ratios[1])
    
    #rebin
    for sf,name in sfDict.iteritems() :
        outDict[sf+'POG'] = SFhisto['Reco'].Clone(sf+'POG')
        for eta in range(1,nEtaBins+1) :
            for pt in range(1,nPtBins+1) :
                ptBin =  tempDict[sf+'POG'+'unreb'].GetYaxis().FindBin(outDict[sf+'POG'].GetYaxis().GetBinCenter(pt))
                if sf!='Trigger' : #trigger defined for abseta
                    etaBin = tempDict[sf+'POG'+'unreb'].GetXaxis().FindBin(outDict[sf+'POG'].GetXaxis().GetBinCenter(eta))
                else :
                    etaBin = tempDict[sf+'POG'+'unreb'].GetXaxis().FindBin(abs(outDict[sf+'POG'].GetXaxis().GetBinCenter(eta)))
                val = tempDict[sf+'POG'+'unreb'].GetBinContent(etaBin,ptBin)
                err = tempDict[sf+'POG'+'unreb'].GetBinError(etaBin,ptBin)
                outDict[sf+'POG'].SetBinContent(eta,pt,val)
                outDict[sf+'POG'].SetBinError(eta,pt,err)
    # outDict['Reco'+'POG'] = outDict['Id'+'POG'].Clone('RecoPOG')
    # outDict['Reco'+'POG'].Multiply(outDict['Iso'+'POG'])
    
    #get varied histos
    for sf,name in sfDict.iteritems() :
        for v in varList :
            if sf=='Trigger' and v!='Stat' : continue
            if sf=='Trigger' : vv = ''
            else : vv = '_'+v.lower()
            tempDict[sf+'BCDEF'+v]=fileDict[sf+'BCDEF'].Get(name+vv)
            tempDict[sf+'GH'+v]=fileDict[sf+'GH'].Get(name+vv)

    #uncertainties
    for sf,name in sfDict.iteritems() :
        for v in varList :
            if sf=='Trigger' and v!='Stat' : continue
            for d in direcList : 
                outDict[sf+'POG'+v+d] = SFhisto['Reco'].Clone(sf+'POG'+v+d)
                for eta in range(1,nEtaBins+1) :
                    for pt in range(1,nPtBins+1) :
                        etaBin = tempDict[sf+'BCDEF'+v].GetXaxis().FindBin(outDict[sf+'POG'].GetXaxis().GetBinCenter(eta))
                        ptBin = tempDict[sf+'BCDEF'+v].GetYaxis().FindBin(outDict[sf+'POG'].GetYaxis().GetBinCenter(pt))
                        valVar = tempDict[sf+'BCDEF'+v].GetBinError(etaBin,ptBin)*era_ratios[0]+tempDict[sf+'BCDEF'+v].GetBinError(etaBin,ptBin)*era_ratios[1]
                        valNom = outDict[sf+'POG'].GetBinContent(eta,pt)
                        if d == 'Up' :
                            val = valNom+valVar
                        if d == 'Down' :
                            val = valNom-valVar
                        outDict[sf+'POG'+v+d].SetBinContent(eta,pt,val)    
    return outDict
            
        
def SFHarvester_WHelicty(fileDict, SFhisto, ratio=True) :
    outDict = {}
    signList = ['Plus', 'Minus']
    parList = ['0','1','2']
    direcList = ['Up','Down']
    
    nEtaBins = SFhisto['Trigger'+'Plus'].GetXaxis().GetNbins()
    nPtBins = SFhisto['Trigger'+'Plus'].GetYaxis().GetNbins()
    
    for s in signList :
        #nominal
        outDict['WHelicity'+s] =  SFhisto['Trigger'+s].Clone('WHelicity'+s)
        outDict['WHelicity'+s].Multiply(SFhisto['Reco'])
        for direc in direcList :
            for par in parList :    
                #variation "stat"
                outDict['WHelicity'+s+'Syst'+par+direc] = SFhisto['Trigger'+s+'Syst'+par+direc].Clone('WHelicity'+s+'Syst'+par+direc)
                outDict['WHelicity'+s+'Syst'+par+direc].Multiply(SFhisto['Reco'])
            
            #variation "syst"
            outDict['WHelicity'+s+'Syst'+'Flat'+direc] = SFhisto['Trigger'+s].Clone('WHelicity'+s)
            for eta in range(1,nEtaBins+1) :
                etaVal = outDict['WHelicity'+s].GetXaxis().GetBinCenter(eta)
                if abs(etaVal)<1.0 :
                    modifier = 0.002
                elif abs(etaVal)<1.5 and abs(etaVal)>1:
                    modifier = 0.004
                else :
                    modifier = 0.014
                
                for pt in range(1,nPtBins+1) :
                    if direc == 'Up' :
                        value = outDict['WHelicity'+s].GetBinContent(eta,pt)*(1+modifier)
                    else :
                        value = outDict['WHelicity'+s].GetBinContent(eta,pt)*(1-modifier)
                    outDict['WHelicity'+s+'Syst'+'Flat'+direc].SetBinContent(eta,pt,value)
    
    if ratio :
        for s in signList :
            parList.append('Flat')
            for par in parList :
                outDict['WHelicity'+s+par+'Ratio'] = outDict['WHelicity'+s+'Syst'+par+'Up'].Clone('WHelicity'+s+par+'Ratio')
                outDict['WHelicity'+s+par+'Ratio'].Divide(outDict['WHelicity'+s])
    
    return outDict

def SFHarvester_OTF(fileDict) :
    outDict = {}
    signList = ['Plus', 'Minus']
    parList = ['0','1','2']
    
    outDict['Trigger'+'Plus'] = fileDict['Trigger'+'Plus'].Get('scaleFactor')
    outDict['Trigger'+'Minus'] = fileDict['Trigger'+'Minus'].Get('scaleFactor')
    outDict['Reco'] = fileDict['Reco'].Get('scaleFactor_etaInterpolated')
    for s in signList :
        for par in parList :
            outDict['Trigger'+s+'Syst'+par]= fileDict['Trigger'+s+'syst'].Get('p'+par)
    
    return outDict

#####################################################################################################################
#
#   AUTHOR: V. Bertacchi (06/2020)
#   THIS DOCUMENTATION IS OUTDATED (10/2/2021)
#
#   USAGE: python SaleFactorHarvester.py (--output_name) (--saveRatio) (--syst)
#   
#   INPUT:
#   expected to run in the directory where is placed data/SF_*.root
#   the .root input name are hardcoded in the main function
#
#   OUTPUT STRUCTURE: 
#   main dir: eta vs pt SF and their Up/Down variation
#   "ratios" dir: eta vs pt (variedSF-nominalSF)/nominalSF
#   "correction" dir: eta vs pt sqrt((RECO)_WHelicity/(ISO/ID)_POG) = ISO SF
#   
#   OUTPUT NAMING CONVENTION: Kind+(Sign)+Syst+("ratio")
#   Kind = "Trigger", "Reco"
#   Sign (Trigger) = "Plus", "Minus"
#   Sign (Reco) = ""
#   Syst = SystName+"Up", SystName+"Down"
#   SystName (Trigger) = "Syst0","Syst1","Syst2" (they are the variation of the parameter of: [0]*erf((x-[1])/[2]) )
#   SystName (Reco) = "Syst", "Stat"
#
#   #Which SF are used?
#   The central values are from WHeliciy (SMP18-012)
#   the variation of trigger are from WHelicity statistical variation of erf
#   the variation of the Reco=ID+ISO are from the ISO SF of the POG 
#   (because in SMP18-012 they overestimated it and they did not use for FR and PR eval., where ISO alone is needed)
#
#
#####################################################################################################################
    

parser = argparse.ArgumentParser("")
parser.add_argument('-o','--output_name', type=str, default='ScaleFactors',help="name of the output root file")
parser.add_argument('-saveRatio', '--saveRatio',type=int, default=True, help="1=Save also sf variation ratios")
parser.add_argument('-syst', '--syst',type=int, default=True, help="1=run the syst. variation")
parser.add_argument('-id', '--id',type=int, default=True, help="1=evaluate ID correction")
parser.add_argument('-pog', '--pog',type=int, default=True, help="1= evaluate the POG SF")
parser.add_argument('-whel', '--whelicity',type=int, default=False, help="1= evaluate the WHelicity SF")
parser.add_argument('-OTF', '--OnTheFly',type=int, default=False, help="1= The SF will be evaluated on the fly in RDF, here only the packer in a single file")



args = parser.parse_args()
RATIO = args.saveRatio  # ratio = (variedSF-nominalSF)/nominalSF
OUTPUT = args.output_name
SYST = args.syst
ID = args.id
POG = args.pog
WHEL = args.whelicity
OTF = args.OnTheFly

era_ratios = [0.548,0.452]

#input files
fileDict = {}

#Whelicity
fileDict['Trigger'+'Plus'] = ROOT.TFile.Open('data/SF_WHelicity_Trigger_Plus.root')
fileDict['Trigger'+'Minus'] = ROOT.TFile.Open('data/SF_WHelicity_Trigger_Minus.root')
fileDict['Reco'] = ROOT.TFile.Open('data/SF_WHelicity_Reco.root')

#Whelicity syst
fileDict['Trigger'+'Plus'+'syst'] = ROOT.TFile.Open('data/SF_WHelicity_Trigger_Plus_syst.root')
fileDict['Trigger'+'Minus'+'syst'] = ROOT.TFile.Open('data/SF_WHelicity_Trigger_Minus_syst.root')

if OTF :
    SFhisto = SFHarvester_OTF(fileDict=fileDict)
    output = ROOT.TFile(OUTPUT+".root", "recreate")
    output.cd()
    for ind, histo in SFhisto.iteritems() :
        histo.SetName(ind)
        histo.Write()
    
else :
    #POG
    fileDict['Iso'+'BCDEF'] = ROOT.TFile.Open('data/SF_POG_Iso_BCDEF.root')
    fileDict['Iso'+'GH'] = ROOT.TFile.Open('data/SF_POG_Iso_GH.root')
    fileDict['Id'+'BCDEF'] = ROOT.TFile.Open('data/SF_POG_Id_BCDEF.root')
    fileDict['Id'+'GH'] = ROOT.TFile.Open('data/SF_POG_Id_GH.root')
    fileDict['Trigger'+'BCDEF'] = ROOT.TFile.Open('data/SF_POG_Trigger_BCDEF.root')
    fileDict['Trigger'+'GH'] = ROOT.TFile.Open('data/SF_POG_Trigger_GH.root')


    #nominal SF
    SFhisto = SFHarvester_nom(fileDict=fileDict)

    if SYST : #variation
        SFhisto.update(SFHarvester_triggerSyst(fileDict=fileDict, SFhisto=SFhisto,ratio=RATIO,WHelicity=WHEL))
        SFhisto.update(SFHarvester_recoSyst(fileDict=fileDict, SFhisto=SFhisto,era_ratios=era_ratios,ratio=RATIO))

    if ID :
        SFhisto.update(SFHarvester_IDCorr(fileDict=fileDict, SFhisto=SFhisto,era_ratios=era_ratios,syst=SYST)) 

    if POG :    
        SFhisto.update(SFHarvester_POG(fileDict=fileDict,SFhisto=SFhisto,era_ratios=era_ratios))

    if WHEL :
        print "WHelicity mode active: sqrt(2) factor on trigger syst!"
        SFhisto.update(SFHarvester_WHelicty(fileDict=fileDict, SFhisto=SFhisto,ratio=RATIO))
        
    #output saving
    output = ROOT.TFile(OUTPUT+".root", "recreate")
    output.cd()
    for ind, histo in SFhisto.iteritems() :
        histo.SetName(ind)
        if 'Ratio' in ind or 'Corr' in ind or 'POG' in ind or 'WHelicity' in ind: continue
        histo.Write()
        
    if RATIO :
        direcRATIO = output.mkdir('ratios')
        direcRATIO.cd()
        for ind, histo in SFhisto.iteritems() :
            if not 'Ratio' in ind : continue
            histo.Write()

    if ID :
        direcID = output.mkdir('Corr_Iso')
        direcID.cd()
        for ind, histo in SFhisto.iteritems() :
            if not 'Corr' in ind : continue
            histo.Write()
    
    if POG :
        direcPOG = output.mkdir('POG')
        direcPOG.cd()
        for ind, histo in SFhisto.iteritems() :
            if not 'POG' in ind : continue
            histo.Write()
            
    if WHEL :
        direcWHelicity = output.mkdir('WHelicity')
        direcWHelicity.cd()
        for ind, histo in SFhisto.iteritems() :
            if not 'WHelicity' in ind : continue
            if 'Ratio' in ind : continue
            histo.Write()