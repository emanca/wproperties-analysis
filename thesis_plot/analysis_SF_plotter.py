import ROOT
import math
import copy
ROOT.gROOT.SetBatch(True)

#get input
inFile = ROOT.TFile.Open("../../analysisOnData/data/ScaleFactors_OnTheFly.root")
signList = ['Plus', 'Minus']    
parList = ['0','1','2']
hDict = {}    
hDict['Trigger'+'Plus'] = inFile.Get('Trigger'+'Plus')
hDict['Trigger'+'Minus'] = inFile.Get('Trigger'+'Minus')
hDict['Reco'] = inFile.Get('Reco')
for s in signList :
    for par in parList :
        hDict['Trigger'+s+'Syst'+par]= inFile.Get('Trigger'+s+'Syst'+par)

for s in signList :
    hDict['tot'+s] = hDict['Trigger'+s].Clone('Tot'+s)
    hDict['tot'+s].Multiply(hDict['Reco'])

#build a disctionary with all the required info for plotting
infoDict = {
    'Trigger'+'Plus'  : {
        'x' : '|#eta^{#mu}|',
        'y' : 'p_{T}^{#mu} [GeV]' ,
        'title' : 'SF_{trig}, #mu^{+}',
        'zlow' : 0.75,
        'zhigh' : 1.14,
        'z' : ''
    },
    'Trigger'+'Plus'+'Syst'+'0' : {
        'x' : '|#eta^{#mu}|',
        'y' : 'p_{T}^{#mu} [GeV]',
        'title' : 'Relative SF variation, parameter 0, #mu^{+}',
        'zlow' : -0.014,
        'zhigh' : 0.014,
        'z' : 'Variation/Nominal'
    }
    
}
infoDict['Trigger'+'Minus'] = copy.deepcopy(infoDict['Trigger'+'Plus'])
infoDict['Trigger'+'Minus']['title'] = 'SF_{trig}, #mu^{-}'

infoDict['Reco'] = copy.deepcopy(infoDict['Trigger'+'Plus'])
infoDict['Reco']['title'] = 'SF_{sel}'

infoDict['tot'+'Plus'] = copy.deepcopy(infoDict['Trigger'+'Plus'])
infoDict['tot'+'Plus']['title'] = 'total SF, #mu^{+}'

infoDict['tot'+'Minus'] = copy.deepcopy(infoDict['Trigger'+'Minus'])
infoDict['tot'+'Minus']['title'] = 'total SF, #mu^{-}'



for s in signList :
    for par in parList : 
        if s=='Plus' and par=='0' : continue 
        if s=='Plus' : signString = '#mu^{+}'
        else : signString = '#mu^{-}'
        infoDict['Trigger'+s+'Syst'+par] = copy.deepcopy(infoDict['Trigger'+'Plus'+'Syst'+'0'])
        infoDict['Trigger'+s+'Syst'+par]['title'] = 'Relative SF variation, parameter '+par+', '+signString

# print(infoDict)

#canvas plotting
cDict = {}
for ind, info in infoDict.items() :
    # print(info)
    cDict[ind] = ROOT.TCanvas(ind,ind,1600,1200)
    cDict[ind].cd()
    cDict[ind].SetGridx()
    cDict[ind].SetGridy()
    hDict[ind].SetTitle(info['title'])
    hDict[ind].GetXaxis().SetTitle(info['x'])
    hDict[ind].GetYaxis().SetTitle(info['y'])
    hDict[ind].GetZaxis().SetTitle(info['z'])
    hDict[ind].GetXaxis().SetTitleOffset(0.8)
    hDict[ind].GetYaxis().SetTitleOffset(0.8)
    hDict[ind].GetZaxis().SetTitleOffset(1.4)
    hDict[ind].SetStats(0)
    hDict[ind].Draw("colz")
    if 'Syst' in ind : 
        cDict[ind].SetRightMargin(0.17)
        hDict[ind].Scale(math.sqrt(2))#inflated to take into account also of selection SF
    hDict[ind].GetZaxis().SetRangeUser(info['zlow'],info['zhigh'])

    
    cDict[ind].Update() 
    palette = hDict[ind].GetListOfFunctions().FindObject("palette")
    if 'Syst' in ind :
        palette.SetX1NDC(0.848)
        palette.SetX2NDC(0.873)
    else : 
        palette.SetX1NDC(0.91)
        palette.SetX2NDC(0.93)
    
    palette.SetY1NDC(0.1)
    palette.SetY2NDC(0.9)
    cDict[ind].Modified()
    cDict[ind].Update()
    
    cDict[ind].SaveAs("Analysis_SF_"+ind+".png")
    cDict[ind].SaveAs("Analysis_SF_"+ind+".pdf")


#output
output = ROOT.TFile("analysis_SF_plotter.root", "recreate")
for i, c in cDict.items() :
    c.Write()
    