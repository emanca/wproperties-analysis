import ROOT
import math
ROOT.gROOT.SetBatch()

inPath = 'OUTPUT_25Apr_noCF_qtMax60_fixLowAcc/PDFscan_58pdf/'

# charges = ["Wplus","Wminus"]
charges = ["Wplus"]
pdfList = set(['LHEPdfWeightHess{}'.format(i+1) for i in range(60)]+['alphaS'])

shiftDict = {}

excludeList = [
    # 'LHEPdfWeightHess22', #skipped, doesn't work! 
    # 'LHEPdfWeightHess3', #skipped, doesn't work!  
    # 'LHEPdfWeightHess16', #skipped, doesn't work! 
    # 'LHEPdfWeightHess59' #skipped, doesn't work! 
]

unitNuisance=50 #mev
print("warning: hardcoded 50 shift on Mw nuisance")

for charge in charges:
    fullList_shift = []
    shift = 0.
    for pdf in pdfList :
        # if pdf in excludeList : 
        #     inPathMod = inPath.replace('58pdf','noBkg') 
        # else : 
        #     inPathMod = inPath
        if pdf in excludeList : continue
        inPathMod = inPath
        fitFile = ROOT.TFile.Open(inPathMod+pdf+'_'+charge+'/fit_'+charge+'_reco.root') 
        fitTree = fitFile.fitresults
        values = []
        for ev in fitTree :
            values.append(eval('ev.mass'))
        values.sort()
        
        print(pdf, unitNuisance*0.5*abs(values[0] - values[-1]))
        
        shift += (0.5*abs(values[0] - values[-1]))**2
        
        fullList_shift.append((0.5*abs(values[0] - values[-1]))**2)#for missing syst
    
    fullList_shift.sort()
    for pdf in excludeList : 
        shift+=fullList_shift[-1]  #take the lasgest PDF, to compensate the missing pdf (conservative)
        
    shiftDict[charge] = unitNuisance*math.sqrt(shift)

print(shiftDict)
        
    
    

