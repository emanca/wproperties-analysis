import ROOT

path = '../test_ACvars_fromRemovedScale_filterTemplate/'
# path = '../test_oneQTYbin_from-ACvars-fromRemovedScale-filterTemplate/'
outList = []
chargeList = ['plus','minus']
cList = ["L", "I", "T", "A", "P", "7", "8", "9", "UL"]
nQt = 9#8
nY = 6

for s in chargeList :
    print(s)
    # inFile = ROOT.TFile.Open(path+'W'+s+'_reco_backup.root')
    inFile = ROOT.TFile.Open(path+'W'+s+'_reco.root')
    outFile = ROOT.TFile('W'+s+'_reco.root',"recreate")
    outFile.cd()
    print("reshaping...")
    for c in cList :
        for qt in range(1,nQt+1) :
            for y in range(1,nY+1) :
                hNom = inFile.Get("helXsecs"+c+"_y_"+str(y)+"_qt_"+str(qt))
                
                # for x in range(1, hNom.GetNbinsX()+1) :
                #     for y in range(1, hNom.GetNbinsX()+1) :
                #         val = abs(hNom.GetBinContent(x,y))
                #         hNom.SetBinContent(x,y,val)
                
                hUp = hNom.Clone(hNom.GetName()+'_massUp')
                hDown = hNom.Clone(hNom.GetName()+'_massDown')
                hUp.Scale(1.1)
                hDown.Scale(1./1.1)
                # shift = hUp.Clone("shift")
                # hUp.Add(shift,10)
                # hDown.Add(shift,-10)

                        
                
                hNom.Write()
                hUp.Write()
                hDown.Write()
                
    #low acc template
    hNom = inFile.Get("LowAcc")
    hUp = hNom.Clone(hNom.GetName()+'_massUp')
    hDown = hNom.Clone(hNom.GetName()+'_massDown')
    hUp.Scale(1.1)
    hDown.Scale(1./1.1)
    shift = hUp.Clone("shift")
    # hUp.Add(shift,10)
    # hDown.Add(shift,-10)
    hNom.Write()
    hUp.Write()
    hDown.Write()
    
    h = inFile.Get("Data")
    h.Write()
    h = inFile.Get("data_obs")
    h.Write()