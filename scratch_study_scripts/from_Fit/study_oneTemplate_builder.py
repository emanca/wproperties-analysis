import ROOT

path = '../test_oneQTYbin_from-ACvars-fromRemovedScale-filterTemplate/'
outList = []
chargeList = ['plus','minus']

for s in chargeList :
    print(s)
    inFile = ROOT.TFile.Open(path+'W'+s+'_reco_backup.root')
    outFile = ROOT.TFile('W'+s+'_reco.root',"recreate")
    outFile.cd()
    print("reshaping...")
    hNom = inFile.Get('data_obs')
    hNom.SetName("helXsecsUL_y_1_qt_1")
    hUp = hNom.Clone(hNom.GetName()+'_massUp')
    hDown = hNom.Clone(hNom.GetName()+'_massDown')
    hUp.Scale(1.1)
    hDown.Scale(1./1.1)
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