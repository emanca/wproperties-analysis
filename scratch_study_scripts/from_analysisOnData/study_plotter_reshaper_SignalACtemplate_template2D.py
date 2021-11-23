import ROOT

# this script produce an alternative version of the shape file. 
# to avoid spikes for each shape, evaluate the alternative tempalte as nominal*(intergral_var/integral_nom)
# only for test and debug purpose
path = '../test_shiftUpDown/'
outList = []
chargeList = ['plus','minus']
for s in chargeList :
    print(s)
    inFile = ROOT.TFile.Open(path+'W'+s+'_reco.root')
    outFile = ROOT.TFile('W'+s+'_reco_reshape.root',"recreate")
    outFile.cd()
    print("reshaping...")
    for key in inFile.GetListOfKeys() :
        # if key.GetName().endswith('Up') or key.GetName().endswith('Down') :
        if key.GetName().endswith('Up'):
            hvar = inFile.Get(key.GetName())
            if 'WtoTau' in key.GetName() : nomName = 'WtoTau'
            elif 'data_obs' in key.GetName() : nomName = 'data_obs'
            elif 'Fake' in key.GetName() : nomName = 'Fake'
            elif 'LowAcc' in key.GetName() : nomName = 'LowAcc'
            elif 'DYJets' in key.GetName() : nomName = 'DYJets'
            else : 
                nomName = key.GetName().split('_')
                nomName = nomName[:-1]
                nomName = '_'.join(nomName)
            hnom = inFile.Get(nomName)
            weight = 1000*hvar.Integral()/hnom.Integral()
            hreshaped = hnom.Clone(hvar.GetName())
            hreshaped.Scale(weight)
            # outList.append(hreshaped)
            hreshaped.Write()
            
            hreshaped_down = hnom.Clone(hvar.GetName().replace('Up','Down'))
            hreshaped_down.Scale(1/weight)
            hreshaped_down.Write()
            
            print(key.GetName(), weight, 1/weight)
        elif key.GetName().endswith('Down'):
            continue 
        else :
            print(key.GetName())
            hnom = inFile.Get(key.GetName())
            # outList.append(hnom)
            hnom.Write()
    # outFile = ROOT.TFile('W'+s+'_reco_reshape.root',"recreate")
    # outFile.cd()
    # print("writing...")
    # for x in outList :
    #     x.Write()
    # outFile.Close()