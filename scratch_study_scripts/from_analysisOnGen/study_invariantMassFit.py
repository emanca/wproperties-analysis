import ROOT
ROOT.gROOT.SetBatch()

mass_low=79
mass_high = 82
FIXMG = False
shift_const = 6.76 #0.00676; //6.76 MeV, gamma^2/(8m)
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

GENLEVEL = True



######################## INV. MASS PROJECTOR, no CUTs #########################
# ROOT.EnableImplicitMT()
# df = ROOT.RDataFrame("Events","/scratchnvme/wmass/WJetsNoCUT_v2/tree_*_*.root")
# histoMass = df.Filter("Wmass_preFSR>75 && Wmass_preFSR<85").Histo1D("Wmass_preFSR")

# outFile = ROOT.TFile("invariantMass.root", "recreate")
# outFile.cd()
# histoMass.Write()
######################## INV. MASS PROJECTOR, no CUTs #########################

if not GENLEVEL :
    inFile =  ROOT.TFile.Open("WToMu_plots.root")
    h2 = inFile.Get('prefit_Signal/Nominal/Wmass_preFSR')

else :
    # inFile = ROOT.TFile.Open("invariantMass.root")
    inFile = ROOT.TFile.Open("genInfo_mass.root")

    

outFile = ROOT.TFile("invariantMassFit_GEN_Down.root", "recreate")


signDict = {
    'minus' : 1,
    'plus' : 2,
}

hDict = {}
fDict = {}
lDict = {}
cDict = {}
lineDict = {}

kList = ['bw','jac','mod']

for s,sign in signDict.items() :
    
    for k in kList : 
        if not GENLEVEL :
            hDict[k+s] = h2.ProjectionX(k+s,sign,sign)
        else :
        #    htemp = inFile.Get("Wmass_preFSR_"+s)
        #    htemp = inFile.Get("angularCoefficients_W"+s+"/Wmass_preFSR")
           htemp = inFile.Get("angularCoefficients_W"+s+"_mass/Wmass_preFSR_massDown")
           hDict[k+s] = htemp.Clone(k+s)
        hDict[k+s].SetStats(0)

    fDict['bw'+s]  = ROOT.TF1("fit_bw","[2]/(pow((pow(x,2)-[0]*[0]),2)+[0]*[0]*[1]*[1])",mass_low,mass_high)
    fDict['jac'+s] = ROOT.TF1("fit_jac","[2]*x/(pow((pow(x,2)-[0]*[0]),2)+[0]*[0]*[1]*[1])",mass_low,mass_high)
    fDict['mod'+s]  = ROOT.TF1("fit_mod","[2]*(x+[3]*(x-[0])*x/[0]+[4]*pow((x-[0]),2)*x/pow([0],2))/(pow((pow(x,2)-[0]*[0]),2)+[0]*[0]*[1]*[1])",mass_low,mass_high)
    
    fDict['bw'+s].SetParameters(80.419,2.085,hDict['bw'+s].GetEntries())
    fDict['jac'+s].SetParameters(80,2,hDict['bw'+s].GetEntries())
    fDict['mod'+s].SetParameters(80.419,2.085,hDict['bw'+s].GetEntries(),-1.5,-1.5)
    if FIXMG :
        fDict['mod'+s].FixParameter(0,80.419)
        fDict['mod'+s].FixParameter(1,2.047)
        
    fDict['bw'+s].SetParNames("M","#Gamma","A")
    fDict['mod'+s].SetParNames("M","#Gamma","A","H","K")
    fDict['jac'+s].SetParNames("M","#Gamma","A")
    
    for k in kList : 
        hDict[k+s].Fit(fDict[k+s],"","",mass_low,mass_high)
        # shift_value_wp = make_pair((fDict['mod'+s].GetParameter(3)+1)*shift_const, (fDict['mod'+s].GetParError(3))*shift_const)
        hDict[k+s+'pull'] = ROOT.TH1F(k+s+'pull',k+s+'pull',hDict[k+s].GetNbinsX(),75,85)
        hDict[k+s+'pull'].SetStats(0)


        for m in range(1,hDict['bw'+s].GetNbinsX()+1) :
            
            hDict[k+s+'pull'].SetBinContent(m,(hDict['bw'+s].GetBinContent(m)-fDict[k+s].Eval(hDict['bw'+s].GetBinCenter(m)))/hDict['bw'+s].GetBinContent(m))

        # print("--------------------COVARIANCE MATRIX INFORMATIONS (W+)------------------")
        # resFit.Print("V")

    fDict['bw'+s].SetLineColor(ROOT.kRed)
    fDict['mod'+s].SetLineColor(ROOT.kGreen)
    fDict['jac'+s].SetLineColor(ROOT.kBlue)
    hDict['bw'+s+'pull'].SetLineColor(ROOT.kRed)
    hDict['mod'+s+'pull'].SetLineColor(ROOT.kGreen)
    hDict['jac'+s+'pull'].SetLineColor(ROOT.kBlue)
    hDict['bw'+s].SetLineColor(ROOT.kBlack)
    hDict['jac'+s].SetLineColor(ROOT.kBlack)
    hDict['mod'+s].SetLineColor(ROOT.kBlack)

    hDict['bw'+s+'pull'].SetLineWidth(1)
    hDict['mod'+s+'pull'].SetLineWidth(1)
    hDict['jac'+s+'pull'].SetLineWidth(1)
    hDict['bw'+s+'pull'].SetMarkerColor(ROOT.kRed)
    hDict['mod'+s+'pull'].SetMarkerColor(ROOT.kGreen)
    hDict['jac'+s+'pull'].SetMarkerColor(ROOT.kBlue)
    hDict['bw'+s+'pull'].SetMarkerStyle(33)
    hDict['mod'+s+'pull'].SetMarkerStyle(33)
    hDict['jac'+s+'pull'].SetMarkerStyle(33)
    hDict['bw'+s].GetXaxis().SetLabelOffset(1)

    hDict['bw'+s].GetYaxis().SetTitle("Events/%.3f GeV"%(hDict['bw'+s].GetBinWidth(10)))
    hDict['bw'+s].SetTitle("")
    hDict['bw'+s].SetMaximum(400000)
    hDict['bw'+s].GetXaxis().SetRangeUser(75,85)
    hDict['bw'+s+'pull'].GetXaxis().SetTitle("Q [GeV]")
    hDict['bw'+s+'pull'].GetYaxis().SetTitle("#frac{N_{MC}-f_{fit}}{N_{MC}}")
    hDict['bw'+s+'pull'].GetXaxis().SetRangeUser(75,85)
    hDict['bw'+s+'pull'].SetTitle("")
    hDict['bw'+s+'pull'].GetYaxis().SetRangeUser(-0.1,0.1)
    hDict['bw'+s+'pull'].GetXaxis().SetTitleSize(0.13)
    hDict['bw'+s+'pull'].GetXaxis().SetLabelSize(0.1)
    hDict['bw'+s+'pull'].GetYaxis().SetTitleSize(0.09)
    hDict['bw'+s+'pull'].GetYaxis().SetTitleOffset(0.48)
    hDict['bw'+s+'pull'].GetYaxis().SetLabelSize(0.08)

    lDict[s] =  ROOT.TLegend(0.30,0.10,0.80,0.4)
    lDict[s].AddEntry(fDict['bw'+s], "#splitline{Breit-Wigner}{#splitline{M=%.4f#pm%.4f}{#splitline{#Gamma=%.3f#pm%.3f}{#chi^{2}/NDF=%.0f/%d}}}"%(fDict['bw'+s].GetParameter(0),fDict['bw'+s].GetParError(0), fDict['bw'+s].GetParameter(1),fDict['bw'+s].GetParError(1),fDict['bw'+s].GetChisquare(),fDict['bw'+s].GetNDF()))
    lDict[s].AddEntry(fDict['jac'+s], "#splitline{Breit-Wigner with Jacobian}{#splitline{M=%.4f#pm%.4f}{#splitline{#Gamma=%.3f#pm%.3f}{#chi^{2}/NDF=%.0f/%d}}}"%(fDict['jac'+s].GetParameter(0),fDict['jac'+s].GetParError(0), fDict['jac'+s].GetParameter(1),fDict['jac'+s].GetParError(1),fDict['jac'+s].GetChisquare(),fDict['jac'+s].GetNDF()))
    lDict[s].AddEntry(fDict['mod'+s], "#splitline{Modified Breit-Wigner}{#splitline{M=%.4f#pm%.4f}{#splitline{#Gamma=%.3f#pm%.3f}{#splitline{H=%.2f#pm%.2f}{#splitline{K=%.0f#pm%.0f}{#chi^{2}/NDF=%.0f/%d}}}}}"%(fDict['mod'+s].GetParameter(0),fDict['mod'+s].GetParError(0), fDict['mod'+s].GetParameter(1),fDict['mod'+s].GetParError(1), fDict['mod'+s].GetParameter(3),fDict['mod'+s].GetParError(3), fDict['mod'+s].GetParameter(4),fDict['mod'+s].GetParError(4),fDict['mod'+s].GetChisquare(),fDict['mod'+s].GetNDF()))
    
    if s=='plus' :
        lDict[s].SetHeader("pp#rightarrow W^{+}+X")
    else :
        lDict[s].SetHeader("pp#rightarrow W^{-}+X")

    #.............. canvas .....................#
    
    lineDict["sx"] =  ROOT.TLine(79,0, 79,400000)
    lineDict["dx"] =  ROOT.TLine(82,0, 82,400000)
    lineDict["sx"].SetLineColor(ROOT.kViolet)
    lineDict["sx"].SetLineWidth(3)
    lineDict["sx"].SetLineStyle(2)
    lineDict["dx"].SetLineColor(ROOT.kViolet)
    lineDict["dx"].SetLineWidth(3)
    lineDict["dx"].SetLineStyle(2)

    lineDict["pull_sx"] =  ROOT.TLine(79,-0.1, 79,0.1)
    lineDict["pull_dx"] =  ROOT.TLine(82,-0.1, 82,0.1)
    lineDict["pull_sx"].SetLineColor(ROOT.kViolet)
    lineDict["pull_sx"].SetLineWidth(3)
    lineDict["pull_sx"].SetLineStyle(2)
    lineDict["pull_dx"].SetLineColor(ROOT.kViolet)
    lineDict["pull_dx"].SetLineWidth(3)
    lineDict["pull_dx"].SetLineStyle(2)


    cDict[s] = ROOT.TCanvas("c_mass"+s,"c_mass"+s,800,700)
    cDict[s].cd()
    
    cDict['h'+s] =  ROOT.TPad("pad_histo"+s,"pad_histo"+s,0,0.3,1,1)
    cDict['pull'+s] =  ROOT.TPad("pad_pull"+s,"pad_pull"+s,0,0,1,0.3)
    cDict['h'+s].SetGridx()
    cDict['h'+s].SetGridy()
    cDict['pull'+s].SetGridx()
    cDict['pull'+s].SetGridy()
    
    cDict['h'+s].SetBottomMargin(0.02)
    cDict['h'+s].Draw()
    cDict['pull'+s].SetTopMargin(0)
    cDict['pull'+s].SetBottomMargin(0.25)
    cDict['pull'+s].Draw()
    
    cDict['h'+s].cd()
    hDict['bw'+s].Draw()
    hDict['jac'+s].Draw("SAME")
    hDict['mod'+s].Draw("SAME")
    fDict['bw'+s].Draw("SAME")
    fDict['jac'+s].Draw("SAME")
    fDict['mod'+s].Draw("SAME")
    lineDict["dx"].Draw("SAME")
    lineDict["sx"].Draw("SAME")
    lDict[s].Draw("SAME")


    cDict['pull'+s].cd()
    hDict['bw'+s+'pull'].Draw()
    hDict['mod'+s+'pull'].Draw("SAME")
    hDict['jac'+s+'pull'].Draw("SAME")
    lineDict["pull_dx"].Draw("SAME")
    lineDict["pull_sx"].Draw("SAME")


    outFile.cd()
    cDict[s].Write()        
        
