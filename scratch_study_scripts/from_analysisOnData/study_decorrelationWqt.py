import ROOT
import math 

ff = ROOT.TFile.Open("WToMu_plots.root")

h3_nom_dec1 = ff.Get("templates_Signal/Nom_decQt0_5/templates") 
h3_var_dec1 = ff.Get("templates_Signal/LHEScaleWeight_decQt0_5/templates_LHEScaleWeight_muR0p5_muF0p5") 
h3_nom_dec2 = ff.Get("templates_Signal/Nom_decQt5_15/templates") 
h3_var_dec2 = ff.Get("templates_Signal/LHEScaleWeight_decQt5_15/templates_LHEScaleWeight_muR0p5_muF0p5") 
h3_nom_dec3 = ff.Get("templates_Signal/Nom_decQt15_32/templates") 
h3_var_dec3 = ff.Get("templates_Signal/LHEScaleWeight_decQt15_32/templates_LHEScaleWeight_muR0p5_muF0p5") 
h3_nom = ff.Get("templates_Signal/Nominal/templates") 
h3_var = ff.Get("templates_Signal/LHEScaleWeight/templates_LHEScaleWeight_muR0p5_muF0p5") 


#nominal
h3_nom_dec = h3_nom_dec1.Clone("h3_nom_dec")
h3_nom_dec.Add(h3_nom_dec2)
h3_nom_dec.Add(h3_nom_dec3)

#variation in region low qt
h3_var_dec = h3_var_dec1.Clone("h3_var_dec")
h3_var_dec.Add(h3_nom_dec2)
h3_var_dec.Add(h3_nom_dec3)

#variation all
h3_var_dec_all = h3_var_dec1.Clone("h3_var_dec_all")
h3_var_dec_all.Add(h3_var_dec2)
h3_var_dec_all.Add(h3_var_dec3)


h1_nom_dec = h3_nom_dec.Project3D("h3_nom_dec_y")
h1_var_dec = h3_var_dec.Project3D("h3_var_dec_y")
h1_nom = h3_nom.Project3D("h3_nom_y")
h1_var = h3_var.Project3D("h3_var_y")
h1_var_dec_all = h3_var_dec_all.Project3D("h3_var_dec_all_y")

# h1_nom_dec1 = h3_nom_dec1.Project3D("h3_nom_dec1_y")
# h1_var_dec1 = h3_var_dec1.Project3D("h3_var_dec1_y")
# h1_nom_dec2 = h3_nom_dec2.Project3D("h3_nom_dec2_y")
# h1_var_dec2 = h3_var_dec2.Project3D("h3_var_dec2_y")
# h1_nom_dec3 = h3_nom_dec3.Project3D("h3_nom_dec3_y")
# h1_var_dec3 = h3_var_dec3.Project3D("h3_var_dec3_y")



out = ROOT.TFile("tester.root", "recreate")
out.cd()
h1_nom_dec.SetName("h1_nom_dec")
h1_var_dec.SetName("h1_var_dec")
h1_nom.SetName("h1_nom")
h1_var.SetName("h1_var")
h1_var_dec_all.SetName("h1_var_dec_all")

h1_nom_dec.Write()
h1_var_dec.Write()
h1_nom.Write()
h1_var.Write()
h1_var_dec_all.Write()

# h1_nom_dec1.SetName("h1_nom_dec1")
# h1_var_dec1.SetName("h1_var_dec1")
# h1_nom_dec2.SetName("h1_nom_dec2")
# h1_var_dec2.SetName("h1_var_dec2")
# h1_nom_dec3.SetName("h1_nom_dec3")
# h1_var_dec3.SetName("h1_var_dec3")

# h1_nom_dec1.Write()
# h1_var_dec1.Write()
# h1_nom_dec2.Write()
# h1_var_dec2.Write()
# h1_nom_dec3.Write()
# h1_var_dec3.Write()