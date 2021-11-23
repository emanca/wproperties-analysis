import ROOT
import math

infile = ROOT.TFile.Open("test_WqtSyst/plot/stackPlots.root")

canvas = infile.Get("prefit_Mu1_eta_plus")
pad = canvas.GetPrimitive("pad_histo_Mu1_eta_plus")
stack = pad.GetPrimitive("stack_Mu1_eta_plus")
stack_list = stack.GetHists()
#histo_W = stack_list.At(5)
h_mc = stack_list.At(0) 
for h in stack_list :
    h_mc.Add(h)
h_data = pad.GetPrimitive("hData_Mu1_eta_plus")
int_data = h_data.Integral()
int_mc = h_mc.Integral()
print("PLUS", int_data, int_mc, int_mc/int_data)

canvas = infile.Get("prefit_Mu1_eta_minus")
pad = canvas.GetPrimitive("pad_histo_Mu1_eta_minus")
stack = pad.GetPrimitive("stack_Mu1_eta_minus")
stack_list = stack.GetHists()
#histo_W = stack_list.At(5)
h_mc = stack_list.At(0) 
for h in stack_list :
    h_mc.Add(h)
h_data = pad.GetPrimitive("hData_Mu1_eta_minus")
int_data = h_data.Integral()
int_mc = h_mc.Integral()
print("MINUS", int_data, int_mc, int_mc/int_data)