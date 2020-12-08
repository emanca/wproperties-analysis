import ROOT
import math
import argparse
import ctypes
import numpy as np
ROOT.gROOT.SetBatch()
np.set_printoptions(linewidth=np.inf,precision=1)

###############################################################################################
#                                      RESULTS
#
#
# With acceptance: python quantileBinBuilder.py --ycut 1 --qtcut 1
# -----------Y quantiles --------------
# edges Plus = [0.  0.4 0.8 1.2 1.6 2.  2.4]
# centers Plus = [0.2 0.5 0.9 1.4 1.9 2.2 0. ]
# edges Minus = [0.  0.4 0.8 1.1 1.5 2.  2.4]
# centers Minus = [0.2 0.5 0.9 1.3 1.7 2.1 0. ]
# -----------Qt quantiles --------------
# edges Plus = [ 0.   4.2  7.4 11.8 18.5 32. ]
# centers Plus = [ 2.6  5.8  9.5 14.8 24.3  0. ]
# edges Minus = [ 0.   4.2  7.5 11.8 18.5 32. ]
# centers Minus = [ 2.6  5.8  9.5 14.9 24.3  0.]
#
#
# Without acceptance: python quantileBinBuilder.py --ycut 0 --qtcut 0 -ny 9 -nq 6
# -----------Y quantiles --------------
# edges Plus = [0.  0.4 0.9 1.3 1.8 2.2 2.7 3.1 3.6 6.1]
# centers Plus = [0.2 0.6 1.1 1.5 2.  2.4 2.9 3.4 4.  0. ]
# edges Minus = [0.  0.4 0.7 1.1 1.5 1.8 2.3 2.7 3.3 6.1]
# centers Minus = [0.2 0.5 0.9 1.2 1.6 2.  2.5 3.  3.7 0. ]
# -----------Qt quantiles --------------
# edges Plus = [  0.    4.    7.   11.   17.4  30.3 200. ]
# centers Plus = [ 2.5  5.4  8.9 13.9 22.8 54.9  0. ]
# edges Minus = [  0.    4.    7.1  11.2  17.6  30.4 200. ]
# centers Minus = [ 2.5  5.5  9.  14.1 23.  54.4  0. ]
###############################################################################################

parser = argparse.ArgumentParser("")

parser.add_argument('-o','--output', type=str, default='quantileBinBuilder',help="name of the output file")
parser.add_argument('-i','--input', type=str, default='../GenInfo/genInfo_fine4quantiles.root',help="name of the fit result root file, after plotter_fitResult")
parser.add_argument('-y','--ycut', type=int, default=True,help="cut at Y= Y acceptance (2.4)")
parser.add_argument('-q','--qtcut', type=int, default=True,help="cut at qt= qt acceptance (32 GeV)")
parser.add_argument('-ny','--ny', type=int, default=6,help="number of y quantiles")
parser.add_argument('-nq','--nqt', type=int, default=5,help="number of qt quantiles")
parser.add_argument('-c','--centers', type=int, default=True,help="evaluate the bin centers, given the quantiles")

args = parser.parse_args()
OUTPUT = args.output
INPUT = args.input
YCUT = args.ycut
QTCUT = args.qtcut
NY = args.ny
NQ = args.nqt
CENT = args.centers

inFile = ROOT.TFile.Open(INPUT)

h2 = {}
h2['Plus'] = inFile.Get('angularCoefficients_Wplus/mapTot')
h2['Minus'] = inFile.Get('angularCoefficients_Wminus/mapTot')

binYcut=-1
binQtcut=-1

if YCUT : binYcut = h2['Plus'].GetXaxis().FindBin(2.39)
if QTCUT : binQtcut = h2['Plus'].GetYaxis().FindBin(31.99)

print("binYcut=",binYcut, "  binQtcut=",binQtcut)
h1 = {}
qDict = {}
cDict= {}


print("-----------Y quantiles --------------")
#definition of quantiles (=[0,1] divided in NQ bins)
xqy = np.zeros((NY))
for i in range(0,NY) : xqy[i] = (i+1)/NY
print("quantiles x=",xqy)

for s, histo2 in h2.items() : 
    histo2.GetXaxis().SetRange(0,binYcut)#y cut 
    h1['y'+s+'int'] = histo2.ProjectionX('h1_y'+s,0,binQtcut) #qt cut
    
    qDict['y'+s+'int'] = np.zeros((NY))
    h1['y'+s+'int'].GetQuantiles(NY,qDict['y'+s+'int'],xqy)

    # for iqt in range(1,histo2.GetNbinsY()+1,100) :
    #     h1['y'+s+str(iqt)] = histo2.ProjectionX('h1_y'+s+str(iqt),iqt,iqt+100)
        
    #     h1['y'+s+str(iqt)].GetQuantiles(NY,yq,xqy)
    #     qDict['y'+s+str(iqt)] = yq.copy()
    
    if CENT :
        cDict['y'+s+'int'] = np.zeros((NY+1))
        qDict['y'+s+'int'] = np.insert(qDict['y'+s+'int'],0,0)
        for b in range(0,np.size(qDict['y'+s+'int'])-1) :
            h1['y'+s+'int'].GetXaxis().SetRangeUser(qDict['y'+s+'int'][b],qDict['y'+s+'int'][b+1])
            cDict['y'+s+'int'][b] = h1['y'+s+'int'].GetMean()
    h1['y'+s+'int'].GetXaxis().SetRange(0,binYcut)#to later on plotting

for s, histo2 in h2.items() :
    print("edges", s,"=", qDict['y'+s+'int'])
    print("centers",s, "=",cDict['y'+s+'int'])



print("-----------Qt quantiles --------------")
#definition of quantiles (=[0,1] divided in NQ bins)
xqqt = np.zeros((NQ))
for i in range(0,NQ) : xqqt[i] = (i+1)/NQ
print("quantiles x=",xqqt)

for s, histo2 in h2.items() :
    histo2.GetYaxis().SetRange(0,binQtcut)#qtcut
    h1['qt'+s+'int'] = histo2.ProjectionY('h1_qt'+s,0,binYcut)#ycut
    
    qDict['qt'+s+'int'] = np.zeros((NQ))
    h1['qt'+s+'int'].GetQuantiles(NQ,qDict['qt'+s+'int'],xqqt)

    # for iy in range(1,histo2.GetNbinsX()+1,4) :
    #     h1['qt'+s+str(iy)] = histo2.ProjectionY('h1_qt'++s+str(iy),iy,iy+4)
        
    #     h1['qt'+s+str(iy)].GetQuantiles(NQ,yq,xqqt)
    #     qDict['qt'+s+str(iy)] = yq.copy()
    
    if CENT :
        cDict['qt'+s+'int'] = np.zeros((NQ+1))
        qDict['qt'+s+'int'] = np.insert(qDict['qt'+s+'int'],0,0)
        for b in range(0,np.size(qDict['qt'+s+'int'])-1) :
            h1['qt'+s+'int'].GetXaxis().SetRangeUser(qDict['qt'+s+'int'][b],qDict['qt'+s+'int'][b+1])
            cDict['qt'+s+'int'][b] = h1['qt'+s+'int'].GetMean()
    h1['qt'+s+'int'].GetXaxis().SetRange(0,binQtcut)#to later on plotting
    
for s, histo2 in h2.items() :
    print("edges",s,"=", qDict['qt'+s+'int'])
    print("centers",s, "=",cDict['qt'+s+'int'])


outFile = ROOT.TFile(OUTPUT+".root", "recreate")
outObj = {}
canvas = {}
leg = {}

# building the output histograms
klist = ['qt','y']
for s, histo2 in h2.items() :
    for k in klist : 
        qDict[k+s+'int'] = np.round(qDict[k+s+'int'],decimals=1)
        cDict[k+s+'int'] = np.round(cDict[k+s+'int'],decimals=1)
                
        outObj[k+s+'int'] = h1[k+s+'int'].Clone(k+s)
        outObj[k+s+'int'].SetTitle(k+s)
        outObj[k+s+'reb'] = h1[k+s+'int'].Rebin(np.size(qDict[k+s+'int'])-1,k+s+'reb', qDict[k+s+'int'])
        outObj[k+s+'reb'].SetTitle(k+s+'reb')
        
        #rescaling to match the shape
        for ib in range(1,outObj[k+s+'reb'].GetNbinsX()+1) :
            val = outObj[k+s+'reb'].GetBinContent(ib)/outObj[k+s+'reb'].GetBinWidth(ib)*outObj[k+s+'int'].GetBinWidth(1)
            outObj[k+s+'reb'].SetBinContent(ib,val)
        
        cDict[k+s+'val'] =  np.zeros(np.size(cDict[k+s+'int'])-1)
        cDict[k+s+'err'] =  np.zeros(np.size(cDict[k+s+'int'])-1)
        cDict[k+s+'errX'] =  np.zeros(np.size(cDict[k+s+'int'])-1)
        for ib in range(1,outObj[k+s+'reb'].GetNbinsX()+1) :
             cDict[k+s+'val'][ib-1] = outObj[k+s+'reb'].GetBinContent(ib)
             cDict[k+s+'err'][ib-1] = outObj[k+s+'reb'].GetBinError(ib)
        outObj[k+s+'g'] = ROOT.TGraphErrors(np.size(cDict[k+s+'int'])-1,cDict[k+s+'int'][:-1], cDict[k+s+'val'], cDict[k+s+'errX'],cDict[k+s+'err'])
        outObj[k+s+'g'].SetName(k+s+'G')
        outObj[k+s+'g'].SetTitle(k+s+'G')

#canvas
for s, histo2 in h2.items() :
    for k in klist : 
        canvas[k+s] = ROOT.TCanvas(k+s,k+s)
        canvas[k+s].SetGridx()
        canvas[k+s].SetGridy()
        leg[k+s] = ROOT.TLegend(0.5,0.75,0.9,0.9)
        
        outObj[k+s+'int'].SetLineColor(ROOT.kGreen+2)
        outObj[k+s+'int'].SetLineWidth(3)
        outObj[k+s+'int'].SetStats(0)
        outObj[k+s+'int'].GetYaxis().SetTitle("Events")
        if k=='qt' :
            outObj[k+s+'int'].GetXaxis().SetTitle("q_{T}^{W} [GeV]")
        if k=='y' :
            outObj[k+s+'int'].GetXaxis().SetTitle("Y_{W}")
            
        outObj[k+s+'reb'].SetLineColor(ROOT.kBlue+1)
        outObj[k+s+'reb'].SetLineWidth(3)
        
        outObj[k+s+'g'].SetMarkerStyle(21)
        outObj[k+s+'g'].SetMarkerColor(ROOT.kRed+3)
        outObj[k+s+'g'].SetLineColor(ROOT.kRed+1)
        outObj[k+s+'g'].SetLineWidth(4)
        
        outObj[k+s+'int'].Draw()
        outObj[k+s+'reb'].Draw("same")
        outObj[k+s+'g'].Draw("same P")
        
        leg[k+s].AddEntry(outObj[k+s+'int'], "original distirbuion")
        leg[k+s].AddEntry(outObj[k+s+'reb'], "rebinned, mid-bin centers")
        leg[k+s].AddEntry(outObj[k+s+'g'], "mean value bin centers")
        leg[k+s].Draw("SAME")
        
outFile.cd()
for i,c in canvas.items() :
    c.Write()