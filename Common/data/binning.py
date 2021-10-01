import ROOT
from math import pi

# ROOT.gInterpreter.ProcessLine('auto qtBins = std::make_tuple<float,float,float,float,float,float,float,float,float,float,float,float,float>(0., 4., 8., 12., 16., 20., 24., 28., 32., 40., 60., 100., 200.);')
# ROOT.gInterpreter.ProcessLine('auto yBins = std::make_tuple<float,float,float,float,float,float,float,float,float>(0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 6.0);')
# qtBins = [0.,2.,4.,6.,8.,10.,12.,16.,20.,26.,36.,60.,100.,200.] #Valerio's binning
qtBins = [0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.,200.]
qtBins_syst = [0., 5., 10.,15.,200.]
# qtBins = [0.,2.14725709,3.75887602,5.0176509,6.35604492,7.97738581,9.8266591,12.08973559,15.1689846,19.09257703,24.54519361,33.19566976,60.,200.]
qtBins_val = [0.,2.53766747,4.21661234,5.70425537,7.53520193,9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.,200.]
yBins = [0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0]
yBins_val =[0., 0.4, 0.8, 1.2 ,1.6, 2.0, 2.4, 2.8, 6.]
# qtBins_val =[0.,  2.,  4., 6.,  8., 10., 12., 16., 20., 26., 36., 60., 200.]
# ptBins = [25.+i for i in range(31)]
ptBins = [25.+i/2. for i in range(61)]
#overflow-Bin
#ptBins.append(200.)
etaBins = [-2.4+i*0.1 for i in range(49)]
chargeBins = [-2. +i*2. for i in range(3)]
mTBins = [0.,40.,1000.]
#mTBins = [0.,1000.]
mTBinsFull = [i*2. for i in range(76)]
metBins = [i*2. for i in range(76)]
isoBins = [0.,0.15,1.]
#isoBins = [0.,100.]
zmassBins = [70. + i*1. for i in range(41)]
pvBins = [9.5 + i*1. for i in range(50)]
cosThetaBins = [round(-1. + 2.*i/100,2) for i in range(101)]
phiBins = [round(-pi + pi*i/100,2) for i in range(101)]

