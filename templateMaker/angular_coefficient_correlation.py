import ROOT
import numpy as np
import matplotlib.pyplot as plt

#fIntoyy = ROOT.TFile.Open('../Fit/FitRes/fit_Wlike_nominal_asimov.root')
fIntoyy = ROOT.TFile.Open('../Fit/FitRes/fit_Wlike_iteration_0_Bnom.root')
corr = fIntoyy.Get('correlation_matrix_channelhelpois')

plot_title = 'correlation plot'

xaxis = corr.GetXaxis()
yaxis = corr.GetYaxis()

xlabels = []
ylabels = []

for i in range(1,xaxis.GetNbins()):
    xlabels.append(xaxis.GetBinLabel(i))
    ylabels.append(yaxis.GetBinLabel(i))

    
n = len(xlabels)


helicities = ['unpolarizedxsec','A0','A1','A2','A3','A4',]

unpolrange = range(40)
A0range    = range(40,80)
A1range    = range(80,120)
A2range    = range(120,160)
A3range    = range(160,200)
A4range    = range(200,240)
fullrange  = range(241)

XRange = A0range
YRange = XRange

arr = np.array([np.array([corr.GetBinContent(xproc,yproc) for xproc in XRange]) for yproc in YRange])

plt.rcParams['font.size'] = 30
plt.figure(figsize=(15,15))

im = plt.imshow(arr.T,origin='lower')
plt.xlabel('Unrolled process bins')
plt.ylabel('Unrolled process bins')
#plt.xticks(ticks= range(len(XRange)) ,labels=[xaxis.GetBinLabel(i) for i in XRange],rotation=45)
#plt.yticks(ticks= range(len(YRange)) ,labels=[yaxis.GetBinLabel(i) for i in YRange])
plt.colorbar(im,fraction=0.046, pad=0.04)
plt.title(plot_title)
plt.clim(-1,1)
plt.tight_layout()
plt.savefig('correlation_plot.png')
plt.show()