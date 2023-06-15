import ROOT
import numpy as np
import matplotlib.pyplot as plt


'''this script makes a plot of the correlation of each of the angular coefficients to mass'''

fIntoyy = ROOT.TFile.Open('../Fit/FitRes/fit_Wlike_nominal_asimov.root')
#fIntoyy = ROOT.TFile.Open('../Fit/FitRes/fit_Wlike_iteration_0_Bnom.root')
corr = fIntoyy.Get('correlation_matrix_channelhelpois')

xaxis = corr.GetXaxis()
yaxis = corr.GetYaxis()

xlabels = []
ylabels = []

for i in range(1,xaxis.GetNbins()):
    xlabels.append(xaxis.GetBinLabel(i))
    ylabels.append(yaxis.GetBinLabel(i))

    
n = len(xlabels)
helicities = ['unpolarizedxsec','A0','A1','A2','A3','A4',]


Mass = np.array([corr.GetBinContent(n+1,xproc) for xproc in range(n)])
plt.rcParams['font.size'] = 30
plt.figure(figsize=(15,10))
plt.plot(Mass,'o')

l = n/len(helicities)


plt.vlines(l,-1,1,'k')
plt.vlines(2*l,-1,1,'k')
plt.vlines(3*l,-1,1,'k')
plt.vlines(4*l,-1,1,'k')
plt.vlines(5*l,-1,1,'k')
plt.text(10,0.9,'unpol',size='x-small')
plt.text(l+20,0.9,'A0',size='x-small')
plt.text(2*l+20,0.9,'A1',size='x-small')
plt.text(3*l+20,0.9,'A2',size='x-small')
plt.text(4*l+20,0.9,'A3',size='x-small')
plt.text(5*l+20,0.9,'A4',size='x-small')
plt.xlabel('Unrolled process bins')
plt.ylabel('Correlation')
plt.title('Correlation to mass')
plt.ylim(-1,1)
plt.tight_layout()
plt.savefig('correlation_to_mass.png')