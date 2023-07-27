import ROOT
import os
import sys
import h5py
import math
import numpy as np
import narf
# import jax
# import jax.numpy as jnp
import matplotlib.pyplot as plt
import mplhep as hep
sys.path.append('data/')
from binning import yBins, qtBins, ptBins, etaBins, mTBins, isoBins, chargeBins
import argparse

# matplotlib stuff
plt.style.use([hep.style.CMS])

parser = argparse.ArgumentParser('')
parser.add_argument('-asimov', '--asimov', default=False, action='store_true', help='plot asimov result')

args = parser.parse_args()
asimov = args.asimov
# file with fit results
f = ROOT.TFile.Open('../Fit/FitRes/fit_Wminus.root')

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
yBins_red = np.array(yBins[:threshold_y+1])
qtBins_red = np.array(qtBins[:threshold_qt+1])
yBinsC = 0.5*(yBins_red[1:]+yBins_red[:-1])
qtBinsC = 0.5*(qtBins_red[1:]+qtBins_red[:-1])

nBins=len(yBinsC)*len(qtBinsC)

# helXsecs = ['helXsec_L', 'helXsec_I', 'helXsec_T', 'helXsec_A', 'helXsec_P','helXsec_UL']
helXsecs = ['helXsec_UL']


# processes =[]
# for hel in helXsecs:
#     for i in range(len(yBinsC)):
#         for j in range(len(qtBinsC)):
#             proc = 'helXsecs' + hel + '_y_{}'.format(i)+'_qt_{}'.format(j)
#             processes.append(proc)

hcov = f.Get('covariance_matrix_channelmu')
cov = narf.root_to_hist(hcov).values()

hel_idx = {}
for hel in helXsecs:
    hel_idx[hel] = []
    for i in range(1,hcov.GetNbinsX()+1):
        if hel in hcov.GetXaxis().GetBinLabel(i):
            print(i,hcov.GetXaxis().GetBinLabel(i))
            hel_idx[hel].append(i)
        # print(i, hcov.GetXaxis().GetBinLabel(i))

impacts = []
impact_names = []
# retrieve impacts per group of POI
for i,hel in enumerate(helXsecs):
    #impact is generalization of per-nuisance impacts above v^T C^-1 v
    #where v is the matrix of poi x nuisance correlations within the group
    #and C is is the subset of the covariance matrix corresponding to the nuisances in the group
    vcovreduced = cov[-1,hel_idx[hel]]
    cov_hel = cov[hel_idx[hel],:]
    cov_hel = cov_hel[:,hel_idx[hel]]
    groupmcov = np.linalg.inv(cov_hel)
    vimpact = np.sqrt(np.matmul(np.matmul(vcovreduced.T,groupmcov),vcovreduced))
    print(vimpact*50.)
    impacts.append(vimpact*50.) # 50 MeV is input variation

# vcovreduced = cov[289-1,0:6*nBins]
# groupmcov = np.linalg.inv(cov[0:6*nBins,0:6*nBins])
# stat = np.sqrt(np.matmul(np.matmul(vcovreduced.T,groupmcov),vcovreduced))*50

# get stat impact
impact_names.extend(helXsecs)

# retrieve impacts per group of nuisance
hgroupimp = f.Get('nuisance_group_impact_nois')
groupimp = narf.root_to_hist(hgroupimp).values()
print(groupimp)
for i in range(groupimp.shape[1]):
    if 'mass' in hgroupimp.GetYaxis().GetBinLabel(i+1):
        continue
    impacts.append(groupimp[0,i]*50.) # 50 MeV is input variation
    impact_names.append(hgroupimp.GetYaxis().GetBinLabel(i+1))

print(np.sqrt(np.sum(np.dot(impacts,impacts))), 'total error summing impacts in quadrature')

plt.rcdefaults()
fig, ax = plt.subplots(figsize=(5, 5))
hep.cms.text('work in progress', loc=0, ax=ax)

y_pos = np.arange(len(impact_names))

ax.hlines(y_pos, -1*np.array(impacts), np.array(impacts), color='coral')
ax.plot(np.array(impacts), y_pos, '|', color='coral')  # Stem ends
ax.plot(-1*np.array(impacts), y_pos, '|', color='coral')  # Stem ends
ax.set_yticks(y_pos)
ax.set_yticklabels(impact_names, fontweight='medium')
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('impact on mass (MeV)', fontweight='medium')

# Add annotation to bars
for i in range(len(impacts)):
    if impact_names[i]=="luminosity":
        plt.text(0.,y_pos[i]-0.065, 
             str(round((impacts[i]), 2))+' MeV',
             fontsize=10,
             color='black', alpha=0.7)
    else:
        plt.text(0.,y_pos[i]-0.065, 
             str(round((impacts[i]), 1))+' MeV',
             fontsize=10,
             color='black', alpha=0.7)

plt.tight_layout()
plt.savefig('impactsOnMass_{}.png'.format("asimov" if asimov else "data"),dpi=300)
# plt.show()
