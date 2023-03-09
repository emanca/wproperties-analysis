import ROOT
import os
import sys
import h5py
import math
import numpy as np
# import jax
# import jax.numpy as jnp
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import mplhep as hep
sys.path.append('data/')

# matplotlib stuff
plt.style.use([hep.style.CMS])
plt.rcdefaults()

fig, ax = plt.subplots()
# hep.cms.text('work in progress', loc=0, ax=ax)

exp_names = ['ALEPH','DELPHI','L3','OPAL','D0','ATLAS','LHCb','CDF', 'Electroweak Fit']
central = [80440,80336,80270,80415, 80367, 80370, 80354, 80433.5, 80350]
stat_errs = [43, 55, 46, 42,  11, 7, 23, 6.4,0.]
tot_errs = [50,67, 55, 52, 26, 19, 32, 9.4,8.]

y_pos = np.arange(len(exp_names))

ax.hlines(y_pos, central-np.array(tot_errs), central+np.array(tot_errs), color='blue', label="total uncertainty")
ax.hlines(y_pos, central-np.array(stat_errs), central+np.array(stat_errs), color='deepskyblue', label="statistical uncertainty")
ax.axvspan(80350-8, 80350+8, alpha=0.1, color='green')
ax.plot(central, y_pos, 'o', color='blue')  # Stem ends
ax.plot(central-np.array(tot_errs), y_pos, '|', color='blue')  # Stem ends
ax.plot(central+np.array(tot_errs), y_pos, '|', color='blue')  # Stem ends

ax.set_yticks(y_pos)
ax.set_yticklabels(exp_names, fontweight='medium')
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('$m_W$ (MeV)', fontweight='medium')
ax.set_xticks(np.arange(80200, 80550, 50))

# # Add annotation to bars
# for i in range(len(impacts)):
#     if impact_names[i]=="luminosity":
#         plt.text(0.,y_pos[i]-0.065, 
#              str(round((impacts[i]), 2))+' MeV',
#              fontsize=10,
#              color='black', alpha=0.7)
#     else:
#         plt.text(0.,y_pos[i]-0.065, 
#              str(round((impacts[i]), 1))+' MeV',
#              fontsize=10,
#              color='black', alpha=0.7)
ax.legend(loc='lower left',bbox_to_anchor=(0.0,1.01), borderaxespad=0, frameon=False)
plt.tight_layout()
plt.savefig('massSummary.pdf')
