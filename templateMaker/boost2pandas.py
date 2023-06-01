

import ROOT
import narf
import pandas as pd
import h5py 
import hist
import hdf5plugin
import math
import boost_histogram as bh
import numpy as np
import matplotlib.pyplot as plt
import re
from collections import OrderedDict
import pdb

def writeFlatInChunks(arr, h5group, outname, maxChunkBytes = 1024**2):    
  arrflat = arr.reshape(-1)
  
  esize = np.dtype(arrflat.dtype).itemsize
  nbytes = arrflat.size*esize

  #special handling for empty datasets, which should not use chunked storage or compression
  if arrflat.size == 0:
    chunksize = 1
    chunks = None
    compression = None
  else:
    chunksize = int(min(arrflat.size,max(1,math.floor(maxChunkBytes/esize))))
    chunks = (chunksize,)
    compression = "gzip"

  h5dset = h5group.create_dataset(outname, arrflat.shape, chunks=chunks, dtype=arrflat.dtype, compression=compression)

  #write in chunks, preserving sparsity if relevant
  for ielem in range(0,arrflat.size,chunksize):
    aout = arrflat[ielem:ielem+chunksize]
    if np.count_nonzero(aout):
      h5dset[ielem:ielem+chunksize] = aout
      
  h5dset.attrs['original_shape'] = np.array(arr.shape,dtype='int64')

  return nbytes

def fillHelGroup(yBinsC,qtBinsC,helXsecs):
    helGroups = OrderedDict()
    for i in range(len(yBinsC)):
        for j in range(len(qtBinsC)):
            s = 'y_{i}_qt_{j}'.format(i=round(yBinsC[i],1),j=round(qtBinsC[j],1))
            helGroups[s] = []
            
            for hel in helXsecs:
                helGroups[s].append('helXsec_'+hel+'_'+s)
            if helGroups[s] == []:
                del helGroups[s]
    return helGroups

def fillHelMetaGroup(yBinsC,qtBinsC,sumGroups):
    helMetaGroups = OrderedDict()
    for i in range(len(yBinsC)):
        s = 'y_{i}'.format(i=round(yBinsC[i],1))
        helMetaGroups[s] = []
        for key in sumGroups:
            if s in key:
                helMetaGroups[s].append(key)
        
        if helMetaGroups[s] == []:
                del helMetaGroups[s]
    
    for j in range(len(qtBinsC)):
        s = 'qt_{j}'.format(j=round(qtBinsC[j],1))
        helMetaGroups[s] = []
        for key in sumGroups:
            if 'qt' in key and key.split('_')[3]==str(round(qtBinsC[j],1)):
                helMetaGroups[s].append(key)
    
        if helMetaGroups[s] == []:
                del helMetaGroups[s]
    return helMetaGroups

def fillSumGroup(yBinsC,qtBinsC,helXsecs,processes):
    sumGroups = OrderedDict()
    for i in range(len(yBinsC)):
        s = 'y_{i}'.format(i=round(yBinsC[i],1))
        for hel in helXsecs:
            for signal in processes:
                sumGroups['helXsec_'+hel+'_'+s] = []
                for j in range(len(qtBinsC)):
                    #if 'helXsecs'+hel+'_'+'y_{i}_qt_{j}'.format(i=i,j=j) in processes:
                    sumGroups['helXsec_'+hel+'_'+s].append('helXsec_'+hel+'_'+s+'_qt_{j}'.format(j=round(qtBinsC[j],1)))
    
    for j in range(len(qtBinsC)):
        s = 'qt_{j}'.format(j=round(qtBinsC[j],1))
        for hel in helXsecs:
            for signal in processes:
                if signal.split('_')[0]+ '_' + signal.split('_')[1]== 'helXsec_'+hel and signal.split('_')[-1] == str(round(qtBinsC[j],1)):
                    sumGroups['helXsec_'+hel+'_'+s] = []
                    for i in range(len(yBinsC)):
                        #print i, signal, 'helXsec_'+hel+'_'+'y_{i}_pt_{j}'.format(i=i,j=j)
                        #print 'append', 'helXsec_'+hel+'_y_{i}_'.format(i=i)+s, 'to', 'helXsec_'+hel+'_'+s
                        sumGroups['helXsec_'+hel+'_'+s].append('helXsec_'+hel+'_y_{i}_'.format(i=round(yBinsC[i],1))+s)
    return sumGroups

def mirrorHisto(nom,var):
    '''
    Parameters
    ==========
    nom: nominal boost histogram
    var: boost histogram corresponding to systematic variation
    Returns
    =======
    Mirrored Histogram: Boost histogram with new two dimensional axis labeled downUpVar. Index '0' corresponds to 
    'down' variation defined as nom/var and index '1' corresponds to 'up' variation defined as var/nom. 
    0/0 division is taken to be 1.
    '''
    downup_axis = common.down_up_axis
    down = hh.divideHists(nom,var)
    up = hh.divideHists(var,nom)
    data = np.stack([down,up],axis=-1)
    new_histo = hist.Hist(*nom.axes,downup_axis, name=var.name, data=data, storage = hist.storage.Weight())
    return new_histo


def setPreconditionVec():
    f=h5py.File('../Fit/FitRes/fit_Wlike_asimov.hdf5.hdf5', 'r')
    hessian = f['hess'][:]
    eig, U = np.linalg.eigh(hessian)
    M1 = np.matmul(np.diag(1./np.sqrt(eig)),U.T)
    # print(M1,np.linalg.inv(np.linalg.inv(M1)))
    preconditioner = M1
    return preconditioner

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~loading boost histograms and cross sections from templates hdf5 file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
f = h5py.File("templatesTest_withSF.hdf5","r")
t = h5py.File('templatesFit.hdf5','r')
results = narf.ioutils.pickle_load_h5py(f["results"])
print(results['ZmumuPostVFP']["output"])
# Hdata_obs = results['dataPostVFP']["output"]["data_obs"].get()

#constants
process = 'ZmumuPostVFP'
V ='Z'
lumi    = results['dataPostVFP']["lumi"]
xsec    = results[process]["dataset"]["xsec"]
weights = results[process]["weight_sum"]
C       = lumi*1000*xsec/weights
Hdata_obs = C*results['ZmumuPostVFP']['output']['signal_nominal'].get()[sum,sum,:,:,:,sum]
# procs = ["lowacc"]
procs = []
systs_groups = {}

#first add nominal boost histogram for signal
H = C*results[process]['output']['signal_nominal'].get()
#Bin information
yBinsC     = H.axes[V+'rap'].centers
qtBinsC    = H.axes[V+'pt'].centers
charges    = H.axes['charge'].centers
eta        = H.axes['mueta'].centers
pt         = H.axes['mupt'].centers
helicities   = list(H.axes['helicities'])
unrolled_dim = len(eta) * len(pt)
qtBins = H.axes[V+'pt'].edges
yBins = H.axes[V+'rap'].edges

#Reshaping the data. 2d format. one row per unrolled pt-eta distribution
unrolled_and_stacked = np.swapaxes(H.to_numpy()[0].reshape(\
                           (len(yBinsC),len(qtBinsC),-1,len(charges),len(helicities)))\
                            ,2,-1).reshape(-1,unrolled_dim)

#TODO add bkg processes to this
sumw = np.concatenate((H[sum,sum,:,:,0,sum].values().ravel(),H[sum,sum,:,:,1,sum].values().ravel()))
sumw2 = np.concatenate((H[sum,sum,:,:,0,sum].variances().ravel(),H[sum,sum,:,:,1,sum].variances().ravel()))

#clean memory
H = None
#Generating multi index 
iterables = [yBinsC, qtBinsC,helicities ,charges]
multi = pd.MultiIndex.from_product(iterables , names = ['rapidity', 'qt' , 'hel','charge'])

'''Building the nominal DataFrame'''
    
#building dataframe
df = pd.Series(list(unrolled_and_stacked),index=multi, name="nominal")
df = pd.DataFrame(df)
print('\nnominal dataframe\n' , df.head())

unrolled_and_stacked = None
#Adding cross section information to our dataframe by creating cross section dataframe and merging
#TODO: pass boost histograms format

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1
T = t['helicity'][:threshold_y,:threshold_qt,:] #cross sections
    
iterables = [yBinsC , qtBinsC , helicities]
multi = pd.MultiIndex.from_product(iterables = [yBinsC , qtBinsC , helicities]
 , names = ['rapidity', 'qt' , 'hel']) #multi index to be joined on
s = pd.Series(T.ravel(), index = multi , name='xsec') #series carrying cross section information

xsec_df = pd.concat([s,s] ,axis=0).reset_index()       #same cross section for both charges, will need double to match dimensions
charges_xsec =  [-1.0]*288 + [1.0]*288
xsec_df['charge'] = charges_xsec

#now the dataframe carries cross section column
df = df.merge(xsec_df ,left_on=['rapidity','qt','hel','charge'], right_on=['rapidity','qt','hel','charge'])
print('\nadded cross-sections\n',df.head())

#setting process as index & cleaning up by removing redundant information
df.set_index(['helXsec_'+df['hel']+'_y_'+df['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+df['qt'].apply(lambda x: round(x,1)).apply(str),df['charge']],inplace=True)
df.drop(columns=['rapidity','qt','charge','hel'],inplace=True)
df.rename_axis(['process','charge'] ,inplace=True)
print('\nre-indexing\n',df.head())

#adding column for helicity group
# df['helgroups'] = df.index.get_level_values(0).map(lambda x: re.search("y.+" , x).group(0))
df['isSignal']  = True

print('\nnominal dataframe\n' , df.head())

#now add other procs
for proc in procs:
    histo=C*results[process]['output']['{}_nominal'.format(proc)].get()
    unrolled = histo.to_numpy()[0].reshape(len(charges),-1)
    histo=None
    #add data
    iterables_proc = [[proc],charges]
    multi_proc = pd.MultiIndex.from_product(iterables_proc , names = ['process','charge'])
    print(charges,iterables_proc,multi_proc)
    df_proc = pd.Series(list(unrolled),index=multi_proc, name="nominal")
    df_proc = pd.DataFrame(df_proc)
    unrolled=None
    df_proc['isSignal'] = False
    df_proc['xsec'] = -1
    df = pd.concat([df,df_proc],axis=0)
    
print('\nreorganizing and adding other procs\n',df.head(),df.tail())

#add systematics

systs_macrogroups = {} # this is a list over groups of systematics
systs_macrogroups['mass']=['mass_var']
# systs_macrogroups['sf']=['effStatTnP_sf_reco','effStatTnP_sf_tracking','effStatTnP_sf_idip','effStatTnP_sf_trigger'] #these correspond to the names of histograms to recall from file

procs = ["signal"]+procs #careful!! this must be the same order as before!
nominal_cols = ['Zrap', 'Zpt', 'mueta', 'mupt', 'charge', 'helicities','downUpVar']
multi = df.index

#loop over systematics:
for proc in procs:
    #get variations
    for syst,nuisances in systs_macrogroups.items():
        #TODO add exception for histograms not found
        syst_dfs = []
        print(syst,nuisances)
        for nuisance in nuisances:
            print('{proc}_{nuisance}'.format(proc=proc,nuisance=nuisance))
            syst_histo = C * results[process]['output']['{proc}_{nuisance}'.format(proc=proc,nuisance=nuisance)].get()
            axes = [axis for axis in syst_histo.axes]
            #select slices in systematics based on "vars"
            syst_arr = syst_histo.to_numpy()[0]
            syst_axes = [axis for axis in syst_histo.axes if axis.name not in nominal_cols]
            nsysts = 1 #this is the total number of systematics after considering all the bins
            for axis in syst_axes:
                nsysts = nsysts*len(axis.centers)
            print("nsysts",nsysts)
            syst_histo = None
            if proc == "signal":
                #Reshaping the data. 2d format. one row per unrolled pt-eta distribution
                syst_arr = syst_arr.reshape((len(yBinsC),len(qtBinsC),-1,len(charges),len(helicities),nsysts,2))#last axis is always up/down
                # rapidity, qt, data, charge, hel, syst, up/down
                syst_arr = np.moveaxis(syst_arr,2,-2)
                # rapidity, qt, charge, hel, syst, data, up/down
                syst_arr = np.swapaxes(syst_arr,2,3).reshape(-1,unrolled_dim*2)
                # rapidity, qt, hel, charge, syst, data, up/down
                names = ['rapidity', 'qt' , 'hel','charge']+[axis.name for axis in syst_axes]
                iterables = [yBinsC, qtBinsC,helicities ,charges] + [axis.centers for axis in syst_axes]
                multi = pd.MultiIndex.from_product(iterables, names = names)
                # print(multi)
                syst_df = pd.Series(list(syst_arr),index=multi,name=syst)
                print(syst_df.head())
                syst_df = pd.DataFrame(syst_df).reset_index()
                idx_strings = ['helXsec_'+syst_df['hel']+'_y_'+syst_df['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+syst_df['qt'].apply(lambda x: round(x,1)).apply(str),syst_df['charge']]
                syst_string = f"{nuisance}_"
                for axis in syst_axes:
                    syst_string+=axis.name.replace(' ','')+'_'+syst_df[axis.name].apply(str)+'_'
                syst_string = syst_string.apply(lambda s: s[:-1] if s.endswith('_') else s)
                idx_strings.append(syst_string)
                syst_df.set_index(idx_strings,inplace=True)
                syst_df.drop(columns=[axis.name for axis in syst_axes],inplace=True)
                syst_df.drop(columns=['rapidity','qt','charge','hel'],inplace=True)
                syst_df.rename_axis(['process','charge','syst'] ,inplace=True)
                syst_dfs.append(syst_df)
            else:
                pass
                #data, charge, syst, up/down
                syst_arr = np.moveaxis(syst_arr,0,-2)
                syst_arr = syst_arr.reshape(-1,unrolled_dim)
        # now merge all the dataframes within the group
        syst_df_merged = pd.concat(syst_dfs, axis=0)
        # get list of systematics and drop index
        syst_list = list(syst_df_merged.query("process == 'helXsec_L_y_0.2_qt_1.5' & charge==1.0").index.get_level_values(2))
        # pdb.set_trace()
        systs_groups[syst]=syst_list
        syst_df_merged= syst_df_merged.droplevel('syst')
        # group by process and charge, and concatenate the arrays
        syst_df_merged = syst_df_merged.groupby(["process", "charge"])[syst].agg(np.concatenate)
        syst_df_merged = syst_df_merged.map(lambda x: x.reshape((-1,unrolled_dim,2)))
        print(syst_df_merged.loc[('helXsec_L_y_0.2_qt_1.5',-1)].shape)
        print(syst_df_merged.head())
        df[syst]=syst_df_merged

print('\nadding systematics\n',df.head(),df.tail())

for syst in systs_macrogroups:
    #now divide by nominal
    df["{}_logk".format(syst)]=df.apply(lambda x: x[syst]/np.expand_dims(x['nominal'],axis=(0,-1)),axis='columns')
    print('\ndivide by nominal\n',df.head(),df.tail())
    #take log
    df["{}_logk".format(syst)]=df["{}_logk".format(syst)].map(lambda x: np.log(x))

    #remove spurious nans
    logkepsilon = math.log(1e-3)
    print('\n before regularization\n',df.head(),df.tail())
    df["{}_logk".format(syst)]=df.apply(lambda x: np.where(np.equal(np.sign(x[syst]*np.expand_dims(x['nominal'],axis=(0,-1))),1),x["{}_logk".format(syst)],logkepsilon*np.ones_like(x[syst])),axis='columns')

    #multiply down times -1
    mask = np.stack([-1*np.ones_like(df[syst].loc[('helXsec_L_y_0.2_qt_1.5',-1)][...,0]),np.ones_like(df[syst].loc[('helXsec_L_y_0.2_qt_1.5',-1)][...,0])],axis=-1)
    print(mask.shape)
    df["{}_logk".format(syst)]=df["{}_logk".format(syst)].apply(lambda x: mask*x)
    
    print('\nfinal df\n',df.head(),df.tail())




'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

#retrieve metadata

procs = list(df.query("charge==1.0").index.get_level_values(0))
signals = list(df.query("charge==1.0 & isSignal==True").index.get_level_values(0))
nproc = len(procs)
nsignals = len(signals)
maskedchans = ['Wlike_minus','Wlike_plus']

#list of groups of signal processes by charge - DON'T NEED THAT
chargegroups = []
chargegroupidxs = []

#list of groups of signal processes by polarization - DON'T NEED THAT
polgroups = []
polgroupidxs = []

#list of groups of signal processes by helicity xsec
helgroups = []
helgroupidxs = []
helGroups = fillHelGroup(yBinsC,qtBinsC,helicities)
for group in helGroups:
  helgroups.append(group)
  helgroupidx = []
  for proc in helGroups[group]:
    helgroupidx.append(procs.index(proc))
  helgroupidxs.append(helgroupidx)

#list of groups of signal processes to be summed
sumgroups = []
sumgroupsegmentids = []
sumgroupidxs = []
sumGroups = fillSumGroup(yBinsC,qtBinsC,helicities,signals)
for igroup,group in enumerate(sumGroups):
  sumgroups.append(group)
  for proc in sumGroups[group]:
    sumgroupsegmentids.append(igroup)
    sumgroupidxs.append(procs.index(proc))
    
#list of groups of signal processes by chargemeta - DON'T NEED THAT
chargemetagroups = []
chargemetagroupidxs = []
  
#list of groups of signal processes by ratiometa - DON'T NEED THAT
ratiometagroups = []
ratiometagroupidxs = []

#list of groups of signal processes by helmeta
helmetagroups = []
helmetagroupidxs = []
helMetaGroups = fillHelMetaGroup(yBinsC,qtBinsC,sumGroups)
for group in helMetaGroups:
  helmetagroups.append(group)
  helmetagroupidx = []
  for proc in helMetaGroups[group]:
    helmetagroupidx.append(sumgroups.index(proc))
  helmetagroupidxs.append(helmetagroupidx)

#list of groups of signal processes for regularization - DON'T NEED THAT
reggroups = []
reggroupidxs = []
  
poly1dreggroups = []
poly1dreggroupfirstorder = []
poly1dreggrouplastorder = []
poly1dreggroupnames = []
poly1dreggroupbincenters = []

poly2dreggroups = []
poly2dreggroupfirstorder = []
poly2dreggrouplastorder = []
poly2dreggroupfullorder = []
poly2dreggroupnames = []
poly2dreggroupbincenters0 = []
poly2dreggroupbincenters1 = []

#list of systematic uncertainties (nuisances)
systsd = OrderedDict()
systs = []
for group, nuisances in systs_groups.items():
    systs.extend(nuisances)
systsnoprofile = []
systsnoconstraint = ['mass_var_mass_var_0.5']

# for syst in systs:
#     if not 'NoProfile' in syst[2]:
#       systsd[syst[0]] = syst
#       systs.append(syst[0])
# for syst in systs:
#     if 'NoProfile' in syst[2]:
#       systsd[syst[0]] = syst
#       systs.append(syst[0])
#       systsnoprofile.append(syst[0])
#     if 'NoConstraint' in syst[2]:
#         systsnoconstraint.append(syst[0])

nsyst = len(systs)
  
#list of groups of systematics (nuisances) and lists of indexes
systgroups = []
systgroupidxs = []
# pdb.set_trace()
for group in systs_groups:
    systgroups.append(group)
    systgroupidx = []
    for syst in systs_groups[group]:
      systgroupidx.append(systs.index(syst))
    systgroupidxs.append(systgroupidx)

#list of groups of systematics to be treated as additional outputs for impacts, etc (aka "nuisances of interest")
noiGroups = {'mass':['mass_var_mass_var_0.5']}
noigroups = []
noigroupidxs = []
for group in noiGroups:
  noigroups.append(group)
  for syst in noiGroups[group]:
    noigroupidxs.append(systs.index(syst))


#write results to hdf5 file

dtype = 'float64'
procSize = nproc*np.dtype(dtype).itemsize
systSize = 2*nsyst*np.dtype(dtype).itemsize
defChunkSize = 4*1024**2
chunkSize = np.amax([defChunkSize,procSize,systSize])

constraintweights = np.ones([nsyst],dtype=dtype)
for syst in systsnoconstraint:
    constraintweights[systs.index(syst)] = 0.

if chunkSize > defChunkSize:
  print("Warning: Maximum chunk size in bytes was increased from %d to %d to align with tensor sizes and allow more efficient reading/writing." % (defChunkSize, chunkSize))

#create HDF5 file (chunk cache set to the chunk size since we can guarantee fully aligned writes
outfilename = "Wlike.hdf5"
f = h5py.File(outfilename, rdcc_nbytes=chunkSize, mode='w')

#save some lists of strings to the file for later use
hprocs = f.create_dataset("hprocs", [len(procs)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hprocs[...] = procs

hsignals = f.create_dataset("hsignals", [len(signals)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hsignals[...] = signals

hsysts = f.create_dataset("hsysts", [len(systs)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hsysts[...] = systs

hsystsnoprofile = f.create_dataset("hsystsnoprofile", [len(systsnoprofile)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hsystsnoprofile[...] = systsnoprofile

hsystsnoconstraint = f.create_dataset("hsystsnoconstraint", [len(systsnoconstraint)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hsystsnoconstraint[...] = systsnoconstraint

hsystgroups = f.create_dataset("hsystgroups", [len(systgroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hsystgroups[...] = systgroups

hsystgroupidxs = f.create_dataset("hsystgroupidxs", [len(systgroupidxs)], dtype=h5py.special_dtype(vlen=np.dtype('int32')), compression="gzip")
hsystgroupidxs[...] = systgroupidxs

hchargegroups = f.create_dataset("hchargegroups", [len(chargegroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hchargegroups[...] = chargegroups

hchargegroupidxs = f.create_dataset("hchargegroupidxs", [len(chargegroups),2], dtype='int32', compression="gzip")
hchargegroupidxs[...] = chargegroupidxs

hpolgroups = f.create_dataset("hpolgroups", [len(polgroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hpolgroups[...] = polgroups

hpolgroupidxs = f.create_dataset("hpolgroupidxs", [len(polgroups),3], dtype='int32', compression="gzip")
hpolgroupidxs[...] = polgroupidxs

hhelgroups = f.create_dataset("hhelgroups", [len(helgroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hhelgroups[...] = helgroups

hhelgroupidxs = f.create_dataset("hhelgroupidxs", [len(helgroups),6], dtype='int32', compression="gzip")
hhelgroupidxs[...] = helgroupidxs

hsumgroups = f.create_dataset("hsumgroups", [len(sumgroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hsumgroups[...] = sumgroups

hsumgroupsegmentids = f.create_dataset("hsumgroupsegmentids", [len(sumgroupsegmentids)], dtype='int32', compression="gzip")
hsumgroupsegmentids[...] = sumgroupsegmentids

hsumgroupidxs = f.create_dataset("hsumgroupidxs", [len(sumgroupidxs)], dtype='int32', compression="gzip")
hsumgroupidxs[...] = sumgroupidxs

hchargemetagroups = f.create_dataset("hchargemetagroups", [len(chargemetagroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hchargemetagroups[...] = chargemetagroups

hchargemetagroupidxs = f.create_dataset("hchargemetagroupidxs", [len(chargemetagroups),2], dtype='int32', compression="gzip")
hchargemetagroupidxs[...] = chargemetagroupidxs

hratiometagroups = f.create_dataset("hratiometagroups", [len(ratiometagroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hratiometagroups[...] = ratiometagroups

hratiometagroupidxs = f.create_dataset("hratiometagroupidxs", [len(ratiometagroups),2], dtype='int32', compression="gzip")
hratiometagroupidxs[...] = ratiometagroupidxs

hhelmetagroups = f.create_dataset("hhelmetagroups", [len(helmetagroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hhelmetagroups[...] = helmetagroups

hhelmetagroupidxs = f.create_dataset("hhelmetagroupidxs", [len(helmetagroups),6], dtype='int32', compression="gzip")
hhelmetagroupidxs[...] = helmetagroupidxs

hreggroups = f.create_dataset("hreggroups", [len(reggroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hreggroups[...] = reggroups

hreggroupidxs = f.create_dataset("hreggroupidxs", [len(reggroupidxs)], dtype=h5py.special_dtype(vlen=np.dtype('int32')), compression="gzip")
hreggroupidxs[...] = reggroupidxs

hpoly1dreggroups = f.create_dataset("hpoly1dreggroups", [len(poly1dreggroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hpoly1dreggroups[...] = poly1dreggroups

hpoly1dreggroupfirstorder = f.create_dataset("hpoly1dreggroupfirstorder", [len(poly1dreggroupfirstorder)], dtype='int32', compression="gzip")
hpoly1dreggroupfirstorder[...] = poly1dreggroupfirstorder

hpoly1dreggrouplastorder = f.create_dataset("hpoly1dreggrouplastorder", [len(poly1dreggrouplastorder)], dtype='int32', compression="gzip")
hpoly1dreggrouplastorder[...] = poly1dreggrouplastorder

hpoly1dreggroupnames = f.create_dataset("hpoly1dreggroupnames", [len(poly1dreggroupnames)], dtype=h5py.special_dtype(vlen="S256"), compression="gzip")
hpoly1dreggroupnames[...] = poly1dreggroupnames

hpoly1dreggroupbincenters = f.create_dataset("hpoly1dreggroupbincenters", [len(poly1dreggroupbincenters)], dtype=h5py.special_dtype(vlen=np.dtype('float64')), compression="gzip")
hpoly1dreggroupbincenters[...] = poly1dreggroupbincenters

hpoly2dreggroups = f.create_dataset("hpoly2dreggroups", [len(poly2dreggroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hpoly2dreggroups[...] = poly2dreggroups

hpoly2dreggroupfirstorder = f.create_dataset("hpoly2dreggroupfirstorder", [len(poly2dreggroupfirstorder),2], dtype='int32', compression="gzip")
hpoly2dreggroupfirstorder[...] = poly2dreggroupfirstorder

hpoly2dreggrouplastorder = f.create_dataset("hpoly2dreggrouplastorder", [len(poly2dreggrouplastorder),2], dtype='int32', compression="gzip")
hpoly2dreggrouplastorder[...] = poly2dreggrouplastorder

hpoly2dreggroupfullorder = f.create_dataset("hpoly2dreggroupfullorder", [len(poly2dreggroupfullorder),2], dtype='int32', compression="gzip")
hpoly2dreggroupfullorder[...] = poly2dreggroupfullorder

hpoly2dreggroupnames = f.create_dataset("hpoly2dreggroupnames", [len(poly2dreggroupnames)], dtype=h5py.special_dtype(vlen="S256"), compression="gzip")
hpoly2dreggroupnames[...] = poly2dreggroupnames

hpoly2dreggroupbincenters0 = f.create_dataset("hpoly2dreggroupbincenters0", [len(poly2dreggroupbincenters0)], dtype=h5py.special_dtype(vlen=np.dtype('float64')), compression="gzip")
hpoly2dreggroupbincenters0[...] = poly2dreggroupbincenters0

hpoly2dreggroupbincenters1 = f.create_dataset("hpoly2dreggroupbincenters1", [len(poly2dreggroupbincenters1)], dtype=h5py.special_dtype(vlen=np.dtype('float64')), compression="gzip")
hpoly2dreggroupbincenters1[...] = poly2dreggroupbincenters1

#Saving Preconditioner
preconditioner = setPreconditionVec()
hpreconditioner = f.create_dataset("hpreconditioner", preconditioner.shape, dtype='float64', compression="gzip")
hpreconditioner[...] = preconditioner


invpreconditioner = np.linalg.inv(preconditioner)
hinvpreconditioner = f.create_dataset("hinvpreconditioner", invpreconditioner.shape, dtype='float64', compression="gzip")
hinvpreconditioner[...] = invpreconditioner

hnoigroups = f.create_dataset("hnoigroups", [len(noigroups)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hnoigroups[...] = noigroups

hnoigroupidxs = f.create_dataset("hnoigroupidxs", [len(noigroupidxs)], dtype='int32', compression="gzip")
hnoigroupidxs[...] = noigroupidxs

hmaskedchans = f.create_dataset("hmaskedchans", [len(maskedchans)], dtype=h5py.special_dtype(vlen=str), compression="gzip")
hmaskedchans[...] = maskedchans

#create h5py datasets with optimized chunk shapes
nbytes = 0

nbytes += writeFlatInChunks(constraintweights, f, "hconstraintweights", maxChunkBytes = chunkSize)
constraintweights = None


# data_obs = np.concatenate((Hdata_obs.to_numpy()[0][...,0].ravel(),Hdata_obs.to_numpy()[0][...,1].ravel()))
data_obs = Hdata_obs.to_numpy()[0][...,0].ravel()
Hdata_obs = None
data_obs = np.random.poisson(lam=data_obs)
data_obs = np.array(data_obs,dtype='float64')

nbytes += writeFlatInChunks(data_obs, f, "hdata_obs", maxChunkBytes = chunkSize)
data_obs = None

#compute poisson parameter for Barlow-Beeston bin-by-bin statistical uncertainties
kstat = np.square(sumw)/sumw2
#numerical protection to avoid poorly defined constraint
kstat = np.where(np.equal(sumw,0.), 1., kstat)

nbytes += writeFlatInChunks(kstat, f, "hkstat", maxChunkBytes = chunkSize)
kstat = None

nbytes += writeFlatInChunks(sumw, f, "hsumw", maxChunkBytes = chunkSize)
sumw = None

nbytes += writeFlatInChunks(sumw2, f, "hsumw2", maxChunkBytes = chunkSize)
sumw2 = None

#n.b data and expected have shape [nbins]
#sumw and sumw2 keep track of total nominal statistical uncertainty per bin and have shape [nbins]

#norm has shape [nbinsfull, nproc] and keeps track of expected normalization

#logk has shape [nbinsfull, nproc, 2, nsyst] and keep track of systematic variations
#per nuisance-parameter, per-process, per-bin
#the second-last dimension, of size 2, indexes "logkavg" and "logkhalfdiff" for asymmetric uncertainties
#where logkavg = 0.5*(logkup + logkdown) and logkhalfdiff = 0.5*(logkup - logkdown)

#n.b, in case of masked channels, nbinsfull includes the masked channels where nbins does not


# retrieve norm
# norm = np.concatenate((np.stack(df.query("charge==-1.")['nominal'].values,axis=-1),np.stack(df.query("charge==1.")['nominal'].values,axis=-1),np.expand_dims(np.stack(df.query("charge==-1.")['xsec'].values,axis=-1),axis=0),np.expand_dims(np.stack(df.query("charge==1.")['xsec'].values,axis=-1),axis=0)),axis=0)
norm = np.concatenate((np.stack(df.query("charge==-1.")['nominal'].values,axis=-1),np.expand_dims(np.stack(df.query("charge==-1.")['xsec'].values,axis=-1),axis=0)),axis=0)
nbytes += writeFlatInChunks(norm, f, "hnorm", maxChunkBytes = chunkSize)

# df = df.drop(columns=['nominal'],inplace=True) #why this doesn't work?
print('\ndrop nominal\n',df.head())
logk_systs = []
for syst in systs_groups:
    print(df.query("charge==-1.")["{}_logk".format(syst)].values[0].shape)
    print(df.query("charge==-1.")['nominal'].values[0].shape)
    # logk_syst = np.moveaxis(np.concatenate((np.stack(df.query("charge==-1.")["{}_logk".format(syst)].values,axis=-2),np.stack(df.query("charge==1.")["{}_logk".format(syst)].values,axis=-2)),axis=1),0,-1)
    logk_syst = np.moveaxis(np.stack(df.query("charge==-1.")["{}_logk".format(syst)].values,axis=-2),0,-1)
    logk_systs.append(logk_syst)
    print(logk_syst.shape)
print(logk_systs[0].shape)
logk_systs = np.concatenate(logk_systs,axis=-1) #concatenate along syst axis
print('logk_systs.shape',logk_systs.shape)

# retrieve logk
logk_up = logk_systs[...,1,:]
logk_down = logk_systs[...,0,:]

print(logk_up.shape,logk_down.shape)

logkavg = 0.5*(logk_up + logk_down)
logkhalfdiff = 0.5*(logk_up - logk_down)

print(logkavg.shape)
#ensure that systematic tensor is sparse where normalization matrix is sparse
# logkavg = np.where(np.equal(np.expand_dims(norm[:-2,:],axis=-1),0.), np.zeros_like(logkavg), logkavg)
# logkhalfdiff = np.where(np.equal(np.expand_dims(norm[:-2,:],axis=-1),0.), np.zeros_like(logkavg), logkhalfdiff)
logkavg = np.where(np.equal(np.expand_dims(norm[:-1,:],axis=-1),0.), np.zeros_like(logkavg), logkavg)
logkhalfdiff = np.where(np.equal(np.expand_dims(norm[:-1,:],axis=-1),0.), np.zeros_like(logkavg), logkhalfdiff)

norm = None

logk = np.stack((logkavg,logkhalfdiff),axis=-2)
print(logk.shape)
# logk = np.concatenate((logk,np.zeros((2,nproc,2,nsyst))),axis=0)
logk = np.concatenate((logk,np.zeros((1,nproc,2,nsyst))),axis=0)
print(logk.shape)
nbytes += writeFlatInChunks(logk, f, "hlogk", maxChunkBytes = chunkSize)
logk = None

print("Total raw bytes in arrays = %d" % nbytes)

