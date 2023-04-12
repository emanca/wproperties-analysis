

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
            s = 'y_{i}_qt_{j}'.format(i=i,j=j)
            helGroups[s] = []
            
            for hel in helXsecs:
                if not '7' in hel and not '8' in hel and not '9' in hel:
                    helGroups[s].append('helXsecs'+hel+'_'+s)
            if helGroups[s] == []:
                del helGroups[s]
    return helGroups

def fillHelMetaGroup(yBinsC,qtBinsC,sumGroups):
    helMetaGroups = OrderedDict()
    for i in range(len(yBinsC)):
        s = 'y_{i}'.format(i=i)
        helMetaGroups[s] = []
        for key in sumGroups:
            if s in key:
                helMetaGroups[s].append(key)
        
        if helMetaGroups[s] == []:
                del helMetaGroups[s]
    
    for j in range(len(qtBinsC)):
        s = 'qt_{j}'.format(j=j)
        helMetaGroups[s] = []
        for key in sumGroups:
            if 'qt' in key and key.split('_')[2]==str(j):
                helMetaGroups[s].append(key)
    
        if helMetaGroups[s] == []:
                del helMetaGroups[s]
    return helMetaGroups

def fillSumGroup(yBinsC,qtBinsC,helXsecs,processes):
    sumGroups = OrderedDict()
    for i in range(len(yBinsC)):
        s = 'y_{i}'.format(i=i)
        for hel in helXsecs:
            for signal in processes:
                sumGroups['helXsecs'+hel+'_'+s] = []
                for j in range(len(qtBinsC)):
                    #if 'helXsecs'+hel+'_'+'y_{i}_qt_{j}'.format(i=i,j=j) in processes:
                    sumGroups['helXsecs'+hel+'_'+s].append('helXsecs'+hel+'_'+s+'_qt_{j}'.format(j=j))
    
    for j in range(len(qtBinsC)):
        s = 'qt_{j}'.format(j=j)
        for hel in helXsecs:
            for signal in processes:
                if signal.split('_')[0] == 'helXsecs'+hel and signal.split('_')[4] == str(j):
                    sumGroups['helXsecs'+hel+'_'+s] = []
                    for i in range(len(yBinsC)):
                        #print i, signal, 'helXsecs'+hel+'_'+'y_{i}_pt_{j}'.format(i=i,j=j)
                        #print 'append', 'helXsecs'+hel+'_y_{i}_'.format(i=i)+s, 'to', 'helXsecs'+hel+'_'+s
                        sumGroups['helXsecs'+hel+'_'+s].append('helXsecs'+hel+'_y_{i}_'.format(i=i)+s)
    return sumGroups

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~loading boost histograms and cross sections from templates hdf5 file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
f = h5py.File("templatesTest2.hdf5","r")
results = narf.ioutils.pickle_load_h5py(f["results"])

data_obs = results['dataPostVFP']["data_obs"].get()
# lumi = results['dataPostVFP']["lumi"]
lumi = 16.8
xsec = results['ZmumuPostVFP']["dataset"]["xsec"]

H = results['ZmumuPostVFP']['output']['signalTemplates_nominal'].get() #boost histogram values
#rescale to match the data luminosity
H = H*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]

t = h5py.File('templatesFit.hdf5','r') # templates file used for getting the cross sections (later will also be boost histogram)
H3 = results['ZmumuPostVFP']['output']['signalTemplates_mass'].get() #boost histogram with systematic mass variations
H3 = H3*lumi*1000*xsec/results["ZmumuPostVFP"]["weight_sum"]

M = H3[:,:,:,:,:,:,[bh.loc('massShift100MeVDown') , bh.loc('massShift100MeVUp')]] #boost histogram for selected mass variations
low_acc      = results['ZmumuPostVFP']['output']['lowacc'].get().to_numpy()[0].reshape(-1,2)
low_acc_mass = results['ZmumuPostVFP']['output']['lowacc_mass'].get()
low_acc_mass_array = low_acc_mass[:,:,:,[bh.loc('massShift100MeVDown') , bh.loc('massShift100MeVUp')]].to_numpy()[0].reshape(-1,2,2)
'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Building the nominal dataframe~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
'''Unpacking the data'''
#first unroll the tensor in eta and pt shape: (6,8,48,60,2,6) -> (6,8,2880,2,6)
unrolled = H.to_numpy()[0].reshape((6,8,-1,2,6))
#next, swap axes such that unrolled eta/pt in in last position (6,8,2880,2,6) -> (6,8,2,6,2880)
#lastly, reshape into 2 dimensional array which will be passed into dataframe (6,8,2,6,2880) -> (576,2880)
#one row corresponds to one unrolled pt/eta distribution (template)
a = np.swapaxes(unrolled , 2 , -1).reshape(-1,2880) 

'''Building the pandas dataframe'''
yBinsC     = H.axes['Zrap'].centers
qtBinsC    = H.axes['Zpt'].centers
charges    = H.axes['charge'].centers
helicities = list(H.axes['helicities'])

#making multi index object
iterables = [yBinsC, qtBinsC,helicities ,charges] #2charges * 6helicities *6y bins * 8qt bins =  576 rows
multi = pd.MultiIndex.from_product(iterables , names = ['rapidity', 'qt' , 'hel','charge'])

#building dataframe
df = pd.DataFrame(a , index = multi)


'''Adding cross section information to our dataframe by creating cross section dataframe and merging'''
qtBins = np.array([0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.,200.]) #these have to be like this for now
yBins = np.array([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0])

threshold_y = np.digitize(2.4,yBins)-1
threshold_qt = np.digitize(60.,qtBins)-1

T = t['helicity'][:threshold_y,:threshold_qt,:] #cross sections
processes = [yBinsC , qtBinsC , helicities]
multi2 = pd.MultiIndex.from_product(processes , names = ['rapidity', 'qt' , 'hel']) #multi index to be joined on
s = pd.Series(T.ravel(), index = multi2 , name='xsec') #series carrying cross section information

xsec_df = pd.concat([s,s] ,axis=0).reset_index()       #same cross section for both charges, will need double to match dimensions
charges =  [-1.0]*288 + [1.0]*288
xsec_df['charge'] = charges


#now the dataframe carries cross section column
df = df.merge(xsec_df ,left_on=['rapidity','qt','hel','charge'], right_on=['rapidity','qt','hel','charge'])


#setting process as index & cleaning up by removing redundant information
df.set_index(['helXsec_'+df['hel']+'_y_'+df['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+df['qt'].apply(lambda x: round(x,1)).apply(str),df['charge']],inplace=True)
df.drop(columns=['rapidity','qt','hel','charge'],inplace=True)
df.rename_axis(['process','charge'] ,inplace=True)

#reorganizing data into a single column labeled 'data' - contains unrolled pt/eta distribution
df['data'] = df.loc[:,0:2879].apply(np.hstack , axis=1) 
df.drop(columns = df.loc[:,0:2879].columns , inplace = True)

#adding column for helicity group
df['helgroups'] = df.index.get_level_values(0).map(lambda x: re.search("y.+" , x).group(0))
df['isSignal']  = True
print('Nominal Dataframe:')
print(df.head())

#low acceptance
df.loc[('low_acc',-1.0),:] = [-1 , low_acc[:,0] ,np.nan, False]
df.loc[('low_acc', 1.0),:] = [-1 , low_acc[:,1] ,np.nan, False]

print(df.loc['low_acc'])

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Building the dataframe for the mass variation systematics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

#reshaping data 
mass_arr = M.to_numpy()[0]
mass_unrolled = mass_arr.reshape((6,8,-1,2,6,2))                     #shape: (6,8,2880,2,6,2)
mass_swapped = np.swapaxes(mass_unrolled , 2,-1).reshape(-1,2880)    #rapidity, qt, mass variation, charge,helicity 

#Building the pandas dataframe indexed by these values
m_yBinsC               = M.axes['Zrap'].centers
m_qtBinsC              = M.axes['Zpt'].centers
m_charges              = M.axes['charge'].centers
m_helicities      = list(M.axes['helicities'])
mass_variations   = list(M.axes['massShift']) #only teo in this case

#multi index object
m_iterables = [m_yBinsC, m_qtBinsC,mass_variations,m_charges ,m_helicities] #order must match mass_swapped array shape
m_multi = pd.MultiIndex.from_product(m_iterables , names = ['rapidity', 'qt' ,'syst' ,'charge','hel'])

#building dataframe and setting index values as process strings
m_df = pd.DataFrame(mass_swapped , index = m_multi)
m_df.reset_index(inplace=True)
m_df.set_index(['helXsec_'+m_df['hel']+'_y_'+m_df['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+m_df['qt'].apply(lambda x: round(x,1)).apply(str),m_df['charge'], m_df['syst']],inplace=True)
m_df.drop(columns=['rapidity','qt','hel','charge','syst'],inplace=True)
m_df.rename_axis(['process','charge','variation'] ,inplace=True)


#reorganizing data into a single column labeled 'data' - contains unrolled pt/eta distribution
m_df['data'] = m_df.loc[:,0:2879].apply(np.hstack , axis=1) 
m_df.drop(columns = m_df.loc[:,0:2879].columns , inplace = True)

print('\nsystematics dataframe')
print(m_df.head())

#retrieve metadata

procs = []
signals = []
nproc = len(processes)
nsignals = len(signals)
maskedchans = []

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
sumGroups = fillSumGroup(yBinsC,qtBinsC,helicities,procs)
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

#list of groups of systematics to be treated as additional outputs for impacts, etc (aka "nuisances of interest")
noigroups = []
noigroupidxs = []
for group in noiGroups:
  noigroups.append(group)
  for syst in noiGroups[group]:
    noigroupidxs.append(systs.index(syst))


#write results to hdf5 file

procSize = nproc*np.dtype(dtype).itemsize
systSize = 2*nsyst*np.dtype(dtype).itemsize
chunkSize = np.amax([options.chunkSize,procSize,systSize])
if chunkSize > options.chunkSize:
  print("Warning: Maximum chunk size in bytes was increased from %d to %d to align with tensor sizes and allow more efficient reading/writing." % (options.chunkSize, chunkSize))

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

# hpreconditioner = f.create_dataset("hpreconditioner", preconditioner.shape, dtype='float64', compression="gzip")
# hpreconditioner[...] = preconditioner

# hinvpreconditioner = f.create_dataset("hinvpreconditioner", invpreconditioner.shape, dtype='float64', compression="gzip")
# hinvpreconditioner[...] = invpreconditioner

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

nbytes += writeFlatInChunks(data_obs, f, "hdata_obs", maxChunkBytes = chunkSize)
data_obs = None

nbytes += writeFlatInChunks(kstat, f, "hkstat", maxChunkBytes = chunkSize)
kstat = None

nbytes += writeFlatInChunks(sumw, f, "hsumw", maxChunkBytes = chunkSize)
sumw = None

nbytes += writeFlatInChunks(sumw2, f, "hsumw2", maxChunkBytes = chunkSize)
sumw2 = None

nbytes += writeFlatInChunks(norm, f, "hnorm", maxChunkBytes = chunkSize)
norm = None

nbytes += writeFlatInChunks(logk, f, "hlogk", maxChunkBytes = chunkSize)
logk = None

print("Total raw bytes in arrays = %d" % nbytes)
print(m_df.head())

m_df.loc[('low_acc',-1.0,'massShift100MeVDown'),:] = [low_acc_mass_array[:,0,0]]
m_df.loc[('low_acc', 1.0,'massShift100MeVDown'),:] = [low_acc_mass_array[:,1,0]]
m_df.loc[('low_acc',-1.0,'massShift100MeVUp'),:]   = [low_acc_mass_array[:,0,1]]
m_df.loc[('low_acc', 1.0,'massShift100MeVUp'),:]   = [low_acc_mass_array[:,1,1]]

print(m_df.loc['low_acc'])
