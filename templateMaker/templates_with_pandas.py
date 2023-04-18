

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
    print(helMetaGroups)
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

def make_nominal_df(nominal_histos):
    
    H = nominal_histos[0] #TODO: pick by name

    #Bin information
    yBinsC     = H.axes[V+'rap'].centers
    qtBinsC    = H.axes[V+'pt'].centers
    charges    = H.axes['charge'].centers
    eta        = H.axes['mueta'].centers
    pt         = H.axes['mupt'].centers
    helicities   = list(H.axes['helicities'])
    unrolled_dim = len(eta) * len(pt)
    
    #Reshaping the data. 2d format. one row per unrolled pt-eta distribution
    unrolled_and_stacked = np.swapaxes(H.to_numpy()[0].reshape(\
                           (len(yBinsC),len(qtBinsC),-1,len(charges),len(helicities)))\
                            ,2,-1).reshape(-1,unrolled_dim)
    
    #Generating multi index 
    iterables = [yBinsC, qtBinsC,helicities ,charges]
    multi = pd.MultiIndex.from_product(iterables , names = ['rapidity', 'qt' , 'hel','charge'])
    
    '''Building the nominal DataFrame'''
    
    #building dataframe
    df = pd.DataFrame(unrolled_and_stacked , index = multi)
    print('\nnominal dataframe\n' , df.head())

    #Adding cross section information to our dataframe by creating cross section dataframe and merging
    #TODO: pass boost histograms format
    qtBins = np.array([0., 3., 6., 9.62315204,12.36966732,16.01207711,21.35210602,29.50001253,60.,200.]) 
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
    print('\nadded cross-sections\n',df.head())
    #setting process as index & cleaning up by removing redundant information
    df.set_index(['helXsec_'+df['hel']+'_y_'+df['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+df['qt'].apply(lambda x: round(x,1)).apply(str),df['charge']],inplace=True)
    df.drop(columns=['rapidity','qt','charge'],inplace=True)
    df.rename_axis(['process','charge'] ,inplace=True)
    print('\nre-indexing\n',df.head())

    #reorganizing data into a single column labeled 'data' - contains unrolled pt/eta distribution
    df['data'] = df.loc[:,0:unrolled_dim-1].apply(np.hstack , axis=1) 
    df.drop(columns = df.loc[:,0:unrolled_dim-1].columns , inplace = True)

    #adding column for helicity group
    df['helgroups'] = df.index.get_level_values(0).map(lambda x: re.search("y.+" , x).group(0))
    df['isSignal']  = True
    
    #now add bkg processes
    for histo in nominal_histos:
      histo = histo.to_numpy()[0].reshape(-1, len(histo.axes['charge'].centers))
      df.loc[('low_acc',-1.0),:] = [np.nan, -1 , histo[:,0] ,np.nan, False]
      df.loc[('low_acc', 1.0),:] = [np.nan, -1 , histo[:,1] ,np.nan, False]
    
    print('\nreorganizing and adding low acceptance\n',df.head(),df.tail())
    
    return df

def addSystVariation(nom_df,syst_dict,results,C):
    '''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~building the systematics DataFrame~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
    
    # syst_dict ={
    # "mass" : {
    #     "vars":["massShift50MeVDown","massShift50MeVUp"],
    #     "procs": ["signalTemplates","lowacc"],
    #     "type": "shapeNoConstraint",
    #     "weight" : 1.
    #   },
    # }
    #loop over systematics:
    for syst,syst_tools in syst_dict.items():
        #get variations
        for proc in syst_tools['procs']:
          syst_histo = C * results['ZmumuPostVFP']['output']['{proc}_{syst}'.format(proc=proc,syst=syst)].get()
          #select slices in systematics based on "vars"
          
          syst_histo = syst_histo[...,[hist.loc(["massShift50MeVDown","massShift50MeVUp"])]] #find a way to make this generic
          print(syst_histo)
    return nom_df
    #building dataframe and setting index values as process strings
    
    m_df = pd.DataFrame(mass_unrolled_and_stacked , index = m_multi)
    print('\nsystematics dataframe\n',m_df.head())
    m_df.reset_index(inplace=True)
    print(m_df.head())
    m_df.set_index(['helXsec_'+m_df['hel']+'_y_'+m_df['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+m_df['qt'].apply(lambda x: round(x,1)).apply(str),m_df['charge'], m_df['syst']],inplace=True)
    print(m_df.head())
    m_df.drop(columns=['rapidity','qt','charge','syst'],inplace=True)
    print(m_df.head())
    m_df.rename_axis(['process','charge','variation_index'] ,inplace=True)
    print(m_df.head())


    #reorganizing data into a single column labeled 'data' - contains unrolled pt/eta distribution
    m_df['data'] = m_df.loc[:,0:unrolled_dim-1].apply(np.hstack , axis=1) 
    m_df.drop(columns = m_df.loc[:,0:unrolled_dim-1].columns , inplace = True)
    print(m_df.head())
    
    #Adding low acceptance for mass shift systematic
    low_acceptance_df = pd.DataFrame(low_acc_mass.swapaxes(0,-1).reshape(-1,unrolled_dim),index=\
           pd.MultiIndex.from_product([['low_acc'],H.axes['charge'].centers,list(M.axes['massShift'])]))
    low_acceptance_df['hel']  = np.nan
    

    m_df = pd.concat([m_df,low_acceptance_df])
    
    m_df['data'] = m_df.loc[:,0:unrolled_dim-1].apply(np.hstack , axis=1) 
    m_df.drop(columns = m_df.loc[:,0:unrolled_dim-1].columns , inplace = True)
    print(m_df.head())

    #logK
    m_df['syst']      = m_df.index.get_level_values(2).map(lambda x: x.split('V')[0] + 'V')
    m_df['variation'] = m_df.index.get_level_values(2).map(lambda x: x.split('V')[1])
    print(m_df.head())
    print(m_df.loc['helXsec_T_y_{}_qt_{}'.format(round(yBinsC[1],1),round(qtBinsC[5],1))].query("variation == 'Down'")['data'].values)
    data = m_df.loc['helXsec_T_y_{}_qt_{}'.format(round(yBinsC[1],1),round(qtBinsC[5],1))].query("variation == 'Down'")['data'].values[1]
    print(data[data!=0.0])

    m_df.reset_index(inplace=True)
    print(m_df.head())
    m_df.drop(columns = 'variation_index' , inplace = True)
    m_df.set_index(['process','charge'] , inplace =True)

    m_down = m_df.query("variation == 'Down'")['data']
    m_up   = m_df.query("variation == 'Up'")['data']
    print(m_down.head())
    nominal= df['data']

    logkepsilon = math.log(1e-3)

    #m_down = m_df.set_index(['process','charge']).query("variation == 'Down'")['data']
    #m_up   = m_df.set_index(['process','charge']).query("variation == 'Up'")['data']

    #nominal= df.sort_index(level=['process','charge'])['data']

    logk_up   = (m_up/nominal).apply(lambda x: np.log(x))
    logk_down = (m_down/nominal).apply(lambda x: -np.log(x))


    #true if product is positive
    truth_up   = (nominal*m_up).apply(lambda x: np.equal(np.sign(x),1))
    truth_down = (nominal*m_down).apply(lambda x: np.equal(np.sign(x),1))


    epsilon_up   =   logk_up.apply(lambda x: logkepsilon*np.ones_like(x))
    epsilon_down = logk_down.apply(lambda x: -logkepsilon*np.ones_like(x))


    logk_up   = pd.Series([np.where(truth_up.values[i]\
              ,(m_up/nominal).apply(lambda x: np.log(x)).values[i]\
              , epsilon_up.values[i])\
                 for i in range(len(logk_up))] , index = logk_up.index)
    logk_down = pd.Series([np.where(truth_down.values[i]\
              ,(m_down/nominal).apply(lambda x: -np.log(x)).values[i]\
              , epsilon_down.values[i])\
                 for i in range(len(logk_up))] , index = logk_down.index)

    logk_up   = pd.DataFrame({'logK':logk_up , 'variation':'Up'})
    logk_down = pd.DataFrame({'logK':logk_down , 'variation':'Down'})

    m_df = pd.concat([logk_up,logk_down])

    print('\nSystematics Dataframe')
    print(m_df.head())

    # print(m_df.loc['helXsec_T_y_{}_qt_{}'.format(round(yBinsC[1],1),round(qtBinsC[5],1))].query("variation == 'Down'")['logK'].values[1])
    # data2 = m_df.loc['helXsec_T_y_{}_qt_{}'.format(round(yBinsC[1],1),round(qtBinsC[5],1))].query("variation == 'Down'")['logK'].values[1]
    # # print(norm_chan[norm_chan>0])
    # # data2 = np.where(np.equal(np.sign(norm_chan*data),1), data2, -logkepsilon*np.ones_like(data2))
    # print(data2[data2<6.9])

    return df,m_df,yBinsC,qtBinsC,helicities

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~loading boost histograms and cross sections from templates hdf5 file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
f = h5py.File("templatesTest2.hdf5","r")
t = h5py.File('templatesFit.hdf5','r')
results = narf.ioutils.pickle_load_h5py(f["results"])
Hdata_obs = results['dataPostVFP']["output"]["data_obs"].get()


#constants
process = 'ZmumuPostVFP'
V ='Z'
two_mass_variations = True
good_idx_mass = [5,15]
lumi    = results['dataPostVFP']["lumi"]
xsec    = results[process]["dataset"]["xsec"]
weights = results[process]["weight_sum"]
C       = lumi*1000*xsec/weights


#nominal boost histogram
H = C * results[process]['output']['signalTemplates_nominal'].get()
#low acceptance -nominal
low_acc      = C * results[process]['output']['lowacc'].get()

nominal_histos = [H,low_acc]

df = make_nominal_df(nominal_histos)

syst_dict ={
    "mass" : {
        "vars":["massShift50MeVDown","massShift50MeVUp"],
        "procs": ["signalTemplates","lowacc"],
        "type": "shapeNoConstraint",
        "weight" : 1.
    },
}

df = addSystVariation(df,syst_dict,results,C)

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

#retrieve metadata

procs = list(df.query("charge==1.0").index.get_level_values(0))
procs2 = list(m_df.query("charge==1.0").index.get_level_values(0))
for p1,p2 in zip(procs,procs2):
    print(p1,p2)
signals = list(df.query("charge==1.0 & isSignal==True").index.get_level_values(0))
nproc = len(procs)
nsignals = len(signals)
maskedchans = ['Wlike_minus']

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

#list of systematic uncertainties (nuisances)
systsd = OrderedDict()
systs = ['mass']
systsnoprofile = []
systsnoconstraint = ['mass']

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
# for group in groups:
#     systgroups.append(group)
#     systgroupidx = []
#     for syst in groups[group]:
#       systgroupidx.append(systs.index(syst))
#     systgroupidxs.append(systgroupidx)

#list of groups of systematics to be treated as additional outputs for impacts, etc (aka "nuisances of interest")
noigroups = []
noigroupidxs = []
# for group in noiGroups:
#   noigroups.append(group)
#   for syst in noiGroups[group]:
#     noigroupidxs.append(systs.index(syst))


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


data_obs = np.concatenate((Hdata_obs.to_numpy()[0][...,0].ravel(),Hdata_obs.to_numpy()[0][...,1].ravel()))
# data_obs = Hdata_obs.to_numpy()[0][...,1].ravel()


nbytes += writeFlatInChunks(data_obs, f, "hdata_obs", maxChunkBytes = chunkSize)
data_obs = None

sumw = np.concatenate((H[sum,sum,:,:,0,sum].values().ravel(),H[sum,sum,:,:,1,sum].values().ravel()))
sumw2 = np.concatenate((H[sum,sum,:,:,0,sum].variances().ravel(),H[sum,sum,:,:,1,sum].variances().ravel()))

# sumw = H[sum,sum,:,:,0,sum].values().ravel()
# sumw2 = H[sum,sum,:,:,0,sum].variances().ravel()


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

norm = np.concatenate((np.stack(df.query("charge==-1.")['data'].values,axis=-1),np.stack(df.query("charge==1.")['data'].values,axis=-1),np.expand_dims(np.stack(df.query("charge==-1.")['xsec'].values,axis=-1),axis=0),np.expand_dims(np.stack(df.query("charge==1.")['xsec'].values,axis=-1),axis=0)),axis=0)
nbytes += writeFlatInChunks(norm, f, "hnorm", maxChunkBytes = chunkSize)
norm = None

# norm = np.concatenate((np.stack(df.query("charge==1.")['data'].values,axis=-1),np.expand_dims(np.stack(df.query("charge==1.")['xsec'].values,axis=-1),axis=0)),axis=0)
# nbytes += writeFlatInChunks(norm, f, "hnorm", maxChunkBytes = chunkSize)

# m_df = m_df[m_df.index.get_level_values(0).map(lambda x: x.find('UL') > 1 or x.find('L') > 1)]
print("logk shape",m_df['logK'].values[0].shape)
logk_ups = np.stack(m_df.query("charge==-1.0  & variation == 'Up'")['logK'].values,axis=-1)
logk_downs = np.stack(m_df.query("charge==-1.0  & variation =='Down'")['logK'].values,axis=-1)

logkavg = 0.5*(logk_ups + logk_downs)
logkhalfdiff = 0.5*(logk_ups - logk_downs)

# #ensure that systematic tensor is sparse where normalization matrix is sparse
# logkavg = np.where(np.equal(norm,0.), 0., logkavg)
# logkhalfdiff = np.where(np.equal(norm,0.), 0., logkhalfdiff)

logk_minus = np.stack((logkavg,logkhalfdiff),axis=-1)

#same for plus

logk_ups = np.stack(m_df.query("charge==1.0 & variation == 'Up'")['logK'].values,axis=-1)
logk_downs = np.stack(m_df.query("charge==1.0 & variation =='Down'")['logK'].values,axis=-1)

logkavg = 0.5*(logk_ups + logk_downs)
logkhalfdiff = 0.5*(logk_ups - logk_downs)

# #ensure that systematic tensor is sparse where normalization matrix is sparse
# logkavg = np.where(np.equal(norm,0.), 0., logkavg)
# logkhalfdiff = np.where(np.equal(norm,0.), 0., logkhalfdiff)

logk_plus = np.stack((logkavg,logkhalfdiff),axis=-1)

print(logk_plus.shape,logk_minus.shape)

logk = np.concatenate((logk_minus,logk_plus,np.zeros((1,nproc,2)),np.zeros((1,nproc,2))), axis=0)

# logk = np.concatenate((logk_minus,np.zeros((1,nproc,2))), axis=0)

# logk = np.zeros([5762,nproc,2,nsyst], dtype)
# print(nproc)
# logk[:,-1,:,0] = np.concatenate((logk_minus,logk_plus,np.zeros((1,289,2)),np.zeros((1,289,2))), axis=0)[:,-1,:]
logk = np.expand_dims(logk,axis=-1) #placeholder

# logk=np.nan_to_num(logk,0.0)

print(logk.shape)
nbytes += writeFlatInChunks(logk, f, "hlogk", maxChunkBytes = chunkSize)
logk = None

print("Total raw bytes in arrays = %d" % nbytes)
