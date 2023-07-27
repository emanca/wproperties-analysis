

import ROOT
import narf
import pandas as pd
import h5py 
import hist
import hdf5plugin
import math
import boost_histogram as bh
from utilities import boostHistHelpers as hh,common
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import re
from collections import OrderedDict
import pdb
import psutil
import argparse

parser = argparse.ArgumentParser('')
parser.add_argument('-bstr_idx', '--bstr_idx', type=int, help='bootstrap index')
parser.add_argument('-seed', '--seed', type=int, help='seed for random gen')

args = parser.parse_args()
bstr_idx = args.bstr_idx
seed = args.seed

# Function to get CPU and memory usage
def get_usage():
    # Get CPU usage as a percentage
    cpu_percent = psutil.cpu_percent()

    # Get memory usage in bytes and convert to megabytes (MB)
    memory_usage = psutil.Process().memory_info().rss / 1024 / 1024

    return cpu_percent, memory_usage

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

def writeSparse(indices, values, dense_shape, h5group, outname, maxChunkBytes = 1024**2):
  outgroup = h5group.create_group(outname)
  
  nbytes = 0
  nbytes += writeFlatInChunks(indices, outgroup, "indices", maxChunkBytes)
  nbytes += writeFlatInChunks(values, outgroup, "values", maxChunkBytes)
  outgroup.attrs['dense_shape'] = np.array(dense_shape, dtype='int64')
  
  return nbytes

def fillHelGroup(yBinsC,qtBinsC,helXsecs, flavour =''):
    helGroups = OrderedDict()
    for i in range(len(yBinsC)):
        for j in range(len(qtBinsC)):
            s = 'y_{i}_qt_{j}'.format(i=round(yBinsC[i],1),j=round(qtBinsC[j],1))
            if not flavour=='': s+='_'+flavour
            helGroups[s] = []
            
            for hel in helXsecs:
                helGroups[s].append('helXsec_'+hel+'_'+s)
            if helGroups[s] == []:
                del helGroups[s]
    return helGroups

def fillHelMetaGroup(yBinsC,qtBinsC,sumGroups,flavour =''):
    helMetaGroups = OrderedDict()
    for i in range(len(yBinsC)):
        s = 'y_{i}'.format(i=round(yBinsC[i],1))
        if not flavour=='': s+='_'+flavour
        helMetaGroups[s] = []
        for key in sumGroups:
            if s in key:
                helMetaGroups[s].append(key)
        
        if helMetaGroups[s] == []:
                del helMetaGroups[s]
    
    for j in range(len(qtBinsC)):
        s = 'qt_{j}'.format(j=round(qtBinsC[j],1))
        if not flavour=='': s+='_'+flavour
        helMetaGroups[s] = []
        for key in sumGroups:
            if 'qt' in key and key.split('_')[3]==str(round(qtBinsC[j],1)):
                helMetaGroups[s].append(key)
    
        if helMetaGroups[s] == []:
                del helMetaGroups[s]
    return helMetaGroups

def fillSumGroup(yBinsC,qtBinsC,helXsecs,processes,flavour =''):
    sumGroups = OrderedDict()
    for i in range(len(yBinsC)):
        s = 'y_{i}'.format(i=round(yBinsC[i],1))
        if not flavour=='': s+='_'+flavour
        for hel in helXsecs:
            for signal in processes:
                sumGroups['helXsec_'+hel+'_'+s] = []
                for j in range(len(qtBinsC)):
                    #if 'helXsecs'+hel+'_'+'y_{i}_qt_{j}'.format(i=i,j=j) in processes:
                    sumGroups['helXsec_'+hel+'_'+s].append('helXsec_'+hel+'_'+s+'_qt_{j}'.format(j=round(qtBinsC[j],1)))
    
    for j in range(len(qtBinsC)):
        s = 'qt_{j}'.format(j=round(qtBinsC[j],1))
        if not flavour=='': s+='_'+flavour
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
    mirr_histo = hist.Hist(*var.axes,downup_axis, name=var.name, data=data, storage = hist.storage.Weight())
    return mirr_histo


def setPreconditionVec():
    f=h5py.File('../Fit/FitRes/fit_Wplus_asimov.hdf5', 'r')
    hessian = f['hess'][:]
    eig, U = np.linalg.eigh(hessian)
    M1 = np.matmul(np.diag(1./np.sqrt(eig)),U.T)
    # print(M1,np.linalg.inv(np.linalg.inv(M1)))
    preconditioner = M1
    return preconditioner

def decorrelateInEta(nominal,rawvars):
    nEtaBins = 48
    j_indices = np.arange(nEtaBins)
    # print(nominal.shape,rawvars.shape)
    # create new histogram with expanded eta axis
    SFaxes = list(rawvars.axes)
    mod_axis = [axis for axis in SFaxes if axis.name=='SF eta'][0]
    idx = SFaxes.index(mod_axis)
    SFaxes[idx] = hist.axis.Regular(48, -2.4, 2.4, underflow=False, overflow=False, name='SF eta')
    dec_histo = hist.Hist(*SFaxes, name=rawvars.name,storage = hist.storage.Double())
    print("here")
    # create data by patching nominal and variations
    diff_shape = dec_histo.shape[len(nominal.shape):]
    # print(diff_shape)
    tmp = np.tile(nominal.values()[...,np.newaxis,np.newaxis,np.newaxis,np.newaxis], diff_shape)
    print("here 2")
    # for i  in range(nEtaBins):
    #     tmp[i,...,i,:,:,:] = rawvars.values()[i,...,0,:,:,:]
    tmp[j_indices,..., j_indices, :, :, :] = rawvars.values()[j_indices, ..., 0, :, :, :]
    # print(tmp.shape)
    print("here 3")
    dec_histo.values()[...] = tmp
    print("here 4")
    # print(rawvars.shape)
    # fig, ax1 = plt.subplots(figsize=(48,10))
    # hep.histplot(dec_histo[0,0,:,:,0, 4, 24, 0, 0, 0].values().ravel()/rawvars[0,0,:,:,0, 4, 0, 0, 0, 0].values().ravel())
    # hep.histplot(dec_histo[0,0,:,:,0, 4, 24, 0, 0, 0].values().ravel()/nominal[0,0,:,:,0, 4].values().ravel())
    # ax1.set_ylim(0.99,1.01)
    # plt.show()

    return dec_histo

def transport_smearing_weights_to_reco(hist_gensmear, nominal_reco, nominal_gensmear):
    print("here 0")
    hist_reco = hist.Hist(
                    *hist_gensmear.axes,
                    storage = hist_gensmear._storage_type()
                )
    print("here 1")
    bin_ratio = hh.divideHists(hist_gensmear, nominal_gensmear)
    print("here 2")
    hist_reco = hh.multiplyHists(nominal_reco, bin_ratio)
    print("here 3")
    return hist_reco

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~loading boost histograms and cross sections from templates hdf5 file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
f = h5py.File("templates_W_halfStat.hdf5","r")
f_bstrp = h5py.File("templates_W_halfStat.hdf5","r")
t = h5py.File('angCoeffWZ.hdf5','r')
results_nom = narf.ioutils.pickle_load_h5py(f["results"])
results = narf.ioutils.pickle_load_h5py(f_bstrp["results"])
# Hdata_obs = results['dataPostVFP']["output"]["data_obs"].get()

Hdata_obs = results_nom['WplusmunuPostVFP']['output']['signal_nominal'].get().project(*["mueta","mupt","charge","passMT","passIso"])
print(Hdata_obs)

#constants
# processes = ['WminusmunuPostVFP','WplusmunuPostVFP']
processes = ['WplusmunuPostVFP']
df = pd.DataFrame()

#TODO add bkg processes to this
# print(hh.findAxes([H],["Vrap","Vpt","helicities"]))
sumw = np.zeros_like(Hdata_obs.values().ravel())
sumw2 = np.zeros_like(Hdata_obs.values().ravel())

for process in processes:
    V ='V'
    lumi    = results['dataPostVFP']["lumi"]
    xsec    = results[process]["dataset"]["xsec"]
    weights = results[process]["weight_sum"]
    C       = lumi*1000*xsec/weights
    weights_nom = results_nom[process]["weight_sum"]
    C_nom       = lumi*1000*xsec/weights_nom
    # procs = ["lowacc"]
    procs = []
    systs_groups = {}

    Hdata_obs = hh.scaleHist(Hdata_obs, C_nom, createNew=False)

    #first add nominal boost histogram for signal
    H = results[process]['output']['signal_nominal'].get()#[...,bstr_idx]
    H = hh.scaleHist(H, C, createNew=False)
    H_nom = results_nom[process]['output']['signal_nominal'].get()
    H_nom = hh.scaleHist(H_nom, C_nom, createNew=False)

    print(H_nom)
    print(H)

    # import mplhep as hep
    # hep.histplot(H.project(*["mueta","mupt","charge","passMT","passIso"])[12,30,-1,True,True,:])
    # plt.show()

    nom_axes = [axis for axis in H.axes]
    nominal_cols = [axis.name for axis in nom_axes]
    unrolled_dim = np.prod(np.array([len(axis.centers) for axis in nom_axes if not axis.name in ["Vrap","Vpt","helicities"]]))
    unrolled_names = [axis.name for axis in nom_axes if not axis.name in ["Vrap","Vpt","helicities"]]
    index_dim = np.prod(np.array([len(axis.centers) for axis in nom_axes if axis.name in ["Vrap","Vpt","helicities"]]))
    index_centers = [axis.centers for axis in nom_axes if axis.name in ["Vrap","Vpt"]]
    helicities   = list(H.axes['helicities'])

    #Bin information
    yBinsC     = np.round(H.axes[V+'rap'].centers,1)
    qtBinsC    = np.round(H.axes[V+'pt'].centers,1)
    charges    = H.axes['charge'].centers
    eta        = H.axes['mueta'].centers
    pt         = H.axes['mupt'].centers

    qtBins = H.axes[V+'pt'].edges
    yBins = H.axes[V+'rap'].edges

    #Reshaping the data. 2d format. one row per unrolled pt-eta-charge(-mt-iso) distribution
    unrolled_and_stacked = H.values().reshape(unrolled_dim,index_dim).T
    print("unrolled_and_stacked",unrolled_and_stacked.shape)

    sumw += H_nom.project(*unrolled_names).values().ravel()
    sumw2 += H_nom.project(*unrolled_names).variances().ravel()

    #clean memory
    H = None
    #Generating multi index 
    multi = pd.MultiIndex.from_product([*index_centers,helicities] , names = ['rapidity', 'qt' , 'hel'])
    '''Building the nominal DataFrame'''
        
    #building dataframe
    df_proc = pd.Series(list(unrolled_and_stacked),index=multi, name="nominal")
    df_proc = pd.DataFrame(df_proc)
    print('\nnominal dataframe\n' , df_proc.head())

    unrolled_and_stacked = None
    #Adding cross section information to our dataframe by creating cross section dataframe and merging
    #TODO: pass boost histograms format

    threshold_y = np.digitize(2.4,yBins)-1
    threshold_qt = np.digitize(60.,qtBins)-1

    hxsecs = narf.ioutils.pickle_load_h5py(t["angCoeffWZ"])
    T = hxsecs[f"hist_coeffs_{process}"].values()[:threshold_y,:threshold_qt,:] #cross sections
    good_idx = [0,1,2,3,4,-1]
    T = T[...,good_idx]

    for i in range(5):
        T[...,i]*=T[...,-1]

    multi = pd.MultiIndex.from_product([*index_centers,helicities], names = ['rapidity', 'qt' , 'hel']) #multi index to be joined on
    s = pd.Series(T.ravel(), index = multi , name='xsec') #series carrying cross section information
    df_proc["xsec"] = s
    print('\nadded cross sections\n',df_proc.head())
    print('\nadded cross sections\n',df_proc.tail())

    df_proc = df_proc.reset_index() #promote all indices to columns
    #setting process as index & cleaning up by removing redundant information
    df_proc.set_index(['helXsec_'+df_proc['hel']+'_y_'+df_proc['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+df_proc['qt'].apply(lambda x: round(x,1)).apply(str)+f'_{process}'],inplace=True)

    df_proc.drop(columns=['rapidity','qt','hel'],inplace=True)
    df_proc.rename_axis(['process'] ,inplace=True)
    print('\nre-indexing\n',df_proc.head())

    # import mplhep as hep
    # hep.hist2dplot(df_proc['nominal'].loc[('helXsec_UL_y_0.2_qt_1.5')].reshape(48,60,2,2,2)[...,-1,1,1])
    # plt.show()

    #adding column for helicity group
    # df_proc['helgroups'] = df_proc.index.get_level_values(0).map(lambda x: re.search("y.+" , x).group(0))
    df_proc['isSignal']  = True
    # df_proc.loc[df_proc.index.str.contains('helXsec_A'),'isSignal']=True
    # df_proc.loc[df_proc.index.str.contains('helXsec_P'),'isSignal']=True
    # df_proc.loc[df_proc.index.str.contains('helXsec_L'),'isSignal']=True
    # df_proc.loc[df_proc.index.str.contains('helXsec_T'),'isSignal']=True
    # df_proc.loc[df_proc.index.str.contains('helXsec_I'),'isSignal']=True
    
    print('\nnominal dataframe\n' , df_proc.head())
    df = pd.concat([df, df_proc], axis=0)

print('\nnominal dataframe concat\n' , df.head())
print('\nnominal dataframe concat\n' , df.tail())

# import mplhep as hep
# plot = df['nominal'].loc[("helXsec_L_y_0.2_qt_1.5",1.0)].reshape(48,60)
# hep.hist2dplot(plot,vmin=-0.005,vmax=0)
# plt.show()
#now add other procs
# for proc in procs:
#     histo=C*results[process]['output']['{}_nominal'.format(proc)].get()
#     unrolled = histo.values().reshape(len(charges),-1)
#     histo=None
#     #add data
#     iterables_proc = [[proc],charges]
#     multi_proc = pd.MultiIndex.from_product(iterables_proc , names = ['process','charge'])
#     print(charges,iterables_proc,multi_proc)
#     df_proc = pd.Series(list(unrolled),index=multi_proc, name="nominal")
#     df_proc = pd.DataFrame(df_proc)
#     unrolled=None
#     df_proc['isSignal'] = False
#     df_proc['xsec'] = -1
#     df = pd.concat([df,df_proc],axis=0)
    
# print('\nreorganizing and adding other procs\n',df.head(),df.tail())

#add systematics

systs_macrogroups = {} # this is a list over groups of systematics
systs_macrogroups['mass']=['mass_var'] #mass and other no-constraint nuisances must be first in the list
# systs_macrogroups['muon_calibration']=['jpsi_var']
# systs_macrogroups['sf']=['effStatTnP_sf_reco','effStatTnP_sf_tracking','effStatTnP_sf_idip','effStatTnP_sf_trigger','effStatTnP_sf_iso','effSystTnP'] #these correspond to the names of histograms to recall from file
# systs_macrogroups['sf']=['effSystTnP'] #these correspond to the names of histograms to recall from file
# systs_macrogroups['prefire']=['muonL1PrefireStat_tensor','muonL1PrefireSyst_tensor','ecalL1Prefire_tensor']

procs = ["signal"]+procs #careful!! this must be the same order as before! -->REVIEW
multi = df.index
nominal_cols.append("downUpVar")

for process in processes:
    lumi    = results['dataPostVFP']["lumi"]
    xsec    = results[process]["dataset"]["xsec"]
    weights = results[process]["weight_sum"]
    C       = lumi*1000*xsec/weights
    #loop over systematics:
    for proc in procs:
        #get variations
        for syst,nuisances in systs_macrogroups.items():
            #TODO add exception for histograms not found
            syst_dfs = []
            print(syst,nuisances)
            for nuisance in nuisances:
                print('{proc}_{nuisance}'.format(proc=proc,nuisance=nuisance))
                syst_histo = results[process]['output']['{proc}_{nuisance}'.format(proc=proc,nuisance=nuisance)].get()#[...,bstr_idx]
                syst_histo = hh.scaleHist(syst_histo, C, createNew=False)
                print(syst_histo)
                print("done")
                axes = [axis for axis in syst_histo.axes]
                # decorrelate in eta if needed
                if 'sf' in syst:
                    nominal = C * results[process]['output']['{proc}_nominal'.format(proc=proc)].get()
                    print("here")
                    syst_histo = mirrorHisto(nominal,syst_histo)
                    # if not 'effSystTnP' in nuisance:
                    #     syst_histo = decorrelateInEta(nominal,syst_histo)
                    nominal = None
                if 'muon_calibration' in syst:
                    print("get {proc}_nominal".format(proc=proc))
                    nominal_reco = C * results[process]['output']['{proc}_nominal'.format(proc=proc)].get()
                    print("get {proc}_nominal_gensmear".format(proc=proc))
                    nominal_gensmear = C * results[process]['output']['{proc}_nominal_gensmear'.format(proc=proc)].get()
                    syst_histo = transport_smearing_weights_to_reco(syst_histo, nominal_reco, nominal_gensmear)
                #select slices in systematics based on "vars"
                
                syst_axes = [axis for axis in syst_histo.axes if axis.name not in nominal_cols]
                nsysts = 1 #this is the total number of systematics after considering all the bins
                for axis in syst_axes:
                    nsysts = nsysts*len(axis.centers)
                print("nsysts",nsysts)
                # syst_arr = np.ascontiguousarray(syst_histo.values(flow=True))
                syst_arr = syst_histo.values()
                print(syst_arr.flags.contiguous)
                print("values")
                syst_histo = None
                if proc == "signal":
                    #Reshaping the data. 2d format. one row per unrolled pt-eta-charge(-mt-iso) distribution
                    syst_arr = np.moveaxis(syst_arr.reshape(unrolled_dim,index_dim,nsysts,2),0,-2).reshape(index_dim*nsysts,unrolled_dim*2)# rapidity, qt, hel, charge, syst, data, up/down
                    # print(syst_arr.shape)
                    # syst_arr = syst_arr.reshape((unrolled_dim,index_dim,nsysts,2),order="A")
                    print(syst_arr.flags.contiguous)
                    print("here 5")
                    names = ['rapidity', 'qt' , 'hel']+[axis.name for axis in syst_axes]
                    iterables = [*index_centers,helicities] + [axis.centers for axis in syst_axes]
                    if syst_axes == []:
                        names.append(f"{nuisance}")
                        iterables.append([0.5])
                    multi = pd.MultiIndex.from_product(iterables, names = names)
                    # print(multi)
                    syst_df = pd.Series(list(syst_arr),index=multi,name=syst)
                    syst_arr = None
                    print("here 6")
                    print(syst_df.head())
                    syst_df = pd.DataFrame(syst_df).reset_index()
                    idx_strings = ['helXsec_'+syst_df['hel']+'_y_'+syst_df['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+syst_df['qt'].apply(lambda x: round(x,1)).apply(str)+f'_{process}']
                    syst_string = f"{nuisance}_"
                    for axis in syst_axes:
                        syst_string+=axis.name.replace(' ','')+'_'+syst_df[axis.name].apply(str)+'_'
                    if syst_axes == []:
                        syst_string+=syst_string+syst_df[syst_string[:-1]].apply(str)+'_'
                    syst_string = syst_string.apply(lambda s: s[:-1] if s.endswith('_') else s)
                    idx_strings.append(syst_string)
                    syst_df.set_index(idx_strings,inplace=True)
                    syst_df.drop(columns=[axis.name for axis in syst_axes],inplace=True)
                    syst_df.drop(columns=['rapidity','qt','hel'],inplace=True)
                    syst_df.rename_axis(['process','syst'] ,inplace=True)
                    syst_dfs.append(syst_df)
                else:
                    pass
                    #data, charge, syst, up/down
                    syst_arr = np.moveaxis(syst_arr,0,-2)
                    syst_arr = syst_arr.reshape(-1,unrolled_dim)
            # now merge all the dataframes within the group
            syst_df_merged = pd.concat(syst_dfs, axis=0)
            print(syst_df_merged.head())
            print(syst_df_merged.tail())
            # get list of systematics and drop index -->REVIEW
            syst_list = list(syst_df_merged.query(f"process == 'helXsec_L_y_0.2_qt_1.5_{process}'").index.get_level_values(1))
            # pdb.set_trace()
            systs_groups[syst]=syst_list
            syst_df_merged= syst_df_merged.droplevel('syst')
            # group by process and charge, and concatenate the arrays
            syst_df_merged = syst_df_merged.groupby(["process"])[syst].agg(np.concatenate)
            syst_df_merged = syst_df_merged.map(lambda x: x.reshape((-1,unrolled_dim,2)))
            df.loc[df.index.str.contains(process),syst]=syst_df_merged

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
    mask = np.stack([-1*np.ones_like(df[syst].loc[('helXsec_L_y_0.2_qt_1.5_WplusmunuPostVFP')][...,0]),np.ones_like(df[syst].loc[('helXsec_L_y_0.2_qt_1.5_WplusmunuPostVFP')][...,0])],axis=-1)
    print(mask.shape)
    df["{}_logk".format(syst)]=df["{}_logk".format(syst)].apply(lambda x: mask*x)
    
    df.drop(columns=[syst],inplace=True)
    print('\nfinal df\n',df.head(),df.tail())

print(df.memory_usage(deep=True))


'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

#retrieve metadata

procs = list(df.index.get_level_values(0))
signals = list(df.query("isSignal==True").index.get_level_values(0))
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
helGroups_plus = fillHelGroup(yBinsC,qtBinsC,helicities, flavour = "WplusmunuPostVFP")
helGroups_minus = fillHelGroup(yBinsC,qtBinsC,helicities, flavour = "WminusmunuPostVFP")

# helGroups_minus.update(helGroups_plus)
# helGroups = helGroups_minus
helGroups = helGroups_plus
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
# sumGroups = fillSumGroup(yBinsC,qtBinsC,helicities,signals)
# for igroup,group in enumerate(sumGroups):
#   sumgroups.append(group)
#   for proc in sumGroups[group]:
#     sumgroupsegmentids.append(igroup)
#     sumgroupidxs.append(procs.index(proc))
    
#list of groups of signal processes by chargemeta - DON'T NEED THAT
chargemetagroups = []
chargemetagroupidxs = []
  
#list of groups of signal processes by ratiometa - DON'T NEED THAT
ratiometagroups = []
ratiometagroupidxs = []

#list of groups of signal processes by helmeta
helmetagroups = []
helmetagroupidxs = []
# helMetaGroups = fillHelMetaGroup(yBinsC,qtBinsC,sumGroups)
# for group in helMetaGroups:
#   helmetagroups.append(group)
#   helmetagroupidx = []
#   for proc in helMetaGroups[group]:
#     helmetagroupidx.append(sumgroups.index(proc))
#   helmetagroupidxs.append(helmetagroupidx)

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
outfilename = f"Wplus_halfStat.hdf5"
# outfilename = f"bstrp_inputs/Wplus_bstrp_{bstr_idx}.hdf5"
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


# data_obs = np.concatenate((Hdata_obs.values()[...,0].ravel(),Hdata_obs.values()[...,1].ravel()))
data_obs = Hdata_obs.values()[...].ravel()
data_obs = np.clip(data_obs, 0, np.inf)

np.random.seed(seed)
data_obs_rand = np.random.poisson(lam=data_obs, size=None)
data_obs_rand = np.array(data_obs_rand,dtype='float64')

Hdata_obs = None

nbytes += writeFlatInChunks(data_obs_rand, f, "hdata_obs", maxChunkBytes = chunkSize)
data_obs = None

#compute poisson parameter for Barlow-Beeston bin-by-bin statistical uncertainties
kstat = np.square(sumw)/(sumw2)
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
norm = np.concatenate((np.stack(df['nominal'].values,axis=-1),np.expand_dims(np.stack(df['xsec'].values,axis=-1),axis=0)),axis=0)
# nbytes += writeFlatInChunks(norm, f, "hnorm", maxChunkBytes = chunkSize)

nonzero = np.nonzero(norm)

norm_sparse_indices = np.argwhere(norm).astype(np.int32)
norm_sparse_values = norm[nonzero].reshape([-1])
norm_sparse_dense_shape = norm.shape
print("norm_sparse_dense_shape",norm_sparse_dense_shape)
print("norm_sparse_indices",norm_sparse_indices.shape)
print("norm_sparse_values",norm_sparse_values.shape)


nbytes += writeSparse(norm_sparse_indices, norm_sparse_values, norm_sparse_dense_shape, f, "hnorm_sparse", maxChunkBytes = chunkSize)
logk_sparse_dense_shape = (norm_sparse_indices.shape[0], 2*nsyst)

# df = df.drop(columns=['nominal'],inplace=True) #why this doesn't work?
print('\ndrop nominal\n',df.head())
logk_systs = []
for syst in systs_groups:
    print(df["{}_logk".format(syst)].values[0].shape)
    print(df['nominal'].values[0].shape)
    # logk_syst = np.moveaxis(np.concatenate((np.stack(df.query("charge==-1.")["{}_logk".format(syst)].values,axis=-2),np.stack(df.query("charge==1.")["{}_logk".format(syst)].values,axis=-2)),axis=1),0,-1)
    logk_syst = np.moveaxis(np.stack(df["{}_logk".format(syst)].values,axis=-2),0,-1)
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

logk = np.stack((logkavg,logkhalfdiff),axis=-2)
# logk = np.concatenate((logk,np.zeros((2,nproc,2,nsyst))),axis=0)
logk = np.concatenate((logk,np.zeros((1,nproc,2,nsyst))),axis=0)

print(logk.shape)

logk = logk.reshape([norm.shape[0],nproc,2*nsyst])

norm = None
print(logk.shape)

logk = logk[nonzero]

print("after nonzero", logk.shape)

nonzero = np.nonzero(logk)

print(np.array(nonzero).shape, np.array(np.transpose(nonzero)).shape)
# nonzero = np.array(nonzero)
logk_sparse_indices = np.argwhere(logk).astype(np.int32)
# logk_sparse_indices = np.transpose(np.array(nonzero).astype(np.int32))
logk_sparse_values = logk[nonzero].reshape([-1])

print("logk_sparse_dense_shape",logk_sparse_dense_shape)
print("logk_sparse_indices",logk_sparse_indices.shape)
print("logk_sparse_values",logk_sparse_values.shape)

nbytes += writeSparse(logk_sparse_indices, logk_sparse_values, logk_sparse_dense_shape, f, "hlogk_sparse", maxChunkBytes = chunkSize)
norm_sparse_indices = None
norm_sparse_values = None
logk_sparse_indices = None
logk_sparse_values = None

# nbytes += writeFlatInChunks(logk, f, "hlogk", maxChunkBytes = chunkSize)
# logk = None

print("Total raw bytes in arrays = %d" % nbytes)

# Call the function to get CPU and memory usage after the script finishes
cpu_usage, memory_usage = get_usage()
print(f"CPU Usage: {cpu_usage}%")
print(f"Memory Usage: {memory_usage} MB")