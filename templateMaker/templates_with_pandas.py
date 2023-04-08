

import ROOT
import narf
import pandas as pd
import h5py 
import hist
import hdf5plugin
import boost_histogram as bh
import numpy as np
import matplotlib.pyplot as plt


'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~loading boost histogram and extracting as numpy array~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
f = h5py.File("templatesTest2.hdf5","r")
results = narf.ioutils.pickle_load_h5py(f["results"])
H = results['ZmumuPostVFP']['output']['signalTemplates_nominal'].get() #boost histogram values
t = h5py.File('templatesFit.hdf5','r') # templates file used for getting the cross sections (later will also be boost histogram)


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
df.set_index('helXsec_'+df['hel']+'_y_'+df['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+df['qt'].apply(lambda x: round(x,1)).apply(str),inplace=True)
df.drop(columns=['rapidity','qt','hel'],inplace=True)
print(df.head())