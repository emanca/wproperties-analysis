

import ROOT
import narf
import pandas as pd
import h5py 
import hist
import hdf5plugin
import boost_histogram as bh
import numpy as np
import matplotlib.pyplot as plt
import re


'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~loading boost histograms and cross sections from templates hdf5 file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
f = h5py.File("templatesTest2.hdf5","r")
results = narf.ioutils.pickle_load_h5py(f["results"])
H = results['ZmumuPostVFP']['output']['signalTemplates_nominal'].get() #boost histogram values
t = h5py.File('templatesFit.hdf5','r') # templates file used for getting the cross sections (later will also be boost histogram)
H3 = results['ZmumuPostVFP']['output']['signalTemplates_mass'].get() #boost histogram with systematic mass variations
M = H3[:,:,:,:,:,:,[bh.loc('massShift100MeVDown') , bh.loc('massShift100MeVUp')]] #boost histogram for selected mass variations


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
df.set_index('helXsec_'+df['hel']+'_y_'+df['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+df['qt'].apply(lambda x: round(x,1)).apply(str),inplace=True)
df.drop(columns=['rapidity','qt','hel'],inplace=True)
df.rename_axis('process' ,inplace=True)

#reorganizing data into a single column labeled 'data' - contains unrolled pt/eta distribution
df['data'] = df.loc[:,0:2879].apply(np.hstack , axis=1) 
df.drop(columns = df.loc[:,0:2879].columns , inplace = True)

#adding column for helicity group
df['helgroups'] = df.index.map(lambda x: re.search("y.+" , x).group(0))
print('Nominal Dataframe:')
print(df.head())





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
m_df.set_index('helXsec_'+m_df['hel']+'_y_'+m_df['rapidity'].apply(lambda x: round(x,1)).apply(str)+'_qt_'+m_df['qt'].apply(lambda x: round(x,1)).apply(str),inplace=True)
m_df.drop(columns=['rapidity','qt','hel'],inplace=True)
m_df.rename_axis('process' ,inplace=True)


#reorganizing data into a single column labeled 'data' - contains unrolled pt/eta distribution
m_df['data'] = m_df.loc[:,0:2879].apply(np.hstack , axis=1) 
m_df.drop(columns = m_df.loc[:,0:2879].columns , inplace = True)

print('\nsystematics dataframe')
print(m_df.head())