import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
import f90nml
# sys.path.insert(0,os.path.join(os.path.dirname(__file__),'dependencies','pyWopwop'))
# import wopwop

#%%
import matplotlib.colors as mcolors

fontName = 'Times New Roman'
fontSize = 16
plt.rc('font',**{'family':'serif','serif':[fontName],'size':fontSize})
plt.rc('mathtext',**{'default':'regular'})
plt.rc('text',**{'usetex':False})
plt.rc('lines',**{'linewidth':2})

#
#%%

cases_directory = os.path.dirname(__file__)
cases = ['single_blade_no_gust','single_blade_gust']

data = {}
for case in cases:
    pred_data = {}
    #   imports reformatted data from wopwop in a dictionary
    with h5py.File(os.path.join(cases_directory,case, 'pressure.h5'), 'r') as dat_file:
        for i,n in dat_file.items():
            pred_data = {**pred_data,**{i:n[()]}}
    data = {**data,**{case:pred_data}}

corrected_ln = (data['single_blade_gust']['function_values'][:,:,:,2]-data['single_blade_no_gust']['function_values'][:,:,:,2]).squeeze()

# pred_data['geometry_values'] = np.flip(pred_data['geometry_values'],axis = 0).squeeze()
# pred_data['function_values'] = np.flip(pred_data['function_values'],axis = 0).squeeze()

#%%


for mic_iter in range(len(corrected_ln)):
    fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
    # ax.plot(pred_data['function_values'][mic_iter,:,0]/pred_data['function_values'][mic_iter,-1,0],pred_data['function_values'][mic_iter,:,1])
    ax.plot(data['single_blade_gust']['function_values'][mic_iter,:,:,0].squeeze(),corrected_ln[mic_iter])
    # ax.plot(pred_data['function_values'][mic_iter,:,0]/pred_data['function_values'][mic_iter,-1,0],pred_data['function_values'][mic_iter,:,-1])
    ax.grid()
    # ax.legend(['Thickness','Loading','Total'])
