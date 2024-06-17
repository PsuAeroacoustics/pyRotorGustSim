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
cases = ['single_blade_arbitrary_gust_.1M','single_blade_arbitrary_gust_.2M','single_blade_arbitrary_gust_.3M','single_blade_arbitrary_gust_.4M','single_blade_arbitrary_gust_.5M','single_blade_arbitrary_gust_.6M']

data = {}
for case in cases:
    pred_data = {}
    #   imports reformatted data from wopwop in a dictionary
    with h5py.File(os.path.join(cases_directory,case, 'pressure.h5'), 'r') as dat_file:
        for i,n in dat_file.items():
            pred_data = {**pred_data,**{i:n[()]}}
    data = {**data,**{case:pred_data}}

# corrected_ln = (data['single_blade_gust']['function_values'][:,:,:,2]-data['single_blade_no_gust']['function_values'][:,:,:,2]).squeeze()

# pred_data['geometry_values'] = np.flip(pred_data['geometry_values'],axis = 0).squeeze()
# pred_data['function_values'] = np.flip(pred_data['function_values'],axis = 0).squeeze()

oaspl = np.zeros((len(cases),len(data[list(data.keys())[0]]['geometry_values'])))
for i,case in enumerate(cases):
    data[case]['geometry_values'] = np.flip(data[case]['geometry_values'],axis = 0).squeeze()
    data[case]['function_values'] = np.flip(data[case]['function_values'],axis = 0).squeeze()
    oaspl[i] = 20*np.log10(np.sqrt(np.mean(data[case]['function_values'][:,:,2]**2,axis = -1))/20e-6)
theta = np.round(np.arctan2(pred_data['geometry_values'][:,0,1],pred_data['geometry_values'][:,0,0])*180/np.pi)%(360)

#%%

M = np.arange(len(data))*.1+.1
for mic_iter in range(len(data[list(data.keys())[0]]['geometry_values'])):

    fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
    plt.subplots_adjust(left = .2,bottom = .15)
    ax.plot(M,oaspl[:,mic_iter],marker = 'o')
    ax.set_xlabel('$M_T$')
    ax.set_ylabel('OASPL, dB')
    ax.set_xticks(M)
    ax.set_ylim([70,110])
    ax.set_title(f'$\\theta = {theta[mic_iter]}$')
    ax.grid()
    plt.savefig(f'oaspl_vs_M_{mic_iter}.png',format = 'png')
    # ax.legend(['Thickness','Loading','Total'])
