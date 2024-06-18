import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'src'))
from help_funcs import *


#%%

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ["Times New Roman"]
plt.rcParams['font.size'] = 16

#
#%%

cases_directory = os.getcwd()
cases = [x for x  in os.listdir(cases_directory) if os.path.isdir(os.path.join(cases_directory,x))]

# data = {}
# for case in cases:
#     pred_data = {}
#     #   imports reformatted data from wopwop in a dictionary
#     with h5py.File(os.path.join(cases_directory,case,'acoustics', 'pressure.h5'), 'r') as dat_file:
#         for i,n in dat_file.items():
#             pred_data = {**pred_data,**{i:n[()]}}
#     data = {**data,**{case:pred_data}}

acs_data = {}
saved_params = {}

for case in cases:
    saved_params.update({case:read_results_from_h5(os.path.join(cases_directory,case))})
    acs_data.update({case:import_results_from_wopwop(saved_params[case]['acs_dir'].decode())})

# corrected_ln = (data['single_blade_gust']['function_values'][:,:,:,2]-data['single_blade_no_gust']['function_values'][:,:,:,2]).squeeze()

# pred_data['geometry_values'] = np.flip(pred_data['geometry_values'],axis = 0).squeeze()
# pred_data['function_values'] = np.flip(pred_data['function_values'],axis = 0).squeeze()

oaspl = np.zeros((len(cases),len(acs_data[list(acs_data.keys())[0]]['geometry_values'])))
for i,case in enumerate(cases):
    acs_data[case]['geometry_values'] = np.flip(acs_data[case]['geometry_values'],axis = 0).squeeze()
    acs_data[case]['function_values'] = np.flip(acs_data[case]['function_values'],axis = 0).squeeze()
    oaspl[i] = 20*np.log10(np.sqrt(np.mean(acs_data[case]['function_values'][:,:,-1]**2,axis = -1))/20e-6)
theta = np.round(np.arctan2(acs_data[cases[0]]['geometry_values'][:,0,1],acs_data[cases[0]]['geometry_values'][:,0,0])*180/np.pi)%(360)

#%%

M = np.array([saved_params[case]['v_gust'].max()/340 for case in cases])
# M = np.array([saved_params[case]['omega']*.19685/340 for case in cases])
plot_ind = M.argsort()
leg_labs = [f'$\\theta = {theta[x]}^\circ$' for x in np.arange(oaspl.shape[-1])][::-1]

for mic_iter in range(oaspl.shape[-1]):

    fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
    plt.subplots_adjust(left = .2,bottom = .15)
    ax.plot(M[plot_ind],oaspl[plot_ind,mic_iter],marker = 'o')
    ax.set_xlabel('$M_g$')
    ax.set_ylabel('OASPL, dB')
    ax.set_ylim([70,110])
    ax.set_title(f'$\\theta = {theta[mic_iter]}^\circ$')
    ax.grid()
    plt.savefig(os.path.join(cases_directory,f'oaspl_vs_Mg_{mic_iter}.png'),format = 'png')
    # ax.legend(['Thickness','Loading','Total'])

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .2,bottom = .15)
ax.plot(M[plot_ind],oaspl[plot_ind,:][:,::-1],marker = 'o')
ax.set_xlabel('$M_g$')
ax.set_ylabel('OASPL, dB')
ax.set_ylim([70,120])
ax.legend(labels =leg_labs)
ax.grid()
plt.savefig(os.path.join(cases_directory,f'oaspl_vs_Mg.png'),format = 'png')
