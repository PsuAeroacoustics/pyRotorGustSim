#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'src'))
from help_funcs import *
import matplotlib.colors as mcolors


#%%

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ["Times New Roman"]
plt.rcParams['font.size'] = 16

def import_case_data(cases_directory):
    acs_data = {}
    saved_params = {}
    cases = [x for x  in os.listdir(os.path.join(cases_directory)) if os.path.isdir(os.path.join(cases_directory,x))]
    for case in cases:
        acs_data.update({case:import_results_from_wopwop(os.path.join(cases_directory,case,'acoustics'))})  
        saved_params.update({case:read_results_from_h5(os.path.join(cases_directory,case))})
    return acs_data,saved_params

#%%
cases_directory = os.getcwd()
baseline_sweep_name = 'Mt_oblique_baseline'
sweeps = ['Mt_parallel_equal_Mg','Mt_parallel_pnt_load']
leglab = ['Steady','Unsteady']

#%%

# sweeps = [baseline_sweep_name]+sweeps

bvispl = {}
M = {}

for sweep in sweeps:
    cases = [x for x  in os.listdir(os.path.join(cases_directory,sweep)) if os.path.isdir(os.path.join(cases_directory,sweep,x))]
    bvispl_temp = np.zeros(len(cases))
    M_temp = np.zeros(len(cases))
    for i,case in enumerate(cases):
        acs_data = import_results_from_wopwop(os.path.join(cases_directory,sweep,case,'acoustics'))
        saved_params=read_results_from_h5(os.path.join(cases_directory,sweep,case))
        # bvispl_temp[i] = np.mean(20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6))
        bvispl_temp[i] = 20*np.log10(np.sqrt(np.mean(acs_data['function_values'][5,:,:,-1]**2,axis  = -1))/20e-6).squeeze()
        M_temp[i] = saved_params['omega']*saved_params['R']/saved_params['sos']
    sort_ind = M_temp.argsort()
    M.update({sweep:M_temp[sort_ind]})
    bvispl.update({sweep:bvispl_temp[sort_ind]})

#%%

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .175,bottom = .15)
for sweep in sweeps:
    ax.plot(M[sweep],bvispl[sweep]-bvispl[sweep][0], marker = 'o')
ax.set_xlabel(r'$M_t$')
ax.set_ylabel(r'$BVI \ SPL, \ dB$')
ax.legend(leglab)
ax.grid()
plt.savefig(os.path.join(cases_directory,f'M_vs_bvispl.png'),format = 'png')

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .175,bottom = .15)
# for sweep in sweeps:
ax.plot(M[sweep],bvispl[sweeps[-1]]-bvispl[sweeps[0]], marker = 'o')
ax.set_xlabel(r'$M_t$')
ax.set_ylabel(r'$ \Delta  \ BVI \ SPL, \ dB$')
ax.legend(leglab)
ax.grid()
plt.savefig(os.path.join(cases_directory,f'M_vs_bvispl_2.png'),format = 'png')


# fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
# plt.subplots_adjust(left = .175,bottom = .15)
# for sweep in sweeps:
#     ax.plot(M[sweep][:-1],np.diff(bvispl[sweep]), marker = 'o')
# ax.set_xlabel(r'$M_t$')
# ax.set_ylabel(r'$ \Delta  \ BVI \ SPL, \ dB$')
# ax.legend(leglab)
# ax.grid()
# plt.savefig(os.path.join(cases_directory,f'M_vs_bvispl.png'),format = 'png')
