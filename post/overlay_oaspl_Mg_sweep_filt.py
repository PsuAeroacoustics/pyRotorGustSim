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

#%%

par_dir = os.getcwd()
cases_dir = [os.path.join(par_dir,'scaled_Mg_sweep'),os.path.join(par_dir,'scaled_Mg_sweep_sdof_dist')]

saved_params = {}
oaspl = {}

for case_dir in cases_dir:
    
    cases = [x for x  in os.listdir(case_dir) if os.path.isdir(os.path.join(case_dir,x))]

    oaspl_temp = []
    saved_params_temp = {}

    for i,case in enumerate(cases):
        saved_params_temp.update({case:read_results_from_h5(os.path.join(case_dir,case))})
        acs_data = import_results_from_wopwop(os.path.join(case_dir,case,'acoustics'))
        
        oaspl_temp.append(20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6))

    saved_params.update({os.path.basename(case_dir):saved_params_temp})
    oaspl.update({os.path.basename(case_dir):np.array(oaspl_temp).squeeze()})

theta = (np.round(np.arctan2(acs_data['geometry_values'][:,:,0,1],acs_data['geometry_values'][:,:,0,0])*180/np.pi)%(360)).squeeze()

#%%
M = np.array([saved_params[list(oaspl.keys())[0]][case]['v_gust'].max()/340 for case in cases])
plot_ind = M.argsort()
leg_labs = [f'$\\theta = {theta[x]}^\circ$' for x in np.arange(oaspl[list(oaspl.keys())[0]].shape[-1])]

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .2,bottom = .15)
ax.plot(M[plot_ind],(oaspl[list(oaspl.keys())[-1]]-oaspl[list(oaspl.keys())[0]])[plot_ind,:],marker = 'o')
ax.set_xlabel(r'$M_g$')
ax.set_ylabel(r'$\Delta \ OASPL, \ dB$')
ax.legend(labels =leg_labs,fontsize = 10)
ax.set_ylim([-6,5])
ax.grid()
plt.savefig(os.path.join(par_dir,f'd_oaspl_{os.path.basename(cases_dir[-1])}.png'),format = 'png')
plt.close()

