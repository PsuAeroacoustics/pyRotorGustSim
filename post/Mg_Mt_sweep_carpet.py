import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import sys
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'src'))
from help_funcs import *

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ["Times New Roman"]
plt.rcParams['font.size'] = 16

#%%

def get_oaspl(cases_dir):
    cases = [x for x  in os.listdir(cases_dir) if os.path.isdir(os.path.join(cases_dir,x))]
    oaspl = []
    for i,case in enumerate(cases):
        acs_data = import_results_from_wopwop(os.path.join(cases_dir,case,'acoustics'))
        oaspl.append(20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6))
    return oaspl

baseline_dir = '/home/dweitsma/codes/github/DanWeitsman/rotor_gust_interaction/cases/nominal/param_sweeps/Mg_Mt_bl'
cases_dir = '/home/dweitsma/codes/github/DanWeitsman/rotor_gust_interaction/cases/nominal/param_sweeps/Mg_Mt'

cases = [x for x  in os.listdir(cases_dir) if os.path.isdir(os.path.join(cases_dir,x))]

oaspl_bl,oaspl_filt = list(map(lambda x: np.array(get_oaspl(x)),[baseline_dir,cases_dir]))

d_oaspl = np.mean(np.mean((oaspl_filt-oaspl_bl),axis = -1),axis = -1)

saved_params = {}
for case in cases:
    saved_params.update({case:read_results_from_h5(os.path.join(cases_dir,case))})
    
gust_width = np.zeros(len(cases))
Mg = np.zeros(len(cases))
Mt = np.zeros(len(cases))
gust_strength = np.zeros(len(cases))

for i,case in enumerate(cases):
    
    max_ind = saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,-1].argmax()
    start_ind = saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,-1][:max_ind+1].argmin()
    end_ind = np.abs(saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,-1][max_ind:]-.1*saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,-1][max_ind]).argmin()
    gust_width[i] = (max_ind-start_ind+end_ind)*saved_params[case]['dpsi']*saved_params[case]['r'][-1]*saved_params[case]['R']/saved_params[case]['c']
    
    # if len(case) ==7:
    #     with open(os.path.join(cases_dir,f'param_{case[-2:]}.json')) as param_file:
    #         gust_strength[i] = json.load(param_file)['gust_params']['strength']
    # else:
    #     with open(os.path.join(cases_dir,f'param_{case[-1:]}.json')) as param_file:
    #         gust_strength[i] = json.load(param_file)['gust_params']['strength']
    Mt[i] = saved_params[case]['omega']*saved_params[case]['R']/saved_params[case]['sos']
    Mg[i] = saved_params[case]['v_gust'][max_ind,-1]/saved_params[case]['sos']

Mg_sort_ind = Mg.argsort()
Mt_sort_ind = Mt[Mg_sort_ind].argsort(kind = 'mergesort')
new_shape = (np.sqrt(len(cases))*np.ones(2)).astype(int)

cmap = plt.cm.Spectral.reversed()

levels = np.linspace(-10,0,21)
levels_c = np.linspace(-10,0,11)

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(right = .85,left = .2,bottom = .15)
dist = ax.contourf(Mg[Mg_sort_ind].reshape(new_shape).T,Mt[Mg_sort_ind][Mt_sort_ind].reshape(new_shape),d_oaspl[Mg_sort_ind][Mt_sort_ind].reshape(new_shape),levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
dist2 = ax.contour(Mg[Mg_sort_ind].reshape(new_shape).T,Mt[Mg_sort_ind][Mt_sort_ind].reshape(new_shape),d_oaspl[Mg_sort_ind][Mt_sort_ind].reshape(new_shape),levels = levels_c,colors = 'k')
plt.clabel(dist2,levels=levels_c)
cbar = fig.colorbar(dist,pad = .05)
cbar.ax.set_ylabel(r'$\Delta \ BVISPL, \ dB \ (re: \ 20 \mu Pa)$')
cbar.set_ticks(levels[::2])
ax.set_ylabel(r'$M_T$')
ax.set_xlabel(r'$M_g$')
ax.set_ylim([0.2,.6])
ax.set_xticks(np.arange(3)*(.3-.1)/2+.1)
plt.savefig(os.path.join(os.path.dirname(cases_dir),f'Mg_Mt_d_oaspl_carpet.png'),format = 'png')
