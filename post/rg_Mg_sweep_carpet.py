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
plt.rcParams['font.size'] = 10

#%%

def get_oaspl(cases_dir):
    cases = [x for x  in os.listdir(cases_dir) if os.path.isdir(os.path.join(cases_dir,x))]
    oaspl = []
    for i,case in enumerate(cases):
        acs_data = import_results_from_wopwop(os.path.join(cases_dir,case,'acoustics'))
        oaspl.append(20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6))
    return oaspl

baseline_dir = '/home/dweitsma/codes/github/DanWeitsman/rotor_gust_interaction/cases/nominal/param_sweeps/rg_Mg_bl'
cases_dir = '/home/dweitsma/codes/github/DanWeitsman/rotor_gust_interaction/cases/nominal/param_sweeps/rg_Mg'
eps = True

cases = [x for x  in os.listdir(cases_dir) if os.path.isdir(os.path.join(cases_dir,x))]

oaspl_bl,oaspl_filt = list(map(lambda x: np.array(get_oaspl(x)),[baseline_dir,cases_dir]))

d_oaspl = np.mean(np.mean((oaspl_filt-oaspl_bl),axis = -1),axis = -1)

saved_params = {}
for case in cases:
    saved_params.update({case:read_results_from_h5(os.path.join(cases_dir,case))})
    
gust_width = np.zeros(len(cases))
Mg = np.zeros(len(cases))
gust_strength = np.zeros(len(cases))

for i,case in enumerate(cases):
    
    max_ind = saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,-1].argmax()
    start_ind = saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,-1][:max_ind+1].argmin()
    end_ind = np.abs(saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,-1][max_ind:]-.645*saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,-1][max_ind]).argmin()
    gust_width[i] = (max_ind-start_ind+end_ind)*saved_params[case]['dpsi']*saved_params[case]['r'][-1]*saved_params[case]['R']/saved_params[case]['c']
    
    if len(case) ==7:
        with open(os.path.join(cases_dir,f'param_{case[-2:]}.json')) as param_file:
            gust_strength[i] = json.load(param_file)['gust_params']['strength']
    else:
        with open(os.path.join(cases_dir,f'param_{case[-1:]}.json')) as param_file:
            gust_strength[i] = json.load(param_file)['gust_params']['strength']

    Mg[i] = saved_params[case]['v_gust'][max_ind,-1]/saved_params[case]['sos']

gust_width_sort_ind = gust_width.argsort()[6:]
Mg_sort_ind = gust_strength[gust_width_sort_ind].argsort(kind = 'mergesort')

cmap = plt.cm.Spectral.reversed()

levels = np.linspace(-10,0,21)
levels_c = np.linspace(-10,0,11)

figsize =((8.27-2)*.49,(11.69-2)*.35)

fig,ax = plt.subplots(1,1, figsize = figsize)
plt.subplots_adjust(left = .2,top = .925,right = .85,bottom = .13)
dist = ax.contourf(gust_width[gust_width_sort_ind].reshape(8,9).T,Mg[gust_width_sort_ind][Mg_sort_ind].reshape(9,8),d_oaspl[gust_width_sort_ind][Mg_sort_ind].reshape(9,8),levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
dist2 = ax.contour(gust_width[gust_width_sort_ind].reshape(8,9).T,Mg[gust_width_sort_ind][Mg_sort_ind].reshape(9,8),d_oaspl[gust_width_sort_ind][Mg_sort_ind].reshape(9,8),levels = levels_c,colors = 'k')
plt.clabel(dist2,levels=levels_c)
cbar = fig.colorbar(dist,pad = .05)
cbar.ax.set_ylabel(r'$\Delta \ BVISPL, \ dB \ (re: \ 20 \mu Pa)$')
cbar.set_ticks(levels[::2])
ax.set_xticks(np.arange(1,5))
ax.set_xlabel(r'$r_g/c$')
ax.set_ylabel(r'$M_g$')
ax.set_ylim([0,.3])
plt.savefig(os.path.join(os.path.dirname(cases_dir),f'rg_Mg_d_oaspl_carpet.png'),format = 'png')
if eps:
    plt.savefig(os.path.join(os.path.dirname(cases_dir),f'rg_Mg_d_oaspl_carpet.eps'),format = 'eps')
plt.close()