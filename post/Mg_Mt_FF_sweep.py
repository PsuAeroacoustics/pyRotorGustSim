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

cases_dir = '/Users/danielweitsman/codes/github/DanWeitsman/unsteady_BEMT/cases/far_field_param_sweep'
cases = [x for x  in os.listdir(cases_dir) if os.path.isdir(os.path.join(cases_dir,x))]

oaspl = []
Mt = np.zeros(len(cases))
Mg = np.zeros(len(cases))
gamma = np.zeros(len(cases))

for i,case in enumerate(cases):
    saved_params = read_results_from_h5(os.path.join(cases_dir,case))
    Mt[i] = saved_params['omega']*saved_params['R']/saved_params['sos']
    Mg[i] = np.round(saved_params['v_gust'][-int(2*np.pi/saved_params['dpsi']):,-1].max()/saved_params['sos'],3)

    with open(os.path.join(cases_dir,f"param_{case.strip('case_')}.json")) as param_file:
        params = json.load(param_file)
        gamma[i] = params['gust_params']['strength']
        
    acs_data = import_results_from_wopwop(os.path.join(cases_dir,case,'acoustics'))
    oaspl.append(20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6))
oaspl = np.array(oaspl)

theta = np.arctan2(acs_data['geometry_values'][:,:,0,1],acs_data['geometry_values'][:,:,0,0])%(2*np.pi)*180/np.pi
phi = np.arctan2(acs_data['geometry_values'][:,:,0,-1],np.linalg.norm((acs_data['geometry_values'][:,:,0,0],acs_data['geometry_values'][:,:,0,1]),axis = 0))*180/np.pi
r = np.linalg.norm(acs_data['geometry_values'][:,:,0],axis = -1)

Lp = 10*np.log10(4*np.pi*r[-1]**2*10**(oaspl[:,-1]/10)).squeeze()
sspreading = (Lp-10*np.log10(4*np.pi*r**2)).T
ff_ind = np.abs(np.abs(oaspl.squeeze()-sspreading)-1).argmin(axis = -1)

Mt_sort_ind = Mt.argsort()
gamma_sort_ind = gamma[Mt_sort_ind].argsort(kind='mergesort')
Mg_sort_ind = Mg[Mt_sort_ind].argsort(kind='mergesort')

new_shape = (9,9)

cmap = plt.cm.Spectral.reversed()

levels = np.linspace(0,4,17)
levels_c = np.linspace(0,4,9)

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .15,top = .925,right = .95,bottom = .13)
dist = ax.tricontourf(Mt.squeeze(),gamma.squeeze(),r[ff_ind].squeeze()/saved_params['R'],levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
dist2 = ax.tricontour(Mt.squeeze(),gamma.squeeze(),r[ff_ind].squeeze()/saved_params['R'],levels = levels_c,colors = 'k')
plt.clabel(dist2,levels=levels_c)
cbar = fig.colorbar(dist,pad = .05)
cbar.ax.set_ylabel(r'$r/R$')
cbar.set_ticks(levels[::2])
ax.set_xlabel(r'$M_T$')
ax.set_ylabel(r'$M_g$')
# ax.set_yticks(np.arange(5)*(.6-.2)/4+.2)
# ax.set_xticks(np.arange(3)*(.3-.1)/2+.05)
# plt.savefig(os.path.join(os.path.dirname(cases_dir),f'Mg_Mt_ff_oaspl_carpet.png'),format = 'png')
