import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import sys
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'src'))
from help_funcs import *
from scipy.interpolate import CubicSpline

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ["Times New Roman"]
plt.rcParams['font.size'] = 10

#%%

def get_oaspl(cases_dir):
    # cases = [x for x  in os.listdir(cases_dir) if os.path.isdir(os.path.join(cases_dir,x))]
    oaspl = []
    for i,case in enumerate(cases):
        acs_data = import_results_from_wopwop(os.path.join(cases_dir,case,'acoustics'))
        oaspl.append(20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6))
    return oaspl

baseline_dir = '/Users/danielweitsman/codes/github/DanWeitsman/unsteady_BEMT/cases/final_designs/sweeps/rg_Mg/baseline_untapered_AR8_unsteady'
cases_dir = '/Users/danielweitsman/codes/github/DanWeitsman/unsteady_BEMT/cases/final_designs/sweeps/rg_Mg/sdof_geom_untapered_AR8_unsteady'
eps = False

cases = [x for x  in os.listdir(cases_dir) if os.path.isdir(os.path.join(cases_dir,x))]

oaspl_bl,oaspl_filt = list(map(lambda x: np.array(get_oaspl(x)),[baseline_dir,cases_dir]))

d_oaspl = (oaspl_filt-oaspl_bl).squeeze().max(axis = 1)

geometry_values = import_results_from_wopwop(os.path.join(cases_dir,cases[0],'acoustics'))['geometry_values']
theta = np.round(np.arctan2(geometry_values[:,:,0,1],geometry_values[:,:,0,0])%(2*np.pi)*180/np.pi,1).squeeze()
phi = np.round(np.arctan2(geometry_values[:,:,0,-1],np.linalg.norm((geometry_values[:,:,0,0],geometry_values[:,:,0,1]),axis = 0))%(2*np.pi)*180/np.pi,1).squeeze()


N = 500
theta_interp = np.arange(N)*(theta[-1]-theta[0])/(N-1)+theta[0]

doaspl_interp = [CubicSpline(x = theta,y = oaspl.squeeze(),axis = -1) for oaspl in [oaspl_bl,oaspl_filt]]
max_oaspl_ind = [i(theta_interp).argmax(axis = -1) for i in doaspl_interp]

dtheta_max = np.diff((theta_interp*np.ones((len(cases),N)))[np.arange(len(cases)),max_oaspl_ind],axis = 0).squeeze()

d_oaspl = doaspl_interp[1](theta_interp)[np.arange(len(cases)),max_oaspl_ind[0]]-doaspl_interp[0](theta_interp)[np.arange(len(cases)),max_oaspl_ind[0]]
# d_oaspl = doaspl_interp[0](theta_interp)[np.arange(len(cases)),max_oaspl_ind[0]]


saved_params_bl = {}
saved_params_filt = {}

for case in cases:
    saved_params_bl.update({case:read_results_from_h5(os.path.join(baseline_dir,case))})
    saved_params_filt.update({case:read_results_from_h5(os.path.join(cases_dir,case))})

gust_width = np.zeros(len(cases))
Mg = np.zeros(len(cases))
# h = np.zeros((len(cases),int(2*np.pi/saved_params[case]['dpsi']),saved_params[cases[0]]['N_elements']))
# v_gust= np.zeros((len(cases),int(2*np.pi/saved_params[case]['dpsi']),saved_params[cases[0]]['N_elements']))
h = []
v_gust = []
loads = []
r_ind_select = int(0.9*saved_params_bl[cases[0]]['N_elements'])  # Select the last 25% of the rotor radius for gust width calculatio

for i,case in enumerate(cases):
    h.append(saved_params_filt[case]['h_gust'][-int(2*np.pi/saved_params_bl[case]['dpsi']):]*saved_params_bl[case]['AR'])
    v_gust.append(saved_params_filt[case]['v_gust'][-int((2*np.pi)/saved_params_bl[case]['dpsi']):])
    loads.append(saved_params_filt[case]['loads'][-int((2*np.pi)/saved_params_bl[case]['dpsi']):])

    max_ind = v_gust[i][:,r_ind_select].argmax(axis = 0)
    start_ind = v_gust[i][:,r_ind_select][:max_ind+1].argmin(axis = 0)+1
    end_ind=np.abs(v_gust[i][max_ind:,r_ind_select]-0.2*v_gust[i][max_ind,r_ind_select]).argmin()
    gust_width[i] = h[i][max_ind,r_ind_select]-h[i][:max_ind+1,r_ind_select][start_ind]


    # if len(case) ==7:
    #     with open(os.path.join(cases_dir,f'param_{case[-2:]}.json')) as param_file:
    #         gust_strength[i] = json.load(param_file)['gust_params']['strength']
    # else:
    #     with open(os.path.join(cases_dir,f'param_{case[-1:]}.json')) as param_file:
    #         gust_strength[i] = json.load(param_file)['gust_params']['strength']

    Mg[i] = v_gust[i][max_ind,r_ind_select]/saved_params_bl[case]['sos']

gust_width_sort_ind = gust_width.argsort()
# Mg_sort_ind = Mg[gust_width_sort_ind].argsort(kind = 'mergesort')

Mg_sort_ind = np.asarray([i[np.argsort(Mg[i])] for i in gust_width_sort_ind.reshape(9,9)]).flatten()

# fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# ax.plot(h[0].T[::5].T,-v_gust[0][::5].T)

fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# plt.subplots_adjust(left = .2,top = .925,right = .85,bottom = .13)
levels = np.linspace(-4,4,17)
levels_c = np.linspace(-4,4,17)

dist = ax.contourf(gust_width[Mg_sort_ind].reshape(9,9),Mg[Mg_sort_ind].reshape(9,9),d_oaspl[Mg_sort_ind].reshape(9,9),levels = levels,cmap = plt.cm.Spectral.reversed(),norm = mcolors.CenteredNorm())

# levels = np.linspace(95,115,41)
# dist = ax.contourf(gust_width[gust_width_sort_ind].reshape(9,9).T,Mg[gust_width_sort_ind][Mg_sort_ind].reshape(9,9),d_oaspl[gust_width_sort_ind][Mg_sort_ind].reshape(9,9),levels = levels,cmap = plt.cm.inferno)

dist2 = ax.contour(gust_width[Mg_sort_ind].reshape(9,9),Mg[Mg_sort_ind].reshape(9,9),d_oaspl[Mg_sort_ind].reshape(9,9),levels = levels,colors = 'k')
plt.clabel(dist2,levels=levels_c)
cbar = fig.colorbar(dist,pad = .05)
cbar.ax.set_ylabel(r'$\Delta \ OASPL_{max}, \ dB \ (re: \ 20 \mu Pa)$')
# cbar.ax.set_ylabel(r'$OASPL_{max}, \ dB \ (re: \ 20 \mu Pa)$')

cbar.set_ticks(levels[::2])
# ax.set_xticks(np.arange(1,5))
ax.set_xlabel(r'$x/c$')
ax.set_ylabel(r'$M_g$')
ax.set_ylim([0,.3])
ax.scatter(gust_width[case_ind_select],Mg[case_ind_select])
# ax.scatter(gust_width[d_oaspl.argmin()],Mg[d_oaspl.argmin()],marker = '^')
# ax.scatter(gust_width[d_oaspl.argmax()],Mg[d_oaspl.argmax()],marker = '*')

plt.savefig(os.path.join(os.path.dirname(cases_dir),f'rg_Mg_d_oaspl_carpet.png'),format = 'png')
if eps:
    plt.savefig(os.path.join(os.path.dirname(cases_dir),f'rg_Mg_d_oaspl_carpet.eps'),format = 'eps')
plt.close()
#%%


fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# plt.subplots_adjust(left = .2,top = .925,right = .85,bottom = .13)
levels = np.linspace(-10,10,41)
dist = ax.contourf(gust_width[Mg_sort_ind].reshape(9,9),Mg[Mg_sort_ind].reshape(9,9),dtheta_max[Mg_sort_ind].reshape(9,9),levels = levels,cmap = plt.cm.Spectral.reversed(),norm = mcolors.CenteredNorm())

# levels = np.linspace(95,100,21)
# dist = ax.contourf(gust_width[Mg_sort_ind].reshape(9,9),Mg[Mg_sort_ind].reshape(9,9),d_oaspl[Mg_sort_ind].reshape(9,9),levels = levels,cmap = plt.cm.inferno)

# dist2 = ax.contour(gust_width[gust_width_sort_ind].reshape(9,9),Mg[gust_width_sort_ind][Mg_sort_ind].reshape(9,9),d_oaspl[gust_width_sort_ind][Mg_sort_ind].reshape(9,9),levels = levels_c,colors = 'k')
# plt.clabel(dist2,levels=levels_c)
cbar = fig.colorbar(dist,pad = .05)
cbar.ax.set_ylabel(r'$\Delta \ \theta \ [deg]$')
cbar.set_ticks(levels[::2])
# ax.set_xticks(np.arange(1,5))
ax.set_xlabel(r'$x/c$')
ax.set_ylabel(r'$M_g$')
ax.set_ylim([0,.3])
ax.scatter(gust_width[case_ind_select],Mg[case_ind_select])

plt.savefig(os.path.join(os.path.dirname(cases_dir),f'rg_Mg_d_theta_carpet.png'),format = 'png')
if eps:
    plt.savefig(os.path.join(os.path.dirname(cases_dir),f'rg_Mg_d_theta_carpet.eps'),format = 'eps')
plt.close()


#%%

# Plot profile at selected gust width and strength 

gust_width_select = np.unique(gust_width)[6]
gust_width_select_ind = np.where(gust_width==gust_width_select)[0]
gust_width_select_ind = gust_width_select_ind[np.argsort(Mg[gust_width_select_ind])]

acs_data = import_results_from_wopwop(os.path.join(cases_dir,cases[gust_width_select_ind],'acoustics'))

fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# ax.plot(h[gust_width_select_ind,:,r_ind_select].T,v_gust[gust_width_select_ind,:,r_ind_select].T/saved_params[case]['sos'])
ax.plot(acs_data['function_values'][...,0].squeeze().T,acs_data['function_values'][...,-1].squeeze().T)
ax.set_xlabel(r'$t \ [sec]$')
ax.set_ylabel(r'$P \ [Pa]$')
# ax.set_ylim([0,.3])
# ax.set_xlim([85,105])
ax.grid()

#%%

gust_width_select = np.unique(gust_width)[4]
gust_width_select_ind = np.where(gust_width==gust_width_select)[0]
gust_width_select_ind = gust_width_select_ind[np.argsort(Mg[gust_width_select_ind])]
case_ind_select = gust_width_select_ind[np.array([1,3,-2])]

default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


leglab = []
fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# ax.plot(h[gust_width_select_ind,:,r_ind_select].T,v_gust[gust_width_select_ind,:,r_ind_select].T/saved_params[case]['sos'])
for i,ind in enumerate(case_ind_select):
    ax.plot((saved_params_filt[cases[ind]]['psi'][-int(2*np.pi/saved_params_filt[cases[ind]]['dpsi']):]*180/np.pi)%360,np.gradient(saved_params_filt[cases[ind]]['loads'][-int(2*np.pi/saved_params_filt[cases[ind]]['dpsi']):,r_ind_select,-1],edge_order=2)/np.diff(saved_params_filt[cases[ind]]['psi'][:2])[0],linestyle = '-',c = default_colors[i],label=f"$w_g = {np.round(gust_width[ind],2)}c, M_g = {np.round(Mg[ind],2)}c$")
    ax.plot((saved_params_filt[cases[ind]]['psi'][-int(2*np.pi/saved_params_filt[cases[ind]]['dpsi']):]*180/np.pi)%360,np.gradient(saved_params_filt[cases[ind]]['filt_loads'][-int(2*np.pi/saved_params_filt[cases[ind]]['dpsi']):,r_ind_select,-1],edge_order=2)/np.diff(saved_params_filt[cases[ind]]['psi'][:2])[0],linestyle = '-.',c = default_colors[i])
    leglab.append(f"$w_g = {np.round(gust_width[ind],2)}c, M_g = {np.round(Mg[ind],2)}c$")

ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel(r'$\partial dFz /\partial \psi\ [N/rad]$')
# ax.set_ylim([0,.3])
ax.set_xlim([85,115])
ax.grid()
ax.legend()


Mg_range = [0.15,0.18]
Mg_select_ind = np.where((Mg>Mg_range[0]) & (Mg<Mg_range[-1]))[0]
Mg_select_ind = Mg_select_ind[np.argsort(gust_width[Mg_select_ind])]
case_ind_select = np.concatenate((case_ind_select, Mg_select_ind[np.array([1,3,-2])]))

leglab = []
fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# ax.plot(h[gust_width_select_ind,:,r_ind_select].T,v_gust[gust_width_select_ind,:,r_ind_select].T/saved_params[case]['sos'])
for i,ind in enumerate(case_ind_select):
    ax.plot((saved_params_filt[cases[ind]]['psi'][-int(2*np.pi/saved_params_filt[cases[ind]]['dpsi']):]*180/np.pi)%360,np.gradient(saved_params_filt[cases[ind]]['loads'][-int(2*np.pi/saved_params_filt[cases[ind]]['dpsi']):,r_ind_select,-1],edge_order=2)/np.diff(saved_params_filt[cases[ind]]['psi'][:2])[0],linestyle = '-',c = default_colors[i],label=f"$w_g = {np.round(gust_width[ind],2)}c, M_g = {np.round(Mg[ind],2)}c$")
    ax.plot((saved_params_filt[cases[ind]]['psi'][-int(2*np.pi/saved_params_filt[cases[ind]]['dpsi']):]*180/np.pi)%360,np.gradient(saved_params_filt[cases[ind]]['filt_loads'][-int(2*np.pi/saved_params_filt[cases[ind]]['dpsi']):,r_ind_select,-1],edge_order=2)/np.diff(saved_params_filt[cases[ind]]['psi'][:2])[0],linestyle = '-.',c = default_colors[i])
    leglab.append(f"$w_g = {np.round(gust_width[ind],2)}c, M_g = {np.round(Mg[ind],2)}c$")

ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel(r'$\partial dFz /\partial \psi\ [N/rad]$')
# ax.set_ylim([0,.3])
ax.set_xlim([85,105])
ax.grid()
ax.legend()










fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# ax.plot(h[gust_width_select_ind,:,r_ind_select].T,v_gust[gust_width_select_ind,:,r_ind_select].T/saved_params[case]['sos'])
for ind in gust_width_select_ind:
    ax.plot(h[ind][:,r_ind_select].T,v_gust[ind][:,r_ind_select].T/saved_params_bl[case]['sos'])
ax.set_xlabel(r'$x/c$')
ax.set_ylabel(r'$M_g$')
ax.set_ylim([0,.3])
ax.set_xlim([0,1])
ax.grid()

Mg_range = [0.14,0.16]
Mg_select_ind = np.where((Mg>Mg_range[0]) & (Mg<Mg_range[-1]))[0]
Mg_select_ind = Mg_select_ind[np.argsort(gust_width[Mg_select_ind])]
leglab  = []
fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# ax.plot(h[Mg_select_ind,:,r_ind_select].T,v_gust[Mg_select_ind,:,r_ind_select].T/saved_params[case]['sos'])

for ind in Mg_select_ind:
    ax.plot(h[ind][:,r_ind_select].T,v_gust[ind][:,r_ind_select].T/saved_params_bl[case]['sos'])
    leglab.append(f"$w_g = {np.round(gust_width[ind],2)}c$")

ax.set_xlabel(r'$x/c$')
ax.set_ylabel(r'$M_g$')
ax.set_ylim([0,.16])
ax.set_xlim([0,1])
ax.grid()
ax.legend(leglab)

fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# ax.plot(h[gust_width_select_ind,:,r_ind_select].T,v_gust[gust_width_select_ind,:,r_ind_select].T/saved_params[case]['sos'])
for ind in Mg_select_ind:
    ax.plot((saved_params_bl[cases[ind]]['psi'][-int(2*np.pi/saved_params_bl[cases[ind]]['dpsi']):]*180/np.pi)%360,saved_params_bl[cases[ind]]['loads'][-int(2*np.pi/saved_params_bl[cases[ind]]['dpsi']):,r_ind_select,-1])
ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel(r'$Fz \ [N]$')
ax.set_xlim([85,105])
# ax.set_xlim([0,1])
ax.grid()
ax.legend(leglab)


fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# ax.plot(h[gust_width_select_ind,:,r_ind_select].T,v_gust[gust_width_select_ind,:,r_ind_select].T/saved_params[case]['sos'])
for ind in Mg_select_ind:
    ax.plot((saved_params_bl[cases[ind]]['psi'][-int(2*np.pi/saved_params_bl[cases[ind]]['dpsi']):]*180/np.pi)%360,np.gradient(saved_params_bl[cases[ind]]['loads'][-int(2*np.pi/saved_params_bl[cases[ind]]['dpsi']):,r_ind_select,-1],edge_order=2)/np.diff(saved_params_bl[cases[ind]]['psi'][:2])[0])
ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel(r'$\partial dFz /\partial \psi\ [N/rad]$')
# ax.set_ylim([0,.3])
ax.set_xlim([85,105])
ax.grid()
ax.legend(leglab)



#%%

gust_width_select = np.unique(gust_width)[6]
gust_width_select_ind = np.where(gust_width==gust_width_select)[0]
gust_width_select_ind = gust_width_select_ind[np.argsort(Mg[gust_width_select_ind])]
case_ind_select = gust_width_select_ind[np.array([1,3,-2])]