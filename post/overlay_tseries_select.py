import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
import f90nml
from scipy.signal import welch
# sys.path.insert(0,os.path.join(os.path.dirname(__file__),'dependencies','pyWopwop'))
# import wopwop
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'src'))
from help_funcs import *

#%%
import matplotlib.colors as mcolors

fontName = 'Times New Roman'
fontSize = 12
plt.rc('font',**{'family':'serif','serif':[fontName],'size':fontSize})
plt.rc('mathtext',**{'default':'regular'})
plt.rc('text',**{'usetex':False})
plt.rc('lines',**{'linewidth':2})

#%%

cases_directory ='/Users/danielweitsman/codes/github/DanWeitsman/unsteady_BEMT/cases/validation/'
cases = ['unsteady_cfd','unsteady_cfd_2']

acs_data = {}
case_data = {}
oaspl = {}
for case in cases:
    acs_data.update({case:import_results_from_wopwop(os.path.join(cases_directory,case,'acoustics'))})
    case_data.update({case:read_results_from_h5(os.path.join(cases_directory,case))})
    oaspl.update({case:np.round(20*np.log10(np.sqrt(np.mean(acs_data[case]['function_values'][:,:,:,-1]**2,axis = -1))/20e-6),1)})

# print(np.abs(oaspl[cases[1]]-oaspl[cases[0]]).max())

theta = np.round(np.arctan2(acs_data[cases[0]]['geometry_values'][:,:,0,1],acs_data[cases[0]]['geometry_values'][:,:,0,0])%(2*np.pi)*180/np.pi,1).squeeze()
phi = np.round(np.arctan2(acs_data[cases[0]]['geometry_values'][:,:,0,-1],np.linalg.norm((acs_data[cases[0]]['geometry_values'][:,:,0,0],acs_data[cases[0]]['geometry_values'][:,:,0,1]),axis = 0))*180/np.pi,1).squeeze()

psi = (((acs_data[cases[0]]['function_values'][0,0,:,0]*(case_data[cases[0]]['omega'])))*180/np.pi)
dpsi = psi[1]-psi[0]

N_obs_theta = acs_data[cases[0]]['function_values'].shape[0]
N_obs_phi = acs_data[cases[0]]['function_values'].shape[1]

phi_select_ind = 0
skip_ind = 1
#%%
for theta_itr in range(int(N_obs_theta/skip_ind)):
        
        fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
        plt.subplots_adjust(left = .2,bottom = .15)
        # ax.plot(acs_data[cases[0]]['function_values'][theta_itr,phi_itr,:,0],acs_data[cases[0]]['function_values'][theta_itr,phi_itr,:,-1])
        # ax.plot(acs_data[cases[1]]['function_values'][theta_itr,phi_itr,:,0],acs_data[cases[1]]['function_values'][theta_itr,phi_itr,:,-1])
        for case in cases:
            ax.plot(psi,acs_data[case]['function_values'][::skip_ind][theta_itr,phi_select_ind,:,-1])
        
        peak_ind =np.abs(acs_data[cases[0]]['function_values'][::skip_ind][theta_itr,phi_select_ind,:,-1]).argmax()
        # ax.plot(psi,acs_data[cases[1]]['function_values'][::skip_ind][theta_itr,phi_select_ind,:,-1])

        ax.set_xlabel(r'$\psi \ [deg]$')
        ax.set_ylabel(r'$Pressure, \ [Pa]$')
        # ax.set_title(rf'$\psi = {theta[::skip_ind][theta_itr]}^\circ$, $\phi = {phi[phi_select_ind]}^\circ$')
        # ax.set_ylim(bottom = 0)
        # ax.set_xlim([0,6e3])
        plt.legend(['Baseline','Treated'])

        ax.set_xlim([psi[int(peak_ind-35/dpsi)],psi[int(peak_ind+35/dpsi)]])
        # ax.set_ylim([-200,100])
        ax.grid()
        # plt.savefig(os.path.join(cases_directory,f'tseries_{theta_itr+phi_itr}.png'),format = 'png')

r_ind_select = np.abs(0.9-case_data[cases[1]]['r']).argmin()
fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# ax.plot(h[gust_width_select_ind,:,r_ind_select].T,v_gust[gust_width_select_ind,:,r_ind_select].T/saved_params[case]['sos'])
ax.plot((case_data[cases[1]]['psi'][-int(2*np.pi/case_data[cases[1]]['dpsi']):]*180/np.pi)%360,case_data[cases[1]]['loads'][-int(2*np.pi/case_data[cases[1]]['dpsi']):,r_ind_select,-1],linestyle = '-')
ax.plot((case_data[cases[1]]['psi'][-int(2*np.pi/case_data[cases[1]]['dpsi']):]*180/np.pi)%360,case_data[cases[1]]['filt_loads'][-int(2*np.pi/case_data[cases[1]]['dpsi']):,r_ind_select,-1],linestyle = '-.')

ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel(r'$dFz [N]$')
# ax.set_ylim([0,.3])
ax.set_xlim([60,150])
ax.grid()
ax.legend(['Baseline','Filtered'])

r_ind_select = np.abs(0.9-case_data[cases[1]]['r']).argmin()
fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
# ax.plot(h[gust_width_select_ind,:,r_ind_select].T,v_gust[gust_width_select_ind,:,r_ind_select].T/saved_params[case]['sos'])
ax.plot((case_data[cases[1]]['psi'][-int(2*np.pi/case_data[cases[1]]['dpsi']):]*180/np.pi)%360,np.gradient(case_data[cases[1]]['loads'][-int(2*np.pi/case_data[cases[1]]['dpsi']):,r_ind_select,-1],edge_order=2)/np.diff(case_data[cases[1]]['psi'][:2])[0],linestyle = '-')
ax.plot((case_data[cases[1]]['psi'][-int(2*np.pi/case_data[cases[1]]['dpsi']):]*180/np.pi)%360,np.gradient(case_data[cases[1]]['filt_loads'][-int(2*np.pi/case_data[cases[1]]['dpsi']):,r_ind_select,-1],edge_order=2)/np.diff(case_data[cases[1]]['psi'][:2])[0],linestyle = '-.')

ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel(r'$\partial dFz /\partial \psi\ [N/rad]$')
# ax.set_ylim([0,.3])
ax.set_xlim([60,150])
ax.grid()
ax.legend(['Baseline','Filtered'])
