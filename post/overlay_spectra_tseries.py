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

cases_directory ='/Users/danielweitsman/codes/github/DanWeitsman/unsteady_BEMT/cases/final_designs'
cases = ['baseline','sdof_geom']

acs_data = {}
case_data = {}
oaspl = {}
for case in cases:
    acs_data.update({case:import_results_from_wopwop(os.path.join(cases_directory,case,'acoustics'))})
    case_data.update({case:read_results_from_h5(os.path.join(cases_directory,case))})
    oaspl.update({case:np.round(20*np.log10(np.sqrt(np.mean(acs_data[case]['function_values'][:,:,:,-1]**2,axis = -1))/20e-6),1)})

theta = np.round(np.arctan2(acs_data[cases[0]]['geometry_values'][:,:,0,1],acs_data[cases[0]]['geometry_values'][:,:,0,0])%(2*np.pi)*180/np.pi,1).squeeze()
phi = np.round(np.arctan2(acs_data[cases[0]]['geometry_values'][:,:,0,-1],np.linalg.norm((acs_data[cases[0]]['geometry_values'][:,:,0,0],acs_data[cases[0]]['geometry_values'][:,:,0,1]),axis = 0))%(2*np.pi)*180/np.pi,1).squeeze()

psi = (((acs_data[cases[0]]['function_values'][0,0,:,0]*(case_data['baseline']['omega'])))*180/np.pi)

N_obs_theta = acs_data[cases[0]]['function_values'].shape[0]
N_obs_phi = acs_data[cases[0]]['function_values'].shape[1]

nperseg = acs_data[cases[0]]['function_values'].shape[-2]
dt = np.diff(acs_data[cases[0]]['function_values'][0,0,:2,0])[0]
fs = dt**-1
df = (nperseg*dt)**-1
f = np.arange(int(nperseg/2)+1)*df

pxx = 10*np.log10(np.asarray([(welch(acs_data[case]['function_values'][:,:,:,-1], fs=fs, window='hann', nperseg=nperseg, noverlap=None, nfft=None, detrend='constant', return_onesided=True, scaling='density', axis=-1, average='mean')[-1]) for case in cases])/20e-6**2)

for theta_itr in range(N_obs_theta):
    for phi_itr in range(N_obs_phi):
        
        fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
        plt.subplots_adjust(left = .2,bottom = .15)
        # ax.plot(acs_data[cases[0]]['function_values'][theta_itr,phi_itr,:,0],acs_data[cases[0]]['function_values'][theta_itr,phi_itr,:,-1])
        # ax.plot(acs_data[cases[1]]['function_values'][theta_itr,phi_itr,:,0],acs_data[cases[1]]['function_values'][theta_itr,phi_itr,:,-1])
        ax.plot(psi,acs_data[cases[0]]['function_values'][theta_itr,phi_itr,:,-1])
        ax.plot(psi,acs_data[cases[1]]['function_values'][theta_itr,phi_itr,:,-1])

        ax.set_xlabel(r'$\psi \ [deg]$')
        ax.set_ylabel(r'$Pressure, \ [Pa]$')
        ax.set_title(rf'$\psi = {theta[theta_itr]}^\circ$, $\phi = {phi[phi_itr]}^\circ$')
        # ax.set_ylim(bottom = 0)
        # ax.set_xlim([0,6e3])
        plt.legend(['Baseline','Treated'])
        ax.grid()
        plt.savefig(os.path.join(cases_directory,f'tseries_{theta_itr+phi_itr}.png'),format = 'png')
        plt.close()

for theta_itr in range(N_obs_theta):
    for phi_itr in range(N_obs_phi):
        
        fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
        plt.subplots_adjust(left = .2,bottom = .15)
        ax.plot(f,pxx[0,theta_itr,phi_itr], marker = 'o')
        ax.plot(f,pxx[1,theta_itr,phi_itr], marker = 'o')
        ax.set_xlabel(r'$Frequency \ [Hz]$')
        ax.set_ylabel(r'$PSD, \ dB/Hz \ (re: \ 20 \mu Pa)$')
        ax.set_title(rf'$\psi = {theta[theta_itr]}^\circ$, $\phi = {phi[phi_itr]}^\circ$')
        ax.set_ylim(bottom = 0)
        ax.set_xlim([0,6e3])
        plt.legend(['Baseline','Treated'])
        ax.grid()
        plt.savefig(os.path.join(cases_directory,f'spectra_{theta_itr+phi_itr}.png'),format = 'png')
        plt.close()


