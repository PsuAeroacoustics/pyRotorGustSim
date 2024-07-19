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
fontSize = 16
plt.rc('font',**{'family':'serif','serif':[fontName],'size':fontSize})
plt.rc('mathtext',**{'default':'regular'})
plt.rc('text',**{'usetex':False})
plt.rc('lines',**{'linewidth':2})

#%%

cases_directory =os.getcwd()
cases = ['unsteady_parallel_lgrid','mdof_geom_lgrid']

#%%
acs_data = {}
oaspl = {}
for case in cases:
    acs_data.update({case:import_results_from_wopwop(os.path.join(cases_directory,case,'acoustics'))})
    oaspl.update({case:20*np.log10(np.sqrt(np.mean(acs_data[case]['function_values'][:,:,:,-1]**2,axis = -1))/20e-6)})

theta = np.round(np.arctan2(acs_data[cases[0]]['geometry_values'][:,:,0,1],acs_data[cases[0]]['geometry_values'][:,:,0,0])%(2*np.pi)*180/np.pi,1).squeeze()
phi = np.round(np.arctan2(acs_data[cases[0]]['geometry_values'][:,:,0,-1],np.linalg.norm((acs_data[cases[0]]['geometry_values'][:,:,0,0],acs_data[cases[0]]['geometry_values'][:,:,0,1]),axis = 0))%(2*np.pi)*180/np.pi,1).squeeze()

N_obs_theta = acs_data[cases[0]]['function_values'].shape[0]
N_obs_phi = acs_data[cases[0]]['function_values'].shape[1]

nperseg = acs_data[cases[0]]['function_values'].shape[-2]
dt = acs_data[cases[0]]['function_values'][0,0,1,0]-acs_data[cases[0]]['function_values'][0,0,0,0]
fs = dt**-1
df = (nperseg*dt)**-1
f = np.arange(nperseg)*df
bpf = np.arange(np.round(nperseg/2+1))

pxx = []
for case in cases:
    pxx.append(welch(acs_data[case]['function_values'][:,:,:,-1], fs=fs, window='hann', nperseg=nperseg, noverlap=None, nfft=None, detrend='constant', return_onesided=True, scaling='density', axis=-1, average='mean')[-1])
pxx = 10*np.log10(np.array(pxx)*df/20e-6**2)

for theta_itr in range(N_obs_theta):
    for phi_itr in range(N_obs_phi):
        
        fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
        plt.subplots_adjust(left = .2,bottom = .15)
        ax.plot(bpf,pxx[:,theta_itr,phi_itr].T, marker = 'o')
        ax.set_xlabel(r'$BPF \ Harmonic$')
        ax.set_ylabel(r'$SPL, \ dB \ (re: \ 20 \mu Pa)$')
        ax.set_title(rf'$\psi = {theta[theta_itr]}^\circ$, $\phi = {phi[phi_itr]}^\circ$')
        ax.set_ylim(bottom = 0)
        ax.set_xlim([0,100])
        plt.legend(['Baseline','Treated'])
        ax.grid()
        plt.savefig(os.path.join(cases_directory,f'spectra_{theta_itr+phi_itr}.png'),format = 'png')
        plt.close()

