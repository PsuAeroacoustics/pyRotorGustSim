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

#

cases_directory ='/Users/danielweitsman/codes/github/DanWeitsman/unsteady_BEMT/cases/final_designs/sweeps/rg_Mg/sdof_geom_untapered_AR10_select_OAR6'
case = 'case_12'

acs_data = {}
saved_params = {}

acs_data = import_results_from_wopwop(os.path.join(cases_directory,case,'acoustics'))


oaspl = 20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6)
theta = np.arctan2(acs_data['geometry_values'][:,:,0,1],acs_data['geometry_values'][:,:,0,0])%(2*np.pi)*180/np.pi
phi = np.arctan2(acs_data['geometry_values'][:,:,0,-1],np.linalg.norm((acs_data['geometry_values'][:,:,0,0],acs_data['geometry_values'][:,:,0,1]),axis = 0))*180/np.pi

np.round(oaspl[:,np.abs(phi[0]-30).argmin()],1).max()
# levels = np.linspace(np.round(oaspl.min()-6),np.round(oaspl.max()+6),50)
levels = np.linspace(85,115,41)
levels_c = np.linspace(85,115,21)

cbar_ticks = np.round(levels)[::8]

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .2,right = .9,bottom = .15)
dist = ax.contourf(theta,phi,oaspl,levels = levels,cmap = 'inferno')
dist2 = ax.contour(theta,phi,oaspl,levels = np.round(levels_c)[::2],colors = 'k')
plt.clabel(dist2,levels=levels_c[::2])
cbar = fig.colorbar(dist,pad = .05)
cbar.ax.set_ylabel(R'$OASPL, \ dB \ (re: \ 20 \ \mu Pa)$')
cbar.set_ticks(cbar_ticks)
ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel(r'$\phi \ [deg]$')
ax.set_xticks(np.round(theta[::10,0]))
plt.savefig(os.path.join(cases_directory,f'oaspl_carpet_{os.path.basename(case)}.png'),format = 'png')

# phi_ind =  12
# theta_ind = 15

# fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
# plt.subplots_adjust(left = .2,bottom = .15)
# ax.plot(acs_data['function_values'][theta_ind,phi_ind,:,0],acs_data['function_values'][theta_ind,phi_ind,:,-1])
# ax.set_xlabel('$Time [s]$')
# ax.set_ylabel('$Pressure [Pa]$')
# ax.grid()
# plt.savefig(os.path.join(cases_directory,f'p_tseries.png'),format = 'png')
