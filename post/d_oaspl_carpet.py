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

cases_directory =os.getcwd()
cases = ['case_0','case_2']

oaspl = {}
for case in cases:

    # saved_params = read_results_from_h5(os.path.join(cases_directory,case))
    acs_data = import_results_from_wopwop(os.path.join(cases_directory,case,'acoustics'))
    oaspl.update({case:20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6)})

theta = np.arctan2(acs_data['geometry_values'][:,:,0,1],acs_data['geometry_values'][:,:,0,0])%(2*np.pi)*180/np.pi
phi = np.arctan2(acs_data['geometry_values'][:,:,0,-1],np.linalg.norm((acs_data['geometry_values'][:,:,0,0],acs_data['geometry_values'][:,:,0,1]),axis = 0))*180/np.pi

cmap = plt.cm.Spectral.reversed()
# levels = np.linspace(np.round(oaspl.min()-6),np.round(oaspl.max()+6),50)
d_oaspl = oaspl[cases[1]]-oaspl[cases[0]]
# levels = np.round(np.linspace(np.round(d_oaspl.min(),1),np.round(d_oaspl.max(),1),30),1)
levels = np.linspace(0,np.round(d_oaspl.max()*.7/2)*2,int(np.round(d_oaspl.max()*.7/2)+1))
levels_c = np.linspace(0,np.round(d_oaspl.max()*.5/2)*2,int(np.round(d_oaspl.max()*.5/2)+1))

# cbar_ticks = levels[::10]

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(right = .85,left = .2,bottom = .15)
dist = ax.contourf(theta,phi,d_oaspl,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
dist2 = ax.contour(theta,phi,d_oaspl,levels = levels_c,colors = 'k')
plt.clabel(dist2,levels=levels_c)
cbar = fig.colorbar(dist,pad = .05)
cbar.ax.set_ylabel(R'$\Delta \ BVISPL, \ dB \ (re: \ 20 \ \mu Pa)$')
cbar.set_ticks(levels[::4])
ax.set_xticks(np.round(theta[::10,0]))
ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel(r'$\phi \ [deg]$')
plt.savefig(os.path.join(cases_directory,f'd_oaspl_carpet_{cases[-1]}.png'),format = 'png')

# phi_ind =  12
# theta_ind = 15

# fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
# plt.subplots_adjust(left = .2,bottom = .15)
# ax.plot(acs_data['function_values'][theta_ind,phi_ind,:,0],acs_data['function_values'][theta_ind,phi_ind,:,-1])
# ax.set_xlabel('$Time [s]$')
# ax.set_ylabel('$Pressure [Pa]$')
# ax.grid()
# plt.savefig(os.path.join(cases_directory,f'p_tseries.png'),format = 'png')
