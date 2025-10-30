import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'src'))
from help_funcs import *
import matplotlib.colors as mcolors
from scipy.interpolate import CubicSpline


#%%

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ["Times New Roman"]
plt.rcParams['font.size'] = 16

#

cases_directory ='/Users/danielweitsman/codes/github/DanWeitsman/unsteady_BEMT/cases/final_designs/sweeps/rg_Mg'
cases = ['baseline_untapered_AR12_select/case_66','sdof_geom_untapered_AR12_select/case_66']

oaspl = {}
for case in cases:

    # saved_params = read_results_from_h5(os.path.join(cases_directory,case))
    acs_data = import_results_from_wopwop(os.path.join(cases_directory,case,'acoustics'))
    oaspl.update({case:np.round(20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6),1)})

theta = (np.arctan2(acs_data['geometry_values'][:,:,0,1],acs_data['geometry_values'][:,:,0,0])*180/np.pi)%360
phi = np.arctan2(acs_data['geometry_values'][:,:,0,-1],np.linalg.norm((acs_data['geometry_values'][:,:,0,0],acs_data['geometry_values'][:,:,0,1]),axis = 0))*180/np.pi

cmap = plt.cm.Spectral.reversed()
# levels = np.linspace(np.round(oaspl.min()-6),np.round(oaspl.max()+6),50)
d_oaspl = oaspl[cases[1]]-oaspl[cases[0]]

N = 500
selected_ind = 7
theta_interp = np.arange(N)*(theta[-1,selected_ind]-theta[0,selected_ind])/(N-1)+theta[0,selected_ind]

doaspl_interp = [CubicSpline(x = theta[:,selected_ind].squeeze(),y = oaspl[case][:,selected_ind].squeeze()) for case in cases]
max_oaspl_ind = [doaspl_interp[i](theta_interp).argmax() for i in range(len(cases))]

fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
plt.subplots_adjust(bottom = .15)
for i in range(len(cases)):
    ax.plot(theta_interp,doaspl_interp[i](theta_interp))
ax.legend(['Baseline','Treated'])
for i in range(len(cases)):
    ax.scatter(theta_interp[max_oaspl_ind[i]],doaspl_interp[i](theta_interp)[max_oaspl_ind[i]])
ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel(r'$OASPL, \ dB \ (re: \ 20 \mu Pa)$')

# levels = np.round(np.linspace(np.round(d_oaspl.min(),1),np.round(d_oaspl.max(),1),30),1)

levels = np.linspace(-3,3,13)
# levels = np.linspace(0,np.round(d_oaspl.max()*.7/2)*2,int(np.round(d_oaspl.max()*.7/2)+1))
# levels_c = np.linspace(0,np.round(d_oaspl.max()*.5/2)*2,int(np.round(d_oaspl.max()*.5/2)+1))

# cbar_ticks = levels[::10]

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(right = .85,left = .2,bottom = .15)
dist = ax.contourf(theta,phi,d_oaspl,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
dist2 = ax.contour(theta,phi,d_oaspl,levels = levels[::2],colors = 'k')
plt.clabel(dist2,levels=levels[::2])
cbar = fig.colorbar(dist,pad = .05)
cbar.ax.set_ylabel(R'$\Delta \ BVISPL, \ dB \ (re: \ 20 \ \mu Pa)$')
cbar.set_ticks(levels[::2])
ax.set_xticks(np.round(theta[::10,0]))
ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel(r'$\phi \ [deg]$')
plt.savefig(os.path.join(cases_directory,f'd_oaspl_carpet_{cases[-1]}.png'),format = 'png')

phi_ind =  7
theta_ind = 15

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .2,bottom = .15)
ax.plot(acs_data['function_values'][theta_ind,phi_ind,:,0],acs_data['function_values'][theta_ind,phi_ind,:,-1])
ax.set_xlabel('$Time [s]$')
ax.set_ylabel('$Pressure [Pa]$')
ax.grid()
# plt.savefig(os.path.join(cases_directory,f'p_tseries.png'),format = 'png')
