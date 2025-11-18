import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
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

#
#%%

cases_directory = os.getcwd()
cases = ['mdof_geom_parallel_compact','mdof_geom_parallel_noncompact']

saved_params = {}
for i,case in enumerate(cases):
    saved_params.update({case:read_results_from_h5(os.path.join(cases_directory,case))})


fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
plt.subplots_adjust(left = 0.15,bottom = .15)
ax.plot(saved_params[cases[0]]['psi']*180/np.pi,saved_params[cases[0]]['loads'][:,int(0.75*saved_params[cases[0]]['loads'].shape[1]),-1])
for k,v in saved_params.items():
    ax.plot(saved_params[k]['psi']*180/np.pi,saved_params[k]['filt_loads'][:,int(0.75*saved_params[k]['filt_loads'].shape[1]),-1])
ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel('$F_z \ [N]$')
ax.legend(['Unfiltered','Compact','Non-compact'])
ax.set_xlim(0,360)
ax.grid()

fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
plt.subplots_adjust(left = 0.15,bottom = .15)
ax.plot(saved_params[cases[0]]['psi']*180/np.pi,np.gradient(saved_params[cases[0]]['loads'][:,int(0.75*saved_params[cases[0]]['loads'].shape[1]),-1]))
for k,v in saved_params.items():
    ax.plot(saved_params[k]['psi']*180/np.pi,np.gradient(saved_params[k]['filt_loads'][:,int(0.75*saved_params[k]['filt_loads'].shape[1]),-1]))
ax.set_xlabel(r'$\psi \ [deg]$')
ax.set_ylabel('$\partial F_z/\partial t \ [N/deg]$')
ax.legend(['Unfiltered','Compact','Non-compact'])
ax.set_xlim(0,360)
ax.grid()
