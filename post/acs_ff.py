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

cases_directory = os.getcwd()
case = 'unsteady_loading_ff'

acs_data = {}
saved_params = {}

acs_data = import_results_from_wopwop(os.path.join(cases_directory,case,'acoustics'))


oaspl = 20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6)
theta = np.arctan2(acs_data['geometry_values'][:,:,0,1],acs_data['geometry_values'][:,:,0,0])%(2*np.pi)*180/np.pi
phi = np.arctan2(acs_data['geometry_values'][:,:,0,-1],np.linalg.norm((acs_data['geometry_values'][:,:,0,0],acs_data['geometry_values'][:,:,0,1]),axis = 0))*180/np.pi
r = np.linalg.norm(acs_data['geometry_values'][:,:,0],axis = -1)

sphere_spread =  20*np.log10(1/r)

W = 10*np.log10(4*np.pi*r[-1]**2*10**(oaspl[-1]/10))
sphere_spread_2 = W-10*np.log10(4*np.pi*r**2)
ff_ind = np.where(np.abs((oaspl-sphere_spread_2)).squeeze()<=1)[0][0]


fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .2)
ax.plot(r.squeeze()/0.19685,oaspl.squeeze())
ax.plot(r.squeeze()/0.19685,sphere_spread_2)
ax.scatter(r[ff_ind]/0.19685,sphere_spread_2[ff_ind])
# ax.plot(r.squeeze()/0.19685,sphere_spread+(oaspl[-1]-sphere_spread[-1]))
ax.set_xscale('log')
ax.set_ylabel(R'$OASPL, \ dB \ (re: \ 20 \ \mu Pa)$')
ax.set_xlabel('$r/R$')
ax.legend(['Predicted','Spherical Spreading'])
ax.set_xlim([1,40])
ax.set_ylim([40,110])
ax.grid()
