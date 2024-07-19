import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'src'))
from help_funcs import *


#%%

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ["Times New Roman"]
plt.rcParams['font.size'] = 16

#%%

cases_directory = os.getcwd()
cases = [x for x  in os.listdir(cases_directory) if os.path.isdir(os.path.join(cases_directory,x))]

acs_data = {}
saved_params = {}

for case in cases:
    saved_params.update({case:read_results_from_h5(os.path.join(cases_directory,case))})

oaspl = np.zeros((len(cases),len(saved_params[cases[0]]['geometry_values'])))
for i,case in enumerate(cases):
    saved_params[case]['geometry_values'] = np.flip(saved_params[case]['geometry_values'],axis = 0).squeeze()
    saved_params[case]['function_values'] = np.flip(saved_params[case]['function_values'],axis = 0).squeeze()
    oaspl[i] = 20*np.log10(np.sqrt(np.mean(saved_params[case]['function_values'][:,:,-1]**2,axis = -1))/20e-6)
theta = np.round(np.arctan2(saved_params[cases[0]]['geometry_values'][:,0,1],saved_params[cases[0]]['geometry_values'][:,0,0])*180/np.pi)%(360)

leg_labs = [f'$\\theta = {theta[x]}^\circ$' for x in np.arange(oaspl.shape[-1])][::-1]

r_ind = -1
gust_width = np.zeros(len(cases))
for i,case in enumerate(cases):
    max_ind = saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,r_ind].argmax()
    start_ind = saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,r_ind][:max_ind+1].argmin()
    end_ind = np.abs(saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,r_ind][max_ind:]-.1*saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,r_ind][max_ind]).argmin()
    gust_width[i] = (max_ind-start_ind+end_ind)*saved_params[case]['dpsi']*saved_params[case]['r'][r_ind]*saved_params[case]['R']/saved_params[case]['c']

plot_ind = gust_width.argsort()
baseline_ind = 2

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .2,bottom = .15)
for case in cases:
    ax.plot(saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,-1])
# ax.plot(saved_params[case]['v_gust'][-int((2*np.pi)/saved_params[case]['dpsi']):,0])
plt.savefig(os.path.join(cases_directory,f'gust.png'),format = 'png')

for mic_iter in range(oaspl.shape[-1]):

    fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
    plt.subplots_adjust(left = .2,bottom = .15)
    ax.plot(gust_width[plot_ind],oaspl[plot_ind,mic_iter]-oaspl[plot_ind[baseline_ind],mic_iter],marker = 'o')
    ax.set_xlabel(r'$r_g/c$')
    ax.set_ylabel(r'$\Delta \ OASPL, \ dB$')
    ax.set_title(rf'$\\theta = {theta[mic_iter]}^\circ$')
    ax.grid()
    plt.savefig(os.path.join(cases_directory,f'oaspl_vs_rg_{mic_iter}.png'),format = 'png')
    # ax.legend(['Thickness','Loading','Total'])

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .2,bottom = .15)
ax.plot(gust_width[plot_ind],(oaspl[plot_ind]-oaspl[plot_ind[baseline_ind]])[:,::-1],marker = 'o')
ax.plot(gust_width[plot_ind[baseline_ind]],0,marker = '^',markersize = 12,color = 'black')
ax.set_xlabel(r'$r_g/c$')
ax.set_ylabel(r'$\Delta \ OASPL, \ dB$')
ax.set_ylim([-4,4])
ax.legend(labels =leg_labs,fontsize = 12)
ax.grid()
plt.savefig(os.path.join(cases_directory,f'oaspl_vs_rg.png'),format = 'png')
