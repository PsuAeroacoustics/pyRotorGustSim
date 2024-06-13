import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
import f90nml
# sys.path.insert(0,os.path.join(os.path.dirname(__file__),'dependencies','pyWopwop'))
# import wopwop

#%%
import matplotlib.colors as mcolors

fontName = 'Times New Roman'
fontSize = 16
plt.rc('font',**{'family':'serif','serif':[fontName],'size':fontSize})
plt.rc('mathtext',**{'default':'regular'})
plt.rc('text',**{'usetex':False})
plt.rc('lines',**{'linewidth':2})

#%%

cases_directory = os.path.dirname(__file__)
cases = 'cases.nam'

#%%
case_info = {}
namelist = f90nml.read(os.path.join(cases_directory,  cases))
for k,v in namelist.items():
    if isinstance(v,dict):
        case_info_temp = {}
        for k_temp,v_temp in v.items():
            case_info_temp = {**case_info_temp, **{k_temp: v_temp}}
        case_info = {**case_info,**case_info_temp}
    else:
        case_info = {**case_info,**{k:v}}

#%%

pred_data = {}
#   imports reformatted data from wopwop in a dictionary
with h5py.File(os.path.join(case_info['globalfoldername'], 'pressure.h5'), 'r') as dat_file:
    for i,n in dat_file.items():
        pred_data = {**pred_data,**{i:n[()]}}

pred_data['geometry_values'] = np.flip(pred_data['geometry_values'],axis = 0).squeeze()
pred_data['function_values'] = np.flip(pred_data['function_values'],axis = 0).squeeze()

theta = np.round(np.arctan2(pred_data['geometry_values'][:,0,1],pred_data['geometry_values'][:,0,0])*180/np.pi)%(360)

dt = pred_data['function_values'][0,1,0]-pred_data['function_values'][0,0,0]

psi = pred_data['function_values'][0,:,0]/(pred_data['function_values'][0,-1,0]/2)*360-360

#%%

# fig,ax = plt.subplots(4,1, figsize = (4.5,4.5))
# for mic_iter in range(len(pred_data['function_values'])):
#     ax[mic_iter].plot(pred_data['function_values'][mic_iter,:,0]/pred_data['function_values'][mic_iter,-1,0],pred_data['function_values'][mic_iter,:,1])
#     ax[mic_iter].plot(pred_data['function_values'][mic_iter,:,0]/pred_data['function_values'][mic_iter,-1,0],pred_data['function_values'][mic_iter,:,2])
#     ax[mic_iter].plot(pred_data['function_values'][mic_iter,:,0]/pred_data['function_values'][mic_iter,-1,0],pred_data['function_values'][mic_iter,:,-1])
# ax[1].legend(['Thickness','Loading','Total','Measured'])

for mic_iter in range(len(pred_data['function_values'])):
    fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
    plt.subplots_adjust(left = .2)
    ax.plot(psi,pred_data['function_values'][mic_iter,:,2])
    ax.set_title(f'$\\theta = {theta[mic_iter]}$')
    ax.set_ylabel('$Pressure \ [Pa]$')
    ax.set_xlabel('$Blade \ Azimuth \ [deg]$')
    # ax.axis([270,320,-120,60])
    min_ind = pred_data['function_values'][mic_iter,:,2].argmin()
    ax.set_xlim([np.round(psi[min_ind]-35/2),np.round(psi[min_ind]+35/2)])
    ax.set_ylim([-120,60])

    ax.set_yticks(np.arange(10)*20-120)
    ax.grid()
    plt.savefig(f'tseries_{mic_iter}.png',format = 'png')

    # ax.plot(pred_data['function_values'][mic_iter,:,0]/pred_data['function_values'][mic_iter,-1,0],pred_data['function_values'][mic_iter,:,-1])
    # ax.legend(['Thickness','Loading','Total'])
