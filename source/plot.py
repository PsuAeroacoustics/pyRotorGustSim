import numpy as np
import matplotlib.pyplot as plt
import h5py
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ["Arial"]
plt.rcParams['font.size'] = 16

#%%

def plot_p_tseries(geom_params,input_params,observer_params,acs_params,saved_params):

    pred_data = {}
#   imports reformatted data from wopwop in a dictionary
    with h5py.File(os.path.join(saved_params['acs_dir'] ,'pressure.h5'), 'r') as dat_file:
        for i,n in dat_file.items():
            pred_data = {**pred_data,**{i:n[()]}}
    pred_data['geometry_values'] = np.flip(pred_data['geometry_values'],axis = 0).squeeze()
    pred_data['function_values'] = np.flip(pred_data['function_values'],axis = 0).squeeze()

    theta = np.round(np.arctan2(pred_data['geometry_values'][:,0,1],pred_data['geometry_values'][:,0,0])*180/np.pi)%(360)
    dt = pred_data['function_values'][0,1,0]-pred_data['function_values'][0,0,0]
    psi = pred_data['function_values'][0,:,0]/(pred_data['function_values'][0,-1,0]/input_params['computational_params']['number_of_revs'])*360-360

    for mic_iter in range(len(pred_data['function_values'])):
        fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
        plt.subplots_adjust(left = .2)
        ax.plot(psi,pred_data['function_values'][mic_iter,:,-1])
        ax.set_title(f'$\\theta = {theta[mic_iter]}$')
        ax.set_ylabel('$Pressure \ [Pa]$')
        ax.set_xlabel('$Blade \ Azimuth \ [deg]$')
        # ax.axis([270,320,-120,60])
        min_ind = pred_data['function_values'][mic_iter,:,-1].argmin()
        ax.set_xlim([np.round(psi[min_ind]-40/2),np.round(psi[min_ind]+40/2)])
        ax.set_ylim([-120,60])
        ax.set_yticks(np.arange(10)*20-120)
        ax.grid()
        plt.savefig(os.path.join(saved_params['case_dir'],f'tseries_{mic_iter}.png'),format = 'png')
