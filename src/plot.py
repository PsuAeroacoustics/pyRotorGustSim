import numpy as np
import matplotlib.pyplot as plt
import os
from help_funcs import *

# import matplotlib.font_manager as fm
# for font in fm.findSystemFonts(fontext='ttf'):
#     print(fm.FontProperties(fname=font).get_name())

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ["Times New Roman"]
plt.rcParams['font.size'] = 16

#%%

def plot_p_tseries(geom_params,input_params,observer_params,acs_params,saved_params):

#   imports reformatted data from wopwop in a dictionary
    pred_data = import_results_from_wopwop(cases_directory=saved_params['acs_dir'])
    pred_data['geometry_values'] = np.flip(pred_data['geometry_values'],axis = 0).squeeze()
    pred_data['function_values'] = np.flip(pred_data['function_values'],axis = 0).squeeze()

    theta = np.round(np.arctan2(pred_data['geometry_values'][:,0,1],pred_data['geometry_values'][:,0,0])*180/np.pi)%(360)
    dt = pred_data['function_values'][0,1,0]-pred_data['function_values'][0,0,0]
    psi = pred_data['function_values'][0,:,0]/(pred_data['function_values'][0,-1,0]/input_params['computational_params']['number_of_revs'])*360-360

    for mic_iter in range(len(pred_data['function_values'])):
        fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
        plt.subplots_adjust(left = .2,bottom = .15)
        ax.plot(psi,pred_data['function_values'][mic_iter,:,-1])
        ax.set_title(f'$\\theta = {theta[mic_iter]}^\circ$')
        ax.set_ylabel('Pressure [Pa]')
        ax.set_xlabel('Blade Azimuth [deg]')
        # ax.axis([270,320,-120,60])
        min_ind = pred_data['function_values'][mic_iter,:,-1].argmin()
        ax.set_xlim([np.round(psi[min_ind]-40/2),np.round(psi[min_ind]+40/2)])
        ax.set_ylim([-120,60])
        ax.set_yticks(np.arange(10)*20-120)
        ax.grid()
        plt.savefig(os.path.join(saved_params['case_dir'],f'tseries_{mic_iter}.png'),format = 'png')

def plot_gust_profile(geom_params,input_params,observer_params,acs_params,saved_params):
    
    h = (np.arange(50+1)*(1.4+.2)/50-.2)/39.37/geom_params['radius']
    n = 2
    v_gust = input_params['gust_params']['strength']/(2*np.pi*geom_params['radius'])*(h/((input_params['gust_params']['core_size']/geom_params['AR'])**(2*n)+(h)**(2*n))**(1/n))
    print(f'max gust velocity = {v_gust.max()*3.281} fps')
    fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
    plt.subplots_adjust(left = .2,bottom = .15)
    ax.plot(h*geom_params['radius']*39.37,v_gust*3.281)
    ax.set_ylabel('V [fps]')
    ax.set_xlabel('Nozzle Width [in]')
    ax.set_xlim([-.2,1.4])
    ax.set_ylim(bottom = 0)
    ax.grid()
    plt.savefig(os.path.join(saved_params['case_dir'],f'gust_profile.png'),format = 'png')

