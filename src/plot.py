import numpy as np
import matplotlib.pyplot as plt
import os
from help_funcs import *
import matplotlib.colors as mcolors

# import matplotlib.font_manager as fm
# for font in fm.findSystemFonts(fontext='ttf'):
#     print(fm.FontProperties(fname=font).get_name())

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ["Times New Roman"]
plt.rcParams['font.size'] = 16

#%%

def plot_p_tseries(geom_params,input_params,res_param,observer_params,acs_params,saved_params):

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

def plot_gust_profile(geom_params,input_params,res_param,observer_params,acs_params,saved_params):
    
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


def plot_load_dist(geom_params,input_params,res_param,observer_params,acs_params,saved_params):
    
    cmap = plt.cm.Spectral.reversed()
    
    fig,ax = plt.subplots(subplot_kw=dict(projection = 'polar'))
    lim = [np.min(saved_params['loads']),np.max(saved_params['loads'])]
    levels = np.linspace(lim[0],lim[1],50)
    cbar_ticks = np.round(levels)[::4]
    # cbar_ticks = np.round(np.arange(50)*lim/50-lim)[::4]
    dist = ax.contourf(saved_params['psi'],saved_params['r'],saved_params['loads'][:,:,-1].T,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
    cbar = fig.colorbar(dist,pad = .1)
    cbar.ax.set_ylabel('$ F_z \ [N]$')
    cbar.ax.set_yticks(cbar_ticks)
    plt.savefig(os.path.join(saved_params['case_dir'],'Fz.png'),format = 'png')

    d_loads = np.gradient(saved_params['loads'][:,:,-1],axis = 0)
    fig,ax = plt.subplots(subplot_kw=dict(projection = 'polar'))
    lim = [np.min(d_loads),np.max(d_loads)]
    levels = np.linspace(lim[0],lim[1],50)
    cbar_ticks = np.round(levels)[::4]
    # cbar_ticks = np.round(np.arange(50)*lim/50-lim)[::4]
    dist = ax.contourf(saved_params['psi'],saved_params['r'],d_loads.T,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
    cbar = fig.colorbar(dist,pad = .1)
    cbar.ax.set_ylabel('$\partial F_z /\partial \psi \ [N/deg]$')
    cbar.ax.set_yticks(cbar_ticks)
    plt.savefig(os.path.join(saved_params['case_dir'],'dFz.png'),format = 'png')

def plot_filt_load_dist(geom_params,input_params,res_param,observer_params,acs_params,saved_params):
    
    cmap = plt.cm.Spectral.reversed()
    
    fig,ax = plt.subplots(subplot_kw=dict(projection = 'polar'))
    lim = [np.min(saved_params['filt_loads']),np.max(saved_params['filt_loads'])]
    levels = np.linspace(lim[0],lim[1],50)
    cbar_ticks = np.round(levels)[::4]
    # cbar_ticks = np.round(np.arange(50)*lim/50-lim)[::4]
    dist = ax.contourf(saved_params['psi'],saved_params['r'],saved_params['filt_loads'][:,:,-1].T,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
    cbar = fig.colorbar(dist,pad = .1)
    cbar.ax.set_ylabel('$ F_z \ [N]$')
    cbar.ax.set_yticks(cbar_ticks)
    plt.savefig(os.path.join(saved_params['case_dir'],'Fz_filt.png'),format = 'png')

    d_loads = np.gradient(saved_params['filt_loads'][:,:,-1],axis = 0)
    fig,ax = plt.subplots(subplot_kw=dict(projection = 'polar'))
    lim = [np.min(d_loads),np.max(d_loads)]
    levels = np.linspace(lim[0],lim[1],50)
    cbar_ticks = np.round(levels)[::4]
    # cbar_ticks = np.round(np.arange(50)*lim/50-lim)[::4]
    dist = ax.contourf(saved_params['psi'],saved_params['r'],d_loads.T,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
    cbar = fig.colorbar(dist,pad = .1)
    cbar.ax.set_ylabel('$\partial F_z /\partial \psi \ [N/deg]$')
    cbar.ax.set_yticks(cbar_ticks)
    plt.savefig(os.path.join(saved_params['case_dir'],'dFz_filt.png'),format = 'png')

def plot_res_params(geom_params,input_params,res_param,observer_params,acs_params,saved_params):

    fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
    ax.scatter(saved_params['r'],saved_params['filt_ind'])
    ax.set_xlabel('r/R')
    plt.savefig(os.path.join(saved_params['case_dir'],'res_dist.png'), dpi=500, bbox_inches="tight", pad_inches=0.0)
    
    if len(saved_params['a'])>1:

        fig,ax = plt.subplots(1,2, figsize = (6.4,4.5))
        plt.subplots_adjust(wspace = 0.45)
        ax[0].plot(saved_params['a'],linestyle = '-.',c = 'black')
        ax[1].plot(saved_params['L'],linestyle = '-.',c = 'black')
        ax[0].stem(saved_params['a'])
        ax[1].stem(saved_params['L'])

        for i in range(2):
            ax[i].set_xlim([0,len(saved_params['a'])])
            ax[i].set_ylim(bottom = 0)
            ax[i].set_xticks(np.arange(len(saved_params['a']))[::2])
            ax[i].set_xlabel('$i$')

        ax[0].set_ylabel('$Radius, \ a_i \ [m]$')
        ax[1].set_ylabel('$Length, \ L_i \ [m]$')

        plt.savefig(os.path.join(saved_params['case_dir'],'res_geom.png'), dpi=500, bbox_inches="tight", pad_inches=0.0)

def plot_res_resp(geom_params,input_params,res_param,observer_params,acs_params,saved_params):

    r_ind = int(.75*saved_params['Z'].shape[-1])

    fig,ax = plt.subplots(2,1, figsize = (6.4,4.5))
    ax[0].tick_params(axis = 'x', labelsize=0)
    ax[0].plot(saved_params['f'],np.real(saved_params['Z'][:,r_ind]))
    ax[0].set_ylabel(r'$Resistance, \ \overline{\theta}$')
    ax[0].set_xlim([500,saved_params['f'][-1]])
    ax[0].set_ylim([0,10])
    ax[0].grid()

    ax[1].plot(saved_params['f'],np.imag(saved_params['Z'][:,r_ind]))
    ax[1].set_ylabel(r'$Reactance, \ \overline{\chi}$')
    ax[1].set_xlim([500,saved_params['f'][-1]])
    ax[-1].set_xlabel('Frequency [Hz]')
    ax[-1].grid()
    ax[-1].set_xlim([500,saved_params['f'][-1]])
    ax[-1].set_ylim([-5, 5])
    plt.savefig(os.path.join(saved_params['case_dir'],'Z.png'),format = 'png')

    R = (saved_params['Z'][:,r_ind]-1)/(saved_params['Z'][:,r_ind]+1)
    fig,ax = plt.subplots(2,1, figsize = (6.4,4.5))
    ax[0].tick_params(axis = 'x', labelsize=0)
    ax[0].plot(saved_params['f'],abs(R))
    ax[0].set_ylabel(r'$Reflection, \ |\mathit{R}|$')
    ax[0].grid()
    ax[0].set_xlim([500,saved_params['f'][-1]])
    ax[0].set_ylim([0, 1])

    ax[-1].plot(saved_params['f'],np.unwrap(np.angle(R)))
    ax[-1].set_ylabel('$Phase, \ \phi \ [rad]$')
    ax[-1].set_xlabel('Frequency [Hz]')
    ax[-1].grid()
    ax[-1].set_xlim([500,saved_params['f'][-1]])
    plt.savefig(os.path.join(saved_params['case_dir'],'R.png'),format = 'png')

    alpha = 1 - abs(R)**2
    fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
    ax.tick_params(axis = 'x', labelsize=0)
    ax.plot(saved_params['f'],alpha)
    ax.set_ylabel(r'$Absorption, \alpha$')
    ax.grid()
    ax.set_xlim([500,saved_params['f'][-1]])
    ax.set_ylim([0, 1])
    plt.savefig(os.path.join(saved_params['case_dir'],'alpha.png'),format = 'png')
