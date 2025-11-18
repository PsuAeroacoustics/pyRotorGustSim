import numpy as np
import matplotlib.pyplot as plt
import os
from help_funcs import *
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon,Rectangle

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
    pred_data['geometry_values'] = np.flip(pred_data['geometry_values'],axis = 0)
    pred_data['function_values'] = np.flip(pred_data['function_values'],axis = 0)

    theta = np.round(np.arctan2(pred_data['geometry_values'][:,:,0,1],pred_data['geometry_values'][:,:,0,0])*180/np.pi)%360
    phi = np.round(np.arctan2(pred_data['geometry_values'][:,:,0,-1],np.linalg.norm((pred_data['geometry_values'][:,:,0,0],pred_data['geometry_values'][:,:,0,1]),axis = 0))*180/np.pi)

    dt = pred_data['function_values'][0,1,0]-pred_data['function_values'][0,0,0]
    psi = (pred_data['function_values'][0,:,0]/((saved_params['omega']/(2*np.pi))**-1)*360)%360

    for theta_iter in range(pred_data['geometry_values'].shape[0]):
        for phi_iter in range(pred_data['geometry_values'].shape[1]):

            fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
            plt.subplots_adjust(left = .22,bottom = .15)

            if acs_params['thicknessNoiseFlag'] or acs_params['totalNoiseFlag']:
                ax.plot(pred_data['function_values'][theta_iter,phi_iter,:,0],pred_data['function_values'][theta_iter,phi_iter,:,1])
                ax.plot(pred_data['function_values'][theta_iter,phi_iter,:,0],pred_data['function_values'][theta_iter,phi_iter,:,2])
                ax.plot(pred_data['function_values'][theta_iter,phi_iter,:,0],pred_data['function_values'][theta_iter,phi_iter,:,-1])
                ax.legend(['Thickness','Loading','Total'])
            else:
                ax.plot(pred_data['function_values'][theta_iter,phi_iter,:,0],pred_data['function_values'][theta_iter,phi_iter,:,-1])


            ax.set(title = rf'$\\theta = {theta[theta_iter,phi_iter]}^\circ, \phi = {phi[theta_iter,phi_iter]}^\circ$',ylabel = r'Pressure \ [Pa]', xlabel =r'Time \ [s]' )
            # min_ind = pred_data['function_values'][mic_iter,:,-1].argmax()
            # ax.set_xlim([0,360])
            # ax.set_xlim([pred_data['function_values'][mic_iter,min_ind,0]-0.5*(saved_params['omega']/(2*np.pi))**-1,pred_data['function_values'][mic_iter,min_ind,0]+0.5*(saved_params['omega']/(2*np.pi))**-1])
            # ax.set_yticks(np.arange(10)*20-120)
            # ax.set_ylim([1.1*pred_data['function_values'][:,:,-1].min(),1.1*pred_data['function_values'][:,:,-1].max()])
            # ax.set_ylim([-25,20])
            ax.grid()
            plt.savefig(os.path.join(saved_params['case_dir'],f'tseries_{(theta_iter+1)*phi_iter}.png'),format = 'png')
            plt.close()

def plot_gust_profile(geom_params,input_params,res_param,observer_params,acs_params,saved_params):
    
    h = (np.arange(50+1)*(1.4+.2)/50-.2)/39.37/geom_params['radius']
    n = 2
    v_gust = input_params['af_params']['radius']/(2*np.pi*geom_params['radius'])*(h/((input_params['af_params']['core_size']/geom_params['AR'])**(2*n)+(h)**(2*n))**(1/n))
    print(f'max gust velocity = {v_gust.max()*3.281} fps')
    fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
    plt.subplots_adjust(left = .2,bottom = .15)
    ax.plot(h*geom_params['radius']*39.37,v_gust*3.281)
    ax.set_ylabel(r'V [fps]')
    ax.set_xlabel(r'Nozzle Width [in]')
    ax.set_xlim([-.2,1.4])
    ax.set_ylim(bottom = 0)
    ax.grid()
    plt.savefig(os.path.join(saved_params['case_dir'],f'gust_profile.png'),format = 'png')
    plt.close()

def plot_load_tseries(geom_params,input_params,res_param,observer_params,acs_params,saved_params):

    fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
    plt.subplots_adjust(left = .15,bottom = .15)
    ax.plot((saved_params['psi'][-int(2*np.pi/saved_params['dpsi']):])%(2*np.pi)*180/np.pi,saved_params['loads'][-int(2*np.pi/saved_params['dpsi']):,int(0.75*saved_params['N_elements']),-1])
    ax.set_ylabel(r'$ F_z \ [N]$')
    ax.set_xlabel(r'$ Azimuthal Angle \ [deg]$')
    ax.set_xlim([0,360])
    ax.grid()
    plt.savefig(os.path.join(saved_params['case_dir'],f'Fz_tseries.png'),format = 'png')
    plt.close()

    fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
    plt.subplots_adjust(left = .15,bottom = .15)
    ax.plot(saved_params['psi'][-int(2*np.pi/saved_params['dpsi']):]%(2*np.pi)*180/np.pi,np.gradient(saved_params['loads'][-int(2*np.pi/saved_params['dpsi']):,int(0.75*saved_params['N_elements']),-1],axis = 0))
    ax.set_ylabel(r'$ dF_z \ [N/deg]$')
    ax.set_xlabel(r'$ Azimuthal Angle \ [deg]$')
    ax.set_xlim([0,360])
    ax.set_ylim([-12,12])
    ax.grid()
    plt.savefig(os.path.join(saved_params['case_dir'],f'dFz_tseries.png'),format = 'png')
    plt.close()

def plot_load_dist(geom_params,input_params,res_param,observer_params,acs_params,saved_params):
    
    cmap = plt.cm.Spectral.reversed()
    
    fig,ax = plt.subplots(subplot_kw=dict(projection = 'polar'))
    lim = [0,15]
    levels = np.linspace(lim[0],lim[1],31)
    cbar_ticks = np.round(levels)[::4]
    # cbar_ticks = np.round(np.arange(50)*lim/50-lim)[::4]
    dist = ax.contourf(saved_params['psi'][:int(2*np.pi/saved_params['dpsi'])],saved_params['r'],saved_params['aoa'][-int(2*np.pi/saved_params['dpsi']):].T*180/np.pi,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
    cbar = fig.colorbar(dist,pad = .1)
    cbar.ax.set_ylabel(r'$\alpha \ [deg]$')
    cbar.ax.set_yticks(cbar_ticks)
    ax.set_rlim([0,saved_params['r'][-1]])
    plt.savefig(os.path.join(saved_params['case_dir'],'aoa.png'),format = 'png')
    plt.close()

    fig,ax = plt.subplots(subplot_kw=dict(projection = 'polar'))
    lim = [np.min(saved_params['loads'][-int(2*np.pi/saved_params['dpsi']):,:,-1]),np.max(saved_params['loads'][-int(2*np.pi/saved_params['dpsi']):,:,-1])]
    levels = np.linspace(lim[0],lim[1],41)
    cbar_ticks = np.round(levels)[::4]
    # cbar_ticks = np.round(np.arange(50)*lim/50-lim)[::4]
    dist = ax.contourf(saved_params['psi'][:int(2*np.pi/saved_params['dpsi'])],saved_params['r'],saved_params['loads'][-int(2*np.pi/saved_params['dpsi']):,:,-1].T,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
    cbar = fig.colorbar(dist,pad = .1)
    cbar.ax.set_ylabel(r'$ F_z \ [N]$')
    cbar.ax.set_yticks(cbar_ticks)
    ax.set_rlim([0,saved_params['r'][-1]])
    plt.savefig(os.path.join(saved_params['case_dir'],'Fz.png'),format = 'png')
    plt.close()

    # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
    # ax.plot(saved_params['loads'][:,-1,-1])

    d_loads = np.gradient(saved_params['loads'][:,:,-1],axis = 0)
    fig,ax = plt.subplots(subplot_kw=dict(projection = 'polar'))
    lim = [-10,10]
    levels = np.linspace(lim[0],lim[1],41)
    cbar_ticks = np.round(levels)[::4]
    # cbar_ticks = np.round(np.arange(50)*lim/50-lim)[::4]
    dist = ax.contourf(saved_params['psi'],saved_params['r'],d_loads.T,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
    cbar = fig.colorbar(dist,pad = .1)
    cbar.ax.set_ylabel(r'$\partial F_z /\partial \psi \ [N/deg]$')
    cbar.ax.set_yticks(cbar_ticks)
    ax.set_rlim([0,saved_params['r'][-1]])
    plt.savefig(os.path.join(saved_params['case_dir'],'dFz.png'),format = 'png')
    plt.close()

def plot_filt_load_dist(geom_params,input_params,res_param,observer_params,acs_params,saved_params):
    
    cmap = plt.cm.Spectral.reversed()
    
    fig,ax = plt.subplots(subplot_kw=dict(projection = 'polar'))
    lim = [np.min(saved_params['filt_loads']),np.max(saved_params['filt_loads'])]
    levels = np.linspace(lim[0],lim[1],41)
    cbar_ticks = np.round(levels)[::5]
    # cbar_ticks = np.round(np.arange(50)*lim/50-lim)[::4]
    dist = ax.contourf(saved_params['psi'],saved_params['r'],saved_params['filt_loads'][:,:,-1].T,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
    cbar = fig.colorbar(dist,pad = .1)
    cbar.ax.set_ylabel(r'$ F_z \ [N]$')
    cbar.ax.set_yticks(cbar_ticks)
    ax.set_rlim([0,saved_params['r'][-1]])
    plt.savefig(os.path.join(saved_params['case_dir'],'Fz_filt.png'),format = 'png')
    plt.close()

    d_loads = np.gradient(saved_params['filt_loads'][:,:,-1],axis = 0)
    fig,ax = plt.subplots(subplot_kw=dict(projection = 'polar'))
    lim = [-10,10]
    levels = np.linspace(lim[0],lim[1],41)
    cbar_ticks = np.round(levels)[::4]
    # cbar_ticks = np.round(np.arange(50)*lim/50-lim)[::4]
    dist = ax.contourf(saved_params['psi'],saved_params['r'],d_loads.T,levels = levels,cmap = cmap,norm = mcolors.CenteredNorm())
    cbar = fig.colorbar(dist,pad = .1)
    cbar.ax.set_ylabel(r'$\partial F_z /\partial \psi \ [N/deg]$')
    cbar.ax.set_yticks(cbar_ticks)
    ax.set_rlim([0,saved_params['r'][-1]])
    plt.savefig(os.path.join(saved_params['case_dir'],'dFz_filt.png'),format = 'png')
    plt.close()

def plot_res_params(geom_params,input_params,res_param,observer_params,acs_params,saved_params):

    dr = np.diff(saved_params['r_elem'][:2])[0]
    c =  saved_params['R']/geom_params['AR']*(geom_params['TR']-(geom_params['TR']-1)*saved_params['r_elem'])
    def get_blade_patch():
        verts = [(saved_params['r_elem'][0], -c[0]/saved_params['R']*3/4), (saved_params['r_elem'][0], c[0]/saved_params['R']/4), (saved_params['r_elem'][-1], c[-1]/saved_params['R']/4),(saved_params['r_elem'][-1], -c[-1]/saved_params['R']*3/4),(saved_params['r_elem'][0], -c[0]/saved_params['R']*3/4)]
        return Polygon(verts,color = 'tab:gray',alpha = .5)
    def get_element_patch(xy,height,color):
        return Rectangle(xy =xy , height = height,width = dr,color = color) 

    figsize =(5.635,1.5)
    fig,ax = plt.subplots(1,1, figsize = figsize)
    plt.subplots_adjust(bottom = .45,left = .05,right = .95,top = .95)
    ax.plot(saved_params['r_elem']*np.ones((2,len(saved_params['r_elem']))),[-c/saved_params['R']*3/4 , c/saved_params['R']/4],c = 'k',linestyle = '-',linewidth = .1)
    ax.add_patch(get_blade_patch())
    for i in range(res_param['N_patches']):
        patch_ind = saved_params['patch_type'][saved_params['filt_ind']]==i
        for ii in range(np.sum(patch_ind)):
            ax.add_patch(get_element_patch((saved_params['r'][saved_params['filt_ind']][patch_ind][ii]-dr/2,saved_params['c'][saved_params['filt_ind']][patch_ind][ii]/saved_params['R']*(.25-res_param['c_extents'][1])),saved_params['c'][saved_params['filt_ind']][patch_ind][ii]/saved_params['R']*np.diff(res_param['c_extents'])[0],color = f'C{i}'))
    ax.set_xlabel(r'r/R')
    ax.set_xlim([0.5,1])
    ax.set_ylim([-.15,.15])
    ax.set_yticks([])
    for i in ['top','right','left']:
        plt.gca().spines[i].set_visible(False) 
    plt.savefig(os.path.join(saved_params['case_dir'],'dist.png'),format = 'png')
    plt.close()


    # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
    # ax.scatter(saved_params['r'],saved_params['filt_ind'])
    # ax.set_xlabel('r/R')
    # plt.savefig(os.path.join(saved_params['case_dir'],'res_dist.png'), dpi=500, bbox_inches="tight", pad_inches=0.0)
    # plt.close()

    # for i in range(res_param['N_patches']):
    #     if len(saved_params[f'patch_{i}']['a'])>1:

    #         fig,ax = plt.subplots(1,2, figsize = (6.4,4.5))
    #         plt.subplots_adjust(wspace = 0.45)
    #         ax[0].plot(saved_params[f'patch_{i}']['a'],linestyle = '-.',c = 'black')
    #         ax[1].plot(saved_params[f'patch_{i}']['L'],linestyle = '-.',c = 'black')
    #         ax[0].stem(saved_params[f'patch_{i}']['a'])
    #         ax[1].stem(saved_params[f'patch_{i}']['L'])

    #         for i in range(2):
    #             ax[i].set_xlim([0,len(saved_params[f'patch_{i}']['a'])])
    #             ax[i].set_ylim(bottom = 0)
    #             ax[i].set_xticks(np.arange(len(saved_params[f'patch_{i}']['a']))[::2])
    #             ax[i].set_xlabel('$i$')

    #         ax[0].set_ylabel('$Radius, \ a_i \ [m]$')
    #         ax[1].set_ylabel('$Length, \ L_i \ [m]$')

    #         plt.savefig(os.path.join(saved_params['case_dir'],f'res_geom_{i}.png'), dpi=500, bbox_inches="tight", pad_inches=0.0)
    #         plt.close()

    # for i in range(res_param['N_patches']):
    #     fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
    #     plt.subplots_adjust(wspace = 0.45,bottom = 0.15)
    #     ax.scatter(saved_params['r'][saved_params['filt_ind']],saved_params[f'patch_{i}']['N']*np.ones(np.sum(saved_params['filt_ind'])))
    #     ax.set_ylabel('Resonators Per Blade Element, N')
    #     ax.set_xlabel('Radius, r/R')
    #     ax.set_ylim(bottom = 0)
    #     plt.savefig(os.path.join(saved_params['case_dir'],f'res_count_{i}.png'), dpi=500, bbox_inches="tight", pad_inches=0.0)
    #     plt.close()

def plot_res_resp(geom_params,input_params,res_param,observer_params,acs_params,saved_params):


    for i in range(res_param['N_patches']):

        r_ind = np.abs(saved_params['r'][i::res_param['N_patches']]-0.75).argmin()
        f_max_ind =  len(saved_params['f'])

        fig,ax = plt.subplots(2,1, figsize = (6.4,4.5))
        plt.subplots_adjust(bottom = 0.15,left = 0.15)
        ax[0].tick_params(axis = 'x', labelsize=0)
        ax[0].plot(saved_params['f'],np.real(saved_params['Z_tot'][:f_max_ind,r_ind]))
        ax[0].set_ylabel(r'$Resistance, \ \overline{\theta}$')
        ax[0].set_xlim([0,5e3])
        ax[0].set_ylim([0,10])
        ax[0].grid()

        ax[1].plot(saved_params['f'],np.imag(saved_params['Z_tot'][:f_max_ind,r_ind]))
        ax[1].set_ylabel(r'$Reactance, \ \overline{\chi}$')
        ax[1].set_xlim([500,saved_params['f'][-1]])
        ax[-1].set_xlabel(r'Frequency [Hz]')
        ax[-1].grid()
        ax[-1].set_xlim([0,5e3])
        ax[-1].set_ylim([-5, 5])
        plt.savefig(os.path.join(saved_params['case_dir'],f'Z_{i}.png'),format = 'png')
        plt.close()

        R = (saved_params['Z_tot'][:f_max_ind,r_ind]-1)/(saved_params['Z_tot'][:f_max_ind,r_ind]+1)
        fig,ax = plt.subplots(2,1, figsize = (6.4,4.5))
        plt.subplots_adjust(bottom = 0.15,left = 0.15)

        ax[0].tick_params(axis = 'x', labelsize=0)
        ax[0].plot(saved_params['f'],abs(R))
        ax[0].set_ylabel(r'$Reflection, \ |\mathit{R}|$')
        ax[0].grid()
        ax[0].set_xlim([0,5e3])
        ax[0].set_ylim([0, 1])

        ax[-1].plot(saved_params['f'],np.unwrap(np.angle(R)))
        ax[-1].set_ylabel(r'$Phase, \ \phi \ [rad]$')
        ax[-1].set_xlabel(r'Frequency [Hz]')
        ax[-1].grid()
        ax[-1].set_xlim([0,5e3])
        plt.savefig(os.path.join(saved_params['case_dir'],f'R_{i}.png'),format = 'png')
        plt.close()

        alpha = 1 - abs(R)**2
        fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        ax.plot(saved_params['f'],alpha)
        ax.set_ylabel(r'$Absorption, \alpha$')
        ax.grid()
        ax.set_xlim([0,5e3])
        ax.set_ylim([0, 1])
        plt.savefig(os.path.join(saved_params['case_dir'],f'alpha_{i}.png'),format = 'png')
        plt.close()
