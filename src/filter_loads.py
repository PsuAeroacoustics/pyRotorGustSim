import numpy as np
from scipy import io
import os
import sys
import h5py
from scipy.optimize import differential_evolution
import aerosandbox as asb
# from xfoil import XFoil
# from xfoil import model

sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'dependencies','resonator'))

from help_funcs import *
from res_funcs import *
from wopwop_input_generator import aperiodic_compact_loading_write
from plot import plot_res_resp

#%%


def filter_loads(geom_params,input_params,res_param,observer_params,acs_params,saved_params,opt = False,noncompact = False):

    def get_oaspl():
        run_wopwop(cases = f"{input_params['case_name']}{os.sep}cases.nam",parallel = False)
        process_wopwop(cases_directory=saved_params['case_dir'],cases = 'cases.nam')
        acs_data = import_results_from_wopwop(cases_directory=saved_params['acs_dir'])
        oaspl = 10*np.log10(np.mean(acs_data['function_values'].squeeze()[:,:,-1]**2,axis = 1)/20e-6**2)
        return oaspl

    def opt_res(x,res_param,geom_params,saved_params):

        res_param.update({'x':list(x)})
        apply_res(res_param,saved_params)

        if 'sigma' not in res_param:
            

            V_res =np.sum([np.sum((np.pi*saved_params[f'patch_{i}']['a']**2*saved_params[f'patch_{i}']['L']*saved_params[f'patch_{i}']['N'].T)[saved_params['patch_type']==i][saved_params['filt_ind'][saved_params['patch_type']==i]]) for i in range(res_param['N_patches'])])
            L= np.array([saved_params[f'patch_{i}']['L'] for i in range(res_param['N_patches'])]).flatten()
            N= np.array([saved_params[f'patch_{i}']['N'] for i in range(res_param['N_patches'])]).flatten()

            if np.all(N>0) and (V_res/V0 <= 0.5) and (V_res/V0 > 0):
                for b_iter in range(geom_params['number_of_blades']):
                    aperiodic_compact_loading_write(os.path.join(saved_params['acs_dir'],f'loading_blade_{b_iter}.dat'),t = saved_params['t'], loads = saved_params['filt_loads'],ascii = False)
                oaspl_filtered = get_oaspl()


                constraint = mu/2*np.max((0,-(1-V_res/V0)))**2+mu/2*np.sum(np.max((np.zeros(len(L)),-(0.25-x[0]/L)),axis = 0))**2
                residual = np.mean(10**(oaspl_filtered/10)/(10**(oaspl_baseline/10)))+constraint

                print(f'Inputs: {x}')
                print(f"Volume fraction: {V_res/V0}")
                print(f'Delta OASPL: {np.mean(oaspl_baseline-oaspl_filtered)} ')
                print(f'Residual: {residual} ')

                x_hist.append(x)
                f_hist.append(residual)
                c_hist.append(constraint)
            else:
                residual = 1

        else:
            for b_iter in range(geom_params['number_of_blades']):
                aperiodic_compact_loading_write(os.path.join(saved_params['acs_dir'],f'loading_blade_{b_iter}.dat'),t = saved_params['t'], loads = saved_params['filt_loads'],ascii = False)
            oaspl_filtered = get_oaspl()

            residual = np.mean(10**(oaspl_filtered/10)/(10**(oaspl_baseline/10)))

            print(f'Inputs: {x} ')
            print(f'Delta OASPL: {np.mean(oaspl_baseline-oaspl_filtered)} ')
            print(f'Residual: {residual} ')

            x_hist.append(x)
            f_hist.append(residual)

        return residual


    def apply_res(res_param,saved_params):
        
        filt_loads = np.zeros(saved_params['loads'].shape)
        # indecies corresponding to the min/max spanwise extents of where the impedance patches are applied
        
        filt_ind = np.zeros(saved_params['N_elements'],dtype=bool)
        filt_ind_extents = (saved_params['r'] >= res_param['r_extents'][0]) & (saved_params['r'] <= res_param['r_extents'][-1])
        filt_ind[filt_ind_extents] = 1


        if 'sigma' in res_param:
            if 'x' in res_param:
                saved_params.update(get_sample_info(filt_ind,A_s,**res_param))

            Z_tot = np.ones((int(saved_params['iterations']/2),np.sum(filt_ind)),dtype=complex)
            Z_tot[:f_ind] = np.expand_dims(init_porous_res(f[:f_ind],**res_param).Z,axis = -1)

        else:
            
            Z_tot = np.ones(((int(saved_params['iterations']/2), saved_params['N_elements'])),dtype=complex)
            filt_resp = np.ones(((saved_params['iterations'], saved_params['N_elements'])),dtype=complex)
            saved_params.update(get_sample_info(filt_ind,A_s,**res_param))
            for i in range(res_param['N_patches']):
                if np.all(saved_params[f'patch_{i}']['N']) > 0:                
                    Z_tot[:f_ind,i::res_param['N_patches']] = smeared_Z(f[:f_ind],A_s,**saved_params[f'patch_{i}'])[:,i::res_param['N_patches']]
                
                # saved_params[f'patch_{i}'].update({'Z':Z_tot[:f_ind]})

                            # Z_tot_2, alpha = smeared_Z(f,N[filt_ind],a,L,res_param['N_res'],bool(res_param['facesheet']),t_fs = res_param['t_fs'],phi_fs = res_param['phi_fs'],N_fs = res_param['N_fs'],A_r = res_param['OAR'],A_s = A_s[filt_ind],M =res_param['M'],SPL = res_param['SPL'])

            #  obtains the filter response using the normal reflection coefficient

            # saved_params.update({'Z':Z_tot,'filt_resp':filt_resp})
            # plot_res_resp(geom_params,input_params,res_param,observer_params,acs_params,saved_params)
        if np.all(np.array([saved_params[f'patch_{i}']['N'] for i in range(res_param['N_patches'])])>0) or 'sigma' in res_param:
            
            Z_tot_nc = np.ones((int(N/2),saved_params['N_elements']),complex)
            Z_tot_nc[:len(Z_tot)] = Z_tot
            filt_resp_nc = get_filt_resp(Z_tot_nc)

            loads_nc = np.zeros((N,saved_params['N_elements'],3))
            loads_nc[:saved_params['iterations']] = saved_params['loads']

            # filt_resp = get_filt_resp(Z_tot)

            if noncompact:
                cp_dist_filt = np.copy(cp_dist)
                xn_filt = apply_filt(cp_dist[:,filt_ind,c_extent_ind[1]:c_extent_ind[0]+1],filt_resp)
                cp_dist_filt[:,filt_ind,c_extent_ind[1]:c_extent_ind[0]+1] = xn_filt
                
                cl_filt = np.trapz(cp_dist_filt,x = af.coordinates[:,0],axis = -1)*np.cos(aoa_cp_dist*np.pi/180)+np.trapz(cp_dist_filt,x = af.coordinates[:,-1],axis = -1)*np.sin(aoa_cp_dist*np.pi/180)
                dCz = cl_filt*np.cos(saved_params['phi_eff'])-saved_params['CD']*np.sin(saved_params['phi_eff'])
                dCx = cl_filt*np.sin(saved_params['phi_eff'])+saved_params['CD']*np.cos(saved_params['phi_eff'])
                dCT = 0.5*saved_params['c']/(np.pi*saved_params['R'])*saved_params['U']**2*dCz
                dCP = 0.5*saved_params['c']/(np.pi*saved_params['R'])*saved_params['r']*saved_params['U']**2*dCx
                filt_loads[:,:,1]= -dCP*input_params['flight_params']['density']*np.pi*saved_params['R']**2*(input_params['flight_params']['omega']*saved_params['R'])**2
                filt_loads[:,:,-1] = dCT*input_params['flight_params']['density']*np.pi*saved_params['R']*(input_params['flight_params']['omega']*saved_params['R'])**2


            else:
                # loads_nc = np.concatenate((saved_params['loads'],np.zeros(saved_params['loads'].shape)),axis = 0)
                # filt_resp_nc= np.concatenate((filt_resp,np.zeros(filt_resp.shape)),axis = 0)
                xn_filt = apply_filt(loads_nc[:,filt_ind,1:],filt_resp_nc[:,filt_ind])

                # saved_params['loads'][int(len(saved_params['loads'])/2):] = 0
                # xn_filt = apply_filt(saved_params['loads'][:,filt_ind,1:],filt_resp[:,filt_ind])
                filt_loads[:,filt_ind,1:]= xn_filt[:saved_params['iterations']]
                filt_loads[:,np.invert(filt_ind),1:] = saved_params['loads'][:,np.invert(filt_ind),1:]

            saved_params.update({'filt_resp':filt_resp,'Z_tot':Z_tot,'filt_loads':filt_loads,'filt_ind':filt_ind,'x':res_param['x']})


    af = asb.Airfoil(geom_params['airfoil'])
    # af = asb.Airfoil('OA209')

    af.coordinates = af.repanel(n_points_per_side = int(geom_params['airfoil_points']/2)).coordinates

    if noncompact:
        aoa_polar = np.arange((20-(-5))/.5+1)*.5-5
        Re = 0.75*geom_params['radius']*input_params['flight_params']['omega']*saved_params['c']/input_params['flight_params']['kinematic_viscosity']
        M = 0.75*input_params['flight_params']['omega']*geom_params['radius']/input_params['flight_params']['sos']
        polar = af.get_aero_from_neuralfoil(alpha=aoa_polar, Re=Re,mach =M,model_size='xlarge')
        aoa_cp_dist = aoa_polar[abs(polar['CL']-0.5*(polar['CL'].max()-np.abs(polar['CL']).min())).argmin()]
        
        xf = XFoil()
        xf.airfoil = model.Airfoil(x = af.coordinates[:,0],y = af.coordinates[:,1])    
        xf.max_iter = 100
        xf.Re = Re
        xf.M = M
        xf.a(aoa_cp_dist, as_dict=False)
        cp = xf.get_cp_distribution()[-1]

        # cl_trapz = lambda cp_: np.trapz(cp_,x = af.coordinates[:,0],axis = -1)*np.cos(aoa_polar*np.pi/180)+np.trapz(cp_,x = af.coordinates[:,-1],axis = -1)*np.sin(aoa_polar*np.pi/180)
        cl = np.trapz(cp,x = af.coordinates[:,0],axis = -1)*np.cos(aoa_cp_dist*np.pi/180)+np.trapz(cp,x = af.coordinates[:,-1],axis = -1)*np.sin(aoa_cp_dist*np.pi/180)

        cp_dist_nominal = cp/cl
        cp_dist = np.expand_dims(saved_params['CL'],axis = -1)*cp_dist_nominal
        # cl = np.trapz(cp_dist,x = af.coordinates[:,0],axis = -1)*np.cos(aoa_cp_dist*np.pi/180)+np.trapz(cp_dist,x = af.coordinates[:,-1],axis = -1)*np.sin(aoa_cp_dist*np.pi/180)

        # import matplotlib.pyplot as plt

        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # plt.subplots_adjust(left = 0.15,bottom = .15)
        # ax.plot(af.coordinates[:,0],af.coordinates[:,1])
        # ax.set_xlabel('x/c')
        # ax.set_xlabel('y/c')
        # ax.set_ylim(-0.5,0.5)
        # ax.grid()

        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # plt.subplots_adjust(left = 0.15,bottom = .15)
        # ax.plot(aoa_polar,polar['CL'])
        # ax.scatter(aoa_cp_dist,polar['CL'][abs(aoa_polar-aoa_cp_dist).argmin()])
        # ax.set_xlabel(r'$\alpha \ [deg]$')
        # ax.set_ylabel('$C_L$')
        # ax.grid()

        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # plt.subplots_adjust(left = 0.15,bottom = .15)
        # ax.plot(af.coordinates[:,0],cp)
        # ax.plot(af.coordinates[:,0],cp_dist_nominal)
        # ax.set_xlabel(r'$x/c$')
        # ax.set_ylabel('$C_P$')
        # ax.legend([r'$\alpha=5.5^\circ$',r'$Scaled \ (C_L=1)$'])
        # ax.set_ylim(2,-3.5)
        # ax.grid()

        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # plt.subplots_adjust(left = 0.15,bottom = .15)
        # ax.plot(saved_params['CL'][:,int(0.75*48)])
        # ax.plot(cl[:,int(0.75*48)])
        # ax.set_xlabel('$\psi$ [deg]')
        # ax.set_ylabel('$C_L$')
        # ax.set_xlim([0,360])
        # ax.legend(['Original','From $C_P$ Distribution'])
        # ax.grid()

        # import matplotlib.pyplot as plt
        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # plt.subplots_adjust(left = 0.15,bottom = .15)
        # # ax.plot(saved_params['v_gust'][:,-1])
        # ax.plot(np.gradient(saved_params['loads'][:,-1,-1]))
        # ax.set_xlabel('$\psi$ [deg]')
        # ax.grid()

        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # plt.subplots_adjust(left = 0.15,bottom = .15)
        # ax.plot(np.gradient(xn_filt[:,-1,-1]))
        # # ax.plot(cl[:,int(0.75*48)])
        # ax.set_xlabel('$\psi$ [deg]')
        # ax.grid()

    t = af.max_thickness()*saved_params['c']
    
    if 'sigma' in res_param:
        res_param['t'] = t 

    # total volume of rotor blade [m^3]
    V0 = np.trapz(np.trapz(np.expand_dims(af.coordinates[:,0],axis = -1)*saved_params['c'],np.expand_dims(af.coordinates[:,-1],axis = -1)*saved_params['c'],axis = 0),x = saved_params['r']*saved_params['R'])
    r_elem = np.arange(saved_params['N_elements']+1)/saved_params['N_elements']*(saved_params['R']-saved_params['e'])+saved_params['e']

    # chordwise extent of the treatment
    c_extent_ind = np.abs(np.array(af.coordinates[af.coordinates[:,1]>0][:,0])-np.expand_dims(res_param['c_extents'],axis = -1)).argmin(axis  = -1)

    # Surface area allocated to the impedance patches
    A_s = np.diff(r_elem)*np.sum(np.linalg.norm(np.diff(af.coordinates[slice(c_extent_ind[1],c_extent_ind[0]+1)],axis = 0),axis = -1))*saved_params['c']
    # frequency resolution [Hz]

        # n_revs = np.ceil(omega/(2*np.pi*res_param['df']))
    # iterations = int(n_revs*2*np.pi/dpsi*2)

    N = int(np.max((np.ceil((saved_params['dt']*res_param['df'])**-1),np.ceil(saved_params['iterations']*2))))

    df = (N*saved_params['dt'])**-1
    # frequency vector [Hz]
    f = np.arange(1,int(saved_params['iterations']/2)+1)*df

    # low pass filter cutoff frequency (dimensions of impedance patch should be less than max acoustic wavelength of interest)
    res_param['f_max'] = 0.75/np.diff(r_elem[:2]).squeeze()*saved_params['sos']/(2*np.pi)
    # f_ind = np.abs(f-res_param['f_max']).argmin()
    f_ind = len(f)

    patch_type = np.zeros(saved_params['N_elements'])
    for i in range(res_param['N_patches']):
        patch_type[i::res_param['N_patches']] = i

    saved_params.update({'f':f[:f_ind],'A_s':A_s,'patch_type':patch_type,'c_extents':res_param['c_extents']})

    if 'sigma' not in res_param:

        mu = 100
        r_min, r_max, L_min, L_max = np.max((.25e-3,res_param['r_min'])),np.min((saved_params['sos']/(2*np.pi*res_param['f_max']),res_param['r_max'])),np.min((saved_params['sos']/(res_param['f_max']*4),res_param['L_min'])),np.min((geom_params['radius']*(1-geom_params['r_c']),res_param['L_max']))
        res_param.update({'r_min':r_min, 'r_max':r_max, 'L_min':L_min, 'L_max':L_max})
        
    if opt:

        oaspl_baseline = get_oaspl()
        saved_params.update({'oaspl_baseline':oaspl_baseline})

        x_hist = []
        f_hist = []
        c_hist = []

        if 'sigma' not in res_param:
            bounds = get_opt_bounds(**res_param)

        else:
            bounds = Bounds(lb = (0,0,0), ub = (1,1,1),keep_feasible=True)

        opt_out = differential_evolution(opt_res,x0 = res_param['x'],bounds = bounds, polish=False,maxiter = int(res_param['maxiter']/(15*len(res_param['x']))),args = (res_param,geom_params,saved_params))
        res_param['x'] = list(opt_out.x)
        saved_params.update({'x_hist':x_hist,'f_hist':f_hist,'c_hist':c_hist})
        print(f"Done! Minimizer: {np.asarray(res_param['x']).squeeze()}")

    apply_res(res_param,saved_params)
    for b_iter in range(geom_params['number_of_blades']):
        aperiodic_compact_loading_write(os.path.join(saved_params['acs_dir'],f'loading_blade_{b_iter}.dat'),t = saved_params['t'], loads = saved_params['filt_loads'],ascii = False)
    
    
