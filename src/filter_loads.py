import numpy as np
from scipy import io
import os
import sys
import h5py
from scipy.optimize import differential_evolution
import aerosandbox as asb

sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'dependencies','resonator'))

from help_funcs import *
from res_funcs import *
from wopwop_input_generator import aperiodic_compact_loading_write
from plot import plot_res_resp
#%%


def filter_loads(geom_params,input_params,res_param,observer_params,acs_params,saved_params,opt = False):

    def get_oaspl():
        run_wopwop(cases = f"{input_params['case_name']}{os.sep}cases.nam",parallel = False)
        process_wopwop(cases_directory=saved_params['case_dir'],cases = 'cases.nam')
        acs_data = import_results_from_wopwop(cases_directory=saved_params['acs_dir'])
        oaspl = 10*np.log10(np.mean(acs_data['function_values'].squeeze()[:,:,-1]**2,axis = 1)/20e-6**2)
        return oaspl

    def opt_res(x,res_param,geom_params,saved_params):

        res_param['x'] = x
        apply_res(res_param,saved_params)

        V_res = np.sum(np.sum(np.pi*saved_params['a']**2*saved_params['L'])*saved_params['N'][saved_params['filt_ind']])

        if not np.all(saved_params['N'] == 0) and (V_res/V0 <= 1):
            for b_iter in range(geom_params['number_of_blades']):
                aperiodic_compact_loading_write(os.path.join(saved_params['acs_dir'],f'loading_blade_{b_iter}.dat'),t = saved_params['t'], loads = saved_params['filt_loads'],ascii = False)
            oaspl_filtered = get_oaspl()

            constraint = mu/2*np.max((0,-(1-V_res/V0)))**2+mu/2*np.max((0,-(0.25-saved_params['x'][0]/saved_params['x'][1])))**2
            residual = np.mean(10**(oaspl_filtered/10)/(10**(oaspl_baseline/10)))+constraint

            print(f'Inputs: {x} ')
            print(f"Volume fraction: {V_res/V0}")
            print(f'Delta OASPL: {np.mean(oaspl_baseline-oaspl_filtered)} ')
            print(f'Residual: {residual} ')

            x_hist.append(x)
            f_hist.append(residual)
            c_hist.append(constraint)
        else:
            residual = 1

        return residual


    def apply_res(res_param,saved_params):
        
        filt_loads = np.zeros(saved_params['loads'].shape)
        # indecies corresponding to the min/max spanwise extents of where the impedance patches are applied
        N_filt_ind = np.abs((np.expand_dims(res_param['filt_extents'],axis = -1)-saved_params['r'])).argmin(axis = -1)

        N,a,L = get_sample_info(res_param['x'],res_param['N_res'],res_param['OAR'],A_s,a_min=r_min,L_min=L_min)
        saved_params.update({'x':res_param['x'],'N':N,'a':a,'L':L,})

        if not np.all(saved_params['N'] == 0):
            
            filt_ind = np.zeros(saved_params['N_elements'],dtype=bool)
            filt_ind[N_filt_ind[0]:N_filt_ind[-1]] = 1  
            if np.any(N== 0):
                filt_ind[N==0] = 0

            Z_tot = np.ones((int(saved_params['iterations']/2),len(N[filt_ind])),dtype=complex)
            Z_tot[:f_ind], alpha = smeared_Z(f[:f_ind],N[filt_ind],a,L,res_param['N_res'],bool(res_param['facesheet']),t_fs = res_param['t_fs'],phi_fs = res_param['phi_fs'],N_fs = res_param['N_fs'],A_r = res_param['OAR'],A_s = A_s[filt_ind],M =res_param['M'],SPL = res_param['SPL'])
            # Z_tot_2, alpha = smeared_Z(f,N[filt_ind],a,L,res_param['N_res'],bool(res_param['facesheet']),t_fs = res_param['t_fs'],phi_fs = res_param['phi_fs'],N_fs = res_param['N_fs'],A_r = res_param['OAR'],A_s = A_s[filt_ind],M =res_param['M'],SPL = res_param['SPL'])

            #  obtains the filter response using the normal reflection coefficient
            filt_resp = get_filt_resp(Z_tot,odd = bool(saved_params['iterations']%2))

            # saved_params.update({'Z':Z_tot,'filt_resp':filt_resp})
            # plot_res_resp(geom_params,input_params,res_param,observer_params,acs_params,saved_params)

            xn_filt = apply_filt(saved_params['loads'][:,filt_ind,1:],filt_resp)
            filt_loads[:,filt_ind,1:]= xn_filt
            filt_loads[:,np.invert(filt_ind),1:] = saved_params['loads'][:,np.invert(filt_ind),1:]
            saved_params.update({'filt_ind':filt_ind,'Z':Z_tot,'filt_resp':filt_resp,'filt_loads':filt_loads})


    af = asb.Airfoil(geom_params['airfoil'])
    t = (af.coordinates[:,-1][af.coordinates[:,-1].argmax()]-af.coordinates[:,-1][-af.coordinates[:,-1].argmax()-1])*saved_params['c']
    # total volume of rotor blade [m^3]
    V0 = (saved_params['R']-saved_params['e'])*np.trapz(af.coordinates[:,0],af.coordinates[:,-1])*saved_params['c']
    r_elem = np.arange(saved_params['N_elements']+1)/saved_params['N_elements']*(saved_params['R']-saved_params['e'])+saved_params['e']

    # Surface area allocated to the impedance patches
    A_s = np.diff(r_elem)*res_param['A_s']*saved_params['c']
    # frequency resolution [Hz]
    df = (saved_params['iterations']*saved_params['dt'])**-1
    # frequency vector [Hz]
    f = np.arange(1,int(saved_params['iterations']/2)+1)*df

    f_ind = np.abs(f-res_param['f_max']).argmin()

    saved_params.update({'f':f})

    mu = 100

    r_min, r_max, L_min, L_max = np.max((.25e-3,res_param['r_min'])),np.min((saved_params['sos']/(2*np.pi*res_param['f_max']),res_param['r_max'])),np.min((saved_params['sos']/(res_param['f_max']*4),res_param['L_min'])),np.max((saved_params['sos']/(150*4),res_param['L_max']))

    if opt:

        oaspl_baseline = get_oaspl()
        saved_params.update({'oaspl_baseline':oaspl_baseline})


        x_hist = []
        f_hist = []
        c_hist = []

        bounds = get_opt_bounds(res_param['x'],r_min, r_max, L_min, L_max)
        opt_out = differential_evolution(opt_res,x0 = res_param['x'],bounds = bounds, polish=False,maxiter = int(res_param['maxiter']/(15*len(res_param['x']))),args = (res_param,geom_params,saved_params))
        res_param['x'] = opt_out.x
        saved_params.update({'x_hist':x_hist,'f_hist':f_hist,'c_hist':c_hist})
        print(f"Done! Minimizer: {res_param['x']}")

    apply_res(res_param,saved_params)
    for b_iter in range(geom_params['number_of_blades']):
        aperiodic_compact_loading_write(os.path.join(saved_params['acs_dir'],f'loading_blade_{b_iter}.dat'),t = saved_params['t'], loads = saved_params['filt_loads'],ascii = False)
    
    
