#!/usr/bin/env python3

import argparse
import os
from help_funcs import *
from scipy.optimize import minimize,Bounds

#%%

def main():
    parser = argparse.ArgumentParser("opt_op_conditions",description='This script wraps a local minimizer around the rotor_gust_interaction script. Optimizes the gust size, velocity and tip speed to determine point that provides the maximum measurable difference for a given resonator treatment.')
    parser.add_argument(
        "-input_geom",
        type= str,
        required=False
    )
    parser.add_argument(
        "-input_param",
        type= str,
        required=False
    )
    parser.add_argument(
        "-observer_param",
        type= str,
        required=False
    )
    parser.add_argument(
        "-acs_param",
        type= str, 
        required=False
    )
    parser.add_argument(
        "-res_param",
        type= str, 
        required=False
    )

    def opt_params(x0):
        # print(f'Input Parameters: M_g: {np.round(x0[0],3)}, r_g: {np.round(x0[1],3)}')
        x0 = dimenstion_params(x0)
        print(f'Input Parameters: M_t: {np.round(x0[2],3)}, M_g: {np.round(x0[0],3)}, r_g: {np.round(x0[1],3)}')
        write_param_files(x0)
        oaspl = []
        for i in range(2):
            run_rotor_gust_interaction(filt = bool(i))
            oaspl.append(get_oaspl())
        d_oaspl = np.mean(np.mean((oaspl[1]-oaspl[0]),axis = -1),axis = -1)
        print(f'Average dOASPL: {np.round(d_oaspl,1)}')
        return d_oaspl

    def write_param_files(x0):
        
        input_params['gust_params']['strength'] = x0[0]
        input_params['gust_params']['core_size']  = x0[1]
        input_params['flight_params']['omega']=x0[2]*input_params['flight_params']['sos']/geom_params['radius']

        with open(os.path.join(cwd,args.input_param),"w") as param_file:
            json.dump(input_params,param_file,indent=2)

    def run_rotor_gust_interaction(filt = False):
        if filt:
            subprocess.run(['rotor_gust_interaction.py','--aero','--acs','--filt','-input_geom', args.input_geom, '-input_param', args.input_param, '-observer_param', args.observer_param, '-acs_param', args.acs_param, '-res_param', args.res_param],check = True)
        else:
            subprocess.run(['rotor_gust_interaction.py','--aero','--acs','-input_geom', args.input_geom, '-input_param', args.input_param, '-observer_param', args.observer_param, '-acs_param', args.acs_param, '-res_param', args.res_param],check = True)

    def get_oaspl():
        acs_data = import_results_from_wopwop(acs_dir)
        oaspl = 20*np.log10(np.sqrt(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -1))/20e-6)
        return oaspl
    

    args = parser.parse_args()
    
    cwd = os.getcwd()

    with open(os.path.join(cwd,args.input_param)) as param_file:
        input_params = json.load(param_file)
    with open(os.path.join(cwd,args.input_geom)) as param_file:
        geom_params = json.load(param_file)


    case_dir = os.path.join(cwd,input_params['case_name'])
    acs_dir = os.path.join(case_dir,'acoustics')

    # sets parameter bounds 
    Mg_bounds = (.75,30)
    rg_bounds = (0.3,3)
    Mt_bounds = (0.2,0.6)
    bounds = np.array((Mg_bounds,rg_bounds,Mt_bounds))

    nondimenstion_params = lambda x: (x-bounds[:,0])/np.diff(bounds).squeeze()
    dimenstion_params = lambda x: x*np.diff(bounds).squeeze()+bounds[:,0]

    
    x0 = [input_params['gust_params']['strength'],input_params['gust_params']['core_size'],input_params['flight_params']['omega']*geom_params['radius']/input_params['flight_params']['sos']]
    
    # x0 = [input_params['gust_params']['strength'],input_params['gust_params']['core_size']]

    res = minimize(opt_params, x0 = x0, args=(), method='L-BFGS-B', jac=None, hess=None, hessp=None, bounds=Bounds(lb = np.zeros(3), ub = np.ones(3),keep_feasible=True), constraints=(), tol=1e-2, callback=None, options={'maxiter':10,'disp':True})
    print('All Done!')

if __name__ == "__main__":

	main()
	print("exiting main.py")