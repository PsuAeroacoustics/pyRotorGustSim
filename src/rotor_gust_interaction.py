#!/usr/bin/env python3

import argparse
import os
from shutil import rmtree
from help_funcs import *
from compute_aero import *
from filter_loads import *
from build_blade_geom import *
from wopwop_input_configure import *
from plot import *

#%%

def main():
    parser = argparse.ArgumentParser("rotor_gust_interaction",description='Simulates a gust interacting with a hovering rotor, only the positive half of the gust profile is considered.')
    parser.add_argument(
		"--aero",
		action='store_true',
		help="Compute aerodynamics",
		default=False,
		required=False
	)
    parser.add_argument(
		"--acs",
		action='store_true',
		help="Compute acoustics",
		default=False,
		required=False
	)
    parser.add_argument(
		"-f","--filt",
		action='store_true',
		help="Filter blade loads with the provided resonator parameters",
		default=False,
		required=False
	)
    parser.add_argument(
		"-o","--opt",
		action='store_true',
		help="Optimize resonator geometry and distribution",
		default=False,
		required=False
	)
    parser.add_argument(
		"-p","--plot",
		action='store_true',
		help="Compute acoustics",
		default=False,
		required=False
	)
    parser.add_argument(
		'-nc',"--noncompact",
		action='store_true',
		help="Include flag to blade loads from estimated surface pressure distribution (scaled nominal distribution)",
		default=False,
		required=False
	)
    
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


    args = parser.parse_args()

    geom_params,input_params,res_param,observer_params,acs_params = read_case_files(args)
    case_dir = os.path.join(os.getcwd(),input_params['case_name'])
    acs_dir = os.path.join(case_dir,'acoustics')

    if args.aero:

        saved_params = {}
        saved_params.update({'case_dir':case_dir,'acs_dir':acs_dir})

        if os.path.exists(case_dir):
            rmtree(case_dir)
        os.mkdir(case_dir)
        os.mkdir(acs_dir)

        compute_aero(geom_params,input_params,res_param,observer_params,acs_params,saved_params,args.filt)
        
        if acs_params['thicknessNoiseFlag'] or acs_params['totalNoiseFlag']:
            build_blade_geom(geom_params,input_params,res_param,observer_params,acs_params,saved_params)

        wopwop_input_configure(geom_params,input_params,res_param,observer_params,acs_params,saved_params)

        if args.filt or args.opt:
            filter_loads(geom_params,input_params,res_param,observer_params,acs_params,saved_params,opt = args.opt,noncompact=args.noncompact)
            # if args.opt:
            #     with open(args.res_param,"w") as res_file:
            #         json.dump(res_param,res_file,indent=2)

    else:
         saved_params  = read_results_from_h5(case_dir)

    if args.acs:
        cpu_count = os.cpu_count()
        if ('nbTheta' in observer_params and (observer_params['nbTheta'] >= int(cpu_count/2) or (observer_params['nbPsi'] >= int(cpu_count/2)))) or ('nbx' in observer_params and (observer_params['nbx'] >= int(cpu_count/2) or (observer_params['nby'] >= int(cpu_count/2)) or (observer_params['nbz'] >= int(cpu_count/2)))):
            parallel = True
        else:
            parallel = False 
        # parallel = False
        run_wopwop(cases = f"{input_params['case_name']}{os.sep}cases.nam",parallel = parallel)
        process_wopwop(cases_directory=case_dir,cases = 'cases.nam')
        
        if args.opt:
            acs_data = import_results_from_wopwop(cases_directory=saved_params['acs_dir'])
            oaspl = 10*np.log10(np.mean(acs_data['function_values'][:,:,:,-1]**2,axis = -2)/20e-6**2)
            saved_params.update({'oaspl':oaspl,'function_values':acs_data['function_values'],'geometry_values':acs_data['geometry_values']})

    if args.plot:
         if args.filt or args.opt:
                if 'sigma' in res_param:
                    list(map(lambda f:f(geom_params,input_params,res_param,observer_params,acs_params,saved_params), [plot_filt_load_dist,plot_res_resp]))
                else:
                    list(map(lambda f:f(geom_params,input_params,res_param,observer_params,acs_params,saved_params), [plot_filt_load_dist,plot_res_resp,plot_res_params]))
         list(map(lambda f:f(geom_params,input_params,res_param,observer_params,acs_params,saved_params), [plot_load_tseries,plot_load_dist]))

    write_results_to_h5(saved_params)
    update_res_params(args.res_param,res_param)

if __name__ == "__main__":
	main()
	print("exiting main.py")
