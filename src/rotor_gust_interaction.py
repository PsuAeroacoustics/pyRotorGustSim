#!/usr/bin/env python3

import argparse
import os
from shutil import rmtree
from help_funcs import *
from compute_aero import *
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
		"--plot",
		action='store_true',
		help="Compute acoustics",
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

    geom_params,input_params,observer_params,acs_params = read_case_files(args)

    saved_params = {}

    case_dir = os.path.join(os.getcwd(),input_params['case_name'])
    acs_dir = os.path.join(case_dir,'acoustics')
    saved_params = {**saved_params,**{'case_dir':case_dir,'acs_dir':acs_dir}}

    if args.aero:

        if os.path.exists(case_dir):
            rmtree(case_dir)
        os.mkdir(case_dir)
        os.mkdir(acs_dir)

        compute_aero(geom_params,input_params,observer_params,acs_params,saved_params)
        wopwop_input_configure(geom_params,input_params,observer_params,acs_params,saved_params)
    
    if args.acs:
        run_wopwop(cases = f"{input_params['case_name']}{os.sep}cases.nam",parallel = False)
        process_wopwop(cases_directory=case_dir,cases = 'cases.nam')
    
    if args.plot:
         list(map(lambda f:f(geom_params,input_params,observer_params,acs_params,saved_params), [plot_p_tseries,plot_gust_profile]))

    write_results_to_h5(saved_params)
    

if __name__ == "__main__":
	main()
	print("exiting main.py")
