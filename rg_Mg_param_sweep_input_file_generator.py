#!/usr/bin/env python3

import numpy as np
import os
import sys
import json
import argparse

#%%

parser = argparse.ArgumentParser("rotor_gust_interaction",description='Simulates a gust interacting with a hovering rotor, only the positive half of the gust profile is considered.')
parser.add_argument(
    "-input_geom",
    type= str,
    required=False,
    default='geom.json'
)
parser.add_argument(
    "-input_param",
    type= str,
    required=False,
    default='param.json'
)
parser.add_argument(
    "--filt",
    action='store_true',
    help="filter?",
    default=False,
    required=False
)
args = parser.parse_args()


with open(args.input_param) as param_file:
    params = json.load(param_file)

with open(args.input_geom) as geom_file:
    geom = json.load(geom_file)

# number of points in the sweep
N = 9
# gust peak location range as a percentage of the rotor chord
rg = np.arange(N)*(0.85-.02)/(N-1)+.02
# gust strength range as a percentage of the reference gust strength as specified in the params fil
Mg = np.arange(N)*(3- 0.5)/(N-1)+0.5

ref_strength_pk_ratio = params['gust_params']['strength']/params['gust_params']['peak_location']

for i,rg_iter in enumerate(rg):
    for ii,Mg_iter in enumerate(Mg):
        
        params['case_name'] = f"case_{i*len(rg)+ii}"
        params['gust_params']['peak_location']  = rg_iter
        params['gust_params']['strength'] = Mg_iter*ref_strength_pk_ratio*rg_iter

        with open(f'param_{i*len(rg)+ii}.json',"w") as param_file:
            json.dump(params,param_file,indent=2)

with open('run.sh',"w") as run_file:
    run_file.write("#!/bin/bash\n")
    for i,rg_iter in enumerate(rg):
        for ii,Mt_iter in enumerate(Mg):
            if args.filt:
                run_file.write(f"rotor_gust_interaction.py --aero --acs --filt -input_geom geom.json -input_param param_{i*len(rg)+ii}.json -observer_param observer_lgrid.json -acs_param acs_param.json -res_param sdof_geom_param.json\n")
            else:
                run_file.write(f"rotor_gust_interaction.py --aero --acs -input_geom geom.json -input_param param_{i*len(rg)+ii}.json -observer_param observer_lgrid.json -acs_param acs_param.json -res_param sdof_geom_param.json\n")
