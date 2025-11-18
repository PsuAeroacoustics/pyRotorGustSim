#!/usr/bin/env python3

import numpy as np
import os
import sys
import json
#%%

cases_directory =os.getcwd()
param_file_name = 'param.json'
geom_file_name = 'geom.json'

N = 9
Mt = np.arange(N)*(.7-.3)/(N-1)+.3
print(Mt)
with open(os.path.join(cases_directory,param_file_name)) as param_file:
    params = json.load(param_file)

with open(os.path.join(cases_directory,geom_file_name)) as geom_file:
    geom = json.load(geom_file)


for i,Mt_iter in enumerate(Mt):
    
    params['case_name'] = f"case_{i}"
    params['gust_params']['strength'] = params['gust_params']['strength']*(params['flight_params']['omega']/(Mt_iter*params['flight_params']['sos']/geom['radius']))
    params['flight_params']['omega'] = Mt_iter*params['flight_params']['sos']/geom['radius']

    with open(os.path.join(cases_directory,f'param_{i}.json'),"w") as param_file:
        json.dump(params,param_file,indent=2)


with open(os.path.join(cases_directory,'run.sh'),"w") as run_file:
    run_file.write("#!/bin/bash\n")
    for i,Mt_iter in enumerate(Mt):
        run_file.write(f"rotor_gust_interaction.py --aero --acs -input_geom geom.json -input_param param_{i}.json -observer_param observer_lgrid.json -acs_param acs_param.json -res_param mdof_dist_param.json\n")
