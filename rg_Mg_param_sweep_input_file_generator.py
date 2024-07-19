import numpy as np
import os
import sys
import json
#%%

cases_directory =os.getcwd()
param_file_name = 'param.json'
geom_file_name = 'geom.json'

N = 9
rg = np.arange(N)*(2-.25)/(N-1)+.25
Mg = np.arange(N)*(15-2)/(N-1)+2

with open(os.path.join(cases_directory,param_file_name)) as param_file:
    params = json.load(param_file)

with open(os.path.join(cases_directory,geom_file_name)) as geom_file:
    geom = json.load(geom_file)


for i,rg_iter in enumerate(rg):
    for ii,Mg_iter in enumerate(Mg):
        
        params['case_name'] = f"case_{i*len(rg)+ii}"
        params['gust_params']['core_size']  = rg_iter
        params['gust_params']['strength'] = Mg_iter

        with open(os.path.join(cases_directory,f'param_{i*len(rg)+ii}.json'),"w") as param_file:
            json.dump(params,param_file,indent=2)


with open(os.path.join(cases_directory,'run.sh'),"w") as run_file:
    run_file.write("#!/bin/bash\n")
    for i,rg_iter in enumerate(rg):
        for ii,Mt_iter in enumerate(Mg):
            run_file.write(f"rotor_gust_interaction.py --aero --acs --filt -input_geom geom.json -input_param param_{i*len(rg)+ii}.json -observer_param observer_lgrid.json -acs_param acs_param.json -res_param mdof_geom_param.json\n")
