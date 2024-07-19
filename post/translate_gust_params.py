#!/usr/bin/env python3

import os
import sys
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'src'))
from help_funcs import *
import numpy as np

#%%

case_dir = os.getcwd()
with open(os.path.join(case_dir,'param.json')) as param_file:
    input_params = json.load(param_file)

saved_params = read_results_from_h5(os.path.join(case_dir,input_params['case_name']))

max_ind = saved_params['v_gust'][-int((2*np.pi)/saved_params['dpsi']):,-1].argmax()
start_ind = saved_params['v_gust'][-int((2*np.pi)/saved_params['dpsi']):,-1][:max_ind+1].argmin()
end_ind = np.abs(saved_params['v_gust'][-int((2*np.pi)/saved_params['dpsi']):,-1][max_ind:]-.1*saved_params['v_gust'][-int((2*np.pi)/saved_params['dpsi']):,-1][max_ind]).argmin()
gust_width = (max_ind-start_ind+end_ind)*saved_params['dpsi']*saved_params['r'][-1]*saved_params['R']/saved_params['c']

Mt = saved_params['omega']*saved_params['R']/saved_params['sos']
Mg = saved_params['v_gust'][max_ind,-1]/saved_params['sos']

print(f'Mt: {np.round(Mt,3)}, Mg: {np.round(Mg,3)}, rg/c: {np.round(gust_width,3)}')