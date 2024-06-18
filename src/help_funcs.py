import json
import subprocess
import os
import sys
import h5py

sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(__file__)),'dependencies','pyWopwop'))
import wopwop

#%%


def read_case_files(args):
    out = []
    if args.input_geom is None:
        print("-input_geom or -case_file must be specified on the command line")
        exit(-1)

    if args.input_geom is not None:
        with open(args.input_geom) as geom_file:
            out.append(json.load(geom_file))

    if args.input_param is not None:
        with open(args.input_param) as param_file:
            out.append(json.load(param_file))

    if args.res_param is not None:
        with open(args.res_param) as res_file:
            out.append(json.load(res_file))

    if args.observer_param is not None:
        with open(args.observer_param) as obs_file:
            out.append(json.load(obs_file))
            
    if args.acs_param is not None:
        with open(args.acs_param) as acs_file:
            out.append(json.load(acs_file))

    return out

def run_wopwop(cases = 'cases.nam',parallel = False):
        
    print(f'Running wopwop...')
    if parallel:
            assert subprocess.run(['mpirun','wopwop3',cases],check = True), 'WOPWOP encountered an error'
    else:
        assert subprocess.run(['wopwop3',cases],check = True), 'WOPWOP encountered an error'

def process_wopwop(cases_directory,cases = 'cases.nam'):
    f1 = lambda a: wopwop.extract_wopwop_quant(case_directory=a, prefix = 'pressure')
    f2 = lambda a: wopwop.extract_wopwop_quant(case_directory=a, prefix = 'spl_spectrum')
    wopwop.apply_to_namelist([f1], cases_directory=cases_directory, cases=cases)

def import_results_from_wopwop(cases_directory):
    pred_data = {}
    #   imports reformatted data from wopwop in a dictionary
    with h5py.File(os.path.join(cases_directory, 'pressure.h5'), 'r') as dat_file:
        for k,v in dat_file.items():
            pred_data.update({k:v[()]})
    return pred_data

def write_results_to_h5(saved_params):
    with h5py.File(os.path.join(saved_params['case_dir'], 'saved_params.h5'), 'w') as f:
        for k,v in saved_params.items():
            f.create_dataset(k, data = v)

def read_results_from_h5(case_dir):
    saved_params ={}
    with h5py.File(os.path.join(case_dir, 'saved_params.h5'), 'r') as f:
        for k,v in f.items():
            saved_params.update({k:v[()]})
    return saved_params

# def process_wopwop(case_path,pressure = True,oaspl = True):
#     if pressure: 
#         with open(os.path.join(case_path,'pressure','pressure.tec')) as f:
#             p_data =np.array(re.split( ",|\n|\t|\\s",f.read())[46:])
#         p_data = p_data[p_data !='']
#         p_data = p_data.reshape(int(len(p_data)/4),4).astype(float)
        
#     if oaspl:
#         with open(os.path.join(case_path,'spl','OASPLdB.tec')) as f:
#             oaspl_data =np.array(re.split( ",|\n|\t|\\s",f.read()))
#         oaspl_data = oaspl_data[oaspl_data !='']
#         oaspl_data = oaspl_data[-3:].astype(float)
    
#     if pressure and oaspl:
#         return p_data, oaspl_data
#     else:
#         if pressure:
#             return p_data
#         else: 
#             return oaspl_data
