import json
import subprocess
import os
import sys
import h5py
import numpy as np
import re

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

def update_res_params(file,res_param):
    with open(file,"w") as res_file:
        json.dump(res_param,res_file,indent=4)

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

def get_wopwop_ascii(case_path,pressure = True,oaspl = True):
    out = []
    if pressure: 
        with open(os.path.join(case_path,'pressure','pressure.tec')) as f:
            p_data =np.array(re.split( ",|\n|\t|\\s",f.read())[46:])
        p_data = p_data[p_data !='']
        out.append(p_data.reshape(int(len(p_data)/4),4).astype(float))
        
    if oaspl:
        with open(os.path.join(case_path,'spl','OASPLdB.tec')) as f:
            oaspl_data =np.array(re.split( ",|\n|\t|\\s",f.read()))
        oaspl_data = oaspl_data[oaspl_data !='']
        out.append(oaspl_data[-3:].astype(float))
    return out

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
            if isinstance(v,dict):
                for k1,v1 in v.items():
                    f.create_dataset(f'{k}/{k1}', data = v1)
            else:
                f.create_dataset(k, data = v)


def write_results_to_h5(saved_params):
    """
    Exports a nested dictionary to an HDF5 file.

    Args:
        data (dict): The nested dictionary to export.
        file_path (str): Path to the HDF5 file to save.
    """
    def dict_to_h5(h5_group, data):
        for key, value in data.items():
            if isinstance(value, dict):
                subgroup = h5_group.create_group(key)
                dict_to_h5(subgroup, value)
            else:
                h5_group.create_dataset(key, data=value)

    with h5py.File(os.path.join(saved_params['case_dir'], 'saved_params.h5'), 'w') as f:
        dict_to_h5(f, saved_params)


def read_results_from_h5(case_dir):
    
    def h5_to_dict(h5_obj):
        """
        Recursively converts an HDF5 file/group into a nested dictionary.

        Args:
            h5_obj (h5py.File or h5py.Group): HDF5 file or group object.

        Returns:
            dict: Nested dictionary representation of the HDF5 structure.
        """
        h5_dict = {}
        for key,value in h5_obj.items():
            if isinstance(value, h5py.Group):
                # Recursively process groups
                h5_dict.update({key:h5_to_dict(value)})
            else:
                if isinstance(value[()], bytes):
                    h5_dict.update({key:value[()].decode()})
                else:
                    h5_dict.update({key:value[()]})
        return h5_dict

    with h5py.File(os.path.join(case_dir, 'saved_params.h5'), 'r') as f:
        saved_params = h5_to_dict(f)
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
