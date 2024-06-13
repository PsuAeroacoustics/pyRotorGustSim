import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
import f90nml
# sys.path.insert(0,os.path.join(os.path.dirname(__file__),'dependencies','pyWopwop'))
# import wopwop

#%%
import matplotlib.colors as mcolors

fontName = 'Times New Roman'
fontSize = 16
plt.rc('font',**{'family':'serif','serif':[fontName],'size':fontSize})
plt.rc('mathtext',**{'default':'regular'})
plt.rc('text',**{'usetex':False})
plt.rc('lines',**{'linewidth':2})

#%%

val_dir = os.path.dirname(__file__)
val_files =['M360_3D_M17.csv','M360_3D_M14.csv','M360_3D_M15.csv','M360_3D_M16.csv']
v_data = [np.loadtxt(os.path.join(val_dir,file),delimiter=",") for file in val_files]
v_data = [m[np.argsort(m[:,0])] for m in v_data]

#%%

cases_directory = os.path.join(os.path.dirname(__file__),'..','..')
cases = 'cases.nam'

#%%
case_info = {}
namelist = f90nml.read(os.path.join(cases_directory,  cases))
for k,v in namelist.items():
    if isinstance(v,dict):
        case_info_temp = {}
        for k_temp,v_temp in v.items():
            case_info_temp = {**case_info_temp, **{k_temp: v_temp}}
        case_info = {**case_info,**case_info_temp}
    else:
        case_info = {**case_info,**{k:v}}

#%%

pred_data = {}
#   imports reformatted data from wopwop in a dictionary
with h5py.File(os.path.join(case_info['globalfoldername'], 'pressure.h5'), 'r') as dat_file:
    for i,n in dat_file.items():
        pred_data = {**pred_data,**{i:n[()]}}

pred_data['geometry_values'] = np.flip(pred_data['geometry_values'],axis = 0).squeeze()
pred_data['function_values'] = np.flip(pred_data['function_values'],axis = 0).squeeze()

#%%

fig,ax = plt.subplots(4,1, figsize = (4.5,4.5))
for mic_iter in range(len(pred_data['function_values'])):
    ax[mic_iter].plot(pred_data['function_values'][mic_iter,:,0]/pred_data['function_values'][mic_iter,-1,0],pred_data['function_values'][mic_iter,:,1])
    ax[mic_iter].plot(pred_data['function_values'][mic_iter,:,0]/pred_data['function_values'][mic_iter,-1,0],pred_data['function_values'][mic_iter,:,2])
    ax[mic_iter].plot(pred_data['function_values'][mic_iter,:,0]/pred_data['function_values'][mic_iter,-1,0],pred_data['function_values'][mic_iter,:,-1])

    # ax[mic_iter].plot(v_data[mic_iter][:,0]/2,v_data[mic_iter][:,-1]-np.mean(v_data[mic_iter][:,-1]))
    
    # ax[mic_iter].set_xlim([0,.5])
    # ax[mic_iter].set_ylim([-2,2])
ax[1].legend(['Thickness','Loading','Total','Measured'])

ph_shift = -.48

for mic_iter in range(len(pred_data['function_values'])):
    fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
    ax.plot(pred_data['function_values'][mic_iter,:,0]/pred_data['function_values'][mic_iter,-1,0],pred_data['function_values'][mic_iter,:,-1]-np.mean(pred_data['function_values'][mic_iter,:,-1]))
    ax.plot(v_data[mic_iter][:,0]/2-ph_shift,v_data[mic_iter][:,-1]-np.mean(v_data[mic_iter][:,-1]))
    ax.grid()
    # ax[mic_iter].set_xlim([0,.5])
    ax.set_ylim([-10,10])
    ax.legend(['Predicted','Measured'])
