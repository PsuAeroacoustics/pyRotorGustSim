import numpy as np
import aerosandbox as asb
from NodeCenteredNorms import *
#%%

def build_blade_geom(geom_params,input_params,res_param,observer_params,acs_params,saved_params):
    
    
    af = asb.Airfoil(geom_params['airfoil'])
    af.coordinates = af.repanel(n_points_per_side = int(input_params['computational_params']['airfoil_elements']/2)).coordinates[...,None]*saved_params['c']
    af.coordinates[:,0] = -af.coordinates[:,0]
    af.coordinates[:,0] = af.coordinates[:,0]+0.25*saved_params['c']

    pnts_per_sections = len(af.coordinates)
    n_sections = saved_params['N_elements']

    blade_nodes = np.zeros((n_sections,pnts_per_sections,3))
    blade_nodes[:,:,1] = af.coordinates[:,0].T
    blade_nodes[:,:,0] = (np.expand_dims(saved_params['r']*saved_params['R'],axis = -1)*np.ones(pnts_per_sections))
    blade_nodes[:,:,-1] = af.coordinates[:,-1].T

    blade_nodes = blade_nodes.reshape((n_sections*pnts_per_sections),3,order = 'F')
    blade_norms = NodeCenteredNorms(blade_nodes,pntsPerXsec =pnts_per_sections ,nXsecs = n_sections)
    blade_nodes = blade_nodes.reshape((pnts_per_sections,n_sections,3),order = 'F')
    blade_norms = blade_norms.reshape((pnts_per_sections,n_sections,3),order = 'F')

    saved_params.update({'blade_nodes':blade_nodes,'blade_norms':blade_norms})

    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # # ax.auto_scale_xyz([-2, 2], [10, 60], [-1, 1])
    # ax.pbaspect = [.09, 1, .05]
    # ax.set(xlabel = 'x',ylabel = 'y',zlabel = 'z')
    # ax.scatter(blade_nodes[:,:,0], blade_nodes[:,:,1], blade_nodes[:,:,2], c='red', linewidths=.2)
    # ax.quiver(blade_nodes[:,:,0], blade_nodes[:,:,1], blade_nodes[:,:,2], blade_norms[:,:,0], blade_norms[:,:,1], blade_norms[:,:,2], length=0.0005)
    # fig.show()
