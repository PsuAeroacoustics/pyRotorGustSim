
from wopwop_input_generator import *

#%%
def wopwop_input_configure(geom_params,input_params,res_param,observer_params,acs_params,saved_params):

    # initializes namelist and environmental object
    nml = []
    
    environment_in = EnvironmentIn(**acs_params)
    
    environment_const = EnvironmentConstants()
    
    # observer_params.update({'Title':'mic array','attachedTo':'aircraft','tMin':saved_params['t'][-int(input_params['computational_params']['number_of_revs']*2*np.pi/saved_params['dpsi'])],'tMax':saved_params['t'][-1]})
    observer_params.update({'Title':'mic array','attachedTo':'aircraft','tMin':saved_params['t'][-int(input_params['computational_params']['number_of_revs']*2*np.pi/saved_params['dpsi'])],'tMax':saved_params['t'][-1]})

    # observer_params.update({'lowPassFrequency':saved_params['omega']/(2*np.pi)*40})
    # observer_params.update({'highPassFrequency':saved_params['omega']/(2*np.pi)*6})

    #   Determines whether to write out a observer file or not based on what is provided in the corresponding JSON file
    if not 'nbTheta' in observer_params or 'nbx' in observer_params or 'xLoc' in observer_params:

        # np.cumsum(np.abs(np.diff(np.sin(2*np.arange(21)*np.pi/20)*30)))+120 - 2sin distribution
        radius = np.array(observer_params['radius'])
        theta = np.array(observer_params['theta'])*np.pi/180
        phi = observer_params['phi']*np.pi/180

        if isinstance(observer_params['radius'],list):
            observer_coordinates = np.array([radius*np.cos(phi)*np.cos(theta),radius*np.cos(phi)*np.sin(theta),radius*np.sin(phi)]).T
        elif isinstance(observer_params['theta'],list):
            observer_coordinates = np.array([radius*np.cos(theta),radius*np.sin(theta),radius*np.sin(phi)*np.ones(len(theta))]).T
        elif isinstance(observer_params['phi'],list):
            observer_coordinates = np.array([radius*np.cos(theta)*np.ones(len(phi)),radius*np.sin(theta)*np.ones(len(phi)),radius*np.sin(phi)]).T

        write_observer_file(os.path.join(saved_params['acs_dir'],f'observer.ascii'),observer=observer_coordinates)
        observer_params.update({'fileName':'observer.ascii'})

    if 'nbTheta' in observer_params:
        observer_params.update({'thetaMin':observer_params['thetaMin']*np.pi/180,'thetaMax':observer_params['thetaMax']*np.pi/180,'psiMin':observer_params['psiMin']*np.pi/180,'psiMax':observer_params['psiMax']*np.pi/180})
    
    observer_in = ObserverIn(**observer_params)

    
    aircraft_container = ContainerIn(Title='aircraft',nbContainer = 1)

    if acs_params['thicknessNoiseFlag'] or acs_params['totalNoiseFlag']:
        rotor_container = ContainerIn(Title='rotor',nbContainer = 2*geom_params['number_of_blades'],nbBase=1)
    else:
        rotor_container = ContainerIn(Title='rotor',nbContainer = geom_params['number_of_blades'],nbBase=1)

    rotor_cb = CB(Title='rotation',Rotation = True,AngleType='KnownFunction',Omega=input_params['flight_params']['omega'],AxisValue=[0,0,1])
    nml.extend([environment_in,environment_const,observer_in,aircraft_container,rotor_container,rotor_cb])

    for b_iter in range(geom_params['number_of_blades']):
        nml.append(ContainerIn(Title=f'blade {b_iter} loading',nbBase=1,patchGeometryFile = f'lifting_line_geometry.dat',patchLoadingFile = f'loading_blade_{b_iter}.dat',dtau = saved_params['t'][-1]/(observer_params['nt']-1)))
        nml.append(CB(Title=f'blade {b_iter} azimuthal offset',AxisValue=[0,0,1],AngleValue=2*np.pi/geom_params['number_of_blades']*b_iter))
        # nml.append(CB(Title=f'align blade {b_iter} with rear of rotor disk',AxisValue=[0,0,1],AngleValue=-np.pi))

        if acs_params['thicknessNoiseFlag'] or acs_params['totalNoiseFlag']:
            nml.append(ContainerIn(Title=f'blade {b_iter} thickness',nbBase=1,patchGeometryFile = f'blade_geometry.dat',dtau = saved_params['t'][-1]/(observer_params['nt']-1)))
            nml.append(CB(Title=f'blade {b_iter} azimuthal offset',AxisValue=[0,0,1],AngleValue=2*np.pi/geom_params['number_of_blades']*b_iter))
            # nml.append(CB(Title=f'align blade {b_iter} with rear of rotor disk',AxisValue=[0,0,1],AngleValue=-np.pi))
            # nml.append(CB(Title=f'blade {b_iter} th0',AxisValue=[1,0,0],AngleValue=saved_params['th0']))

    write_nml_file(os.path.join(saved_params['acs_dir'],f"{input_params['case_name']}.nam"),nml)
    case = caseName(globalFolderName = saved_params['acs_dir'],caseNameFile=f"{input_params['case_name']}.nam")
    write_nml_file(os.path.join(saved_params['case_dir'],'cases.nam'),[case])
    constant_compact_geometry_write(os.path.join(saved_params['acs_dir'],f'lifting_line_geometry.dat'),nodes=saved_params['lifting_line_nodes'],norms=saved_params['lifting_line_norms'],ascii=False)
    
    if acs_params['thicknessNoiseFlag'] or acs_params['totalNoiseFlag']:
        constant_compact_geometry_write(os.path.join(saved_params['acs_dir'],f'blade_geometry.dat'),nodes=saved_params['blade_nodes'],norms=saved_params['blade_norms'],ascii=False)

    for b_iter in range(geom_params['number_of_blades']):
        aperiodic_compact_loading_write(os.path.join(saved_params['acs_dir'],f'loading_blade_{b_iter}.dat'),t = saved_params['t'], loads = saved_params['loads'],ascii = False)

